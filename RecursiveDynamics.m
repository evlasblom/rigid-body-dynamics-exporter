% RECURSIVEDYNAMICS Symbolically computes the model dynamics.
%   RECURSIVEDYNAMICS(Model,...)
%   uses a model structure created by CreateModel to symbolically compute
%   the dynamic properties of that model.
%
%   name:       name for exported data. If not specified, uses Model.name 
%               string
%   export:     export requested data to .mat-files and .m-files.
%               boolean, default true
%   exportloc:  export location
%               string, default '' (cd)
%   exportdata: specify the functions to be created. The options are:
%                   - 'M' mass matrix
%                   - 'C' coriolis matrix
%                   - 'G' gravity vector
%               cell, default {'M','C','G'}
%   exportlang: specify the exported language. The options are:
%                   - 'Matlab'
%               string, default 'Matlab'
%   removezero: removes all zero expressions in the exported function. It
%               might make the exported file less easy to read, but the
%               numerical implementation will be faster.
%               boolean, default true
%   simplify:   simplify symbolic equations. When true, all atomic elements
%               are simplified. Derivation will take longer, but the final
%               numerical implementation will be faster.
%               boolean, default true
%   par_edit:   make parameters easily editable in the exported functions.
%               This is done by introducing a simple parameter structure
%               that will  function as an interface to modify the
%               parameters. This is done by supplying this structure to the 
%               exported functions (and the visualization) as for example:
%                   Mass = SomeModel_M(x,parameters);
%                   Animation3D(SomeModel,t,x,'parameters',parameters);
%               boolean, default false
%   access:     make all important internal variables global such that they
%               can be accessed outside to use them for other calculations.
%               boolean, default false
%
% Made by Erik Vlasblom
% Last modified: 21-05-2015
function RecursiveDynamics(inModel,varargin)

% ------------------ OPTIONS
options = struct(...
    'name','',...
    'export',true,...
    'exportloc',[],...
    'exportdata',cell(1,1),...
    'exportlang','Matlab',...
    'removezero',true,...
    'simplify',true,...
    'par_edit',false,...
    'access',false);
options.exportdata = {'M','C','G'};
onames = fieldnames(options);
if ~isempty(varargin)
    nvarargin = length(varargin);
    if mod(nvarargin,2) == 0
        for i = 1:2:nvarargin
            if ~isempty(find(strcmp(onames,varargin{i}),1))
                options.(varargin{i}) = varargin{i+1};
            else
                error(['Input ',varargin{i},' is not a valid option']);
            end
        end
    else
        error('Please specify all additional input')
    end
end

% set name
if isempty(options.name)
    options.name = inModel.name;
end

% set export functions
if strcmp(options.exportlang,'Matlab')
    NAME = @MTLB_name;
    ATOM = @MTLB_atom;
    GETEL = @MTLB_get_element;
    ADD = @MTLB_addition;
    SUB = @MTLB_subtraction;
    MULT = @MTLB_multiplication;
    TRNS = @MTLB_transpose;
    INVS = @MTLB_inverse;
    RBAR = @MTLB_rbar;
    OMBAR = @MTLB_ombar;
    ADJ = @MTLB_adjoint;
    IADJ = @MTLB_inverseadjoint;
    AD = @MTLB_ad;
    CORIO = @MTLB_coriolis_element;    
    WR_ATOM = @MTLB_write_atom;    
    WR_FUNCTION = @MTLB_export_all;
end

% simplify or just pass trough symbolic variables
if options.simplify; S = @simplify;
else S = @(a) a;
end

% create symbolic parameters such that they can be edited later
if options.par_edit
    [symModel, parameters] = EditParameters(inModel);
    Model = symModel;
else
    [~, parameters] = EditParameters(inModel);
    Model = inModel;
end

% give acces to variables by defining them globally
if options.access
    clearvars -global
    global q dq ddq;
    global H Hm unittwist 
    global dAd Jl J Ja dJ
    global Im I dMM Fz 
    global MM CC GG
    global list exportlist
end

% Labels:
% Relative to some other local frame: H.l21 and twists.l43
l = @(i,j) ['l_',num2str(i),'_',num2str(j)];
% Relative to global frame: H.g20
g = @(i,j) ['g_',num2str(i),'_',num2str(j)];
% Remaining body properties: mass.b1, inertia.b3, J.b6, ...
b = @(i) ['b_',num2str(i)];
% Partial derivatives: dJ.b6 / dq_6
pd = @(i,j) ['pd_',num2str(i),'_',num2str(j)];
% Elements: CC.e_i_j
el = @(i,j) ['el_',num2str(i),'_',num2str(j)];


% ------------------ INITIALIZATION

% Parameters
n = Model.dof;
gravity = Model.gravity;

% Degrees of freedom
q = sym('q',[n 1]);
dq = sym('dq',[n 1]);
ddq = sym('ddq',[n 1]);

% Connectivity: as done in [4] and [5].
p = Model.connectivity.lambda;

% ZEROS
H.(g(0,0)) = ATOM('H',[0 0],HomogeneousTransform([0,0,0]));
J.(b(0)) = zero_atom(ATOM('J',[0 0],zeros(6,n)));
for k = 1:n
    dJ.(pd(0,k)) = zero_atom(ATOM('dJ',[k 0],zeros(6,n))); 
end
for k = 1:n
    dMM.(pd(0,k)) = zero_atom(ATOM('dMM',[k 0],zeros(n,n)));
end
MM.(b(0)) = zero_atom(ATOM('MM',0,zeros(n,n)));
GG.(b(0)) = zero_atom(ATOM('GG',0,zeros(n,1)));

list{1} = H.(g(0,0));
list{2} = J.(b(0));
for k = 1:n list{end+1} = dJ.(pd(0,k)); end
for k = 1:n list{end+1} = dMM.(pd(0,k)); end
list{end+1} = MM.(b(0));
list{end+1} = GG.(b(0));


% ------------------ RECURSIVE DYNAMICS 
% The dynamics are derived using a geometric approach. It is based on
% methods used in [1] and [2].
tic
disp('RECURSIVE DYNAMICS')

B = Model.rigidbody;
for i = 1:n
    fprintf(['  - body ',num2str(i),'\n']);
    
    % --- KINEMATICS 
    % Local transformations (recursion)    
    hloc = CreateTransform(B(i).joint.type,B(i).joint.axis,B(i).joint.offset,q(i));
    H.(l(i,p{i})) = ATOM('H',[i,p{i}],hloc);
    list{end+1} = H.(l(i,p{i}));
    
    hm = HomogeneousTransform(B(i).inertial.origin);
    Hm.(b(i)) = ATOM('Hm',i,hm);
    list{end+1} = Hm.(b(i));

    % Local unit twists, expressed in body coordinate frame  (recursion)
    utw = simplify(... % always simplify
        UnTilde(TwistFromHomogeneous(H.(l(i,p{i})).exp,q(i),'b')));
    unittwist.(l(i,p{i})) = ATOM('unittwist',[i,p{i}],utw);
   
    % Global transformation (recursion)
    H.(g(i,0)) = MULT('H',[i 0],H.(g(p{i},0)),H.(l(i,p{i})));
    list{end+1} = H.(g(i,0));
 
    % Jacobians, expressed in body coordinate frame (recursion) 
    Jtmp = zeros(6,n); Jtmp(:,i) = unittwist.(l(i,p{i})).exp;
    Jl.(b(i)) = ATOM('Jl',i,Jtmp);
    J.(b(i)) = ADD('J',i,...
        MULT([],[],IADJ(H.(l(i,p{i}))),J.(b(p{i}))),...
        Jl.(b(i)));
    Ja.(b(i)) = MULT('Ja',i,... 
        ADJ(RBAR(H.(g(i,0)))),IADJ(Hm.(b(i))),J.(b(i)));
    list{end+1} = Jl.(b(i));
    list{end+1} = J.(b(i));
    list{end+1} = Ja.(b(i));
    
    % --- DYNAMICS
    % Mass matrix (summation)
    Im.(b(i)) = ATOM('Im',i,...
        diag([B(i).inertial.mass*ones(1,3),B(i).inertial.i]));
    I.(b(i)) = MULT('I',i,...
        TRNS(IADJ(Hm.(b(i)))),Im.(b(i)),IADJ(Hm.(b(i))));
    MM.(b(i)) = ADD('MM',i,...
        MM.(b(i-1)),...
        MULT([],[],TRNS(J.(b(i))),I.(b(i)),J.(b(i))));
    list{end+1} = Im.(b(i));
    list{end+1} = I.(b(i));
    list{end+1} = MM.(b(i));
    
    for k = 1:n
        % Adjoint derivatives
        if i == k
            dd = S(... % optional simplify
                -Ad(unittwist.(l(k,p{k})).exp)*InverseAdjoint(H.(l(i,p{i})).exp));
            dAd.(b(i)) = ATOM('dAd',i,dd);
            list{end+1} = dAd.(b(i));
        end
        
        % Jacobian derivative, expressed in body coordinate frame (recursion)
        if i == k
            dJ.(pd(i,k)) = ADD('dJ',[i k],...
                MULT([],[],dAd.(b(i)),J.(b(p{i}))),...
                MULT([],[],IADJ(H.(l(i,p{i}))),dJ.(pd(p{i},k))));
        else
            dJ.(pd(i,k)) = MULT('dJ',[i k],...
                IADJ(H.(l(i,p{i}))),dJ.(pd(p{i},k)));
        end
        list{end+1} = dJ.(pd(i,k));
        
       % Mass matrix derivatives (summation)
       dMM.(pd(i,k)) = ADD('dMM',[i k],...
           dMM.(pd(i-1,k)),...
           MULT([],[],TRNS(dJ.(pd(i,k))),I.(b(i)),J.(b(i))),...
           MULT([],[],TRNS(J.(b(i))),I.(b(i)),dJ.(pd(i,k))));
       list{end+1} = dMM.(pd(i,k));
    end
    
    % Gravity wrench (summation) 
    Fz.(b(i)) = ATOM('Fz',i,...
        [B(i).inertial.mass*-gravity 0 0 0].');
    GG.(b(i)) = ADD('GG',i,...
        GG.(b(i-1)),MULT([],[],TRNS(Ja.(b(i))),Fz.(b(i))));
    list{end+1} = Fz.(b(i));
    list{end+1} = GG.(b(i));

    fprintf('  ')
    toc
    fprintf('\n');
end

% Coriolis elements
fprintf('  - coriolis matrix \n');
CCname = 'CC';
dMlist = cell(1,n);
for k = 1:n; dMlist{k} = dMM.(pd(n,k)); end
for i = 1:n
    for j = 1:n
        CC.(el(i,j)) = CORIO(CCname,[i j],dMlist,dq,n); 
        list{end+1} = CC.(el(i,j));
    end
end
fprintf('  ')
toc
fprintf('\n');

% --- OUTPUT
for ii = 1:length(options.exportdata)
    data = options.exportdata{ii};
    switch data
        case 'M'
            list{end+1} = struct('caller',data,'exp',MM.(b(n)).caller,'zero',false);
        case 'C'
            list{end+1} = struct('caller',data,'exp',CCname,'zero',false);
        case 'G'
            list{end+1} = struct('caller',data,'exp',GG.(b(n)).caller,'zero',false);
    end
end

if options.removezero
    exportlist = cell(0);
    curr = length(list);
    for ii = 1:curr
        if ~list{ii}.zero
            exportlist{end+1} = list{ii};
        end
    end
else
    exportlist = list;
end

% ------------------ EXPORT DATA 
% --- Export settings
if ~isempty(options.exportloc); 
    if ~strcmp(exploc(end),'\'); options.exportloc = [options.exportloc,'\']; end;
end

% --- Export parameters
save([options.exportloc,options.name,'_parameters.mat'],'parameters');

% --- Export function
if options.export
    tic
    disp(['EXPORT TO ',upper(options.exportlang),' FUNCTION'])
    
    WR_FUNCTION(options,exportlist,[q; dq],parameters);
    
    fprintf('  ')
    toc
    fprintf('\n');
end


J = rmfield(J,b(0));
for k = 1:n; dJ = rmfield(dJ,pd(0,k)); end
for k = 1:n; dMM = rmfield(dMM,pd(0,k)); end
MM = rmfield(MM,b(0));
GG = rmfield(GG,b(0));


end


% ------------------ REFERENCES------------------------------------------

% [1] Stramigioli, Stefano, and Herman Bruyninckx. 
%     Geometry and screw theory for robotics.
%     Tutorial during ICRA 2001 (2001).

% [2] Murray, Richard M., et al. 
%     A mathematical introduction to robotic manipulation. 
%     CRC press, 1994.

% [3] Spong, Mark W., Seth Hutchinson, and Mathukumalli Vidyasagar. 
%     Robot modeling and control. Vol. 3. 
%     New York: Wiley, 2006.

% [4] Featherstone, Roy. 
%     A Beginner's Guide to 6-D Vectors (Part 1).
%     Robotics & Automation Magazine, IEEE 17.3 (2010): 83-94.

% [5] Featherstone, Roy. 
%     A Beginner's Guide to 6-D Vectors (Part 2)[Tutorial].
%     Robotics & Automation Magazine, IEEE 17.4 (2010): 88-99.
