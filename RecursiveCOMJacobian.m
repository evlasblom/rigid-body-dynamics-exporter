% RECURSIVECOMJACOBIAN Symbolically computes the model dynamics.
%   RECURSIVECOMJACOBIAN(Model,...)
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
%                   - 'Jc' COM Jacobian
%               cell, default {'Jc'}
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
% Last modified: 21-06-2015
function RecursiveCOMJacobian(inModel,varargin)

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
options.exportdata = {'Jc'};
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
    global dAd Jl J Ja Jcomt Jcom dJ
    global Im I dMM Fz totalmass mass
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
Jcomt.(b(0)) = zero_atom(ATOM('Jcom',[0 0],zeros(6,n)));
Jcom.(b(0)) = zero_atom(ATOM('Jcom',[0 0],zeros(3,n)));
Jcom.(g(0,0)) = zero_atom(ATOM('Jcom',[0 0],zeros(3,n)));
totalmass.(b(0)) = zero_atom(ATOM('totalmass',[0],0));

list{1} = H.(g(0,0));
list{2} = J.(b(0));
list{3} = Jcomt.(b(0));
list{4} = Jcom.(b(0));
list{5} = Jcom.(g(0,0));
list{6} = totalmass.(b(0));



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
    list{end+1} = unittwist.(l(i,p{i}));
   
    % Global transformation (recursion)
    H.(g(i,0)) = MULT('H',[i 0],H.(g(p{i},0)),H.(l(i,p{i})));
    list{end+1} = H.(g(i,0));
 
    % Jacobians, expressed in body coordinate frame (recursion) 
    Jtmp = zeros(6,n); Jtmp(:,i) = unittwist.(l(i,p{i})).exp;
    Jl.(b(i)) = ATOM('Jl',i,Jtmp);
    J.(b(i)) = ADD('J',i,...
        MULT([],[],IADJ(H.(l(i,p{i}))),J.(b(p{i}))),...
        Jl.(b(i)));
    list{end+1} = Jl.(b(i));
    list{end+1} = J.(b(i));
    
    % Center of Mass Jacobian
    mass.(b(i)) = ATOM('mass',i,B(i).inertial.mass);
    Jcomt.(b(i)) = MULT('Jcomt',i,...
        ADJ(RBAR(H.(g(i,0)))),IADJ(Hm.(b(i))),J.(b(i)));
    Jcom.(b(i)) = GETEL('Jcomi',{[1:3],[1:n]},Jcomt.(b(i)));
    Jcom.(g(i,0)) = ADD('Jcom',i,...
        Jcom.(g(p{i},0)),MULT([],[],mass.(b(i)),Jcom.(b(i))));
    totalmass.(b(i)) = ADD('totalmass',i,...
        totalmass.(b(p{i})),mass.(b(i)));
    list{end+1} = mass.(b(i));
    list{end+1} = Jcomt.(b(i));
    list{end+1} = Jcom.(b(i));
    list{end+1} = Jcom.(g(i,0));    
    list{end+1} = totalmass.(b(i));

    fprintf('  ')
    toc
    fprintf('\n');
end

% Center of Mass Jacobian
Jcom.final = MULT('Jcom_final',n,...
    INVS(totalmass.(b(i))),Jcom.(g(i,0)));
list{end+1} = Jcom.final;


% --- OUTPUT
for ii = 1:length(options.exportdata)
    data = options.exportdata{ii};
    switch data
        case 'Jc'
            list{end+1} = struct('caller',data,'exp',Jcom.final.caller,'zero',false);
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
Jcomt = rmfield(Jcomt,b(0));
Jcom = rmfield(Jcom,b(0));
Jcom = rmfield(Jcom,g(0,0));
totalmass = rmfield(totalmass,b(0));


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
