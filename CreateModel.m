% CREATEMODEL Creates a model structure from parameters.
%   Model = CREATEMODEL(lambda,joints,axes,d,r,mass,inertia,...) 
%   creates a model structure with the following inputs:
%
% 	 lambda:    connectivity vector     
%               lamda{i} gives parent of body i
%    joints:    type of joint
%               joints{i} gives a instance of the jointtype class
%    axes:      rotation or translation axis 
%               axes{i} gives an 'x', 'y' or 'z' string
% 	 d:         offsets between the joints as seen from parent body frame
%               d(i,:) gives the x y and z coordinates
% 	 r:         center of mass offsets as seen from the body frame
%               r(i,:) gives the x y and z coordinates
% 	 mass:      mass of the rigid bodies
%               m(i) gives body i's mass
% 	 inertia:   inertias of the rigid bodies along principle axes
%               inertia(i,:) gives the x y and z components
%
%   The final input is optional:
%
%   name:       model name
%               string
%   export:     export the model to a .mat-file or not 
%               boolean, default true
%   exportloc:  export location
%               string
%   modfnc:     custom modification of the model structure
%               handle to Model = modfnc(Model)
%   base:       draw fixed base, size depends on first joint offset
%               boolean, default false
%
%   For example:
%
%   Model = CreateModel(...
%     {0},jointtype.Revolute,'z',[0 0 0],[0 0 0.5],1,[0.1 0.1 0.1],...
%     'name',UltraBot,'modfnc',@add_visual);
%
%   Creates a simple one degree of freedom robot, named Ultrabot, with a
%   single revolute joint of which the joint frame is located at the global
%   frame (d = [0,0,0]). The function @add_visual could look like:
%
%   function Model = add_visual(Model)
%       Model.rigidbody(1).visual.vertices = load_vertices();
%       Model.rigidbody(1).visual.faces = load_faces();
%   end
%
% Made by Erik Vlasblom
% Last modified: 03-02-2015
function varargout = CreateModel(lambda,joints,axes,d,r,mass,inertia, varargin)

% -------- OPTIONS ---------------------------------
options = struct(...
    'name','Some_Robot',...
    'export',true,...
    'exportloc',[],...
    'modfnc',[],...
    'base',false);
onames = fieldnames(options);
if ~isempty(varargin)
    nvarargin = length(varargin);
    if mod(nvarargin,2) == 0
        for ii = 1:2:nvarargin
            if ~isempty(find(strcmp(onames,varargin{ii}),1))
                options.(varargin{ii}) = varargin{ii+1};
            else
                error(['Input ',varargin{ii},' is not a valid option']);
            end
        end
    else
        error('Please specify all additional input')
    end
end

tic
disp('CREATE MODEL')

% --- Initialize
n = length(lambda);
Model.name = options.name;
Model.dof = n;
Model.gravity = [0 0 -9.81];


% --- Connectivity
kappa = cell(size(lambda));
mu = cell(size(lambda));
for ii = 1:n; kappa{ii} = GetPathToRoot(lambda,ii); end;
for ii = 1:n; mu{ii} = GetChildren(lambda,ii); end;
G = GetAdjacency(lambda);

Model.connectivity.lambda = lambda;
Model.connectivity.kappa = kappa;
Model.connectivity.mu = mu;
Model.connectivity.G = G;


% --- Rigid bodies
Model.rigidbody = struct;
for ii = 1:n      
    Model.rigidbody(ii).inertial.mass = mass(ii);
    Model.rigidbody(ii).inertial.i = inertia(ii,:);
    Model.rigidbody(ii).inertial.origin = r(ii,:);
    %Model.rigidbody(ii).inertial.transform = Hm;
    
    Model.rigidbody(ii).joint.type = joints{ii};
    Model.rigidbody(ii).joint.axis = axes{ii};
    Model.rigidbody(ii).joint.offset = d(ii,:);
    %Model.rigidbody(ii).joint.transform = H;
    
    Model.rigidbody(ii).visual.bodysizemax = [];
    Model.rigidbody(ii).visual.bodysizemin = [];
    Model.rigidbody(ii).visual.vertices = [];
    Model.rigidbody(ii).visual.faces = [];
    
    Model.rigidbody(ii).collision = [];
end
if options.base
    Model.base.visual.bodysizemax = [];
    Model.base.visual.bodysizemin = [];
    Model.base.visual.vertices = [];
    Model.base.visual.faces = [];
    Model.base.collision = [];
end

% --- Visual
Model.visual.bodywidth = 0.02;

% Modification function
if ~isempty(options.modfnc); Model = options.modfnc(Model); end

% --- Export
exploc = options.exportloc;
if ~isempty(exploc); 
    if ~strcmp(exploc(end),'\'); exploc = [options.exportloc,'\']; end;
end
if options.export 
    tmp.(Model.name) = Model;
    save([exploc,options.name,'.mat'],'-struct','tmp')
end


% --- Output
if nargout == 0; assignin('base',Model.name,Model);
elseif nargout == 1; varargout{1} = Model;
else error('Wrong amount of outputs')
end

fprintf('  ')
toc
fprintf('\n');


end