
clear all
clc

% Floating humanoid robot TUlip undergoing an external force.

% NOTE: for complex structures like this, it is advised to use the function
% rbde2spatial to change model definitions and use Featherstone's toolbox
% called spatial_v2 (see Upperbody example). It performs a lot faster as
% the number of degrees of freedom increases. 
% Later versions of this toolbox will have increased computational speed.

%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'TUlipFloating';

% --- Parameters
% connectivity
lambda = {0 1 2 3 4 5 6 7 8 9 10 11 6 13 14 15 16 17};

% joint types
joints = {jointtype.Prismatic, ... % floating base
    jointtype.Prismatic, ...
    jointtype.Prismatic, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ... % left leg
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ... % right leg
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute};

% joint axes
axes = {'x','y','z','z','y','x', ... % floating base
    'z','x','y','y','y','x',... % left leg
    'z','x','y','y','y','x'}; % right leg

% joint offsets, center of mass, inertial parameters
[d,r,mass,inertia] = TUlipFloating_parameters();

% --- Create model structure
% for the dynamics, the length of the last link does not matter
% only the com position counts, but for visualization, you want
% to see nice stuff, hence the modfnc. The size of the other links
% is automatically determined.
TUlipFloating = CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@TUlipFloating_edit);

t0 = 0;
x0 = zeros(length(lambda))';
axv = [-1 1 -1 1 -1.5 0.5];
Animation3D(TUlipFloating,t0,x0,'showframe',true,'showcom',true,'axis',axv)



%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------

RecursiveDynamics(TUlipFloating,...
    'simplify',true,'par_edit',false,'export',true);

point = [0 0 0]; 
RecursiveConstraints(TUlipFloating,...
    'simplify',true,'par_edit',false,'export',true,...
    'body',6,'point',point); 


%% ------------------ SIMULATE  -----------------------------------------

disp('SIMULATION')

load([modelname,'.mat']); n = TUlipFloating.dof;

% --- Simulate
x0 = [[0, 0, 0], [0, 0, 0, ...           % trunk
    10, 30, -30, 15, 0, 0, ...           % left leg
    0, 0, 0, 0, 0, 0]*pi/180 ...         % right leg
       [0, 0, 0], [0, 0, 0, ...          % velocity trunk
       0, 0, 0, 0, 0, 0, ...             % velocity left leg
       0, -10, 0, 0, 0, 0]*pi/180];      % velocity right leg
   
foot_depth = 0.05;
H18_0 = H_Transform(TUlipFloating,x0,18,0);
x0(1:3) = -Transform([0 0 -foot_depth],H18_0); % place right leg at origin
   
t_end = 5;
h = 1/250; 

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;
u = zeros(length(t),n);     
p = zeros(length(t),3)';
pdes = zeros(length(t),3)';
e = zeros(length(t),3)';
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = TUlipFloating_M_C_G(x(ii,:));              % eom matrices
    [~,J6,~] = TUlipFloatingCon6_D_dD_dconv(x(ii,:));    % trunk Jacobian
    if t(ii) <= 0.1
        Fext = [250 250 500]';
    else
        Fext = [0 0 0]';
    end
    
    % Forward Dynamics
    qn = x(ii,1:n)'; 
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-0.1*qd+u(ii,:)'+J6.'*Fext);
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    x(ii+1,:) = [q_next', qd_next'];
end
fprintf('  ')
toc
fprintf('\n');

%% ------------------ VISUALIZE  -----------------------------------------

options.ylabel = {'x (m)','y (m)','z (m)'};
options.title = {'Position of the Trunk'};
PlotStates(t,x,[1 2 3],options)

axv = [-1 1 -1 1 -1.5 0.5];
Animation3D(TUlipFloating,t,x,...
    'showcom',false,...
    'showworld',true,...
    'follow',3,...
    'axis',axv)
