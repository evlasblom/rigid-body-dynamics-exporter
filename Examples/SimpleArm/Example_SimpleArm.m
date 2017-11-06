
clear all
clc

% Simulation of a simple robot. By changing the derivation settings, the
% parameters can be modified in the exported files. This means that as long
% as the structure stays the same, but when parameters change, it is not
% necessary to redo the derivation.

%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'SimpleArm';

% --- Parameters
% connectivity
lambda = {0 1};

% joint types
joints = {jointtype.Revolute, ...
    jointtype.Revolute};

% joint axes
axes = {'y','y'};
    
% joint offsets
d = [0,     0,      0.1; 
    0,      0,      0.2];

% center of mass
r = [0,     0,      0.1; 
    0,      0,      0.1];

% inertial properties
mass = [0.5, 0.5];
inertia = [0.02, 0.02, 0.003; 
    0.02, 0.02, 0.003];

% --- Create model structure
% for the dynamics, the length of the last link does not matter
% only the com position counts, but for visualization, you want
% to see nice stuff, hence the modfnc. The size of the other links
% is automatically determined.
SimpleArm = CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@SimpleArm_edit);

t0 = 0;
x0 = zeros(length(lambda))';

axv = [-0.5 0.5 -0.5 0.5 -0.5 0.5];
Animation3D(SimpleArm,t0,x0,'showframe',true,'showcom',true,'axis',axv)


%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------


RecursiveDynamics(SimpleArm,'par_edit',true)


%% ------------------ SIMULATE  -----------------------------------------


disp('SIMULATION')

modelname = 'SimpleArm';
load([modelname,'.mat']); n = SimpleArm.dof;
load([modelname,'_parameters.mat']); p = parameters;

% --- Simulation
x0 = [[120 0]*pi/180 ...
       [0, 0]*pi/180];
t_end = 2;
h = 1/100; 

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;       
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = SimpleArm_M_C_G(x(ii,:),p); 
    
    % Forward Dynamics
    qn = x(ii,1:n)'; 
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-0.025*qd);
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    x(ii+1,:) = [q_next', qd_next'];
end
fprintf('  ')
toc
fprintf('\n');


% --- Modify parameters
pmod = p;
pmod.m1   = p.m1*2;
pmod.d2_z = p.d2_z*2;


% --- Modified simulation
x0 = [[120 0]*pi/180 ...
       [0, 0]*pi/180];
t_end = 2;
h = 1/100; 

tic
t = 0:h:t_end;
xm = zeros(length(t),length(x0)); xm(1,:) = x0;       
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = SimpleArm_M_C_G(xm(ii,:),pmod); 
    
    % Forward Dynamics
    qn = xm(ii,1:n)'; 
    qd = xm(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-0.025*qd);
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    xm(ii+1,:) = [q_next', qd_next'];
end
fprintf('  ')
toc
fprintf('\n');


%% ------------------ VISUALIZE  -----------------------------------------

options.ylabel = {'\phi_1 (rad)','\phi_2 (rad)'};
options.legend = {'Normal','Modified'};
options.title = {'Simulation with Modified Parameters'};
PlotStates(t,x,t,xm,[1 2],options)

axv = [-0.5 0.5 -0.5 0.5 -0.5 0.5];

Animation3D(SimpleArm,t,x,...
    'showtime',true,'showcom',true,'axis',axv)

Animation3D(SimpleArm,t,xm,...
    'showtime',true,'showcom',true,'axis',axv,'parameters',pmod)






