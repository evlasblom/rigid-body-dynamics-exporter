
clear all
clc

% Derivation of dynamics that can be used to do gravity compensation in
% V-REP using the LessSimpleArm_VREP.m file.


%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'LessSimpleArm';

% --- Parameters
% connectivity
lambda = {0 1 2};

% joint types
joints = {jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute};

% joint axes
axes = {'y','y','x'};

% joint offsets
d = [0,     0,      0.1; 
    0,      0,      0.2;
    0,      0,      0.2];

% center of mass
r = [0,     0,      0.1; 
    0,      0,      0.1
    0,      0,      0.1];

% inertial properties
mass = [0.5, 0.5, 0.5];
inertia = [0.02, 0.02, 0.003; 
    0.02, 0.02, 0.003;
    0.02, 0.02, 0.003];

% --- Create model structure
% for the dynamics, the length of the last link does not matter
% only the com position counts, but for visualization, you want
% to see nice stuff, hence the modfnc. The size of the other links
% is automatically determined.
LessSimpleArm = CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@SimpleArm_edit);

t0 = 0;
x0 = zeros(length(lambda))';
Animation3D(LessSimpleArm,t0,x0,'showframe',true,'showcom',true)


%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------


RecursiveDynamics(LessSimpleArm)


%% ------------------ SIMULATE  -----------------------------------------


disp('SIMULATION')

modelname = 'LessSimpleArm';
load([modelname,'.mat']); n = LessSimpleArm.dof;
load([modelname,'_parameters.mat']); p = parameters;

% --- Simulate
x0 = [[120 -30 20]*pi/180 ...
       [0, 0, 0]*pi/180];
t_end = 2;
h = 1/100; 

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;       
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = LessSimpleArm_M_C_G(x(ii,:)); 
    
    % Forward Dynamics
    qn = x(ii,1:n)'; 
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-0.01*qd);
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    x(ii+1,:) = [q_next', qd_next'];
end
fprintf('  ')
toc
fprintf('\n');


%% ------------------ VISUALIZE  -----------------------------------------

PlotStates(t,x,[1 2 3])

Animation3D(LessSimpleArm,t,x,...
    'showtime',true,'showcom',true)









