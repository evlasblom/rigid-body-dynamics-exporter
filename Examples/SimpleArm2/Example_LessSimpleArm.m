
clear all
clc

% Derivation of dynamics that are used for state estmation with an
% Unscented Kalman Filter

addpath('UKF')

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
load([modelname,'.mat']); n = LessSimpleArm.dof; L = 2*n;
load([modelname,'_parameters.mat']); p = parameters;

% --- Simulate
x0 = [[120 -30 20]*pi/180 ...
       [0, 0, 0]*pi/180];
t_end = 2;
h = 1/100; % also defined in LessSimpleArm_dynamics

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;       
for ii = 1:length(t)  
    x(ii+1,:) = LessSimpleArm_dynamics(x(ii,:)');
end
fprintf('  ')
toc
fprintf('\n');


% --- Estimate
ut = ut_parameters('general_symmetric',length(x0));

h = @(x) x;
f = @LessSimpleArm_dynamics;
Q = diag(1e-6*ones(1,2*n)); % process noise
R = diag(1e-1*ones(1,2*n)); % measurement noise

tic
x_est = zeros(size(x));
x_est(1,:) = (x0' + (Q*randn(2*n,1))); % first guess estimate
P = 4*Q; % first guess covariance
z = zeros(size(x)); % measurement
for ii = 1:length(t)
    z(ii+1,:) = x(ii+1,:) + (R*randn(2*n,1))';
    [~, x_est(ii+1,:), P, ~] = ukf(f,h,Q,R,P,x_est(ii,:)',z(ii+1,:)',ut,L);
end
fprintf('  ')
toc
fprintf('\n');


%% ------------------ VISUALIZE  -----------------------------------------

options.title = {'State Estimation'};
options.ylabel = {'\phi_1 (rad)','\phi_2 (rad)','\phi_3 (rad)'};
options.legend = {'Real','Measured','Estimated'};
PlotStates(t,x(1:end-1,:),t,z(1:end-1,:),t,x_est(1:end-1,:),[1 2 3],options)

Animation3D(LessSimpleArm,t,x,...
    'showtime',true,'showcom',true)









