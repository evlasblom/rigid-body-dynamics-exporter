
clear all
clc

% Simulation and control of a standard robot arm using exported matrices,
% including constraints. The arm is swung against its joint angle limit
% after which the end-effector is dropped on the floor, halting the motion
% due to friction.

%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'SCARA';

% --- Parameters
% connectivity
lambda = {0 1 2 3};

% joint types
joints = {jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Prismatic};

% joint axes
axes = {'z','z','z','z'};

% joint offsets
d = [0,     0,      0.2; % location joint 1
    0,      0.1,   0;    % distance joint 1 - joint 2
    0,      0.1,   0;    % distance joint 2 - joint 3
    0,      0,    0]; % distance joint 3 - joint 4

% center of mass
r = [0,     0.05,   0;  % distance from joint 1
    0,      0.05    0;  % distance from joint 2
    0,      0,   0;  % distance from joint 3
    0,      0,   0.025]; % distance from joint 4

% inertial properties
mass = [0.5, 0.5, 0.1, 0.1];
inertia = [0.1, 0.1, 0.1; 
    0.01, 0.01, 0.01; 
    0.01, 0.01, 0.01; 
    0.01, 0.01, 0.01];

% --- Create model structure
% for the dynamics, the length of the last link does not matter
% only the com position counts, but for visualization, you want
% to see nice stuff, hence the modfnc. The size of the other links
% is automatically determined.
SCARA = CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@SCARA_edit);

t0 = 0;
x0 = zeros(length(lambda))';
Animation3D(SCARA,t0,x0,'showcom',true,'showframe',true,...
    'axis',[-0.1 0.3 -0.1 0.3 0 0.4])


%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------


RecursiveDynamics(SCARA)

point = [0 0 -0.05]';
RecursiveConstraints(SCARA,'body',4,'point',point)


%% ------------------ SIMULATE  -----------------------------------------

disp('SIMULATION')

load([modelname,'.mat']); n = SCARA.dof;

con = true;         % constraints on/off
lam = 0.001;        % damped pseudo-inverse coefficient (not used now)
tol = 1e-12;        % constraint tolerance
beta = 0.01;        % joint damping
e = 0.1;            % restitution coefficient 
FDYN = @SCARA_M_C_G;                        % Dynamics function
FCON = @SCARA_constraints;                  % Constraints function
FCONDET = @SCARA_constraintdetection;       % Constraint detection

% --- Simulate
d2r = @(d) d*pi/180;
x0 = [d2r([10, 30, 45]), 0 ...
       d2r([-40, -40, 0]), 0];
t_end = 5;
h = 1/100; 

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;
xref = x;
u = zeros(length(t),n);     
p = zeros(length(t),3)';
pdes = zeros(length(t),3)';
err = zeros(length(t),3)';
for ii = 1:length(t)-1   
    % Control
    Kp = [10 0 10 5];
    Kd = [1 0 1 0.5];
    if t(ii) <= 2; xref(ii,1) = d2r(-5);
    elseif t(ii) > 2; xref(ii,1) = d2r(15);
    end
    xref(ii,2) = x(ii,2);
    xref(ii,3) = 0;
    if t(ii) <= 3; xref(ii,4) = 0; 
    elseif t(ii) > 3; xref(ii,4) = - 0.15;
    end
    u(ii,1) = -Kp(1)*(x(ii,1)-xref(ii,1)) - Kd(1)*(x(ii,5)-0);
    u(ii,2) = -Kp(2)*(x(ii,2)-xref(ii,2)) - Kd(2)*(x(ii,6)-0);
    u(ii,3) = -Kp(3)*(x(ii,3)-xref(ii,3)) - Kd(3)*(x(ii,7)-0);
    u(ii,4) = -Kp(4)*(x(ii,4)-xref(ii,4)) - Kd(4)*(x(ii,8)-0) + 9.81*0.1;
    
    % Numerical Integration
    x(ii+1,:) = CNI(SCARA,FDYN,FCON,FCONDET,...
        t(ii),x(ii,:),u(ii,:)',h,beta,con,e,lam,tol);
    
end
fprintf('  ')
toc
fprintf('\n');

%% ------------------ VISUALIZE  -----------------------------------------

clear options
options.ylabel = {'\theta_1 (rad)','\theta_2 (rad)','\theta_3 (rad)','\theta_4 (rad)'};
options.legend = {'Measured','Reference'};
options.title = {'Constrained Degrees of Freedom'};
options.events = [2 3];
PlotStates(t,x,t,xref,[1 2 3 4],options)

Animation3D(SCARA,t,x,...
    'showtime',true,...
    'showcom',true,...
    'axis',[-0.25 0.25 -0.25 0.25 -0.1 0.4],...
    'export',false)


