
clear all
clc

% Simulation and control of a standard robot arm using exported matrices,
% including an end-effector Jacobian. Results are saved to a video.

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

point = [0 0 0]';
RecursiveConstraints(SCARA,'body',4,'point',point)


%% ------------------ SIMULATE  -----------------------------------------

disp('SIMULATION')

load([modelname,'.mat']); n = SCARA.dof;

% --- Simulate
x0 = [[10, 30, 45]*pi/180, 0 ...
       [0, 0, 0]*pi/180, 0];
t_end = 9;
h = 1/100; 

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;
u = zeros(length(t),n);     
p = zeros(length(t),3)';
pdes = zeros(length(t),3)';
e = zeros(length(t),3)';
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = SCARA_M_C_G(x(ii,:));             % eom matrices
    [~,J,~] = SCARACon4_D_dD_dconv(x(ii,:));    % only jacobian
    Jp = J(1:3,:);                              % only position jacobian
    
    % Gravity compensation
    u(ii,:) = G;
    
    % End-effector control
    % - position, reference, error
    if t(ii) <= 3; pdes(:,ii) = [0.05 0.05 0.05]';
    elseif t(ii) <= 6; pdes(:,ii) = [-0.05 0.15 0.15]';
    elseif t(ii) <=9; pdes(:,ii) = [-0.05 -0.15 0.05]';
    end
    H4_0 = H_Transform(SCARA,x(ii,:),4,0);
    p(:,ii) = Transform(point,H4_0);
    e(:,ii) = pdes(:,ii) - p(:,ii);  
    % - joint space control
    lam = 0.25;
    Ji = Jp'*inv(Jp*Jp' + lam^2*eye(3)); % damped pseudo inverse
    dphi = Ji*e(:,ii);
    u(ii,1) = 80*dphi(1) - 2*(x(ii,5)-0) + u(ii,1);
    u(ii,2) = 80*dphi(2) - 1*(x(ii,6)-0) + u(ii,2);
    u(ii,3) = 80*dphi(3) - 1*(x(ii,7)-0) + u(ii,3);
    u(ii,4) = 20*dphi(4) - 2*(x(ii,8)-0) + u(ii,4);
    
    % Forward Dynamics
    qn = x(ii,1:n)'; 
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-0.005*qd+u(ii,:)');
    
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
options.legend = {'Position','Reference'};
options.title = {'Controlling the End-Effector'};
options.events = [3 6 9];
PlotStates(t,p',t,pdes',[1 2 3],options)

Animation3D(SCARA,t,x,...
    'showtime',true,...
    'showcom',true,...
    'axis',[-0.2 0.2 -0.2 0.2 0 0.4],...
    'export',false)


