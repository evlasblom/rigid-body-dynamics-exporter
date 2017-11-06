
clear all
clc

% Simulation with the nonlinear and linearized dynamics of a template model
% for the humanoid robot TUlip.

%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'TriplePendulum';

% --- Parameters
% connectivity
lambda = {0 1 2};

% joint types
joints = {jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute};

% joint axes
axes = {'y','y','y'};

% joint offsets, center of mass, inertial properties
[dout,rout,mout,iout] = TUlip_parameters();
    
d = [dout(1,:);   % foot depth
    dout(3,:);    % lowerleg
    dout(4,:)];   % upperleg

r = [rout(2,1),    0,      rout(2,3);   % lowerleg
    rout(3,1),     0,      rout(3,3); % upperleg
    rout(6,1),0,rout(6,3)+dout(6,3)]; % trunk

mass = [2*mout(2),... % lowerleg
    2*mout(3), ...     % upperleg
    2*mout(4)+2*mout(5)+mout(6)];        % trunk

inertia = [2*iout(2,:);  % lowerleg
    2*iout(3,:);         % upperleg
    iout(6,:)];            % trunk

% --- Create model structure
% for the dynamics, the length of the last link does not matter
% only the com position counts, but for visualization, you want
% to see nice stuff, hence the modfnc. The size of the other links
% is automatically determined.
TriplePendulum = CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@TriplePendulum_edit);

t0 = 0;
x0 = zeros(length(lambda))';
Animation3D(TriplePendulum,t0,x0,'showframe',true,'showcom',true)


%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------

RecursiveDynamics(TriplePendulum,'par_edit',true)

RecursiveLinearizedDynamics(TriplePendulum,'par_edit',true)

RecursiveCOMJacobian(TriplePendulum,'par_edit',true)


%% ------------------ SIMULATE  -----------------------------------------


disp('SIMULATION')

modelname = 'TriplePendulum';
load([modelname,'.mat']); n = TriplePendulum.dof;
load([modelname,'_parameters.mat']); par = parameters;

% When the center of mass is located in the middle of the segments, the
% nonlinear and linearized dynamics behave similar around the equilibrium
% point with sufficient damping. Though, commenting the code below (thereby
% using the original center of mass positions) shows that linearized
% dynamics do not always match the nonlinear system very well.
par.r1_x = 0;
par.r2_x = 0;
par.r3_x = 0;


% --- Simulation Nonlinear System
x0 = [[170 5 5]*pi/180 ...
       [0 0 0]*pi/180];
t_end = 5;
h = 1/100; 
damp = 0.25;

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;
r = zeros(length(t),3); r(1,:) = COMCalculator(TriplePendulum,x0,par)';     
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = TriplePendulum_M_C_G(x(ii,:),par); 
    
    % Forward Dynamics
    qn = x(ii,1:n)'; 
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-damp*qd);
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    x(ii+1,:) = [q_next', qd_next'];
    
    % Post-Processing
    r(ii+1,:) = COMCalculator(TriplePendulum,x(ii+1,:),par)';
end
fprintf('  ')
toc
fprintf('\n');


% --- Simulation Linearized System
x0 = [[170 5 5]*pi/180 ...
       [0 0 0]*pi/180];
t_end = 5;
h = 1/100; 
damp = 0.25;

tic
t = 0:h:t_end;
xlin = [[180 0 0]*pi/180, [0 0 0]];
rlin = COMCalculator(TriplePendulum,xlin,par)';
[Ml,Kl] = TriplePendulumLin_M_K(xlin,par); 
Jc = TriplePendulum_Jc(xlin,par);

xl = zeros(length(t),length(x0)); xl(1,:) = x0-xlin;
rl = zeros(length(t),3); rl(1,:) = Jc*xl(1,1:n)';   
for ii = 1:length(t)-1    
    % Forward Dynamics
    qn = xl(ii,1:n)'; 
    qd = xl(ii,n+1:end)';
    qdd = inv(Ml)*(-Kl*qn-damp*qd);
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    xl(ii+1,:) = [q_next', qd_next'];
    
     % Post-Processing
    rl(ii+1,:) = Jc*xl(ii+1,1:n)';
end
xl = xl+xlin(ones(size(xl,1),1),:);
rl = rl+rlin(ones(size(rl,1),1),:);
fprintf('  ')
toc
fprintf('\n');


%% --- Visualize
plotoptions.ylabel = {'\theta_1 (rad)','\theta_2 (rad)','\theta_3 (rad)'};
plotoptions.legend = {'Nonlinear','Linear'};
plotoptions.title = {'Joint Angles'};
PlotStates(t,x,t,xl,[1 2 3],plotoptions)

plotoptions.ylabel = {'x (m)','y (m)','z (m)'};
plotoptions.legend = {'Nonlinear','Linear'};
plotoptions.title = {'Center of Mass'};
PlotStates(t,r,t,rl,[1 2 3],plotoptions)

axv = [-1 1 -1 1 -1 1]*1;

Animation3D(TriplePendulum,t,x,...
    'timescale',true,...
    'showtime',true,'showcom',true,'showframe',true,'axis',axv)

Animation3D(TriplePendulum,t,xl,...
    'timescale',true,...
    'showtime',true,'showcom',true,'showframe',true,'axis',axv)








