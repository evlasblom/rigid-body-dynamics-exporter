
clear all
clc

% Humanoid robot TUlip using the Jacobian of the center of mass to move
% it's center of mass up and down. It's a simple example on a complex
% system so the final movement looks somewhat... weird.

% NOTE: for complex structures like this, it is advised to use the function
% rbde2spatial to change model definitions and use Featherstone's toolbox
% called spatial_v2 (see Upperbody example). It performs a lot faster as
% the number of degrees of freedom increases. 
% Later versions of this toolbox will have increased computational speed.


%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'TUlip';

% --- Parameters
% connectivity
lambda = {0 1 2 3 4 5 6 7 8 9 10 11};

% joint types
joints = {jointtype.Revolute, ... % left leg
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
axes = {'x','y','y','y','x','z',... % left leg
    'z','x','y','y','y','x'}; % right leg

% joint offsets, center of mass, inertial properties
[d,r,mass,inertia] = TUlip_parameters();

% --- Create model structure
% for the dynamics, the length of the last link does not matter
% only the com position counts, but for visualization, you want
% to see nice stuff, hence the modfnc. The size of the other links
% is automatically determined.
CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@TUlip_edit,'base',true);

t0 = 0;
x0 = zeros(length(lambda))';
axv = [-1 1 -1 1 -0.5 1.5];
Animation3D(TUlip,t0,x0,'showframe',true,'showcom',true,'axis',axv)


%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------

RecursiveDynamics(TUlip);

RecursiveCOMJacobian(TUlip);

point = [0.0850, 0, 0.3890]; 
RecursiveConstraints(TUlip,'body',6,'point',point);


%% ------------------ SIMULATE  -----------------------------------------

disp('SIMULATION')

load([modelname,'.mat']); n = TUlip.dof;

% --- Simulate
x0 = [[-10 5 -10 5 -5 0 0 -10 -10 20 -5 0]*pi/180, ...
    [0 0 0 0 0 0 0 0 0 0 0 0]*pi/180];
t_end = 4;
h = 1/500; 
damp = 0.25;

tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;
u = zeros(length(t),n);     
r = zeros(length(t),3); r(1,:) = COMCalculator(TUlip,x0)';
rdes = zeros(length(t),3);
e = zeros(length(t),3);
de = zeros(length(t),3);
for ii = 1:length(t)-1   
    % Dynamics
    [M,C,G] = TUlip_M_C_G(x(ii,:));             % eom matrices
    Jc = TUlip_Jc(x(ii,:));                     % com jacobian
    
    % Gravity compensation
    u(ii,:) = G;
    
    % Center of mass control
    % - position, reference, error
    if t(ii) <= 2; rdes(ii,:) = [0 0 0.68];
    elseif t(ii) <= 4; rdes(ii,:) = [0 0 0.72];
    end
    r(ii,:) = COMCalculator(TUlip,x(ii,:))';
    e(ii,:) = rdes(ii,:) - r(ii,:); 
    if ii>1 de(ii,:) = (e(ii,:) - e(ii-1,:))./h;
    end
    
    % - center of mass control
    lam = 0.25;
    Ji = Jc'*inv(Jc*Jc' + lam^2*eye(3)); % damped pseudo inverse
    dphi = Ji(1:6,:)*e(ii,:)';
    Kp = [2 1 1 1 2 0.5]*220;
    Kd = [4 1 1 1 4 1]*10;
    u(ii,1) = Kp(1)*dphi(1) - Kd(1)*x(ii,1+n) + u(ii,1);
    u(ii,2) = Kp(2)*dphi(2) - Kd(2)*x(ii,2+n) + u(ii,2);
    u(ii,3) = Kp(3)*dphi(3) - Kd(3)*x(ii,3+n) + u(ii,3);
    u(ii,4) = Kp(4)*dphi(4) - Kd(4)*x(ii,4+n) + u(ii,4);
    u(ii,5) = Kp(5)*dphi(5) - Kd(5)*x(ii,5+n) + u(ii,5);
    u(ii,6) = Kp(6)*dphi(6) - Kd(6)*x(ii,6+n) + u(ii,6);
    
    % Forward Dynamics
    qn = x(ii,1:n)'; 
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-damp*qd+u(ii,:)');
    
    % Euler Step
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    x(ii+1,:) = [q_next', qd_next'];
end
fprintf('  ')
toc
fprintf('\n');

r(end,:) = COMCalculator(TUlip,x(end,:))';
rdes(end,:) = rdes(end-1,:);

%% ------------------ VISUALIZE  -----------------------------------------

clear options
options.ylabel = {'r_x (m)','r_y (m)','r_z (m)'};
options.legend = {'Position','Reference'};
options.title = {'Controlling the Center of Mass'};
PlotStates(t,r,t,rdes,[1 2 3],options)

axv = [-1 1 -1 1 -0.5 1.5];
Animation3D(TUlip,t,x,...
    'timescale',true,...
    'showtime',true,...
    'showcom',true,...
    'axis',axv,...
    'export',false)

