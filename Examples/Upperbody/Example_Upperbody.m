
% Comparison with the spatial_v2 toolbox by Featherstone.

clear all
clc

%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'Upperbody';

% --- Parameters
lambda = {0 1 2 3 1 5 6};

joints = {jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute};

axes = {'y','x','z','z','x','z','z'};

d = [0,     0,      0;
    0,      0.05,   0.1;
    0,      0,      0;
    0,      0.1,    0;
    0,     -0.05,   0.1;
    0,      0,      0;
    0,     -0.1,    0];

r = [0,     0,      0.02;
    0,      0,      0;
    0,      0.05,   0;
    0,      0.05,   0;
    0,      0,      0;
    0,     -0.05,   0;
    0,     -0.05,   0];

mass = [1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];

inertia = [0.1, 0.1, 0.1;
    0.01, 0.01, 0.01;
    0.01, 0.01, 0.01;
    0.01, 0.01, 0.01;
    0.01, 0.01, 0.01;
    0.01, 0.01, 0.01;
    0.01, 0.01, 0.01];

% --- Create model structure
Upperbody = CreateModel(lambda,joints,axes,d,r,mass,inertia,...
    'name',modelname,'modfnc',@Upperbody_edit);

Upperbody_spatial = rbde2spatial(Upperbody);

t0 = 0;
x0 = zeros(1,length(lambda));
axv = [-0.25 0.25 -0.25 0.25 -0.25 0.25];

Animation3D(Upperbody,t0,x0,'showcom',true,'axis',axv)

% visualization using spatial_v2 toolbox
try
    showmotion(Upperbody_spatial,[0 1],[x0' x0'])
catch
    warning('Featherstone''s spatial_v2 toolbox not found.');
end


%% ------------------ DERIVE AND EXPORT EQUATIONS OF MOTION --------------


RecursiveDynamics(Upperbody)



%% ------------------ SIMULATE  -----------------------------------------
% We can now modify the parameters and change the dynamics without having
% to derive them again because of the paredit = true option.

disp('SIMULATION')

load([modelname,'.mat']); n = Upperbody.dof;
load([modelname,'_parameters.mat']); p = parameters;

x0 = [[10, 30, 45, 20, -30, 45, 20]*pi/180 ...
    [0, 0, 0, 0, 0, 0, 0]*pi/180]; x0(5:7) = -x0(5:7);
t_end = 10;
h = 1/100;

% --- Simulation
tic
t = 0:h:t_end;
x = zeros(length(t),length(x0)); x(1,:) = x0;
for ii = 1:length(t)-1
    [M,C,G] = Upperbody_M_C_G(x(ii,:));
    qn = x(ii,1:n)';
    qd = x(ii,n+1:end)';
    qdd = inv(M)*(-C*qd-G-0.05*qd);
    
    q_next = qn + qd*h;
    qd_next = qd + qdd*h;
    x(ii+1,:) = [q_next', qd_next'];
end
fprintf('  ')
toc
fprintf('\n');

% --- Simulation with spatial toolbox
try
    tic
    t = 0:h:t_end;
    xs = zeros(length(t),length(x0)); xs(1,:) = x0;
    for ii = 1:length(t)-1
        [M,c] = HandC(Upperbody_spatial,xs(ii,1:n),xs(ii,n+1:end));
        qn = xs(ii,1:n)';
        qd = xs(ii,n+1:end)';
        qdd = inv(M)*(-c-0.05*qd);
        
        q_next = qn + qd*h;
        qd_next = qd + qdd*h;
        xs(ii+1,:) = [q_next', qd_next'];
    end
    fprintf('  ')
    toc
    fprintf('\n');
catch
    % do nothing
end



%% ------------------ VISUALIZE  -----------------------------------------

try
    options.ylabel = {'\phi_1 (rad)','\phi_2 (rad)','\phi_3 (rad)','\phi_4 (rad)'};
    options.legend = {'Exported','Spatial Toolbox'};
    options.title = {'Comparison to Spatial Toolbox'};
    PlotStates(t,x,t,xs,[1 2 3 4],options)
catch
    options.ylabel = {'\phi_1 (rad)','\phi_2 (rad)','\phi_3 (rad)','\phi_4 (rad)'};
    options.legend = {'Exported'};
    options.title = {'Comparison to Spatial Toolbox'};
    PlotStates(t,x,[1 2 3 4],options)
end

axv = [-0.25 0.25 -0.25 0.25 -0.25 0.25];

try  
    Animation3D(Upperbody,t,x,...
        'showtime',true,'showcom',true,'axis',axv)
    
    % visualization using spatial_v2 toolbox:
    showmotion(Upperbody_spatial,t,xs(:,1:n)')
catch
    Animation3D(Upperbody,t,x,...
    'showtime',true,'showcom',true,'axis',axv)
end


