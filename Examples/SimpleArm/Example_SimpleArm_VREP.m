clear all
pause(0.5);
clc

% Add Rigid Body Dynamics Exporter to startup.m

% Add all subfolders
working_dir = cd;
addpath(genpath(working_dir));


%% --------- VREP SETUP ----------------------

% Make sure a synchronous connection with the following port is set up 
% in the file remoteApiConnections.txt. Then start Vrep and run this file.
port = 19997;
% ip = '192.168.0.93';
ip = '127.0.0.1';

% Specify paths to matlab api files
path = 'C:\Program Files (x86)\V-REP3\V-REP_PRO_EDU';
libpath = { fullfile(path, 'programming', 'remoteApi')
    fullfile(path, 'programming', 'remoteApiBindings', 'matlab', 'matlab')
    fullfile(path, 'programming', 'remoteApiBindings', 'lib', 'lib')
    }; addpath( libpath{:} );

disp('Starting vrep connection.')
vrep = remApi('remoteApi','extApi.h');
client = vrep.simxStart(ip, port, true, false, 2000, 5);
if client < 0; 
    error('Connecting failed.'); 
end

% Some constants
vrep_steptime = 50; % copy from vrep [ms]
oneshot_wait = vrep.simx_opmode_oneshot_wait;
oneshot = vrep.simx_opmode_oneshot;
stream = vrep.simx_opmode_streaming;
buffer = vrep.simx_opmode_buffer;

% Start scene
vrep.simxAddStatusbarMessage(client,'Connection started.',oneshot_wait);
% vrep.simxLoadScene(client,'SimpleArm_Matlab.ttt',1,oneshot_wait);

% Get handles
[~, h{1}] = vrep.simxGetObjectHandle(client,'joint1#', oneshot_wait);
[~, h{2}] = vrep.simxGetObjectHandle(client,'joint2#', oneshot_wait);
% [~,h,intData,floatData,stringData] = vrep.simxGetObjectGroupData(...
%     client,vrep.sim_object_joint_type,15,oneshot_wait)



%% --------- USER INPUT ----------------------

% Load models
load('SimpleArm.mat'); model = SimpleArm; n = model.dof;
load('SimpleArm_parameters.mat'); p = parameters;

% Options
t_end = 2; % end time
vis = 1; % visualisation on / off

options.h = vrep_steptime/1000; % step time [s]
options.f = 1/options.h; % simulation frequency [Hz]
options.lam = 0.01; % damping coefficient for least squares solution

% Initial positions
x0 = [[120 0]*pi/180,... % joint position
    0 0]; % joint speeds are 0
for ii = 1:length(h);
    vrep.simxSetJointPosition(client,h{ii},x0(ii),oneshot_wait);
end


%% --------- REAL WORLD SIMULATION ----------------------

% 1. Pre-allocate
t = 0:options.h:t_end;
x = zeros(length(t),length(x0));
u = zeros(length(t),n);
ref = cell(n,1);
for ii = 1:n; ref{ii} = ['ref',num2str(ii)]; end


% 2. Initial values
x(1,:) = x0;

% 3. Simulation via Vrep
tic
vrep.simxSynchronous(client,true);
vrep.simxStartSimulation(client,oneshot);
fprintf('\n-- Ground truth simulation \n')

% Vrep state
for jj = 1:n
    % Joint angle: ok
    [~, x(1,jj)] = vrep.simxGetJointPosition(client, h{jj},oneshot_wait);
    % JOINT ANG VEL?
end

for ii = 1:length(t)-1           
    % Control TRANSFORM TO GLOBAL!
    [M,C,G] = SimpleArm_M_C_G(x(ii,:));
%     u(ii,:) = G;
    
    % Vrep commands: ok
    signal = vrep.simxPackFloats(u(ii,:));
    vrep.simxSetStringSignal(client,'u',signal,oneshot_wait);

    % Simulation step: ok
    vrep.simxSynchronousTrigger(client);
    
    % Vrep state
    for jj = 1:n
        % Joint angle: ok
        [~, x(ii+1,jj)] = vrep.simxGetJointPosition(client, h{jj},oneshot_wait);
        % JOINT ANG VEL?
    end
    
end

vrep.simxPauseSimulation(client,oneshot);
toc




%% --------- VISUALIZATION --------------------------

if vis == 1
    fprintf('\n-- Visualizing \n')
    close all
    
    % ---- STATES ----
    % black, green, blue, red
    PlotStates(t,x,[1 2])
    PlotStates(t,u,[1 2])
            
   
    % ---- ANIMATION ----
    % animation of real simulation
    xview = [-0.4 0.4 -0.4 0.4 -0.4 0.4];
    Animation3D(model,t,x,'showtime',true,'axis',xview)
end



%% --------- VREP CLOSE  ---------------
if exist('vrep','var')
    if client >= -1
        disp('Closing vrep connection')
        vrep.simxAddStatusbarMessage(client,'Connection terminated.',oneshot_wait);
        vrep.simxFinish(client);
        pause(0.5);
    end
end

