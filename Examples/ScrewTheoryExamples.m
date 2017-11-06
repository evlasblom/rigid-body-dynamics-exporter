
clear all
clc

% Manual derivation of some properties using different functions from the
% folder /utilities.

%% ------------------ CREATE ROBOT MODEL ---------------------------------

modelname = 'Upperbody';

% --- Parameters
% connectivity
lambda = {0 1 2 3 1 5 6}; p = lambda; n = length(lambda);
kappa = cell(size(lambda));
mu = cell(size(lambda));
for i = 1:n; kappa{i} = GetPathToRoot(lambda,i); end;
for i = 1:n; mu{i} = GetChildren(lambda,i); end;

% joint types
joints = {jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute, ...
    jointtype.Revolute};

% joint axes
axes = {'y','x','z','z','x','z','z'};

% joint offsets
d = [0,     0,      0; 
    0,      0.05,   0.1; 
    0,      0,      0; 
    0,      0.1,    0; 
    0,     -0.05,   0.1; 
    0,      0,      0; 
    0,     -0.1,    0];

% center of mass
r = [0,     0,      0.02; 
    0,      0,      0; 
    0,      0.05,   0; 
    0,      0.05,   0; 
    0,      0,      0; 
    0,     -0.05,   0; 
    0,     -0.05,   0];

% inertial properties
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
    'name',modelname,'export',false);

% modify visualization
bw = 0.02;
Upperbody.rigidbody(4).visual.bodysizemin = -[bw/2,bw/2,bw/2];
Upperbody.rigidbody(4).visual.bodysizemax = [bw/2,0.1,bw/2];
Upperbody.rigidbody(7).visual.bodysizemin = -[bw/2,0.1,bw/2];
Upperbody.rigidbody(7).visual.bodysizemax = [bw/2,bw/2,bw/2];
Upperbody.visual.bodywidth = bw;

% visualize
t0 = 0;
x0 = zeros(n)';
axv = [-0.25 0.25 -0.25 0.25 -0.25 0.25];
Animation3D(Upperbody,t0,x0,'showcom',true,'showframe',true,...
    'axis',axv)


%% ------------------ RECURSIVE DERIVATONS -------------------------------
% References are used to dynamically access structure elements.
% Relative to some other local frame: H and twists
l = @(i,j) ['l',num2str(i),'_',num2str(j)];
% Relative to global frame: H
g = @(i,j) ['g',num2str(i),'_',num2str(j)];
% Remaining body properties: mass, inertia, J, ...
b = @(i) ['b',num2str(i)];

% --- DEFINITIONS
q = sym('q',[n 1]);
dq = sym('dq',[n 1]);
ddq = sym('ddq',[n 1]);

% Local transformations:
H.(g(0,0)) = HomogeneousTransform([0,0,0]);
H.(l(1,p{1})) = HomogeneousTransform(q(1),'y',d(1,:));
H.(l(2,p{2})) = HomogeneousTransform(q(2),'x',d(2,:));
H.(l(3,p{3})) = HomogeneousTransform(q(3),'z',d(3,:));
H.(l(4,p{4})) = HomogeneousTransform(q(4),'z',d(4,:));
H.(l(5,p{5})) = HomogeneousTransform(-q(5),'x',d(5,:));
H.(l(6,p{6})) = HomogeneousTransform(-q(6),'z',d(6,:));
H.(l(7,p{7})) = HomogeneousTransform(-q(7),'z',d(7,:));

% Local Center of Mass transformations:
for i = 1:n; Hm.(b(i)) = HomogeneousTransform(r(i,:)); end;

% --- RECURSIVE EQUATIONS
disp('RECURSIVE EQUATIONS')
tic
J.(b(0)) = zeros(6,n);
GG.(b(0)) = zeros(n,1);
for i = 1:n
    fprintf(['  - body ',num2str(i),'\n']);
 
    % Global transformation
    H.(g(i,0)) = simplify(H.(g(p{i},0))*H.(l(i,p{i})));
    Rbar.(g(i,0)) = h2rbar(H.(g(i,0)));
    
    % Local unit twists, expressed in body coordinate frame
    unittwist.(l(i,p{i})) = simplify(...
        UnTilde(TwistFromHomogeneous(H.(l(i,p{i})),q(i),'b')));
 
    % Jacobian, expressed in body coordinate frame 
    fprintf('  - - jacobian \n');
    Jl = zeros(6,n); Jl(:,i) = unittwist.(l(i,p{i}));
    J.(b(i)) = simplify(...
        InverseAdjoint(H.(l(i,p{i})))*J.(b(p{i})) + Jl);   
    
    % Jacobian, transformed to analytically derived version
    Jta.(b(i)) = simplify(...
        Adjoint(Rbar.(g(i,0)))*J.(b(i)));

    fprintf('  ')
    toc
    fprintf('\n');
end

J = rmfield(J,b(0));
GG = rmfield(GG,b(0));
fprintf('\n\n')


%% ------------------ ANALYTIC DERIVATIONS -------------------------------


% --- ANALYTIC EQUATIONS
disp('ANALYTIC EQUATIONS')
tic
V.(b(0)) = 0;
dV.(b(0)) = zeros(n,1);
for i = 1:n
    fprintf(['  - body ',num2str(i),'\n']);
    
    % Jacobian, analytically derived
    % (should be related to geometric one)
    Ja.(b(i)) = jacobian(h2p(H.(g(i,0))),q);
    
    fprintf('  ')
    toc
    fprintf('\n');
end

V = rmfield(V,b(0));
dV = rmfield(dV,b(0));
fprintf('\n\n')



%% ------------------ COMPARISON RECURSIVE AND ANALYTIC DERIVATIONS ------

tol = 1e-12;
qr = rand(n,1);

% ---  COMPARISON RECURSIVE AND ANALYTIC DERIVATIONS 
disp('COMPARISON RECURSIVE AND ANALYTIC DERIVATIONS ')
tic
for i = 1:n
    fprintf(['  - body ',num2str(i),'\n']);
    
    % compute difference
    tmp = Jta.(b(i)); Jtai = tmp(1:3,:);
    Jai = Ja.(b(i)); 
    diff = subs(Jtai - Jai,q,qr);
    
    % difference
    if sum(sum(diff)) < tol
        fprintf('  Analytic Jacobian: OK \n')
    else
        fprintf('  Analytic Jacobian: WRONG \n')
    end
    
    fprintf('  ')
    toc
    fprintf('\n');
end

fprintf('\n\n')


%% ------------------ REFERENCES------------------------------------------

% [1] Stramigioli, Stefano, and Herman Bruyninckx. 
%     Geometry and screw theory for robotics.
%     Tutorial during ICRA 2001 (2001).

% [2] Murray, Richard M., et al. 
%     A mathematical introduction to robotic manipulation. 
%     CRC press, 1994.

% [3] Spong, Mark W., Seth Hutchinson, and Mathukumalli Vidyasagar. 
%     Robot modeling and control. Vol. 3. 
%     New York: Wiley, 2006.

% [4] Featherstone, Roy. 
%     A Beginner's Guide to 6-D Vectors (Part 1).
%     Robotics & Automation Magazine, IEEE 17.3 (2010): 83-94.

% [5] Featherstone, Roy. 
%     A Beginner's Guide to 6-D Vectors (Part 2)[Tutorial].
%     Robotics & Automation Magazine, IEEE 17.4 (2010): 88-99.







