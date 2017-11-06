% Constrained Numerical Integration
function xo = CNI(model,FDYN,FCON,FCONDET,t,x,u,h,beta,con,e,lam,tol)

x = x';

persistent model_int n ncon 
persistent act stab loc imp

% 0. ----- Settings
if isempty(model_int) || t == 0
    model_int = model;
    n = model_int.dof;
    ncon = size(FCON(x));
    act = zeros(ncon); % active constraints (for FD)
    stab = zeros(ncon); % require stabilization (for CS)
    loc = zeros(ncon); % stabilization position (for CS)
    imp = zeros(ncon); % require impact (for CI)
end

% 1. ----- Forward Dynamics (FD)
dx = ForwardDynamics(...
        model_int,...
        FDYN,...
        FCON,...
        beta,x,u,act,n,lam);


% 2. ----- Numerical Integration (NI)
xn = x + (dx.*h);


% 3. ----- Constraint Detection (CD)
if con
    [~, act, stab, loc, imp, ~] = ...
        ConstraintDetection(model_int,...
        FDYN,...
        FCON,...
        FCONDET,...
        beta,xn,u,act,stab,loc,imp,n,tol,lam);
end


% 4. ----- Constraint Stabilization (CS)
if con
    if sum(stab) > 0
        xn = GaussNewton(...
            FCON,...
            xn,stab,loc,tol,lam);
    end
    stab = zeros(ncon);
end


% 5. ----- Contact Impulse (CI)
if con
    if sum(imp) > 0
        xn = ContactImpulse(...
            model_int,...
            FDYN,...
            FCON,...
            e,xn,imp,n,lam);
    end
    imp = zeros(ncon);
end


xo = xn';

end