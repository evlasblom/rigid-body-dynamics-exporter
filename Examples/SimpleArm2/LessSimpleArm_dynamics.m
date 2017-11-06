function xo = LessSimpleArm_dynamics(xi)

h = 1/100; 

% Dynamics
[M,C,G] = LessSimpleArm_M_C_G(xi);

% Forward Dynamics
n = length(xi)/2;
qn = xi(1:n);
qd = xi(n+1:end);
qdd = inv(M)*(-C*qd-G-0.01*qd);

% Euler Step
q_next = qn + qd*h;
qd_next = qd + qdd*h;
xo = [q_next', qd_next'];

end