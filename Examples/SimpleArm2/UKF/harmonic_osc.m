function xn = harmonic_osc(x)

h = 0.01; % stepsize
m = 0.1; % mass
c = 0.2; % spring damping
k = 5; % spring stiffness
f = 0; % control signal

% state space eq
dx = [x(2);
    1/m*(-m*9.81-c*x(2)-k*x(1)+f)];

% numerical integration
xn = x + dx*h;

end