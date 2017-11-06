

%% ---------------- REAL SYSTEM
h = 0.01;
t = 0:h:5;
x0 = [0;2];
L = length(x0);
x = zeros(L,length(t));


for ii = 1:length(t)
    x(:,ii+1) = harmonic_osc(x(:,ii));
end

figure
ax1 = subplot(2,1,1);
plot(ax1,t,x(1,1:end-1))
legend('position')
ax2 = subplot(2,1,2);
plot(ax2,t,x(2,1:end-1))
legend('velocity')
xlabel('Time [s]')



%% ------------ ESTIMATION
ut = ut_parameters('general_symmetric',length(x0));

what_meas = 'position';
f = @harmonic_osc;
h = @(x) x;
Q = diag([1e-6, 1e-6]); % process noise
R = diag([1e-1, 1e-1]); % measurement noise

x_est = zeros(size(x));
x_est(:,1) = x0 + (Q*randn(2,1)); % first guess estimate
P = 4*Q; % first guess covariance
z = zeros(size(x)); % measurement
for ii = 1:length(t)
    z(:,ii+1) = x(:,ii+1) + (R*randn(2,1));
    switch what_meas
        case 'position'
            z_in = z(1,ii+1);
            z(2,ii+1) = 0;
        case 'velocity'
            z_in = z(2,ii+1);
            z(1,ii+1) = 0;
        case 'both'
            z_in = z(:,ii+1);
    end
    [~, x_est(:,ii+1), P, ~] = ukf(f,h,Q,R,P,x_est(:,ii),z_in,ut,L);
end

figure
ax1 = subplot(2,1,1);
plot(ax1,t,x(1,1:end-1),t,z(1,1:end-1),t,x_est(1,1:end-1))
legend('position','measured','estimated')
ax2 = subplot(2,1,2);
plot(ax2,t,x(2,1:end-1),t,z(2,1:end-1),t,x_est(2,1:end-1))
legend('velocity','measured','estimated')
xlabel('Time [s]')



