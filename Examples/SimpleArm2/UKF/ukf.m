%UKF   Unsctended Kalman Filter
%   [x_p, x, P, K] = ukf(h,t,Q,R,P,x,z,ut,L) estimates the current
%   state and covariance using the previous state x and covariance P. The
%   function handle f should make a predicction for the next state where as
%   h is the measurement function.
%   Q and R are the process and measurement noise respectively and z is the
%   measurement signal.
%   For ut, use ut = ut_parameters();
%
%   This filter assumes additive noise and uses the scaled
%   unsctented transformation. 
%
% Erik Vlasblom
% 10-04-2013

% ============== UNSCTENTED KALMAN FILTER ===================

function [x_priori, x_posteriori, P_posteriori, K] = ukf(f,h,Q,R,P,x,z,ut,L)

if (size(x,2)>size(x,1))
    x = x';
end

Wm = ut.Wm;
Wc = ut.Wc;

% 1. Sigma points
X = sigmapoints(x,P,L,ut);

% 2. Prediction step
X_priori = []; x_priori = zeros(size(x));
for i = 1:size(X,2)
    Xtmp = f(X(:,i));
    X_priori(:,i) = Xtmp;
    x_priori = x_priori + Wm(i)*X_priori(:,i);
end
P_priori = covariance(Wc,X_priori,x_priori,X_priori,x_priori) + Q;

% 3. Correction step
Y = []; y = zeros(size(z));
X_y = sigmapoints(x_priori, P_priori, L, ut);
for i = 1:size(X,2)
    Y(:,i) = h(X_y(:,i));
    y = y + Wm(i)*Y(:,i);
end
P_yy = covariance(Wc,Y,y,Y,y) + R;
P_xy = covariance(Wc,X_priori,x_priori,Y,y);
K = P_xy/(P_yy);
x_posteriori = x_priori + K*(z-y);
P_posteriori = P_priori - K*P_yy*K';


    
% ----------- subfunctions ---------------------------------

    function X = sigmapoints(x,P,L,ut)
        % x     mean
        % P     covariance
        % par   structure with state length and measurement length
        % ut    structure with scaling parameters
        X1 = x(:,ones(1,L));
        X2 = sqrt(ut.scale)*chol(P)';
        X = [x X1+X2 X1-X2];
    end

    function P = covariance(W,A,a,B,b)
        % W     weights
        % A     sigma set
        % a     mean of sigma set
        P = (A-a(:,ones(size(A,2),1)))*diag(W)*(B-b(:,ones(size(B,2),1)))';
    end


end
