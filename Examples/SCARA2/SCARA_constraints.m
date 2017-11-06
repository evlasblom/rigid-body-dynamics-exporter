function [D,dD,dconv] = SCARA_constraints(q)

% upperbody angle
% D1 <= -tol;   ---> constrain
D1 = q(2);
dD1 = [0 1 0 0];
dconv1 = 0;

% arm location
[D2,dD2,dconv2] = SCARACon4_D_dD_dconv(q);

% combine
D = [D1;
    D2];

dD = [dD1;
    dD2];

dconv = [dconv1;
    dconv2];

end