%TWISTEXPONENT calculates the exponent of a twist.
%   H = TWISTEXPONENT(T,theta) computes the exponent of the twist T
%   giving the transformation matrix as H = exp(T*theta), with T being
%   the matrix form. Though, the input here can also be a vector.
%
%   References(2), eq. 2.32.
%   References(2), eq. 2.36.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013

function expT = TwistExponent(T,theta)
if size(T,1) == 4 && size(T,2) == 4
    t = UnTilde(T);
    v = t(1:3);
    w = t(4:6);
    if (sum(w) == 0) 
        % in case of a pure translation
        expW = eye(3);
        p = v*theta;
    else 
        % in case of a rotation and possible rotation induced translation
        expW = Rodrigues(T(1:3,1:3),theta);
        p = (eye(3)-expW)*(cross(w,v))+w*w'*v*theta;
    end
    expT = [expW p; zeros(1,3), 1];
elseif size(T,1) == 3 && size(T,2) == 3
    expT = Rodrigues(T,theta);
elseif (size(T,1) == 6 || size(T,1) == 3)  && size(T,2) == 1
    expT = TwistExponent(Tilde(T),theta);
elseif (size(T,2) == 6 || size(T,2) == 3)  && size(T,1) == 1
    expT = TwistExponent(Tilde(T'),theta);
else
    error('Input should be a twist or a column vector of twist coordinates')
end
end