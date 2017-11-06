%RODRIGUES calculates the exponent of the rotational velocities.
%   R = RODRIGUES(W,theta) computes the exponent of the rotational
%   velocities giving the rotation matrix as R = exp(W*theta), with W being
%   the skew symmetric form. Though, the input here can also be a vector.
%
%   References(2), eq. 2.14.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013

function expA = Rodrigues(A,theta)
if size(A,1) == 3 && size(A,2) == 3
    a = UnTilde(A);
    norma = norm(a);
    if ~isa(norma,'sym')
        % in case there is a zero rotational twist
        if (norma <= 1e-10)
            expA = eye(3);
        % speed up calculation in case norm(a) = 1
        elseif (norma == 1)
            expA = eye(3) + A*sin(theta) + A*A*(1-cos(theta));
        end
    % the following is for the general case (i.e. including norm(a) ~= 1)
    else
    expA = eye(3) + A./norma*sin(norma*theta) + ...
        A*A./(norma^2)*(1 - cos(norma*theta)); 
    end
elseif (size(A,1) == 1 && size(A,2) == 3) || (size(A,2) == 1 && size(A,1) == 3)
    expA = Rodrigues(Tilde(A),theta);
else
    error('Matrix must be 3-by-3')
end
end