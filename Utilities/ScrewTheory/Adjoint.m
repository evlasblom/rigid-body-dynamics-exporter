%ADJOINT Computes the adjoint matrix.
%   Ad = ADJOINT(H) uses the homogeneous transformation H to compute the
%   Adjoint matrix.
%
%   Ad = ADJOINT(H,-1) uses the homogeneous transformation H to compute the
%   inverse of the Adjoint matrix.
%
%   Made by Erik Vlasblom
%   Last modified: 10-10-2014
function Ad = Adjoint(H,varargin)
[n,m] = size(H);
if (n == 4 && m == 4) || (n == 3 && m == 3)
    R = H(1:3,1:3);
    if n == 4; p = H(1:3,4);
    elseif n == 3; p = [0;0;0];
    end;
    if nargin == 1
        Ad = [R Tilde(p)*R; zeros(3,3) R];
    elseif nargin>1
        if varargin{1} == -1
            Ad = [R.' -R.'*Tilde(p); zeros(3,3) R.'];
        else
            error('Wrong argument specified');
        end
    end
else
    error('The input should be a homogeneous transformation')
end
end