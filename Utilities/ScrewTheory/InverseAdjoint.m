%INVERSEADJOINT Computes the inverse adjoint matrix.
%   Ad = INVERSEADJOINT(H) uses the homogeneous transformation H to compute the
%   inverse of the Adjoint matrix.
%
%   Ad = INVERSEADJOINT(H,-1) uses the homogeneous transformation H to compute the
%   Adjoint matrix.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function Ad = InverseAdjoint(H,varargin)
if size(H,1) == 4 && size(H,2) == 4
    R = H(1:3,1:3);
    p = H(1:3,4);
    if nargin == 1
        Ad = [R.' -R.'*Tilde(p); zeros(3,3) R.'];
    elseif nargin>1
        if varargin{1} == -1
            Ad = [R Tilde(p)*R; zeros(3,3) R];
        else
            error('Wrong argument specified');
        end
    end
else
    error('The input should be a homogeneous transformation')
end
end