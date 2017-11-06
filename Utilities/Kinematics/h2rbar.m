%H2RBAR Gives the rotation matrix as a homogeneous transform.
%   R = H2RBAR(H) returns the rotation matrix of the homogeneous
%   transformation matrix H, as a homogeneous transformation matrix. That
%   is, the translation part is set to zero.
%
%   Made by Erik Vlasblom
%   Last modified: 31-10-2015
function H = h2rbar(H)
if size(H,1) == 4 && size(H,2) == 4
    H(1:3,4) = [0 0 0]';
else
    error('The input is not a homogeneous transform');
end
end