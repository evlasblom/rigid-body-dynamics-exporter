%H2PBAR Gives the translation as a homogeneous transform.
%   R = H2PBAR(H) returns the translation vector of the homogeneous
%   transformation matrix H, as a homogeneous transformation matrix. That
%   is, the rotation matrix is set to identity.
%
%   Made by Erik Vlasblom
%   Last modified: 31-10-2015
function H = h2pbar(H)
if size(H,1) == 4 && size(H,2) == 4
    H(1:3,1:3) = eye(3);
else
    error('The input is not a homogeneous transform');
end
end