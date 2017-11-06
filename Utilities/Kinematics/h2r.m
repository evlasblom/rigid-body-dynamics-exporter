%H2R Gives a rotation matrix from a transform.
%   R = H2R(H) returns the rotation matrix of the homogeneous
%   transformation matrix H.
%
%   Made by Erik Vlasblom
%   Last modified: 29-05-2013
function R = h2r(H)
if size(H,1) == 4 && size(H,2) == 4
    R = H(1:3,1:3);
else
    error('The input is not a homogeneous transform');
end
end