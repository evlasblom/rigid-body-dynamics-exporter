%H2P Gives the translation of a transform.
%   P = H2P(H) returns the translation vector of the homogeneous
%   transformation matrix H.
%
%   Made by Erik Vlasblom
%   Last modified: 29-05-2013
function p = h2p(H)
if size(H,1) == 4 && size(H,2) == 4
    p = H(1:3,4);
else
    error('The input is not a homogeneous transform');
end
end