%TILDE Returns the twist coordinates from a twist.
%   T = TILDE(t) returns a 4-by-4 twist matrix if the six twist coordinates
%   are supplied or it returns a 3-by-3 skew symmetric matrix if three 
%   parameters are given.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function x_tilde = Tilde(x)
if size(x,2) == 6
    x = x.';
end
if length(x) == 6
    x_tilde = [Skew(x(4:6)) x(1:3); zeros(1,4)];
elseif length(x) == 3
    x_tilde = Skew(x);
else
    error('Input should be either 3 or 6 elements long')
end
end