%UNTILDE Returns the twist coordinates from a twist.
%   t = UNTILDE(T) returns the six twist parameters from a 4-by-4 Twist or
%   it returns the three parameters when a 3-by-3 skew symmetric matrix is
%   supplied as T.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function unTilde = UnTilde(Tilde)
if size(Tilde,1) == 4 && size(Tilde,2) == 4
    unTilde = [Tilde(1:3,4);UnSkew(Tilde(1:3,1:3))];
elseif size(Tilde,1) == 3 && size(Tilde,2) == 3
    unTilde = UnSkew(Tilde);
end
end