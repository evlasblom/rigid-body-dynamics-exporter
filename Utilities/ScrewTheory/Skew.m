%SKEW Turns a vector into a skew symmetric matrix.
%   X = SKEW(x) returns a skew symmetric matrix of three parameters are
%   supplied in a vector.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function x_tilde = Skew(x)
x_tilde = [ 0    -x(3)  x(2);
            x(3)  0    -x(1);
           -x(2)  x(1)  0];
end