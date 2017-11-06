%UNSKEW Turns a skew symmetric matrix into a vector.
%   x = UNSKEW(X) returns the three parameters when of the skew symmetric
%   matrix in a vector.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function x = UnSkew(x_tilde)
x = [x_tilde(3,2); x_tilde(1,3); x_tilde(2,1)];
end