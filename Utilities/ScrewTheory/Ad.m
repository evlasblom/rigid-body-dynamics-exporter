%AD Computes the adjoint matrix.
%   ad = AD(T) uses the matrix representation of a twist to compute
%   the adjoint matrix.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function ad = Ad(T)
if size(T,1) == 4 && size(T,2) == 4
    w = T(1:3,1:3);
    v = Tilde(T(1:3,4));
elseif size(T,1) == 6 && size(T,2) == 1;
    v = Skew(T(1:3));
    w = Skew(T(4:6));
else
    error('The input should be a twist represented by a matrix')
end
ad = [w v; zeros(3) w];
end