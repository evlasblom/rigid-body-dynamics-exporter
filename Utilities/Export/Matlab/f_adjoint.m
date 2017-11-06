function Ad = f_adjoint(H)

R = H(1:3,1:3);
d = H(1:3,4);
Ad = [R f_skew(d)*R; zeros(3,3) R];

end