function iAd = f_inverseadjoint(H)

R = H(1:3,1:3);
d = H(1:3,4);
iAd = [R.' -R.'*f_skew(d); zeros(3,3) R.'];

end