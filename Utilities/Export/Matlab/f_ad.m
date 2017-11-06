function ad = f_ad(T)

v = f_skew(T(1:3));
om = f_skew(T(4:6));
ad = [om, v; zeros(3,3) om];

end