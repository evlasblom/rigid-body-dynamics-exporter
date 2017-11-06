function [M,C,G] = Upperbody_M_C_G(q)


q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5); q6 = q(6); q7 = q(7); dq1 = q(8); dq2 = q(9); dq3 = q(10); dq4 = q(11); dq5 = q(12); dq6 = q(13); dq7 = q(14);



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [cos(q1),0,sin(q1),0;
	0,1,0,0;
	-sin(q1),0,cos(q1),0;
	0,0,0,1];

Hm1 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0.02;
	0,0,0,1];

H1_0 = H0_0*H1_0;

Jl1 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	1,0,0,0,0,0,0;
	0,0,0,0,0,0,0];

J1 = (Jl1);

Ja1 = f_adjoint(f_rbar(H1_0))*f_inverseadjoint(Hm1)*J1;

Im1 = [1,0,0,0,0,0;
	0,1,0,0,0,0;
	0,0,1,0,0,0;
	0,0,0,0.1,0,0;
	0,0,0,0,0.1,0;
	0,0,0,0,0,0.1];

I1 = (f_inverseadjoint(Hm1)).'*Im1*f_inverseadjoint(Hm1);

MM1 = ((J1).'*I1*J1);

dAd1 = [-sin(q1),0,-cos(q1),0,0,0;
	0,0,0,0,0,0;
	cos(q1),0,-sin(q1),0,0,0;
	0,0,0,-sin(q1),0,-cos(q1);
	0,0,0,0,0,0;
	0,0,0,cos(q1),0,-sin(q1)];

Fz1 = [0;
	0;
	9.81;
	0;
	0;
	0];

GG1 = ((Ja1).'*Fz1);

H2_1 = [1,0,0,0;
	0,cos(q2),-sin(q2),1/20;
	0,sin(q2),cos(q2),1/10;
	0,0,0,1];

Hm2 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H2_0 = H1_0*H2_1;

Jl2 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,1,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0];

J2 = (f_inverseadjoint(H2_1)*J1+Jl2);

Ja2 = f_adjoint(f_rbar(H2_0))*f_inverseadjoint(Hm2)*J2;

Im2 = [0.1,0,0,0,0,0;
	0,0.1,0,0,0,0;
	0,0,0.1,0,0,0;
	0,0,0,0.01,0,0;
	0,0,0,0,0.01,0;
	0,0,0,0,0,0.01];

I2 = (f_inverseadjoint(Hm2)).'*Im2*f_inverseadjoint(Hm2);

MM2 = (MM1+(J2).'*I2*J2);

dAd2 = [0,0,0,0,0,0;
	0,-sin(q2),cos(q2),cos(q2)/20 + sin(q2)/10,0,0;
	0,-cos(q2),-sin(q2),cos(q2)/10 - sin(q2)/20,0,0;
	0,0,0,0,0,0;
	0,0,0,0,-sin(q2),cos(q2);
	0,0,0,0,-cos(q2),-sin(q2)];

dJ2_2 = (dAd2*J1);

dMM2_2 = ((dJ2_2).'*I2*J2+(J2).'*I2*dJ2_2);

Fz2 = [0;
	0;
	0.981;
	0;
	0;
	0];

GG2 = (GG1+(Ja2).'*Fz2);

H3_2 = [cos(q3),-sin(q3),0,0;
	sin(q3),cos(q3),0,0;
	0,0,1,0;
	0,0,0,1];

Hm3 = [1,0,0,0;
	0,1,0,0.05;
	0,0,1,0;
	0,0,0,1];

H3_0 = H2_0*H3_2;

Jl3 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,1,0,0,0,0];

J3 = (f_inverseadjoint(H3_2)*J2+Jl3);

Ja3 = f_adjoint(f_rbar(H3_0))*f_inverseadjoint(Hm3)*J3;

Im3 = [0.1,0,0,0,0,0;
	0,0.1,0,0,0,0;
	0,0,0.1,0,0,0;
	0,0,0,0.01,0,0;
	0,0,0,0,0.01,0;
	0,0,0,0,0,0.01];

I3 = (f_inverseadjoint(Hm3)).'*Im3*f_inverseadjoint(Hm3);

MM3 = (MM2+(J3).'*I3*J3);

dJ3_2 = f_inverseadjoint(H3_2)*dJ2_2;

dMM3_2 = (dMM2_2+(dJ3_2).'*I3*J3+(J3).'*I3*dJ3_2);

dAd3 = [-sin(q3),cos(q3),0,0,0,0;
	-cos(q3),-sin(q3),0,0,0,0;
	0,0,0,0,0,0;
	0,0,0,-sin(q3),cos(q3),0;
	0,0,0,-cos(q3),-sin(q3),0;
	0,0,0,0,0,0];

dJ3_3 = (dAd3*J2);

dMM3_3 = ((dJ3_3).'*I3*J3+(J3).'*I3*dJ3_3);

Fz3 = [0;
	0;
	0.981;
	0;
	0;
	0];

GG3 = (GG2+(Ja3).'*Fz3);

H4_3 = [cos(q4),-sin(q4),0,0;
	sin(q4),cos(q4),0,1/10;
	0,0,1,0;
	0,0,0,1];

Hm4 = [1,0,0,0;
	0,1,0,0.05;
	0,0,1,0;
	0,0,0,1];

H4_0 = H3_0*H4_3;

Jl4 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,1,0,0,0];

J4 = (f_inverseadjoint(H4_3)*J3+Jl4);

Ja4 = f_adjoint(f_rbar(H4_0))*f_inverseadjoint(Hm4)*J4;

Im4 = [0.1,0,0,0,0,0;
	0,0.1,0,0,0,0;
	0,0,0.1,0,0,0;
	0,0,0,0.01,0,0;
	0,0,0,0,0.01,0;
	0,0,0,0,0,0.01];

I4 = (f_inverseadjoint(Hm4)).'*Im4*f_inverseadjoint(Hm4);

MM4 = (MM3+(J4).'*I4*J4);

dJ4_2 = f_inverseadjoint(H4_3)*dJ3_2;

dMM4_2 = (dMM3_2+(dJ4_2).'*I4*J4+(J4).'*I4*dJ4_2);

dJ4_3 = f_inverseadjoint(H4_3)*dJ3_3;

dMM4_3 = (dMM3_3+(dJ4_3).'*I4*J4+(J4).'*I4*dJ4_3);

dAd4 = [-sin(q4),cos(q4),0,0,0,sin(q4)/10;
	-cos(q4),-sin(q4),0,0,0,cos(q4)/10;
	0,0,0,0,0,0;
	0,0,0,-sin(q4),cos(q4),0;
	0,0,0,-cos(q4),-sin(q4),0;
	0,0,0,0,0,0];

dJ4_4 = (dAd4*J3);

dMM4_4 = ((dJ4_4).'*I4*J4+(J4).'*I4*dJ4_4);

Fz4 = [0;
	0;
	0.981;
	0;
	0;
	0];

GG4 = (GG3+(Ja4).'*Fz4);

H5_1 = [1,0,0,0;
	0,cos(q5),-sin(q5),-1/20;
	0,sin(q5),cos(q5),1/10;
	0,0,0,1];

Hm5 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H5_0 = H1_0*H5_1;

Jl5 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,1,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0];

J5 = (f_inverseadjoint(H5_1)*J1+Jl5);

Ja5 = f_adjoint(f_rbar(H5_0))*f_inverseadjoint(Hm5)*J5;

Im5 = [0.1,0,0,0,0,0;
	0,0.1,0,0,0,0;
	0,0,0.1,0,0,0;
	0,0,0,0.01,0,0;
	0,0,0,0,0.01,0;
	0,0,0,0,0,0.01];

I5 = (f_inverseadjoint(Hm5)).'*Im5*f_inverseadjoint(Hm5);

MM5 = (MM4+(J5).'*I5*J5);

dMM5_2 = (dMM4_2);

dMM5_3 = (dMM4_3);

dMM5_4 = (dMM4_4);

dAd5 = [0,0,0,0,0,0;
	0,-sin(q5),cos(q5),sin(q5)/10 - cos(q5)/20,0,0;
	0,-cos(q5),-sin(q5),cos(q5)/10 + sin(q5)/20,0,0;
	0,0,0,0,0,0;
	0,0,0,0,-sin(q5),cos(q5);
	0,0,0,0,-cos(q5),-sin(q5)];

dJ5_5 = (dAd5*J1);

dMM5_5 = ((dJ5_5).'*I5*J5+(J5).'*I5*dJ5_5);

Fz5 = [0;
	0;
	0.981;
	0;
	0;
	0];

GG5 = (GG4+(Ja5).'*Fz5);

H6_5 = [cos(q6),-sin(q6),0,0;
	sin(q6),cos(q6),0,0;
	0,0,1,0;
	0,0,0,1];

Hm6 = [1,0,0,0;
	0,1,0,-0.05;
	0,0,1,0;
	0,0,0,1];

H6_0 = H5_0*H6_5;

Jl6 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,1,0];

J6 = (f_inverseadjoint(H6_5)*J5+Jl6);

Ja6 = f_adjoint(f_rbar(H6_0))*f_inverseadjoint(Hm6)*J6;

Im6 = [0.1,0,0,0,0,0;
	0,0.1,0,0,0,0;
	0,0,0.1,0,0,0;
	0,0,0,0.01,0,0;
	0,0,0,0,0.01,0;
	0,0,0,0,0,0.01];

I6 = (f_inverseadjoint(Hm6)).'*Im6*f_inverseadjoint(Hm6);

MM6 = (MM5+(J6).'*I6*J6);

dMM6_2 = (dMM5_2);

dMM6_3 = (dMM5_3);

dMM6_4 = (dMM5_4);

dJ6_5 = f_inverseadjoint(H6_5)*dJ5_5;

dMM6_5 = (dMM5_5+(dJ6_5).'*I6*J6+(J6).'*I6*dJ6_5);

dAd6 = [-sin(q6),cos(q6),0,0,0,0;
	-cos(q6),-sin(q6),0,0,0,0;
	0,0,0,0,0,0;
	0,0,0,-sin(q6),cos(q6),0;
	0,0,0,-cos(q6),-sin(q6),0;
	0,0,0,0,0,0];

dJ6_6 = (dAd6*J5);

dMM6_6 = ((dJ6_6).'*I6*J6+(J6).'*I6*dJ6_6);

Fz6 = [0;
	0;
	0.981;
	0;
	0;
	0];

GG6 = (GG5+(Ja6).'*Fz6);

H7_6 = [cos(q7),-sin(q7),0,0;
	sin(q7),cos(q7),0,-1/10;
	0,0,1,0;
	0,0,0,1];

Hm7 = [1,0,0,0;
	0,1,0,-0.05;
	0,0,1,0;
	0,0,0,1];

H7_0 = H6_0*H7_6;

Jl7 = [0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,0;
	0,0,0,0,0,0,1];

J7 = (f_inverseadjoint(H7_6)*J6+Jl7);

Ja7 = f_adjoint(f_rbar(H7_0))*f_inverseadjoint(Hm7)*J7;

Im7 = [0.1,0,0,0,0,0;
	0,0.1,0,0,0,0;
	0,0,0.1,0,0,0;
	0,0,0,0.01,0,0;
	0,0,0,0,0.01,0;
	0,0,0,0,0,0.01];

I7 = (f_inverseadjoint(Hm7)).'*Im7*f_inverseadjoint(Hm7);

MM7 = (MM6+(J7).'*I7*J7);

dMM7_2 = (dMM6_2);

dMM7_3 = (dMM6_3);

dMM7_4 = (dMM6_4);

dJ7_5 = f_inverseadjoint(H7_6)*dJ6_5;

dMM7_5 = (dMM6_5+(dJ7_5).'*I7*J7+(J7).'*I7*dJ7_5);

dJ7_6 = f_inverseadjoint(H7_6)*dJ6_6;

dMM7_6 = (dMM6_6+(dJ7_6).'*I7*J7+(J7).'*I7*dJ7_6);

dAd7 = [-sin(q7),cos(q7),0,0,0,-sin(q7)/10;
	-cos(q7),-sin(q7),0,0,0,-cos(q7)/10;
	0,0,0,0,0,0;
	0,0,0,-sin(q7),cos(q7),0;
	0,0,0,-cos(q7),-sin(q7),0;
	0,0,0,0,0,0];

dJ7_7 = (dAd7*J6);

dMM7_7 = ((dJ7_7).'*I7*J7+(J7).'*I7*dJ7_7);

Fz7 = [0;
	0;
	0.981;
	0;
	0;
	0];

GG7 = (GG6+(Ja7).'*Fz7);

CC(1,1) = 0 + 0.5*(+dMM7_2(1,1))*dq2 + 0.5*(+dMM7_3(1,1))*dq3 + 0.5*(+dMM7_4(1,1))*dq4 + 0.5*(+dMM7_5(1,1))*dq5 + 0.5*(+dMM7_6(1,1))*dq6 + 0.5*(+dMM7_7(1,1))*dq7 ;

CC(1,2) = 0.5*(+dMM7_2(1,1))*dq1 + 0.5*(+dMM7_2(1,2)+dMM7_2(1,2))*dq2 + 0.5*(+dMM7_3(1,2)+dMM7_2(1,3))*dq3 + 0.5*(+dMM7_4(1,2)+dMM7_2(1,4))*dq4 + 0.5*(+dMM7_5(1,2)+dMM7_2(1,5))*dq5 + 0.5*(+dMM7_6(1,2)+dMM7_2(1,6))*dq6 + 0.5*(+dMM7_7(1,2)+dMM7_2(1,7))*dq7 ;

CC(1,3) = 0.5*(+dMM7_3(1,1))*dq1 + 0.5*(+dMM7_2(1,3)+dMM7_3(1,2))*dq2 + 0.5*(+dMM7_3(1,3)+dMM7_3(1,3))*dq3 + 0.5*(+dMM7_4(1,3)+dMM7_3(1,4))*dq4 + 0.5*(+dMM7_5(1,3)+dMM7_3(1,5))*dq5 + 0.5*(+dMM7_6(1,3)+dMM7_3(1,6))*dq6 + 0.5*(+dMM7_7(1,3)+dMM7_3(1,7))*dq7 ;

CC(1,4) = 0.5*(+dMM7_4(1,1))*dq1 + 0.5*(+dMM7_2(1,4)+dMM7_4(1,2))*dq2 + 0.5*(+dMM7_3(1,4)+dMM7_4(1,3))*dq3 + 0.5*(+dMM7_4(1,4)+dMM7_4(1,4))*dq4 + 0.5*(+dMM7_5(1,4)+dMM7_4(1,5))*dq5 + 0.5*(+dMM7_6(1,4)+dMM7_4(1,6))*dq6 + 0.5*(+dMM7_7(1,4)+dMM7_4(1,7))*dq7 ;

CC(1,5) = 0.5*(+dMM7_5(1,1))*dq1 + 0.5*(+dMM7_2(1,5)+dMM7_5(1,2))*dq2 + 0.5*(+dMM7_3(1,5)+dMM7_5(1,3))*dq3 + 0.5*(+dMM7_4(1,5)+dMM7_5(1,4))*dq4 + 0.5*(+dMM7_5(1,5)+dMM7_5(1,5))*dq5 + 0.5*(+dMM7_6(1,5)+dMM7_5(1,6))*dq6 + 0.5*(+dMM7_7(1,5)+dMM7_5(1,7))*dq7 ;

CC(1,6) = 0.5*(+dMM7_6(1,1))*dq1 + 0.5*(+dMM7_2(1,6)+dMM7_6(1,2))*dq2 + 0.5*(+dMM7_3(1,6)+dMM7_6(1,3))*dq3 + 0.5*(+dMM7_4(1,6)+dMM7_6(1,4))*dq4 + 0.5*(+dMM7_5(1,6)+dMM7_6(1,5))*dq5 + 0.5*(+dMM7_6(1,6)+dMM7_6(1,6))*dq6 + 0.5*(+dMM7_7(1,6)+dMM7_6(1,7))*dq7 ;

CC(1,7) = 0.5*(+dMM7_7(1,1))*dq1 + 0.5*(+dMM7_2(1,7)+dMM7_7(1,2))*dq2 + 0.5*(+dMM7_3(1,7)+dMM7_7(1,3))*dq3 + 0.5*(+dMM7_4(1,7)+dMM7_7(1,4))*dq4 + 0.5*(+dMM7_5(1,7)+dMM7_7(1,5))*dq5 + 0.5*(+dMM7_6(1,7)+dMM7_7(1,6))*dq6 + 0.5*(+dMM7_7(1,7)+dMM7_7(1,7))*dq7 ;

CC(2,1) = 0.5*(-dMM7_2(1,1))*dq1 + 0.5*(+dMM7_2(2,1)-dMM7_2(2,1))*dq2 + 0.5*(+dMM7_3(2,1)-dMM7_2(3,1))*dq3 + 0.5*(+dMM7_4(2,1)-dMM7_2(4,1))*dq4 + 0.5*(+dMM7_5(2,1)-dMM7_2(5,1))*dq5 + 0.5*(+dMM7_6(2,1)-dMM7_2(6,1))*dq6 + 0.5*(+dMM7_7(2,1)-dMM7_2(7,1))*dq7 ;

CC(2,2) = 0.5*(+dMM7_2(2,1)-dMM7_2(1,2))*dq1 + 0.5*(+dMM7_2(2,2)+dMM7_2(2,2)-dMM7_2(2,2))*dq2 + 0.5*(+dMM7_3(2,2)+dMM7_2(2,3)-dMM7_2(3,2))*dq3 + 0.5*(+dMM7_4(2,2)+dMM7_2(2,4)-dMM7_2(4,2))*dq4 + 0.5*(+dMM7_5(2,2)+dMM7_2(2,5)-dMM7_2(5,2))*dq5 + 0.5*(+dMM7_6(2,2)+dMM7_2(2,6)-dMM7_2(6,2))*dq6 + 0.5*(+dMM7_7(2,2)+dMM7_2(2,7)-dMM7_2(7,2))*dq7 ;

CC(2,3) = 0.5*(+dMM7_3(2,1)-dMM7_2(1,3))*dq1 + 0.5*(+dMM7_2(2,3)+dMM7_3(2,2)-dMM7_2(2,3))*dq2 + 0.5*(+dMM7_3(2,3)+dMM7_3(2,3)-dMM7_2(3,3))*dq3 + 0.5*(+dMM7_4(2,3)+dMM7_3(2,4)-dMM7_2(4,3))*dq4 + 0.5*(+dMM7_5(2,3)+dMM7_3(2,5)-dMM7_2(5,3))*dq5 + 0.5*(+dMM7_6(2,3)+dMM7_3(2,6)-dMM7_2(6,3))*dq6 + 0.5*(+dMM7_7(2,3)+dMM7_3(2,7)-dMM7_2(7,3))*dq7 ;

CC(2,4) = 0.5*(+dMM7_4(2,1)-dMM7_2(1,4))*dq1 + 0.5*(+dMM7_2(2,4)+dMM7_4(2,2)-dMM7_2(2,4))*dq2 + 0.5*(+dMM7_3(2,4)+dMM7_4(2,3)-dMM7_2(3,4))*dq3 + 0.5*(+dMM7_4(2,4)+dMM7_4(2,4)-dMM7_2(4,4))*dq4 + 0.5*(+dMM7_5(2,4)+dMM7_4(2,5)-dMM7_2(5,4))*dq5 + 0.5*(+dMM7_6(2,4)+dMM7_4(2,6)-dMM7_2(6,4))*dq6 + 0.5*(+dMM7_7(2,4)+dMM7_4(2,7)-dMM7_2(7,4))*dq7 ;

CC(2,5) = 0.5*(+dMM7_5(2,1)-dMM7_2(1,5))*dq1 + 0.5*(+dMM7_2(2,5)+dMM7_5(2,2)-dMM7_2(2,5))*dq2 + 0.5*(+dMM7_3(2,5)+dMM7_5(2,3)-dMM7_2(3,5))*dq3 + 0.5*(+dMM7_4(2,5)+dMM7_5(2,4)-dMM7_2(4,5))*dq4 + 0.5*(+dMM7_5(2,5)+dMM7_5(2,5)-dMM7_2(5,5))*dq5 + 0.5*(+dMM7_6(2,5)+dMM7_5(2,6)-dMM7_2(6,5))*dq6 + 0.5*(+dMM7_7(2,5)+dMM7_5(2,7)-dMM7_2(7,5))*dq7 ;

CC(2,6) = 0.5*(+dMM7_6(2,1)-dMM7_2(1,6))*dq1 + 0.5*(+dMM7_2(2,6)+dMM7_6(2,2)-dMM7_2(2,6))*dq2 + 0.5*(+dMM7_3(2,6)+dMM7_6(2,3)-dMM7_2(3,6))*dq3 + 0.5*(+dMM7_4(2,6)+dMM7_6(2,4)-dMM7_2(4,6))*dq4 + 0.5*(+dMM7_5(2,6)+dMM7_6(2,5)-dMM7_2(5,6))*dq5 + 0.5*(+dMM7_6(2,6)+dMM7_6(2,6)-dMM7_2(6,6))*dq6 + 0.5*(+dMM7_7(2,6)+dMM7_6(2,7)-dMM7_2(7,6))*dq7 ;

CC(2,7) = 0.5*(+dMM7_7(2,1)-dMM7_2(1,7))*dq1 + 0.5*(+dMM7_2(2,7)+dMM7_7(2,2)-dMM7_2(2,7))*dq2 + 0.5*(+dMM7_3(2,7)+dMM7_7(2,3)-dMM7_2(3,7))*dq3 + 0.5*(+dMM7_4(2,7)+dMM7_7(2,4)-dMM7_2(4,7))*dq4 + 0.5*(+dMM7_5(2,7)+dMM7_7(2,5)-dMM7_2(5,7))*dq5 + 0.5*(+dMM7_6(2,7)+dMM7_7(2,6)-dMM7_2(6,7))*dq6 + 0.5*(+dMM7_7(2,7)+dMM7_7(2,7)-dMM7_2(7,7))*dq7 ;

CC(3,1) = 0.5*(-dMM7_3(1,1))*dq1 + 0.5*(+dMM7_2(3,1)-dMM7_3(2,1))*dq2 + 0.5*(+dMM7_3(3,1)-dMM7_3(3,1))*dq3 + 0.5*(+dMM7_4(3,1)-dMM7_3(4,1))*dq4 + 0.5*(+dMM7_5(3,1)-dMM7_3(5,1))*dq5 + 0.5*(+dMM7_6(3,1)-dMM7_3(6,1))*dq6 + 0.5*(+dMM7_7(3,1)-dMM7_3(7,1))*dq7 ;

CC(3,2) = 0.5*(+dMM7_2(3,1)-dMM7_3(1,2))*dq1 + 0.5*(+dMM7_2(3,2)+dMM7_2(3,2)-dMM7_3(2,2))*dq2 + 0.5*(+dMM7_3(3,2)+dMM7_2(3,3)-dMM7_3(3,2))*dq3 + 0.5*(+dMM7_4(3,2)+dMM7_2(3,4)-dMM7_3(4,2))*dq4 + 0.5*(+dMM7_5(3,2)+dMM7_2(3,5)-dMM7_3(5,2))*dq5 + 0.5*(+dMM7_6(3,2)+dMM7_2(3,6)-dMM7_3(6,2))*dq6 + 0.5*(+dMM7_7(3,2)+dMM7_2(3,7)-dMM7_3(7,2))*dq7 ;

CC(3,3) = 0.5*(+dMM7_3(3,1)-dMM7_3(1,3))*dq1 + 0.5*(+dMM7_2(3,3)+dMM7_3(3,2)-dMM7_3(2,3))*dq2 + 0.5*(+dMM7_3(3,3)+dMM7_3(3,3)-dMM7_3(3,3))*dq3 + 0.5*(+dMM7_4(3,3)+dMM7_3(3,4)-dMM7_3(4,3))*dq4 + 0.5*(+dMM7_5(3,3)+dMM7_3(3,5)-dMM7_3(5,3))*dq5 + 0.5*(+dMM7_6(3,3)+dMM7_3(3,6)-dMM7_3(6,3))*dq6 + 0.5*(+dMM7_7(3,3)+dMM7_3(3,7)-dMM7_3(7,3))*dq7 ;

CC(3,4) = 0.5*(+dMM7_4(3,1)-dMM7_3(1,4))*dq1 + 0.5*(+dMM7_2(3,4)+dMM7_4(3,2)-dMM7_3(2,4))*dq2 + 0.5*(+dMM7_3(3,4)+dMM7_4(3,3)-dMM7_3(3,4))*dq3 + 0.5*(+dMM7_4(3,4)+dMM7_4(3,4)-dMM7_3(4,4))*dq4 + 0.5*(+dMM7_5(3,4)+dMM7_4(3,5)-dMM7_3(5,4))*dq5 + 0.5*(+dMM7_6(3,4)+dMM7_4(3,6)-dMM7_3(6,4))*dq6 + 0.5*(+dMM7_7(3,4)+dMM7_4(3,7)-dMM7_3(7,4))*dq7 ;

CC(3,5) = 0.5*(+dMM7_5(3,1)-dMM7_3(1,5))*dq1 + 0.5*(+dMM7_2(3,5)+dMM7_5(3,2)-dMM7_3(2,5))*dq2 + 0.5*(+dMM7_3(3,5)+dMM7_5(3,3)-dMM7_3(3,5))*dq3 + 0.5*(+dMM7_4(3,5)+dMM7_5(3,4)-dMM7_3(4,5))*dq4 + 0.5*(+dMM7_5(3,5)+dMM7_5(3,5)-dMM7_3(5,5))*dq5 + 0.5*(+dMM7_6(3,5)+dMM7_5(3,6)-dMM7_3(6,5))*dq6 + 0.5*(+dMM7_7(3,5)+dMM7_5(3,7)-dMM7_3(7,5))*dq7 ;

CC(3,6) = 0.5*(+dMM7_6(3,1)-dMM7_3(1,6))*dq1 + 0.5*(+dMM7_2(3,6)+dMM7_6(3,2)-dMM7_3(2,6))*dq2 + 0.5*(+dMM7_3(3,6)+dMM7_6(3,3)-dMM7_3(3,6))*dq3 + 0.5*(+dMM7_4(3,6)+dMM7_6(3,4)-dMM7_3(4,6))*dq4 + 0.5*(+dMM7_5(3,6)+dMM7_6(3,5)-dMM7_3(5,6))*dq5 + 0.5*(+dMM7_6(3,6)+dMM7_6(3,6)-dMM7_3(6,6))*dq6 + 0.5*(+dMM7_7(3,6)+dMM7_6(3,7)-dMM7_3(7,6))*dq7 ;

CC(3,7) = 0.5*(+dMM7_7(3,1)-dMM7_3(1,7))*dq1 + 0.5*(+dMM7_2(3,7)+dMM7_7(3,2)-dMM7_3(2,7))*dq2 + 0.5*(+dMM7_3(3,7)+dMM7_7(3,3)-dMM7_3(3,7))*dq3 + 0.5*(+dMM7_4(3,7)+dMM7_7(3,4)-dMM7_3(4,7))*dq4 + 0.5*(+dMM7_5(3,7)+dMM7_7(3,5)-dMM7_3(5,7))*dq5 + 0.5*(+dMM7_6(3,7)+dMM7_7(3,6)-dMM7_3(6,7))*dq6 + 0.5*(+dMM7_7(3,7)+dMM7_7(3,7)-dMM7_3(7,7))*dq7 ;

CC(4,1) = 0.5*(-dMM7_4(1,1))*dq1 + 0.5*(+dMM7_2(4,1)-dMM7_4(2,1))*dq2 + 0.5*(+dMM7_3(4,1)-dMM7_4(3,1))*dq3 + 0.5*(+dMM7_4(4,1)-dMM7_4(4,1))*dq4 + 0.5*(+dMM7_5(4,1)-dMM7_4(5,1))*dq5 + 0.5*(+dMM7_6(4,1)-dMM7_4(6,1))*dq6 + 0.5*(+dMM7_7(4,1)-dMM7_4(7,1))*dq7 ;

CC(4,2) = 0.5*(+dMM7_2(4,1)-dMM7_4(1,2))*dq1 + 0.5*(+dMM7_2(4,2)+dMM7_2(4,2)-dMM7_4(2,2))*dq2 + 0.5*(+dMM7_3(4,2)+dMM7_2(4,3)-dMM7_4(3,2))*dq3 + 0.5*(+dMM7_4(4,2)+dMM7_2(4,4)-dMM7_4(4,2))*dq4 + 0.5*(+dMM7_5(4,2)+dMM7_2(4,5)-dMM7_4(5,2))*dq5 + 0.5*(+dMM7_6(4,2)+dMM7_2(4,6)-dMM7_4(6,2))*dq6 + 0.5*(+dMM7_7(4,2)+dMM7_2(4,7)-dMM7_4(7,2))*dq7 ;

CC(4,3) = 0.5*(+dMM7_3(4,1)-dMM7_4(1,3))*dq1 + 0.5*(+dMM7_2(4,3)+dMM7_3(4,2)-dMM7_4(2,3))*dq2 + 0.5*(+dMM7_3(4,3)+dMM7_3(4,3)-dMM7_4(3,3))*dq3 + 0.5*(+dMM7_4(4,3)+dMM7_3(4,4)-dMM7_4(4,3))*dq4 + 0.5*(+dMM7_5(4,3)+dMM7_3(4,5)-dMM7_4(5,3))*dq5 + 0.5*(+dMM7_6(4,3)+dMM7_3(4,6)-dMM7_4(6,3))*dq6 + 0.5*(+dMM7_7(4,3)+dMM7_3(4,7)-dMM7_4(7,3))*dq7 ;

CC(4,4) = 0.5*(+dMM7_4(4,1)-dMM7_4(1,4))*dq1 + 0.5*(+dMM7_2(4,4)+dMM7_4(4,2)-dMM7_4(2,4))*dq2 + 0.5*(+dMM7_3(4,4)+dMM7_4(4,3)-dMM7_4(3,4))*dq3 + 0.5*(+dMM7_4(4,4)+dMM7_4(4,4)-dMM7_4(4,4))*dq4 + 0.5*(+dMM7_5(4,4)+dMM7_4(4,5)-dMM7_4(5,4))*dq5 + 0.5*(+dMM7_6(4,4)+dMM7_4(4,6)-dMM7_4(6,4))*dq6 + 0.5*(+dMM7_7(4,4)+dMM7_4(4,7)-dMM7_4(7,4))*dq7 ;

CC(4,5) = 0.5*(+dMM7_5(4,1)-dMM7_4(1,5))*dq1 + 0.5*(+dMM7_2(4,5)+dMM7_5(4,2)-dMM7_4(2,5))*dq2 + 0.5*(+dMM7_3(4,5)+dMM7_5(4,3)-dMM7_4(3,5))*dq3 + 0.5*(+dMM7_4(4,5)+dMM7_5(4,4)-dMM7_4(4,5))*dq4 + 0.5*(+dMM7_5(4,5)+dMM7_5(4,5)-dMM7_4(5,5))*dq5 + 0.5*(+dMM7_6(4,5)+dMM7_5(4,6)-dMM7_4(6,5))*dq6 + 0.5*(+dMM7_7(4,5)+dMM7_5(4,7)-dMM7_4(7,5))*dq7 ;

CC(4,6) = 0.5*(+dMM7_6(4,1)-dMM7_4(1,6))*dq1 + 0.5*(+dMM7_2(4,6)+dMM7_6(4,2)-dMM7_4(2,6))*dq2 + 0.5*(+dMM7_3(4,6)+dMM7_6(4,3)-dMM7_4(3,6))*dq3 + 0.5*(+dMM7_4(4,6)+dMM7_6(4,4)-dMM7_4(4,6))*dq4 + 0.5*(+dMM7_5(4,6)+dMM7_6(4,5)-dMM7_4(5,6))*dq5 + 0.5*(+dMM7_6(4,6)+dMM7_6(4,6)-dMM7_4(6,6))*dq6 + 0.5*(+dMM7_7(4,6)+dMM7_6(4,7)-dMM7_4(7,6))*dq7 ;

CC(4,7) = 0.5*(+dMM7_7(4,1)-dMM7_4(1,7))*dq1 + 0.5*(+dMM7_2(4,7)+dMM7_7(4,2)-dMM7_4(2,7))*dq2 + 0.5*(+dMM7_3(4,7)+dMM7_7(4,3)-dMM7_4(3,7))*dq3 + 0.5*(+dMM7_4(4,7)+dMM7_7(4,4)-dMM7_4(4,7))*dq4 + 0.5*(+dMM7_5(4,7)+dMM7_7(4,5)-dMM7_4(5,7))*dq5 + 0.5*(+dMM7_6(4,7)+dMM7_7(4,6)-dMM7_4(6,7))*dq6 + 0.5*(+dMM7_7(4,7)+dMM7_7(4,7)-dMM7_4(7,7))*dq7 ;

CC(5,1) = 0.5*(-dMM7_5(1,1))*dq1 + 0.5*(+dMM7_2(5,1)-dMM7_5(2,1))*dq2 + 0.5*(+dMM7_3(5,1)-dMM7_5(3,1))*dq3 + 0.5*(+dMM7_4(5,1)-dMM7_5(4,1))*dq4 + 0.5*(+dMM7_5(5,1)-dMM7_5(5,1))*dq5 + 0.5*(+dMM7_6(5,1)-dMM7_5(6,1))*dq6 + 0.5*(+dMM7_7(5,1)-dMM7_5(7,1))*dq7 ;

CC(5,2) = 0.5*(+dMM7_2(5,1)-dMM7_5(1,2))*dq1 + 0.5*(+dMM7_2(5,2)+dMM7_2(5,2)-dMM7_5(2,2))*dq2 + 0.5*(+dMM7_3(5,2)+dMM7_2(5,3)-dMM7_5(3,2))*dq3 + 0.5*(+dMM7_4(5,2)+dMM7_2(5,4)-dMM7_5(4,2))*dq4 + 0.5*(+dMM7_5(5,2)+dMM7_2(5,5)-dMM7_5(5,2))*dq5 + 0.5*(+dMM7_6(5,2)+dMM7_2(5,6)-dMM7_5(6,2))*dq6 + 0.5*(+dMM7_7(5,2)+dMM7_2(5,7)-dMM7_5(7,2))*dq7 ;

CC(5,3) = 0.5*(+dMM7_3(5,1)-dMM7_5(1,3))*dq1 + 0.5*(+dMM7_2(5,3)+dMM7_3(5,2)-dMM7_5(2,3))*dq2 + 0.5*(+dMM7_3(5,3)+dMM7_3(5,3)-dMM7_5(3,3))*dq3 + 0.5*(+dMM7_4(5,3)+dMM7_3(5,4)-dMM7_5(4,3))*dq4 + 0.5*(+dMM7_5(5,3)+dMM7_3(5,5)-dMM7_5(5,3))*dq5 + 0.5*(+dMM7_6(5,3)+dMM7_3(5,6)-dMM7_5(6,3))*dq6 + 0.5*(+dMM7_7(5,3)+dMM7_3(5,7)-dMM7_5(7,3))*dq7 ;

CC(5,4) = 0.5*(+dMM7_4(5,1)-dMM7_5(1,4))*dq1 + 0.5*(+dMM7_2(5,4)+dMM7_4(5,2)-dMM7_5(2,4))*dq2 + 0.5*(+dMM7_3(5,4)+dMM7_4(5,3)-dMM7_5(3,4))*dq3 + 0.5*(+dMM7_4(5,4)+dMM7_4(5,4)-dMM7_5(4,4))*dq4 + 0.5*(+dMM7_5(5,4)+dMM7_4(5,5)-dMM7_5(5,4))*dq5 + 0.5*(+dMM7_6(5,4)+dMM7_4(5,6)-dMM7_5(6,4))*dq6 + 0.5*(+dMM7_7(5,4)+dMM7_4(5,7)-dMM7_5(7,4))*dq7 ;

CC(5,5) = 0.5*(+dMM7_5(5,1)-dMM7_5(1,5))*dq1 + 0.5*(+dMM7_2(5,5)+dMM7_5(5,2)-dMM7_5(2,5))*dq2 + 0.5*(+dMM7_3(5,5)+dMM7_5(5,3)-dMM7_5(3,5))*dq3 + 0.5*(+dMM7_4(5,5)+dMM7_5(5,4)-dMM7_5(4,5))*dq4 + 0.5*(+dMM7_5(5,5)+dMM7_5(5,5)-dMM7_5(5,5))*dq5 + 0.5*(+dMM7_6(5,5)+dMM7_5(5,6)-dMM7_5(6,5))*dq6 + 0.5*(+dMM7_7(5,5)+dMM7_5(5,7)-dMM7_5(7,5))*dq7 ;

CC(5,6) = 0.5*(+dMM7_6(5,1)-dMM7_5(1,6))*dq1 + 0.5*(+dMM7_2(5,6)+dMM7_6(5,2)-dMM7_5(2,6))*dq2 + 0.5*(+dMM7_3(5,6)+dMM7_6(5,3)-dMM7_5(3,6))*dq3 + 0.5*(+dMM7_4(5,6)+dMM7_6(5,4)-dMM7_5(4,6))*dq4 + 0.5*(+dMM7_5(5,6)+dMM7_6(5,5)-dMM7_5(5,6))*dq5 + 0.5*(+dMM7_6(5,6)+dMM7_6(5,6)-dMM7_5(6,6))*dq6 + 0.5*(+dMM7_7(5,6)+dMM7_6(5,7)-dMM7_5(7,6))*dq7 ;

CC(5,7) = 0.5*(+dMM7_7(5,1)-dMM7_5(1,7))*dq1 + 0.5*(+dMM7_2(5,7)+dMM7_7(5,2)-dMM7_5(2,7))*dq2 + 0.5*(+dMM7_3(5,7)+dMM7_7(5,3)-dMM7_5(3,7))*dq3 + 0.5*(+dMM7_4(5,7)+dMM7_7(5,4)-dMM7_5(4,7))*dq4 + 0.5*(+dMM7_5(5,7)+dMM7_7(5,5)-dMM7_5(5,7))*dq5 + 0.5*(+dMM7_6(5,7)+dMM7_7(5,6)-dMM7_5(6,7))*dq6 + 0.5*(+dMM7_7(5,7)+dMM7_7(5,7)-dMM7_5(7,7))*dq7 ;

CC(6,1) = 0.5*(-dMM7_6(1,1))*dq1 + 0.5*(+dMM7_2(6,1)-dMM7_6(2,1))*dq2 + 0.5*(+dMM7_3(6,1)-dMM7_6(3,1))*dq3 + 0.5*(+dMM7_4(6,1)-dMM7_6(4,1))*dq4 + 0.5*(+dMM7_5(6,1)-dMM7_6(5,1))*dq5 + 0.5*(+dMM7_6(6,1)-dMM7_6(6,1))*dq6 + 0.5*(+dMM7_7(6,1)-dMM7_6(7,1))*dq7 ;

CC(6,2) = 0.5*(+dMM7_2(6,1)-dMM7_6(1,2))*dq1 + 0.5*(+dMM7_2(6,2)+dMM7_2(6,2)-dMM7_6(2,2))*dq2 + 0.5*(+dMM7_3(6,2)+dMM7_2(6,3)-dMM7_6(3,2))*dq3 + 0.5*(+dMM7_4(6,2)+dMM7_2(6,4)-dMM7_6(4,2))*dq4 + 0.5*(+dMM7_5(6,2)+dMM7_2(6,5)-dMM7_6(5,2))*dq5 + 0.5*(+dMM7_6(6,2)+dMM7_2(6,6)-dMM7_6(6,2))*dq6 + 0.5*(+dMM7_7(6,2)+dMM7_2(6,7)-dMM7_6(7,2))*dq7 ;

CC(6,3) = 0.5*(+dMM7_3(6,1)-dMM7_6(1,3))*dq1 + 0.5*(+dMM7_2(6,3)+dMM7_3(6,2)-dMM7_6(2,3))*dq2 + 0.5*(+dMM7_3(6,3)+dMM7_3(6,3)-dMM7_6(3,3))*dq3 + 0.5*(+dMM7_4(6,3)+dMM7_3(6,4)-dMM7_6(4,3))*dq4 + 0.5*(+dMM7_5(6,3)+dMM7_3(6,5)-dMM7_6(5,3))*dq5 + 0.5*(+dMM7_6(6,3)+dMM7_3(6,6)-dMM7_6(6,3))*dq6 + 0.5*(+dMM7_7(6,3)+dMM7_3(6,7)-dMM7_6(7,3))*dq7 ;

CC(6,4) = 0.5*(+dMM7_4(6,1)-dMM7_6(1,4))*dq1 + 0.5*(+dMM7_2(6,4)+dMM7_4(6,2)-dMM7_6(2,4))*dq2 + 0.5*(+dMM7_3(6,4)+dMM7_4(6,3)-dMM7_6(3,4))*dq3 + 0.5*(+dMM7_4(6,4)+dMM7_4(6,4)-dMM7_6(4,4))*dq4 + 0.5*(+dMM7_5(6,4)+dMM7_4(6,5)-dMM7_6(5,4))*dq5 + 0.5*(+dMM7_6(6,4)+dMM7_4(6,6)-dMM7_6(6,4))*dq6 + 0.5*(+dMM7_7(6,4)+dMM7_4(6,7)-dMM7_6(7,4))*dq7 ;

CC(6,5) = 0.5*(+dMM7_5(6,1)-dMM7_6(1,5))*dq1 + 0.5*(+dMM7_2(6,5)+dMM7_5(6,2)-dMM7_6(2,5))*dq2 + 0.5*(+dMM7_3(6,5)+dMM7_5(6,3)-dMM7_6(3,5))*dq3 + 0.5*(+dMM7_4(6,5)+dMM7_5(6,4)-dMM7_6(4,5))*dq4 + 0.5*(+dMM7_5(6,5)+dMM7_5(6,5)-dMM7_6(5,5))*dq5 + 0.5*(+dMM7_6(6,5)+dMM7_5(6,6)-dMM7_6(6,5))*dq6 + 0.5*(+dMM7_7(6,5)+dMM7_5(6,7)-dMM7_6(7,5))*dq7 ;

CC(6,6) = 0.5*(+dMM7_6(6,1)-dMM7_6(1,6))*dq1 + 0.5*(+dMM7_2(6,6)+dMM7_6(6,2)-dMM7_6(2,6))*dq2 + 0.5*(+dMM7_3(6,6)+dMM7_6(6,3)-dMM7_6(3,6))*dq3 + 0.5*(+dMM7_4(6,6)+dMM7_6(6,4)-dMM7_6(4,6))*dq4 + 0.5*(+dMM7_5(6,6)+dMM7_6(6,5)-dMM7_6(5,6))*dq5 + 0.5*(+dMM7_6(6,6)+dMM7_6(6,6)-dMM7_6(6,6))*dq6 + 0.5*(+dMM7_7(6,6)+dMM7_6(6,7)-dMM7_6(7,6))*dq7 ;

CC(6,7) = 0.5*(+dMM7_7(6,1)-dMM7_6(1,7))*dq1 + 0.5*(+dMM7_2(6,7)+dMM7_7(6,2)-dMM7_6(2,7))*dq2 + 0.5*(+dMM7_3(6,7)+dMM7_7(6,3)-dMM7_6(3,7))*dq3 + 0.5*(+dMM7_4(6,7)+dMM7_7(6,4)-dMM7_6(4,7))*dq4 + 0.5*(+dMM7_5(6,7)+dMM7_7(6,5)-dMM7_6(5,7))*dq5 + 0.5*(+dMM7_6(6,7)+dMM7_7(6,6)-dMM7_6(6,7))*dq6 + 0.5*(+dMM7_7(6,7)+dMM7_7(6,7)-dMM7_6(7,7))*dq7 ;

CC(7,1) = 0.5*(-dMM7_7(1,1))*dq1 + 0.5*(+dMM7_2(7,1)-dMM7_7(2,1))*dq2 + 0.5*(+dMM7_3(7,1)-dMM7_7(3,1))*dq3 + 0.5*(+dMM7_4(7,1)-dMM7_7(4,1))*dq4 + 0.5*(+dMM7_5(7,1)-dMM7_7(5,1))*dq5 + 0.5*(+dMM7_6(7,1)-dMM7_7(6,1))*dq6 + 0.5*(+dMM7_7(7,1)-dMM7_7(7,1))*dq7 ;

CC(7,2) = 0.5*(+dMM7_2(7,1)-dMM7_7(1,2))*dq1 + 0.5*(+dMM7_2(7,2)+dMM7_2(7,2)-dMM7_7(2,2))*dq2 + 0.5*(+dMM7_3(7,2)+dMM7_2(7,3)-dMM7_7(3,2))*dq3 + 0.5*(+dMM7_4(7,2)+dMM7_2(7,4)-dMM7_7(4,2))*dq4 + 0.5*(+dMM7_5(7,2)+dMM7_2(7,5)-dMM7_7(5,2))*dq5 + 0.5*(+dMM7_6(7,2)+dMM7_2(7,6)-dMM7_7(6,2))*dq6 + 0.5*(+dMM7_7(7,2)+dMM7_2(7,7)-dMM7_7(7,2))*dq7 ;

CC(7,3) = 0.5*(+dMM7_3(7,1)-dMM7_7(1,3))*dq1 + 0.5*(+dMM7_2(7,3)+dMM7_3(7,2)-dMM7_7(2,3))*dq2 + 0.5*(+dMM7_3(7,3)+dMM7_3(7,3)-dMM7_7(3,3))*dq3 + 0.5*(+dMM7_4(7,3)+dMM7_3(7,4)-dMM7_7(4,3))*dq4 + 0.5*(+dMM7_5(7,3)+dMM7_3(7,5)-dMM7_7(5,3))*dq5 + 0.5*(+dMM7_6(7,3)+dMM7_3(7,6)-dMM7_7(6,3))*dq6 + 0.5*(+dMM7_7(7,3)+dMM7_3(7,7)-dMM7_7(7,3))*dq7 ;

CC(7,4) = 0.5*(+dMM7_4(7,1)-dMM7_7(1,4))*dq1 + 0.5*(+dMM7_2(7,4)+dMM7_4(7,2)-dMM7_7(2,4))*dq2 + 0.5*(+dMM7_3(7,4)+dMM7_4(7,3)-dMM7_7(3,4))*dq3 + 0.5*(+dMM7_4(7,4)+dMM7_4(7,4)-dMM7_7(4,4))*dq4 + 0.5*(+dMM7_5(7,4)+dMM7_4(7,5)-dMM7_7(5,4))*dq5 + 0.5*(+dMM7_6(7,4)+dMM7_4(7,6)-dMM7_7(6,4))*dq6 + 0.5*(+dMM7_7(7,4)+dMM7_4(7,7)-dMM7_7(7,4))*dq7 ;

CC(7,5) = 0.5*(+dMM7_5(7,1)-dMM7_7(1,5))*dq1 + 0.5*(+dMM7_2(7,5)+dMM7_5(7,2)-dMM7_7(2,5))*dq2 + 0.5*(+dMM7_3(7,5)+dMM7_5(7,3)-dMM7_7(3,5))*dq3 + 0.5*(+dMM7_4(7,5)+dMM7_5(7,4)-dMM7_7(4,5))*dq4 + 0.5*(+dMM7_5(7,5)+dMM7_5(7,5)-dMM7_7(5,5))*dq5 + 0.5*(+dMM7_6(7,5)+dMM7_5(7,6)-dMM7_7(6,5))*dq6 + 0.5*(+dMM7_7(7,5)+dMM7_5(7,7)-dMM7_7(7,5))*dq7 ;

CC(7,6) = 0.5*(+dMM7_6(7,1)-dMM7_7(1,6))*dq1 + 0.5*(+dMM7_2(7,6)+dMM7_6(7,2)-dMM7_7(2,6))*dq2 + 0.5*(+dMM7_3(7,6)+dMM7_6(7,3)-dMM7_7(3,6))*dq3 + 0.5*(+dMM7_4(7,6)+dMM7_6(7,4)-dMM7_7(4,6))*dq4 + 0.5*(+dMM7_5(7,6)+dMM7_6(7,5)-dMM7_7(5,6))*dq5 + 0.5*(+dMM7_6(7,6)+dMM7_6(7,6)-dMM7_7(6,6))*dq6 + 0.5*(+dMM7_7(7,6)+dMM7_6(7,7)-dMM7_7(7,6))*dq7 ;

CC(7,7) = 0.5*(+dMM7_7(7,1)-dMM7_7(1,7))*dq1 + 0.5*(+dMM7_2(7,7)+dMM7_7(7,2)-dMM7_7(2,7))*dq2 + 0.5*(+dMM7_3(7,7)+dMM7_7(7,3)-dMM7_7(3,7))*dq3 + 0.5*(+dMM7_4(7,7)+dMM7_7(7,4)-dMM7_7(4,7))*dq4 + 0.5*(+dMM7_5(7,7)+dMM7_7(7,5)-dMM7_7(5,7))*dq5 + 0.5*(+dMM7_6(7,7)+dMM7_7(7,6)-dMM7_7(6,7))*dq6 + 0.5*(+dMM7_7(7,7)+dMM7_7(7,7)-dMM7_7(7,7))*dq7 ;

M = MM7;

C = CC;

G = GG7;




% ---------- sub functions

	function Ad = f_adjoint(H)
		
		R = H(1:3,1:3);
		d = H(1:3,4);
		Ad = [R f_skew(d)*R; zeros(3,3) R];
		
	end

	function iAd = f_inverseadjoint(H)
		
		R = H(1:3,1:3);
		d = H(1:3,4);
		iAd = [R.' -R.'*f_skew(d); zeros(3,3) R.'];
		
	end

	function ad = f_ad(T)
		
		v = f_skew(T(1:3));
		om = f_skew(T(4:6));
		ad = [om, v; zeros(3,3) om];
		
	end

	function x_tilde = f_skew(x)
		
		x_tilde = [ 0    -x(3)  x(2);
		            x(3)  0    -x(1);
		           -x(2)  x(1)  0];
		       
	end

	function H = f_rbar(H)
		
		H(1:3,4) = [0; 0; 0];
		
	end

	function T = f_ombar(T)
		
		T(1:3) = [0; 0; 0];
		
	end


end