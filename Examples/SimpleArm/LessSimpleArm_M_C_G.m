function [M,C,G] = LessSimpleArm_M_C_G(q)


q1 = q(1); q2 = q(2); q3 = q(3); dq1 = q(4); dq2 = q(5); dq3 = q(6);



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [cos(q1),0,sin(q1),0;
	0,1,0,0;
	-sin(q1),0,cos(q1),1/10;
	0,0,0,1];

Hm1 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0.1;
	0,0,0,1];

H1_0 = H0_0*H1_0;

Jl1 = [0,0,0;
	0,0,0;
	0,0,0;
	0,0,0;
	1,0,0;
	0,0,0];

J1 = (Jl1);

Ja1 = f_adjoint(f_rbar(H1_0))*f_inverseadjoint(Hm1)*J1;

Im1 = [0.5,0,0,0,0,0;
	0,0.5,0,0,0,0;
	0,0,0.5,0,0,0;
	0,0,0,0.02,0,0;
	0,0,0,0,0.02,0;
	0,0,0,0,0,0.003];

I1 = (f_inverseadjoint(Hm1)).'*Im1*f_inverseadjoint(Hm1);

MM1 = ((J1).'*I1*J1);

dAd1 = [-sin(q1),0,-cos(q1),0,-sin(q1)/10,0;
	0,0,0,0,0,0;
	cos(q1),0,-sin(q1),0,cos(q1)/10,0;
	0,0,0,-sin(q1),0,-cos(q1);
	0,0,0,0,0,0;
	0,0,0,cos(q1),0,-sin(q1)];

Fz1 = [0;
	0;
	4.905;
	0;
	0;
	0];

GG1 = ((Ja1).'*Fz1);

H2_1 = [cos(q2),0,sin(q2),0;
	0,1,0,0;
	-sin(q2),0,cos(q2),1/5;
	0,0,0,1];

Hm2 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0.1;
	0,0,0,1];

H2_0 = H1_0*H2_1;

Jl2 = [0,0,0;
	0,0,0;
	0,0,0;
	0,0,0;
	0,1,0;
	0,0,0];

J2 = (f_inverseadjoint(H2_1)*J1+Jl2);

Ja2 = f_adjoint(f_rbar(H2_0))*f_inverseadjoint(Hm2)*J2;

Im2 = [0.5,0,0,0,0,0;
	0,0.5,0,0,0,0;
	0,0,0.5,0,0,0;
	0,0,0,0.02,0,0;
	0,0,0,0,0.02,0;
	0,0,0,0,0,0.003];

I2 = (f_inverseadjoint(Hm2)).'*Im2*f_inverseadjoint(Hm2);

MM2 = (MM1+(J2).'*I2*J2);

dAd2 = [-sin(q2),0,-cos(q2),0,-sin(q2)/5,0;
	0,0,0,0,0,0;
	cos(q2),0,-sin(q2),0,cos(q2)/5,0;
	0,0,0,-sin(q2),0,-cos(q2);
	0,0,0,0,0,0;
	0,0,0,cos(q2),0,-sin(q2)];

dJ2_2 = (dAd2*J1);

dMM2_2 = ((dJ2_2).'*I2*J2+(J2).'*I2*dJ2_2);

Fz2 = [0;
	0;
	4.905;
	0;
	0;
	0];

GG2 = (GG1+(Ja2).'*Fz2);

H3_2 = [1,0,0,0;
	0,cos(q3),-sin(q3),0;
	0,sin(q3),cos(q3),1/5;
	0,0,0,1];

Hm3 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0.1;
	0,0,0,1];

H3_0 = H2_0*H3_2;

Jl3 = [0,0,0;
	0,0,0;
	0,0,0;
	0,0,1;
	0,0,0;
	0,0,0];

J3 = (f_inverseadjoint(H3_2)*J2+Jl3);

Ja3 = f_adjoint(f_rbar(H3_0))*f_inverseadjoint(Hm3)*J3;

Im3 = [0.5,0,0,0,0,0;
	0,0.5,0,0,0,0;
	0,0,0.5,0,0,0;
	0,0,0,0.02,0,0;
	0,0,0,0,0.02,0;
	0,0,0,0,0,0.003];

I3 = (f_inverseadjoint(Hm3)).'*Im3*f_inverseadjoint(Hm3);

MM3 = (MM2+(J3).'*I3*J3);

dJ3_2 = f_inverseadjoint(H3_2)*dJ2_2;

dMM3_2 = (dMM2_2+(dJ3_2).'*I3*J3+(J3).'*I3*dJ3_2);

dAd3 = [0,0,0,0,0,0;
	0,-sin(q3),cos(q3),sin(q3)/5,0,0;
	0,-cos(q3),-sin(q3),cos(q3)/5,0,0;
	0,0,0,0,0,0;
	0,0,0,0,-sin(q3),cos(q3);
	0,0,0,0,-cos(q3),-sin(q3)];

dJ3_3 = (dAd3*J2);

dMM3_3 = ((dJ3_3).'*I3*J3+(J3).'*I3*dJ3_3);

Fz3 = [0;
	0;
	4.905;
	0;
	0;
	0];

GG3 = (GG2+(Ja3).'*Fz3);

CC(1,1) = 0 + 0.5*(+dMM3_2(1,1))*dq2 + 0.5*(+dMM3_3(1,1))*dq3 ;

CC(1,2) = 0.5*(+dMM3_2(1,1))*dq1 + 0.5*(+dMM3_2(1,2)+dMM3_2(1,2))*dq2 + 0.5*(+dMM3_3(1,2)+dMM3_2(1,3))*dq3 ;

CC(1,3) = 0.5*(+dMM3_3(1,1))*dq1 + 0.5*(+dMM3_2(1,3)+dMM3_3(1,2))*dq2 + 0.5*(+dMM3_3(1,3)+dMM3_3(1,3))*dq3 ;

CC(2,1) = 0.5*(-dMM3_2(1,1))*dq1 + 0.5*(+dMM3_2(2,1)-dMM3_2(2,1))*dq2 + 0.5*(+dMM3_3(2,1)-dMM3_2(3,1))*dq3 ;

CC(2,2) = 0.5*(+dMM3_2(2,1)-dMM3_2(1,2))*dq1 + 0.5*(+dMM3_2(2,2)+dMM3_2(2,2)-dMM3_2(2,2))*dq2 + 0.5*(+dMM3_3(2,2)+dMM3_2(2,3)-dMM3_2(3,2))*dq3 ;

CC(2,3) = 0.5*(+dMM3_3(2,1)-dMM3_2(1,3))*dq1 + 0.5*(+dMM3_2(2,3)+dMM3_3(2,2)-dMM3_2(2,3))*dq2 + 0.5*(+dMM3_3(2,3)+dMM3_3(2,3)-dMM3_2(3,3))*dq3 ;

CC(3,1) = 0.5*(-dMM3_3(1,1))*dq1 + 0.5*(+dMM3_2(3,1)-dMM3_3(2,1))*dq2 + 0.5*(+dMM3_3(3,1)-dMM3_3(3,1))*dq3 ;

CC(3,2) = 0.5*(+dMM3_2(3,1)-dMM3_3(1,2))*dq1 + 0.5*(+dMM3_2(3,2)+dMM3_2(3,2)-dMM3_3(2,2))*dq2 + 0.5*(+dMM3_3(3,2)+dMM3_2(3,3)-dMM3_3(3,2))*dq3 ;

CC(3,3) = 0.5*(+dMM3_3(3,1)-dMM3_3(1,3))*dq1 + 0.5*(+dMM3_2(3,3)+dMM3_3(3,2)-dMM3_3(2,3))*dq2 + 0.5*(+dMM3_3(3,3)+dMM3_3(3,3)-dMM3_3(3,3))*dq3 ;

M = MM3;

C = CC;

G = GG3;




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