function [M,C,G] = SimpleArm_M_C_G(q,p)


q1 = q(1); q2 = q(2); dq1 = q(3); dq2 = q(4);

d1_z = p.d1_z; r1_z = p.r1_z; m1 = p.m1; i1_x = p.i1_x; i1_y = p.i1_y; i1_z = p.i1_z; d2_z = p.d2_z; r2_z = p.r2_z; m2 = p.m2; i2_x = p.i2_x; i2_y = p.i2_y; i2_z = p.i2_z;



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [cos(q1),0,sin(q1),0;
	0,1,0,0;
	-sin(q1),0,cos(q1),d1_z;
	0,0,0,1];

Hm1 = [1,0,0,0;
	0,1,0,0;
	0,0,1,r1_z;
	0,0,0,1];

H1_0 = H0_0*H1_0;

Jl1 = [0,0;
	0,0;
	0,0;
	0,0;
	1,0;
	0,0];

J1 = (Jl1);

Ja1 = f_adjoint(f_rbar(H1_0))*f_inverseadjoint(Hm1)*J1;

Im1 = [m1,0,0,0,0,0;
	0,m1,0,0,0,0;
	0,0,m1,0,0,0;
	0,0,0,i1_x,0,0;
	0,0,0,0,i1_y,0;
	0,0,0,0,0,i1_z];

I1 = (f_inverseadjoint(Hm1)).'*Im1*f_inverseadjoint(Hm1);

MM1 = ((J1).'*I1*J1);

dAd1 = [-sin(q1),0,-cos(q1),0,-d1_z*sin(q1),0;
	0,0,0,0,0,0;
	cos(q1),0,-sin(q1),0,d1_z*cos(q1),0;
	0,0,0,-sin(q1),0,-cos(q1);
	0,0,0,0,0,0;
	0,0,0,cos(q1),0,-sin(q1)];

Fz1 = [0;
	0;
	(981*m1)/100;
	0;
	0;
	0];

GG1 = ((Ja1).'*Fz1);

H2_1 = [cos(q2),0,sin(q2),0;
	0,1,0,0;
	-sin(q2),0,cos(q2),d2_z;
	0,0,0,1];

Hm2 = [1,0,0,0;
	0,1,0,0;
	0,0,1,r2_z;
	0,0,0,1];

H2_0 = H1_0*H2_1;

Jl2 = [0,0;
	0,0;
	0,0;
	0,0;
	0,1;
	0,0];

J2 = (f_inverseadjoint(H2_1)*J1+Jl2);

Ja2 = f_adjoint(f_rbar(H2_0))*f_inverseadjoint(Hm2)*J2;

Im2 = [m2,0,0,0,0,0;
	0,m2,0,0,0,0;
	0,0,m2,0,0,0;
	0,0,0,i2_x,0,0;
	0,0,0,0,i2_y,0;
	0,0,0,0,0,i2_z];

I2 = (f_inverseadjoint(Hm2)).'*Im2*f_inverseadjoint(Hm2);

MM2 = (MM1+(J2).'*I2*J2);

dAd2 = [-sin(q2),0,-cos(q2),0,-d2_z*sin(q2),0;
	0,0,0,0,0,0;
	cos(q2),0,-sin(q2),0,d2_z*cos(q2),0;
	0,0,0,-sin(q2),0,-cos(q2);
	0,0,0,0,0,0;
	0,0,0,cos(q2),0,-sin(q2)];

dJ2_2 = (dAd2*J1);

dMM2_2 = ((dJ2_2).'*I2*J2+(J2).'*I2*dJ2_2);

Fz2 = [0;
	0;
	(981*m2)/100;
	0;
	0;
	0];

GG2 = (GG1+(Ja2).'*Fz2);

CC(1,1) = 0 + 0.5*(+dMM2_2(1,1))*dq2 ;

CC(1,2) = 0.5*(+dMM2_2(1,1))*dq1 + 0.5*(+dMM2_2(1,2)+dMM2_2(1,2))*dq2 ;

CC(2,1) = 0.5*(-dMM2_2(1,1))*dq1 + 0.5*(+dMM2_2(2,1)-dMM2_2(2,1))*dq2 ;

CC(2,2) = 0.5*(+dMM2_2(2,1)-dMM2_2(1,2))*dq1 + 0.5*(+dMM2_2(2,2)+dMM2_2(2,2)-dMM2_2(2,2))*dq2 ;

M = MM2;

C = CC;

G = GG2;




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