function [M,K] = TriplePendulumLin_M_K(q,p)


q1 = q(1); q2 = q(2); q3 = q(3);

d1_z = p.d1_z; r1_x = p.r1_x; r1_z = p.r1_z; m1 = p.m1; i1_x = p.i1_x; i1_y = p.i1_y; i1_z = p.i1_z; d2_z = p.d2_z; r2_x = p.r2_x; r2_z = p.r2_z; m2 = p.m2; i2_x = p.i2_x; i2_y = p.i2_y; i2_z = p.i2_z; d3_z = p.d3_z; r3_x = p.r3_x; r3_z = p.r3_z; m3 = p.m3; i3_x = p.i3_x; i3_y = p.i3_y; i3_z = p.i3_z;



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

dGcol0_0 = [0;
	0;
	0];

H1_0 = [cos(q1),0,sin(q1),0;
	0,1,0,0;
	-sin(q1),0,cos(q1),d1_z;
	0,0,0,1];

unittwist1_0 = [0;
	0;
	0;
	0;
	1;
	0];

H1_0 = H0_0*H1_0;

Jl1 = [0,0,0;
	0,0,0;
	0,0,0;
	0,0,0;
	1,0,0;
	0,0,0];

J1 = (Jl1);

AdHm1 = [1,0,0,0,r1_z,0;
	0,1,0,-r1_z,0,r1_x;
	0,0,1,0,-r1_x,0;
	0,0,0,1,0,0;
	0,0,0,0,1,0;
	0,0,0,0,0,1];

I1 = [m1,0,0,0,m1*r1_z,0;
	0,m1,0,-m1*r1_z,0,m1*r1_x;
	0,0,m1,0,-m1*r1_x,0;
	0,-m1*r1_z,0,i1_x + m1*r1_z^2,0,-m1*r1_x*r1_z;
	m1*r1_z,0,-m1*r1_x,0,i1_y + m1*r1_x^2 + m1*r1_z^2,0;
	0,m1*r1_x,0,-m1*r1_x*r1_z,0,i1_z + m1*r1_x^2];

MM1 = ((J1).'*I1*J1);

dTddq1_1 = (unittwist1_0);

dAd1 = [-sin(q1),0,-cos(q1),0,-d1_z*sin(q1),0;
	0,0,0,0,0,0;
	cos(q1),0,-sin(q1),0,d1_z*cos(q1),0;
	0,0,0,-sin(q1),0,-cos(q1);
	0,0,0,0,0,0;
	0,0,0,cos(q1),0,-sin(q1)];

dJa1_1 = (f_adjoint(f_rbar(H1_0))*f_ad(f_ombar(dTddq1_1))*AdHm1*J1);

Fz1 = [0;
	0;
	(981*m1)/100;
	0;
	0;
	0];

dGcol1_1 = (dJa1_1).'*Fz1;

dG1 = [dGcol1_1,dGcol0_0,dGcol0_0];

KK1 = (dG1);

H2_1 = [cos(q2),0,sin(q2),0;
	0,1,0,0;
	-sin(q2),0,cos(q2),d2_z;
	0,0,0,1];

unittwist2_1 = [0;
	0;
	0;
	0;
	1;
	0];

H2_0 = H1_0*H2_1;

Jl2 = [0,0,0;
	0,0,0;
	0,0,0;
	0,0,0;
	0,1,0;
	0,0,0];

J2 = (f_inverseadjoint(H2_1)*J1+Jl2);

AdHm2 = [1,0,0,0,r2_z,0;
	0,1,0,-r2_z,0,r2_x;
	0,0,1,0,-r2_x,0;
	0,0,0,1,0,0;
	0,0,0,0,1,0;
	0,0,0,0,0,1];

I2 = [m2,0,0,0,m2*r2_z,0;
	0,m2,0,-m2*r2_z,0,m2*r2_x;
	0,0,m2,0,-m2*r2_x,0;
	0,-m2*r2_z,0,i2_x + m2*r2_z^2,0,-m2*r2_x*r2_z;
	m2*r2_z,0,-m2*r2_x,0,i2_y + m2*r2_x^2 + m2*r2_z^2,0;
	0,m2*r2_x,0,-m2*r2_x*r2_z,0,i2_z + m2*r2_x^2];

MM2 = (MM1+(J2).'*I2*J2);

dTddq2_1 = f_inverseadjoint(H2_1)*dTddq1_1;

dJa2_1 = (f_adjoint(f_rbar(H2_0))*f_ad(f_ombar(dTddq2_1))*AdHm2*J2);

dTddq2_2 = (unittwist2_1);

dAd2 = [-sin(q2),0,-cos(q2),0,-d2_z*sin(q2),0;
	0,0,0,0,0,0;
	cos(q2),0,-sin(q2),0,d2_z*cos(q2),0;
	0,0,0,-sin(q2),0,-cos(q2);
	0,0,0,0,0,0;
	0,0,0,cos(q2),0,-sin(q2)];

dJ2_2 = (dAd2*J1);

dJa2_2 = (f_adjoint(f_rbar(H2_0))*AdHm2*dJ2_2+f_adjoint(f_rbar(H2_0))*f_ad(f_ombar(dTddq2_2))*AdHm2*J2);

Fz2 = [0;
	0;
	(981*m2)/100;
	0;
	0;
	0];

dGcol2_1 = (dJa2_1).'*Fz2;

dGcol2_2 = (dJa2_2).'*Fz2;

dG2 = [dGcol2_1,dGcol2_2,dGcol0_0];

KK2 = (KK1+dG2);

H3_2 = [cos(q3),0,sin(q3),0;
	0,1,0,0;
	-sin(q3),0,cos(q3),d3_z;
	0,0,0,1];

unittwist3_2 = [0;
	0;
	0;
	0;
	1;
	0];

H3_0 = H2_0*H3_2;

Jl3 = [0,0,0;
	0,0,0;
	0,0,0;
	0,0,0;
	0,0,1;
	0,0,0];

J3 = (f_inverseadjoint(H3_2)*J2+Jl3);

AdHm3 = [1,0,0,0,r3_z,0;
	0,1,0,-r3_z,0,r3_x;
	0,0,1,0,-r3_x,0;
	0,0,0,1,0,0;
	0,0,0,0,1,0;
	0,0,0,0,0,1];

I3 = [m3,0,0,0,m3*r3_z,0;
	0,m3,0,-m3*r3_z,0,m3*r3_x;
	0,0,m3,0,-m3*r3_x,0;
	0,-m3*r3_z,0,i3_x + m3*r3_z^2,0,-m3*r3_x*r3_z;
	m3*r3_z,0,-m3*r3_x,0,i3_y + m3*r3_x^2 + m3*r3_z^2,0;
	0,m3*r3_x,0,-m3*r3_x*r3_z,0,i3_z + m3*r3_x^2];

MM3 = (MM2+(J3).'*I3*J3);

dTddq3_1 = f_inverseadjoint(H3_2)*dTddq2_1;

dJa3_1 = (f_adjoint(f_rbar(H3_0))*f_ad(f_ombar(dTddq3_1))*AdHm3*J3);

dTddq3_2 = f_inverseadjoint(H3_2)*dTddq2_2;

dJ3_2 = f_inverseadjoint(H3_2)*dJ2_2;

dJa3_2 = (f_adjoint(f_rbar(H3_0))*AdHm3*dJ3_2+f_adjoint(f_rbar(H3_0))*f_ad(f_ombar(dTddq3_2))*AdHm3*J3);

dTddq3_3 = (unittwist3_2);

dAd3 = [-sin(q3),0,-cos(q3),0,-d3_z*sin(q3),0;
	0,0,0,0,0,0;
	cos(q3),0,-sin(q3),0,d3_z*cos(q3),0;
	0,0,0,-sin(q3),0,-cos(q3);
	0,0,0,0,0,0;
	0,0,0,cos(q3),0,-sin(q3)];

dJ3_3 = (dAd3*J2);

dJa3_3 = (f_adjoint(f_rbar(H3_0))*AdHm3*dJ3_3+f_adjoint(f_rbar(H3_0))*f_ad(f_ombar(dTddq3_3))*AdHm3*J3);

Fz3 = [0;
	0;
	(981*m3)/100;
	0;
	0;
	0];

dGcol3_1 = (dJa3_1).'*Fz3;

dGcol3_2 = (dJa3_2).'*Fz3;

dGcol3_3 = (dJa3_3).'*Fz3;

dG3 = [dGcol3_1,dGcol3_2,dGcol3_3];

KK3 = (KK2+dG3);

M = MM3;

K = KK3;




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