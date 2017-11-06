function [D,dD,dconv] = TUlipCon6_D_dD_dconv(q)


q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5); q6 = q(6); q7 = q(7); q8 = q(8); q9 = q(9); q10 = q(10); q11 = q(11); q12 = q(12); dq1 = q(13); dq2 = q(14); dq3 = q(15); dq4 = q(16); dq5 = q(17); dq6 = q(18); dq7 = q(19); dq8 = q(20); dq9 = q(21); dq10 = q(22); dq11 = q(23); dq12 = q(24);



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [1,0,0,0;
	0,cos(q1),-sin(q1),0;
	0,sin(q1),cos(q1),1/20;
	0,0,0,1];

unittwist1_0 = [0;
	0;
	0;
	1;
	0;
	0];

H1_0 = H0_0*H1_0;

Jl1 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	1,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J1 = (Jl1);

dQ1 = [dq1];

twist1_0 = unittwist1_0*dQ1;

twist1_0 = (twist1_0);

H2_1 = [cos(q2),0,sin(q2),0;
	0,1,0,0;
	-sin(q2),0,cos(q2),0;
	0,0,0,1];

unittwist2_1 = [0;
	0;
	0;
	0;
	1;
	0];

H2_0 = H1_0*H2_1;

Jl2 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,1,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J2 = (f_inverseadjoint(H2_1)*J1+Jl2);

dQ2 = [dq2];

twist2_1 = unittwist2_1*dQ2;

twist2_0 = (f_inverseadjoint(H2_1)*twist1_0+twist2_1);

dJ2 = (-f_ad(twist2_1)*f_inverseadjoint(H2_1)*J1);

H3_2 = [cos(q3),0,sin(q3),0;
	0,1,0,0;
	-sin(q3),0,cos(q3),8/25;
	0,0,0,1];

unittwist3_2 = [0;
	0;
	0;
	0;
	1;
	0];

H3_0 = H2_0*H3_2;

Jl3 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,1,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J3 = (f_inverseadjoint(H3_2)*J2+Jl3);

dQ3 = [dq3];

twist3_2 = unittwist3_2*dQ3;

twist3_0 = (f_inverseadjoint(H3_2)*twist2_0+twist3_2);

dJ3 = (f_inverseadjoint(H3_2)*dJ2-f_ad(twist3_2)*f_inverseadjoint(H3_2)*J2);

H4_3 = [cos(q4),0,sin(q4),0;
	0,1,0,0;
	-sin(q4),0,cos(q4),11/40;
	0,0,0,1];

unittwist4_3 = [0;
	0;
	0;
	0;
	1;
	0];

H4_0 = H3_0*H4_3;

Jl4 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,1,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J4 = (f_inverseadjoint(H4_3)*J3+Jl4);

dQ4 = [dq4];

twist4_3 = unittwist4_3*dQ4;

twist4_0 = (f_inverseadjoint(H4_3)*twist3_0+twist4_3);

dJ4 = (f_inverseadjoint(H4_3)*dJ3-f_ad(twist4_3)*f_inverseadjoint(H4_3)*J3);

H5_4 = [1,0,0,0;
	0,cos(q5),-sin(q5),-231/4000;
	0,sin(q5),cos(q5),1/25;
	0,0,0,1];

unittwist5_4 = [0;
	0;
	0;
	1;
	0;
	0];

H5_0 = H4_0*H5_4;

Jl5 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,1,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J5 = (f_inverseadjoint(H5_4)*J4+Jl5);

dQ5 = [dq5];

twist5_4 = unittwist5_4*dQ5;

twist5_0 = (f_inverseadjoint(H5_4)*twist4_0+twist5_4);

dJ5 = (f_inverseadjoint(H5_4)*dJ4-f_ad(twist5_4)*f_inverseadjoint(H5_4)*J4);

H6_5 = [cos(q6),-sin(q6),0,0;
	sin(q6),cos(q6),0,0;
	0,0,1,141/2000;
	0,0,0,1];

unittwist6_5 = [0;
	0;
	0;
	0;
	0;
	1];

H6_0 = H5_0*H6_5;

Jl6 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,1,0,0,0,0,0,0];

J6 = (f_inverseadjoint(H6_5)*J5+Jl6);

dQ6 = [dq6];

twist6_5 = unittwist6_5*dQ6;

twist6_0 = (f_inverseadjoint(H6_5)*twist5_0+twist6_5);

dJ6 = (f_inverseadjoint(H6_5)*dJ5-f_ad(twist6_5)*f_inverseadjoint(H6_5)*J5);

Hp6 = [1,0,0,0.085;
	0,1,0,0;
	0,0,1,0.389;
	0,0,0,1];

Hp6_0 = H6_0*Hp6;

Hp1_2_3_0_4 = Hp6_0([1  2  3],[4]);

Jp6 = f_adjoint(f_rbar(H6_0))*f_inverseadjoint(Hp6)*J6;

Jp1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jp6([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

dQ = [dq1;
	dq2;
	dq3;
	dq4;
	dq5;
	dq6;
	dq7;
	dq8;
	dq9;
	dq10;
	dq11;
	dq12];

dJp6 = (f_adjoint(f_rbar(H6_0))*f_inverseadjoint(Hp6)*dJ6+f_adjoint(f_rbar(H6_0))*f_ad(f_ombar(twist6_0))*f_inverseadjoint(Hp6)*J6);

dJp1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = dJp6([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

ap6_0 = dJp1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12*dQ;

D = Hp1_2_3_0_4;

dD = Jp1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12;

dconv = ap6_0;




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