function [D,dD,dconv] = SCARACon4_D_dD_dconv(q)


q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); dq1 = q(5); dq2 = q(6); dq3 = q(7); dq4 = q(8);



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [cos(q1),-sin(q1),0,0;
	sin(q1),cos(q1),0,0;
	0,0,1,1/5;
	0,0,0,1];

unittwist1_0 = [0;
	0;
	0;
	0;
	0;
	1];

H1_0 = H0_0*H1_0;

Jl1 = [0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	1,0,0,0];

J1 = (Jl1);

dQ1 = [dq1];

twist1_0 = unittwist1_0*dQ1;

twist1_0 = (twist1_0);

H2_1 = [cos(q2),-sin(q2),0,0;
	sin(q2),cos(q2),0,1/10;
	0,0,1,0;
	0,0,0,1];

unittwist2_1 = [0;
	0;
	0;
	0;
	0;
	1];

H2_0 = H1_0*H2_1;

Jl2 = [0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,1,0,0];

J2 = (f_inverseadjoint(H2_1)*J1+Jl2);

dQ2 = [dq2];

twist2_1 = unittwist2_1*dQ2;

twist2_0 = (f_inverseadjoint(H2_1)*twist1_0+twist2_1);

dJ2 = (-f_ad(twist2_1)*f_inverseadjoint(H2_1)*J1);

H3_2 = [cos(q3),-sin(q3),0,0;
	sin(q3),cos(q3),0,1/10;
	0,0,1,0;
	0,0,0,1];

unittwist3_2 = [0;
	0;
	0;
	0;
	0;
	1];

H3_0 = H2_0*H3_2;

Jl3 = [0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0;
	0,0,1,0];

J3 = (f_inverseadjoint(H3_2)*J2+Jl3);

dQ3 = [dq3];

twist3_2 = unittwist3_2*dQ3;

twist3_0 = (f_inverseadjoint(H3_2)*twist2_0+twist3_2);

dJ3 = (f_inverseadjoint(H3_2)*dJ2-f_ad(twist3_2)*f_inverseadjoint(H3_2)*J2);

H4_3 = [1,0,0,0;
	0,1,0,0;
	0,0,1,q4;
	0,0,0,1];

unittwist4_3 = [0;
	0;
	1;
	0;
	0;
	0];

H4_0 = H3_0*H4_3;

Jl4 = [0,0,0,0;
	0,0,0,0;
	0,0,0,1;
	0,0,0,0;
	0,0,0,0;
	0,0,0,0];

J4 = (f_inverseadjoint(H4_3)*J3+Jl4);

dQ4 = [dq4];

twist4_3 = unittwist4_3*dQ4;

twist4_0 = (f_inverseadjoint(H4_3)*twist3_0+twist4_3);

dJ4 = (f_inverseadjoint(H4_3)*dJ3-f_ad(twist4_3)*f_inverseadjoint(H4_3)*J3);

Hp4 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

Hp4_0 = H4_0*Hp4;

Hp1_2_3_0_4 = Hp4_0([1  2  3],[4]);

Jp4 = f_adjoint(f_rbar(H4_0))*f_inverseadjoint(Hp4)*J4;

Jp1_2_3_0_1_2_3_4 = Jp4([1  2  3],[1  2  3  4]);

dQ = [dq1;
	dq2;
	dq3;
	dq4];

dJp4 = (f_adjoint(f_rbar(H4_0))*f_inverseadjoint(Hp4)*dJ4+f_adjoint(f_rbar(H4_0))*f_ad(f_ombar(twist4_0))*f_inverseadjoint(Hp4)*J4);

dJp1_2_3_0_1_2_3_4 = dJp4([1  2  3],[1  2  3  4]);

ap4_0 = dJp1_2_3_0_1_2_3_4*dQ;

D = Hp1_2_3_0_4;

dD = Jp1_2_3_0_1_2_3_4;

dconv = ap4_0;




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