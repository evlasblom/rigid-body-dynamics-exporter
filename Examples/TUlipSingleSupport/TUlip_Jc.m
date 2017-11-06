function [Jc] = TUlip_Jc(q)


q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5); q6 = q(6); q7 = q(7); q8 = q(8); q9 = q(9); q10 = q(10); q11 = q(11); q12 = q(12); dq1 = q(13); dq2 = q(14); dq3 = q(15); dq4 = q(16); dq5 = q(17); dq6 = q(18); dq7 = q(19); dq8 = q(20); dq9 = q(21); dq10 = q(22); dq11 = q(23); dq12 = q(24);



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [1,0,0,0;
	0,cos(q1),-sin(q1),0;
	0,sin(q1),cos(q1),1/20;
	0,0,0,1];

Hm1 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
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

mass1 = [0.075];

Jcomt1 = f_adjoint(f_rbar(H1_0))*f_inverseadjoint(Hm1)*J1;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt1([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom1 = (mass1*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass1 = (mass1);

H2_1 = [cos(q2),0,sin(q2),0;
	0,1,0,0;
	-sin(q2),0,cos(q2),0;
	0,0,0,1];

Hm2 = [1,0,0,0.0375;
	0,1,0,-0.00059;
	0,0,1,0.17;
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

mass2 = [1.023];

Jcomt2 = f_adjoint(f_rbar(H2_0))*f_inverseadjoint(Hm2)*J2;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt2([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom2 = (Jcom1+mass2*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass2 = (totalmass1+mass2);

H3_2 = [cos(q3),0,sin(q3),0;
	0,1,0,0;
	-sin(q3),0,cos(q3),8/25;
	0,0,0,1];

Hm3 = [1,0,0,0.005;
	0,1,0,0.0102;
	0,0,1,0.17772;
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

mass3 = [2.14];

Jcomt3 = f_adjoint(f_rbar(H3_0))*f_inverseadjoint(Hm3)*J3;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt3([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom3 = (Jcom2+mass3*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass3 = (totalmass2+mass3);

H4_3 = [cos(q4),0,sin(q4),0;
	0,1,0,0;
	-sin(q4),0,cos(q4),11/40;
	0,0,0,1];

Hm4 = [1,0,0,0;
	0,1,0,-0.05775;
	0,0,1,0.04;
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

mass4 = [0.075];

Jcomt4 = f_adjoint(f_rbar(H4_0))*f_inverseadjoint(Hm4)*J4;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt4([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom4 = (Jcom3+mass4*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass4 = (totalmass3+mass4);

H5_4 = [1,0,0,0;
	0,cos(q5),-sin(q5),-231/4000;
	0,sin(q5),cos(q5),1/25;
	0,0,0,1];

Hm5 = [1,0,0,-0.009;
	0,1,0,0;
	0,0,1,0.0775;
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

mass5 = [0.61];

Jcomt5 = f_adjoint(f_rbar(H5_0))*f_inverseadjoint(Hm5)*J5;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt5([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom5 = (Jcom4+mass5*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass5 = (totalmass4+mass5);

H6_5 = [cos(q6),-sin(q6),0,0;
	sin(q6),cos(q6),0,0;
	0,0,1,141/2000;
	0,0,0,1];

Hm6 = [1,0,0,-0.025;
	0,1,0,-0.07688;
	0,0,1,0.19;
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

mass6 = [8.154];

Jcomt6 = f_adjoint(f_rbar(H6_0))*f_inverseadjoint(Hm6)*J6;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt6([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom6 = (Jcom5+mass6*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass6 = (totalmass5+mass6);

H7_6 = [cos(q7),-sin(q7),0,0;
	sin(q7),cos(q7),0,-3797/25000;
	0,0,1,0;
	0,0,0,1];

Hm7 = [1,0,0,-0.009;
	0,1,0,0;
	0,0,1,0.007;
	0,0,0,1];

unittwist7_6 = [0;
	0;
	0;
	0;
	0;
	1];

H7_0 = H6_0*H7_6;

Jl7 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,1,0,0,0,0,0];

J7 = (f_inverseadjoint(H7_6)*J6+Jl7);

mass7 = [0.61];

Jcomt7 = f_adjoint(f_rbar(H7_0))*f_inverseadjoint(Hm7)*J7;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt7([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom7 = (Jcom6+mass7*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass7 = (totalmass6+mass7);

H8_7 = [1,0,0,0;
	0,cos(q8),-sin(q8),0;
	0,sin(q8),cos(q8),-141/2000;
	0,0,0,1];

Hm8 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

unittwist8_7 = [0;
	0;
	0;
	1;
	0;
	0];

H8_0 = H7_0*H8_7;

Jl8 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,1,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J8 = (f_inverseadjoint(H8_7)*J7+Jl8);

mass8 = [0.075];

Jcomt8 = f_adjoint(f_rbar(H8_0))*f_inverseadjoint(Hm8)*J8;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt8([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom8 = (Jcom7+mass8*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass8 = (totalmass7+mass8);

H9_8 = [cos(q9),0,sin(q9),0;
	0,1,0,-231/4000;
	-sin(q9),0,cos(q9),-1/25;
	0,0,0,1];

Hm9 = [1,0,0,0.005;
	0,1,0,-0.0102;
	0,0,1,-0.09728;
	0,0,0,1];

unittwist9_8 = [0;
	0;
	0;
	0;
	1;
	0];

H9_0 = H8_0*H9_8;

Jl9 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,1,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J9 = (f_inverseadjoint(H9_8)*J8+Jl9);

mass9 = [2.14];

Jcomt9 = f_adjoint(f_rbar(H9_0))*f_inverseadjoint(Hm9)*J9;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt9([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom9 = (Jcom8+mass9*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass9 = (totalmass8+mass9);

H10_9 = [cos(q10),0,sin(q10),0;
	0,1,0,0;
	-sin(q10),0,cos(q10),-11/40;
	0,0,0,1];

Hm10 = [1,0,0,0.0375;
	0,1,0,0.00059;
	0,0,1,-0.15;
	0,0,0,1];

unittwist10_9 = [0;
	0;
	0;
	0;
	1;
	0];

H10_0 = H9_0*H10_9;

Jl10 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,1,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J10 = (f_inverseadjoint(H10_9)*J9+Jl10);

mass10 = [1.023];

Jcomt10 = f_adjoint(f_rbar(H10_0))*f_inverseadjoint(Hm10)*J10;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt10([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom10 = (Jcom9+mass10*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass10 = (totalmass9+mass10);

H11_10 = [cos(q11),0,sin(q11),0;
	0,1,0,0;
	-sin(q11),0,cos(q11),-8/25;
	0,0,0,1];

Hm11 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

unittwist11_10 = [0;
	0;
	0;
	0;
	1;
	0];

H11_0 = H10_0*H11_10;

Jl11 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,1,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J11 = (f_inverseadjoint(H11_10)*J10+Jl11);

mass11 = [0.075];

Jcomt11 = f_adjoint(f_rbar(H11_0))*f_inverseadjoint(Hm11)*J11;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt11([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom11 = (Jcom10+mass11*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass11 = (totalmass10+mass11);

H12_11 = [1,0,0,0;
	0,cos(q12),-sin(q12),0;
	0,sin(q12),cos(q12),0;
	0,0,0,1];

Hm12 = [1,0,0,0.035;
	0,1,0,-0.01;
	0,0,1,-0.01;
	0,0,0,1];

unittwist12_11 = [0;
	0;
	0;
	1;
	0;
	0];

H12_0 = H11_0*H12_11;

Jl12 = [0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,1;
	0,0,0,0,0,0,0,0,0,0,0,0;
	0,0,0,0,0,0,0,0,0,0,0,0];

J12 = (f_inverseadjoint(H12_11)*J11+Jl12);

mass12 = [0.37];

Jcomt12 = f_adjoint(f_rbar(H12_0))*f_inverseadjoint(Hm12)*J12;

Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12 = Jcomt12([1  2  3],[1   2   3   4   5   6   7   8   9  10  11  12]);

Jcom12 = (Jcom11+mass12*Jcomi1_2_3_0_1_2_3_4_5_6_7_8_9_10_11_12);

totalmass12 = (totalmass11+mass12);

Jcom_final12 = inv(totalmass12)*Jcom12;

Jc = Jcom_final12;




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