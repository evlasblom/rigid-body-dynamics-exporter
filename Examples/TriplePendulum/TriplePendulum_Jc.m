function [Jc] = TriplePendulum_Jc(q,p)


q1 = q(1); q2 = q(2); q3 = q(3); dq1 = q(4); dq2 = q(5); dq3 = q(6);

d1_z = p.d1_z; r1_x = p.r1_x; r1_z = p.r1_z; m1 = p.m1; i1_x = p.i1_x; i1_y = p.i1_y; i1_z = p.i1_z; d2_z = p.d2_z; r2_x = p.r2_x; r2_z = p.r2_z; m2 = p.m2; i2_x = p.i2_x; i2_y = p.i2_y; i2_z = p.i2_z; d3_z = p.d3_z; r3_x = p.r3_x; r3_z = p.r3_z; m3 = p.m3; i3_x = p.i3_x; i3_y = p.i3_y; i3_z = p.i3_z;



H0_0 = [1,0,0,0;
	0,1,0,0;
	0,0,1,0;
	0,0,0,1];

H1_0 = [cos(q1),0,sin(q1),0;
	0,1,0,0;
	-sin(q1),0,cos(q1),d1_z;
	0,0,0,1];

Hm1 = [1,0,0,r1_x;
	0,1,0,0;
	0,0,1,r1_z;
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

mass1 = [m1];

Jcomt1 = f_adjoint(f_rbar(H1_0))*f_inverseadjoint(Hm1)*J1;

Jcomi1_2_3_0_1_2_3 = Jcomt1([1  2  3],[1  2  3]);

Jcom1 = (mass1*Jcomi1_2_3_0_1_2_3);

totalmass1 = (mass1);

H2_1 = [cos(q2),0,sin(q2),0;
	0,1,0,0;
	-sin(q2),0,cos(q2),d2_z;
	0,0,0,1];

Hm2 = [1,0,0,r2_x;
	0,1,0,0;
	0,0,1,r2_z;
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

mass2 = [m2];

Jcomt2 = f_adjoint(f_rbar(H2_0))*f_inverseadjoint(Hm2)*J2;

Jcomi1_2_3_0_1_2_3 = Jcomt2([1  2  3],[1  2  3]);

Jcom2 = (Jcom1+mass2*Jcomi1_2_3_0_1_2_3);

totalmass2 = (totalmass1+mass2);

H3_2 = [cos(q3),0,sin(q3),0;
	0,1,0,0;
	-sin(q3),0,cos(q3),d3_z;
	0,0,0,1];

Hm3 = [1,0,0,r3_x;
	0,1,0,0;
	0,0,1,r3_z;
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

mass3 = [m3];

Jcomt3 = f_adjoint(f_rbar(H3_0))*f_inverseadjoint(Hm3)*J3;

Jcomi1_2_3_0_1_2_3 = Jcomt3([1  2  3],[1  2  3]);

Jcom3 = (Jcom2+mass3*Jcomi1_2_3_0_1_2_3);

totalmass3 = (totalmass2+mass3);

Jcom_final3 = inv(totalmass3)*Jcom3;

Jc = Jcom_final3;




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