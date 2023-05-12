
// Written by C. Nguyen and P. Fernandez

void flux_AV2_ns3d(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param, double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{

// Artificial viscosity fluxes. Latest version of the model.
		
	double rampFactor = app.rampFactor;
	double porder = (double) app.porder;
	double alpha = 100.0;
	double bbeta = 1.0e-2;
	double k_h = 1.5;
	double eps_v = 1.0e-8;
	double epsilon0 = 0.01;

	double gam = param[0];
	double Minf = param[4];

	for (int i = 0; i <ng; i++) {
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];
		double u10 = udg[9*ng+i];
		double u11 = udg[10*ng+i];
		double u12 = udg[11*ng+i];
		double u13 = udg[12*ng+i];
		double u14 = udg[13*ng+i];
		double u15 = udg[14*ng+i];
		double u16 = udg[15*ng+i];
		double u17 = udg[16*ng+i];
		double u18 = udg[17*ng+i];
		double u19 = udg[18*ng+i];
		double u20 = udg[19*ng+i];

		double h = pg[3*ng+i];

		double r = u1;
		double ru = u2;
		double rv = u3;
		double rw = u4;
		double rE = u5;
		double rx = u6;
		double rux = u7;
		double rvx = u8;
		double rwx = u9;
		double rEx = u10;
		double ry = u11;
		double ruy = u12;
		double rvy = u13;
		double rwy = u14;
		double rEy = u15;
		double rz = u16;
		double ruz = u17;
		double rvz = u18;
		double rwz = u19;
		double rEz = u20;

		double u = ru/r;
		double v = rv/r;
		double w = rw/r;
		double ux = (rux - u*rx)/r;
		double vy = (rvy - v*ry)/r;
		double wz = (rwz - w*rz)/r;

		double t2 = 1.0/3.141592653589793;
		double t3 = 1.0/r;
		double t12 = ru*rx*t3;
		double t4 = rux-t12;
		double t5 = t3*t4;
		double t13 = rv*ry*t3;
		double t6 = rvy-t13;
		double t7 = t3*t6;
		double t14 = rw*rz*t3;
		double t8 = rwz-t14;
		double t9 = t3*t8;
		double t10 = t5+t7+t9;
		double t11 = atan(alpha);
		double t15 = alpha*h*t10;
		double t16 = atan(t15);
		double t17 = t2*t16;
		double t18 = t17+1.0/2.0;
		double t19 = h*t10*t18;
		double t21 = t2*t11;
		double t20 = t19-t21+1.0/2.0;
		double t22 = 1.0/(r*r);
		double t23 = gam*(1.0/2.0);
		double t24 = t23-1.0/2.0;
		double t25 = ru*ru;
		double t26 = rv*rv;
		double t27 = rw*rw;
		f[0*ng+i] += epsilon0*rampFactor*rx*t20;
		f[1*ng+i] += epsilon0*rampFactor*rux*t20;
		f[2*ng+i] += epsilon0*rampFactor*rvx*t20;
		f[3*ng+i] += epsilon0*rampFactor*rwx*t20;
		f[4*ng+i] += epsilon0*rampFactor*t20*(gam*rEx-t24*(ru*t3*t4*2.0+rx*t22*t25+rx*t22*t26+rx*t22*t27+rv*t3*(rvx-rv*rx*t3)*2.0+rw*t3*(rwx-rw*rx*t3)*2.0));
		f[5*ng+i] += epsilon0*rampFactor*ry*t20;
		f[6*ng+i] += epsilon0*rampFactor*ruy*t20;
		f[7*ng+i] += epsilon0*rampFactor*rvy*t20;
		f[8*ng+i] += epsilon0*rampFactor*rwy*t20;
		f[9*ng+i] += epsilon0*rampFactor*t20*(gam*rEy-t24*(rv*t3*t6*2.0+ry*t22*t25+ry*t22*t26+ry*t22*t27+ru*t3*(ruy-ru*ry*t3)*2.0+rw*t3*(rwy-rw*ry*t3)*2.0));
		f[10*ng+i] += epsilon0*rampFactor*rz*t20;
		f[11*ng+i] += epsilon0*rampFactor*ruz*t20;
		f[12*ng+i] += epsilon0*rampFactor*rvz*t20;
		f[13*ng+i] += epsilon0*rampFactor*rwz*t20;
		f[14*ng+i] += epsilon0*rampFactor*t20*(gam*rEz-t24*(rw*t3*t8*2.0+rz*t22*t25+rz*t22*t26+rz*t22*t27+ru*t3*(ruz-ru*rz*t3)*2.0+rv*t3*(rvz-rv*rz*t3)*2.0));

	}

	if (computeJacobian == 1) {

		for (int i = 0; i <ng; i++) {
			double u1 = udg[0*ng+i];
			double u2 = udg[1*ng+i];
			double u3 = udg[2*ng+i];
			double u4 = udg[3*ng+i];
			double u5 = udg[4*ng+i];
			double u6 = udg[5*ng+i];
			double u7 = udg[6*ng+i];
			double u8 = udg[7*ng+i];
			double u9 = udg[8*ng+i];
			double u10 = udg[9*ng+i];
			double u11 = udg[10*ng+i];
			double u12 = udg[11*ng+i];
			double u13 = udg[12*ng+i];
			double u14 = udg[13*ng+i];
			double u15 = udg[14*ng+i];
			double u16 = udg[15*ng+i];
			double u17 = udg[16*ng+i];
			double u18 = udg[17*ng+i];
			double u19 = udg[18*ng+i];
			double u20 = udg[19*ng+i];

			double h = pg[3*ng+i];

			double r = u1;
			double ru = u2;
			double rv = u3;
			double rw = u4;
			double rE = u5;
			double rx = u6;
			double rux = u7;
			double rvx = u8;
			double rwx = u9;
			double rEx = u10;
			double ry = u11;
			double ruy = u12;
			double rvy = u13;
			double rwy = u14;
			double rEy = u15;
			double rz = u16;
			double ruz = u17;
			double rvz = u18;
			double rwz = u19;
			double rEz = u20;

			double u = ru/r;
			double v = rv/r;
			double w = rw/r;
			double ux = (rux - u*rx)/r;
			double vy = (rvy - v*ry)/r;
			double wz = (rwz - w*rz)/r;

			double t2 = 1.0/r;
			double t10 = ru*rx*t2;
			double t3 = rux-t10;
			double t4 = 1.0/(r*r);
			double t12 = rv*ry*t2;
			double t5 = rvy-t12;
			double t14 = rw*rz*t2;
			double t6 = rwz-t14;
			double t7 = 1.0/(r*r*r);
			double t8 = 1.0/3.141592653589793;
			double t9 = h*h;
			double t11 = t2*t3;
			double t13 = t2*t5;
			double t15 = t2*t6;
			double t16 = t11+t13+t15;
			double t17 = t3*t4;
			double t18 = t4*t5;
			double t19 = t4*t6;
			double t25 = ru*rx*t7;
			double t26 = rv*ry*t7;
			double t27 = rw*rz*t7;
			double t20 = t17+t18+t19-t25-t26-t27;
			double t21 = alpha*h*t16;
			double t22 = atan(t21);
			double t23 = t8*t22;
			double t24 = t23+1.0/2.0;
			double t28 = h*t20*t24;
			double t29 = alpha*alpha;
			double t30 = t16*t16;
			double t31 = t9*t29*t30;
			double t32 = t31+1.0;
			double t33 = 1.0/t32;
			double t34 = alpha*t8*t9*t16*t20*t33;
			double t35 = t28+t34;
			double t36 = gam*(1.0/2.0);
			double t37 = t36-1.0/2.0;
			double t59 = rv*rx*t2;
			double t38 = rvx-t59;
			double t61 = rw*rx*t2;
			double t39 = rwx-t61;
			double t40 = ru*ru;
			double t41 = rv*rv;
			double t42 = rw*rw;
			double t69 = ru*ry*t2;
			double t43 = ruy-t69;
			double t72 = rw*ry*t2;
			double t44 = rwy-t72;
			double t45 = atan(alpha);
			double t46 = h*t16*t24;
			double t50 = t8*t45;
			double t47 = t46-t50+1.0/2.0;
			double t80 = ru*rz*t2;
			double t48 = ruz-t80;
			double t82 = rv*rz*t2;
			double t49 = rvz-t82;
			double t51 = h*rx*t4*t24;
			double t52 = alpha*rx*t4*t8*t9*t16*t33;
			double t53 = t51+t52;
			double t54 = gam*rEx;
			double t55 = rx*t4*t40;
			double t56 = rx*t4*t41;
			double t57 = rx*t4*t42;
			double t58 = ru*t2*t3*2.0;
			double t60 = rv*t2*t38*2.0;
			double t62 = rw*t2*t39*2.0;
			double t63 = t55+t56+t57+t58+t60+t62;
			double t90 = t37*t63;
			double t64 = t54-t90;
			double t65 = gam*rEy;
			double t66 = ry*t4*t40;
			double t67 = ry*t4*t41;
			double t68 = ry*t4*t42;
			double t70 = ru*t2*t43*2.0;
			double t71 = rv*t2*t5*2.0;
			double t73 = rw*t2*t44*2.0;
			double t74 = t66+t67+t68+t70+t71+t73;
			double t91 = t37*t74;
			double t75 = t65-t91;
			double t76 = gam*rEz;
			double t77 = rz*t4*t40;
			double t78 = rz*t4*t41;
			double t79 = rz*t4*t42;
			double t81 = ru*t2*t48*2.0;
			double t83 = rv*t2*t49*2.0;
			double t84 = rw*t2*t6*2.0;
			double t85 = t77+t78+t79+t81+t83+t84;
			double t92 = t37*t85;
			double t86 = t76-t92;
			double t87 = h*ry*t4*t24;
			double t88 = alpha*ry*t4*t8*t9*t16*t33;
			double t89 = t87+t88;
			double t93 = h*rz*t4*t24;
			double t94 = alpha*rz*t4*t8*t9*t16*t33;
			double t95 = t93+t94;
			double t96 = h*ru*t4*t24;
			double t97 = alpha*ru*t4*t8*t9*t16*t33;
			double t98 = t96+t97;
			double t99 = epsilon0*rampFactor*t47;
			double t100 = h*t2*t24;
			double t101 = alpha*t2*t8*t9*t16*t33;
			double t102 = t100+t101;
			double t103 = h*rv*t4*t24;
			double t104 = alpha*rv*t4*t8*t9*t16*t33;
			double t105 = t103+t104;
			double t106 = t4*t40;
			double t107 = t4*t41;
			double t108 = t4*t42;
			double t109 = t106+t107+t108;
			double t110 = epsilon0*rampFactor*t37*t47*t109;
			double t111 = epsilon0*rampFactor*rx*t102;
			double t112 = epsilon0*rampFactor*rux*t102;
			double t113 = epsilon0*rampFactor*rvx*t102;
			double t114 = epsilon0*rampFactor*rwx*t102;
			double t115 = epsilon0*rampFactor*t64*t102;
			double t116 = epsilon0*rampFactor*ry*t102;
			double t117 = epsilon0*rampFactor*ruy*t102;
			double t118 = epsilon0*rampFactor*rvy*t102;
			double t119 = epsilon0*rampFactor*rwy*t102;
			double t120 = epsilon0*rampFactor*t75*t102;
			double t121 = epsilon0*rampFactor*rz*t102;
			double t122 = epsilon0*rampFactor*ruz*t102;
			double t123 = epsilon0*rampFactor*rvz*t102;
			double t124 = epsilon0*rampFactor*rwz*t102;
			double t125 = epsilon0*rampFactor*t86*t102;
			double t126 = epsilon0*gam*rampFactor*t47;
			double t127 = h*rw*t4*t24;
			double t128 = alpha*rw*t4*t8*t9*t16*t33;
			double t129 = t127+t128;
			f_udg[0*ng+i] += -epsilon0*rampFactor*rx*t35;
			f_udg[1*ng+i] += -epsilon0*rampFactor*rux*t35;
			f_udg[2*ng+i] += -epsilon0*rampFactor*rvx*t35;
			f_udg[3*ng+i] += -epsilon0*rampFactor*rwx*t35;
			f_udg[4*ng+i] += -epsilon0*rampFactor*t35*t64+epsilon0*rampFactor*t37*t47*(ru*t3*t4*2.0+rv*t4*t38*2.0+rw*t4*t39*2.0);
			f_udg[5*ng+i] += -epsilon0*rampFactor*ry*t35;
			f_udg[6*ng+i] += -epsilon0*rampFactor*ruy*t35;
			f_udg[7*ng+i] += -epsilon0*rampFactor*rvy*t35;
			f_udg[8*ng+i] += -epsilon0*rampFactor*rwy*t35;
			f_udg[9*ng+i] += -epsilon0*rampFactor*t35*t75+epsilon0*rampFactor*t37*t47*(ru*t4*t43*2.0+rv*t4*t5*2.0+rw*t4*t44*2.0);
			f_udg[10*ng+i] += -epsilon0*rampFactor*rz*t35;
			f_udg[11*ng+i] += -epsilon0*rampFactor*ruz*t35;
			f_udg[12*ng+i] += -epsilon0*rampFactor*rvz*t35;
			f_udg[13*ng+i] += -epsilon0*rampFactor*rwz*t35;
			f_udg[14*ng+i] += -epsilon0*rampFactor*t35*t86+epsilon0*rampFactor*t37*t47*(ru*t4*t48*2.0+rv*t4*t49*2.0+rw*t4*t6*2.0);
			f_udg[15*ng+i] += -epsilon0*rampFactor*rx*t53;
			f_udg[16*ng+i] += -epsilon0*rampFactor*rux*t53;
			f_udg[17*ng+i] += -epsilon0*rampFactor*rvx*t53;
			f_udg[18*ng+i] += -epsilon0*rampFactor*rwx*t53;
			f_udg[19*ng+i] += -epsilon0*rampFactor*t53*t64-epsilon0*rampFactor*t2*t3*t37*t47*2.0;
			f_udg[20*ng+i] += -epsilon0*rampFactor*ry*t53;
			f_udg[21*ng+i] += -epsilon0*rampFactor*ruy*t53;
			f_udg[22*ng+i] += -epsilon0*rampFactor*rvy*t53;
			f_udg[23*ng+i] += -epsilon0*rampFactor*rwy*t53;
			f_udg[24*ng+i] += -epsilon0*rampFactor*t53*t75-epsilon0*rampFactor*t2*t37*t43*t47*2.0;
			f_udg[25*ng+i] += -epsilon0*rampFactor*rz*t53;
			f_udg[26*ng+i] += -epsilon0*rampFactor*ruz*t53;
			f_udg[27*ng+i] += -epsilon0*rampFactor*rvz*t53;
			f_udg[28*ng+i] += -epsilon0*rampFactor*rwz*t53;
			f_udg[29*ng+i] += -epsilon0*rampFactor*t53*t86-epsilon0*rampFactor*t2*t37*t47*t48*2.0;
			f_udg[30*ng+i] += -epsilon0*rampFactor*rx*t89;
			f_udg[31*ng+i] += -epsilon0*rampFactor*rux*t89;
			f_udg[32*ng+i] += -epsilon0*rampFactor*rvx*t89;
			f_udg[33*ng+i] += -epsilon0*rampFactor*rwx*t89;
			f_udg[34*ng+i] += -epsilon0*rampFactor*t64*t89-epsilon0*rampFactor*t2*t37*t38*t47*2.0;
			f_udg[35*ng+i] += -epsilon0*rampFactor*ry*t89;
			f_udg[36*ng+i] += -epsilon0*rampFactor*ruy*t89;
			f_udg[37*ng+i] += -epsilon0*rampFactor*rvy*t89;
			f_udg[38*ng+i] += -epsilon0*rampFactor*rwy*t89;
			f_udg[39*ng+i] += -epsilon0*rampFactor*t75*t89-epsilon0*rampFactor*t2*t5*t37*t47*2.0;
			f_udg[40*ng+i] += -epsilon0*rampFactor*rz*t89;
			f_udg[41*ng+i] += -epsilon0*rampFactor*ruz*t89;
			f_udg[42*ng+i] += -epsilon0*rampFactor*rvz*t89;
			f_udg[43*ng+i] += -epsilon0*rampFactor*rwz*t89;
			f_udg[44*ng+i] += -epsilon0*rampFactor*t86*t89-epsilon0*rampFactor*t2*t37*t47*t49*2.0;
			f_udg[45*ng+i] += -epsilon0*rampFactor*rx*t95;
			f_udg[46*ng+i] += -epsilon0*rampFactor*rux*t95;
			f_udg[47*ng+i] += -epsilon0*rampFactor*rvx*t95;
			f_udg[48*ng+i] += -epsilon0*rampFactor*rwx*t95;
			f_udg[49*ng+i] += -epsilon0*rampFactor*t64*t95-epsilon0*rampFactor*t2*t37*t39*t47*2.0;
			f_udg[50*ng+i] += -epsilon0*rampFactor*ry*t95;
			f_udg[51*ng+i] += -epsilon0*rampFactor*ruy*t95;
			f_udg[52*ng+i] += -epsilon0*rampFactor*rvy*t95;
			f_udg[53*ng+i] += -epsilon0*rampFactor*rwy*t95;
			f_udg[54*ng+i] += -epsilon0*rampFactor*t75*t95-epsilon0*rampFactor*t2*t37*t44*t47*2.0;
			f_udg[55*ng+i] += -epsilon0*rampFactor*rz*t95;
			f_udg[56*ng+i] += -epsilon0*rampFactor*ruz*t95;
			f_udg[57*ng+i] += -epsilon0*rampFactor*rvz*t95;
			f_udg[58*ng+i] += -epsilon0*rampFactor*rwz*t95;
			f_udg[59*ng+i] += -epsilon0*rampFactor*t86*t95-epsilon0*rampFactor*t2*t6*t37*t47*2.0;
			f_udg[60*ng+i] += 0.0;
			f_udg[61*ng+i] += 0.0;
			f_udg[62*ng+i] += 0.0;
			f_udg[63*ng+i] += 0.0;
			f_udg[64*ng+i] += 0.0;
			f_udg[65*ng+i] += 0.0;
			f_udg[66*ng+i] += 0.0;
			f_udg[67*ng+i] += 0.0;
			f_udg[68*ng+i] += 0.0;
			f_udg[69*ng+i] += 0.0;
			f_udg[70*ng+i] += 0.0;
			f_udg[71*ng+i] += 0.0;
			f_udg[72*ng+i] += 0.0;
			f_udg[73*ng+i] += 0.0;
			f_udg[74*ng+i] += 0.0;
			f_udg[75*ng+i] += t99-epsilon0*rampFactor*rx*t98;
			f_udg[76*ng+i] += -epsilon0*rampFactor*rux*t98;
			f_udg[77*ng+i] += -epsilon0*rampFactor*rvx*t98;
			f_udg[78*ng+i] += -epsilon0*rampFactor*rwx*t98;
			f_udg[79*ng+i] += t110-epsilon0*rampFactor*t64*t98;
			f_udg[80*ng+i] += -epsilon0*rampFactor*ry*t98;
			f_udg[81*ng+i] += -epsilon0*rampFactor*ruy*t98;
			f_udg[82*ng+i] += -epsilon0*rampFactor*rvy*t98;
			f_udg[83*ng+i] += -epsilon0*rampFactor*rwy*t98;
			f_udg[84*ng+i] += -epsilon0*rampFactor*t75*t98;
			f_udg[85*ng+i] += -epsilon0*rampFactor*rz*t98;
			f_udg[86*ng+i] += -epsilon0*rampFactor*ruz*t98;
			f_udg[87*ng+i] += -epsilon0*rampFactor*rvz*t98;
			f_udg[88*ng+i] += -epsilon0*rampFactor*rwz*t98;
			f_udg[89*ng+i] += -epsilon0*rampFactor*t86*t98;
			f_udg[90*ng+i] += t111;
			f_udg[91*ng+i] += t99+t112;
			f_udg[92*ng+i] += t113;
			f_udg[93*ng+i] += t114;
			f_udg[94*ng+i] += t115-epsilon0*rampFactor*ru*t2*t37*t47*2.0;
			f_udg[95*ng+i] += t116;
			f_udg[96*ng+i] += t117;
			f_udg[97*ng+i] += t118;
			f_udg[98*ng+i] += t119;
			f_udg[99*ng+i] += t120;
			f_udg[100*ng+i] += t121;
			f_udg[101*ng+i] += t122;
			f_udg[102*ng+i] += t123;
			f_udg[103*ng+i] += t124;
			f_udg[104*ng+i] += t125;
			f_udg[105*ng+i] += 0.0;
			f_udg[106*ng+i] += 0.0;
			f_udg[107*ng+i] += t99;
			f_udg[108*ng+i] += 0.0;
			f_udg[109*ng+i] += epsilon0*rampFactor*rv*t2*t37*t47*-2.0;
			f_udg[110*ng+i] += 0.0;
			f_udg[111*ng+i] += 0.0;
			f_udg[112*ng+i] += 0.0;
			f_udg[113*ng+i] += 0.0;
			f_udg[114*ng+i] += 0.0;
			f_udg[115*ng+i] += 0.0;
			f_udg[116*ng+i] += 0.0;
			f_udg[117*ng+i] += 0.0;
			f_udg[118*ng+i] += 0.0;
			f_udg[119*ng+i] += 0.0;
			f_udg[120*ng+i] += 0.0;
			f_udg[121*ng+i] += 0.0;
			f_udg[122*ng+i] += 0.0;
			f_udg[123*ng+i] += t99;
			f_udg[124*ng+i] += epsilon0*rampFactor*rw*t2*t37*t47*-2.0;
			f_udg[125*ng+i] += 0.0;
			f_udg[126*ng+i] += 0.0;
			f_udg[127*ng+i] += 0.0;
			f_udg[128*ng+i] += 0.0;
			f_udg[129*ng+i] += 0.0;
			f_udg[130*ng+i] += 0.0;
			f_udg[131*ng+i] += 0.0;
			f_udg[132*ng+i] += 0.0;
			f_udg[133*ng+i] += 0.0;
			f_udg[134*ng+i] += 0.0;
			f_udg[135*ng+i] += 0.0;
			f_udg[136*ng+i] += 0.0;
			f_udg[137*ng+i] += 0.0;
			f_udg[138*ng+i] += 0.0;
			f_udg[139*ng+i] += t126;
			f_udg[140*ng+i] += 0.0;
			f_udg[141*ng+i] += 0.0;
			f_udg[142*ng+i] += 0.0;
			f_udg[143*ng+i] += 0.0;
			f_udg[144*ng+i] += 0.0;
			f_udg[145*ng+i] += 0.0;
			f_udg[146*ng+i] += 0.0;
			f_udg[147*ng+i] += 0.0;
			f_udg[148*ng+i] += 0.0;
			f_udg[149*ng+i] += 0.0;
			f_udg[150*ng+i] += -epsilon0*rampFactor*rx*t105;
			f_udg[151*ng+i] += -epsilon0*rampFactor*rux*t105;
			f_udg[152*ng+i] += -epsilon0*rampFactor*rvx*t105;
			f_udg[153*ng+i] += -epsilon0*rampFactor*rwx*t105;
			f_udg[154*ng+i] += -epsilon0*rampFactor*t64*t105;
			f_udg[155*ng+i] += t99-epsilon0*rampFactor*ry*t105;
			f_udg[156*ng+i] += -epsilon0*rampFactor*ruy*t105;
			f_udg[157*ng+i] += -epsilon0*rampFactor*rvy*t105;
			f_udg[158*ng+i] += -epsilon0*rampFactor*rwy*t105;
			f_udg[159*ng+i] += t110-epsilon0*rampFactor*t75*t105;
			f_udg[160*ng+i] += -epsilon0*rampFactor*rz*t105;
			f_udg[161*ng+i] += -epsilon0*rampFactor*ruz*t105;
			f_udg[162*ng+i] += -epsilon0*rampFactor*rvz*t105;
			f_udg[163*ng+i] += -epsilon0*rampFactor*rwz*t105;
			f_udg[164*ng+i] += -epsilon0*rampFactor*t86*t105;
			f_udg[165*ng+i] += 0.0;
			f_udg[166*ng+i] += 0.0;
			f_udg[167*ng+i] += 0.0;
			f_udg[168*ng+i] += 0.0;
			f_udg[169*ng+i] += 0.0;
			f_udg[170*ng+i] += 0.0;
			f_udg[171*ng+i] += t99;
			f_udg[172*ng+i] += 0.0;
			f_udg[173*ng+i] += 0.0;
			f_udg[174*ng+i] += epsilon0*rampFactor*ru*t2*t37*t47*-2.0;
			f_udg[175*ng+i] += 0.0;
			f_udg[176*ng+i] += 0.0;
			f_udg[177*ng+i] += 0.0;
			f_udg[178*ng+i] += 0.0;
			f_udg[179*ng+i] += 0.0;
			f_udg[180*ng+i] += t111;
			f_udg[181*ng+i] += t112;
			f_udg[182*ng+i] += t113;
			f_udg[183*ng+i] += t114;
			f_udg[184*ng+i] += t115;
			f_udg[185*ng+i] += t116;
			f_udg[186*ng+i] += t117;
			f_udg[187*ng+i] += t99+t118;
			f_udg[188*ng+i] += t119;
			f_udg[189*ng+i] += t120-epsilon0*rampFactor*rv*t2*t37*t47*2.0;
			f_udg[190*ng+i] += t121;
			f_udg[191*ng+i] += t122;
			f_udg[192*ng+i] += t123;
			f_udg[193*ng+i] += t124;
			f_udg[194*ng+i] += t125;
			f_udg[195*ng+i] += 0.0;
			f_udg[196*ng+i] += 0.0;
			f_udg[197*ng+i] += 0.0;
			f_udg[198*ng+i] += 0.0;
			f_udg[199*ng+i] += 0.0;
			f_udg[200*ng+i] += 0.0;
			f_udg[201*ng+i] += 0.0;
			f_udg[202*ng+i] += 0.0;
			f_udg[203*ng+i] += t99;
			f_udg[204*ng+i] += epsilon0*rampFactor*rw*t2*t37*t47*-2.0;
			f_udg[205*ng+i] += 0.0;
			f_udg[206*ng+i] += 0.0;
			f_udg[207*ng+i] += 0.0;
			f_udg[208*ng+i] += 0.0;
			f_udg[209*ng+i] += 0.0;
			f_udg[210*ng+i] += 0.0;
			f_udg[211*ng+i] += 0.0;
			f_udg[212*ng+i] += 0.0;
			f_udg[213*ng+i] += 0.0;
			f_udg[214*ng+i] += 0.0;
			f_udg[215*ng+i] += 0.0;
			f_udg[216*ng+i] += 0.0;
			f_udg[217*ng+i] += 0.0;
			f_udg[218*ng+i] += 0.0;
			f_udg[219*ng+i] += t126;
			f_udg[220*ng+i] += 0.0;
			f_udg[221*ng+i] += 0.0;
			f_udg[222*ng+i] += 0.0;
			f_udg[223*ng+i] += 0.0;
			f_udg[224*ng+i] += 0.0;
			f_udg[225*ng+i] += -epsilon0*rampFactor*rx*t129;
			f_udg[226*ng+i] += -epsilon0*rampFactor*rux*t129;
			f_udg[227*ng+i] += -epsilon0*rampFactor*rvx*t129;
			f_udg[228*ng+i] += -epsilon0*rampFactor*rwx*t129;
			f_udg[229*ng+i] += -epsilon0*rampFactor*t64*t129;
			f_udg[230*ng+i] += -epsilon0*rampFactor*ry*t129;
			f_udg[231*ng+i] += -epsilon0*rampFactor*ruy*t129;
			f_udg[232*ng+i] += -epsilon0*rampFactor*rvy*t129;
			f_udg[233*ng+i] += -epsilon0*rampFactor*rwy*t129;
			f_udg[234*ng+i] += -epsilon0*rampFactor*t75*t129;
			f_udg[235*ng+i] += t99-epsilon0*rampFactor*rz*t129;
			f_udg[236*ng+i] += -epsilon0*rampFactor*ruz*t129;
			f_udg[237*ng+i] += -epsilon0*rampFactor*rvz*t129;
			f_udg[238*ng+i] += -epsilon0*rampFactor*rwz*t129;
			f_udg[239*ng+i] += t110-epsilon0*rampFactor*t86*t129;
			f_udg[240*ng+i] += 0.0;
			f_udg[241*ng+i] += 0.0;
			f_udg[242*ng+i] += 0.0;
			f_udg[243*ng+i] += 0.0;
			f_udg[244*ng+i] += 0.0;
			f_udg[245*ng+i] += 0.0;
			f_udg[246*ng+i] += 0.0;
			f_udg[247*ng+i] += 0.0;
			f_udg[248*ng+i] += 0.0;
			f_udg[249*ng+i] += 0.0;
			f_udg[250*ng+i] += 0.0;
			f_udg[251*ng+i] += t99;
			f_udg[252*ng+i] += 0.0;
			f_udg[253*ng+i] += 0.0;
			f_udg[254*ng+i] += epsilon0*rampFactor*ru*t2*t37*t47*-2.0;
			f_udg[255*ng+i] += 0.0;
			f_udg[256*ng+i] += 0.0;
			f_udg[257*ng+i] += 0.0;
			f_udg[258*ng+i] += 0.0;
			f_udg[259*ng+i] += 0.0;
			f_udg[260*ng+i] += 0.0;
			f_udg[261*ng+i] += 0.0;
			f_udg[262*ng+i] += 0.0;
			f_udg[263*ng+i] += 0.0;
			f_udg[264*ng+i] += 0.0;
			f_udg[265*ng+i] += 0.0;
			f_udg[266*ng+i] += 0.0;
			f_udg[267*ng+i] += t99;
			f_udg[268*ng+i] += 0.0;
			f_udg[269*ng+i] += epsilon0*rampFactor*rv*t2*t37*t47*-2.0;
			f_udg[270*ng+i] += t111;
			f_udg[271*ng+i] += t112;
			f_udg[272*ng+i] += t113;
			f_udg[273*ng+i] += t114;
			f_udg[274*ng+i] += t115;
			f_udg[275*ng+i] += t116;
			f_udg[276*ng+i] += t117;
			f_udg[277*ng+i] += t118;
			f_udg[278*ng+i] += t119;
			f_udg[279*ng+i] += t120;
			f_udg[280*ng+i] += t121;
			f_udg[281*ng+i] += t122;
			f_udg[282*ng+i] += t123;
			f_udg[283*ng+i] += t99+t124;
			f_udg[284*ng+i] += t125-epsilon0*rampFactor*rw*t2*t37*t47*2.0;
			f_udg[285*ng+i] += 0.0;
			f_udg[286*ng+i] += 0.0;
			f_udg[287*ng+i] += 0.0;
			f_udg[288*ng+i] += 0.0;
			f_udg[289*ng+i] += 0.0;
			f_udg[290*ng+i] += 0.0;
			f_udg[291*ng+i] += 0.0;
			f_udg[292*ng+i] += 0.0;
			f_udg[293*ng+i] += 0.0;
			f_udg[294*ng+i] += 0.0;
			f_udg[295*ng+i] += 0.0;
			f_udg[296*ng+i] += 0.0;
			f_udg[297*ng+i] += 0.0;
			f_udg[298*ng+i] += 0.0;
			f_udg[299*ng+i] += t126;
		}

	}
}

