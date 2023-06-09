
#define PI (3.141592653589793)

// Written by: C. Nguyen & P. Fernandez

void source_ransSA2dNEW(double *s, double *s_udg, double *pg, double *udg, appstruct &app, double *param,
					    double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];

	/* Allocate dynamic memory */
	double * mu  = new double [ng];
	double * mu_udg;
	double * s_mu;
	if (computeJacobian == 1 && app.viscosityModel != 0) {
		mu_udg = new double [ng * nc];
		s_mu = new double [ng * ncu];
	}

	// Closure parameters of SA model:
	double rMax = 5.0;     // rMax in [5,10] is typically used.
	double cv1 = 7.1;
	double sigma = 2.0/3.0;
	double cb1 = 0.1355;
	double cb2 = 0.622;
	double kappa = 0.41;
	double cw1 = cb1/(kappa*kappa) + (1+cb2)/sigma;
	double cw2 = 0.3;
	double cw3 = 2.0;
	double ct1 = 1.0;      // From David's paper
	double ct2 = 2.0;      // From David's paper
	double ct3 = 1.2;      // From NASA Langley webpage (ct3 = 1.1 from David's paper)
	double ct4 = 0.5;      // From NASA Langley webpage (ct4 = 2.0 from David's paper)

	// SA regularization parameters (Ref. Hemant Chaurasia PhD Thesis):
	double b = 100.0;
	double c = 1.0/2.0 - atan(b)/PI;

	// No laminar suppression and forced transition are considered:
	double deltaU = 0.0;
	double ft1 = 0.0;             // ft1 = ct1*gt*exp(-ct2*(wt*wt/(deltaU*deltaU))*(walDistance*wallDistance+gt*gt*dt*dt))
	double ft2 = 0.0;             // ft2 = ct3*exp(-ct4*psi [or chi])

	getViscosity(mu, mu_udg, pg, udg, param, app.viscosityModel, ng, nc, ncu, nd, computeJacobian);

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];

		double wallDistance = max(pg[3*ng+i], 1.0e-8);

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

		double t4 = 1.0/u1;
		double t37 = t4*u5*u6;
		double t2 = -t37+u10;
		double t3 = 1.0/(u1*u1);
		double t38 = t4*u5*u11;
		double t5 = -t38+u15;
		double t6 = 1.0/mu[i];
		double t18 = t4*u3*u6;
		double t19 = t18-u8;
		double t20 = t4*t19;
		double t21 = t4*u2*u11;
		double t22 = t4*(t21-u12);
		double t7 = -t20+t22;
		double t8 = 1.0/3.141592653589793;
		double t9 = 1.0/param3;
		double t10 = b*t6*u5;
		double t11 = atan(t10);
		double t12 = t8*t11;
		double t13 = t12+1.0/2.0;
		double t14 = t6*t13*u5;
		double t15 = c+t14;
		double t16 = t15*t15;
		double t17 = cv1*cv1;
		double t23 = t7*t7;
		double t24 = t23+1.0E-20;
		double t25 = sqrt(t24);
		double t26 = 1.0/(kappa*kappa);
		double t27 = 1.0/(wallDistance*wallDistance);
		double t28 = t16*t16;
		double t29 = cv1*t17+t15*t16;
		double t30 = 1.0/t29;
		double t31 = t28*t30;
		double t32 = t31+1.0;
		double t33 = 1.0/t32;
		double t34 = t15*t33;
		double t35 = t34-1.0;
		double t36 = 1.0/sigma;
		double t39 = cw3*cw3;
		double t40 = t39*t39;
		double t41 = t25*(1.0/1.0E1);
		double t42 = 1.0/sqrt(t24);
		double t43 = mu[i]*t4*t9*t15*t26*t27*t35*t42;
		double t44 = t43-9.0/1.0E1;
		double t45 = b*t44;
		double t46 = atan(t45);
		double t47 = t8*t46;
		double t48 = t47-1.0/2.0;
		double t49 = t25*(9.0/1.0E1);
		double t53 = mu[i]*t4*t9*t15*t26*t27*t35;
		double t50 = t49-t53;
		double t51 = c*t25;
		double t54 = t48*t50;
		double t52 = t41+t51-t54;
		double t55 = 1.0/t52;
		double t60 = mu[i]*t4*t9*t15*t26*t27*t55;
		double t56 = rMax-t60;
		double t61 = b*t56;
		double t62 = atan(t61);
		double t63 = t8*t62;
		double t64 = t63+1.0/2.0;
		double t65 = t56*t64;
		double t57 = c-rMax+t65;
		double t58 = t57*t57;
		double t59 = t58*t58;
		double t69 = c-rMax+t65+t58*t59;
		double t70 = cw2*t69;
		double t66 = c-rMax+t65-t70;
		double t67 = t66*t66;
		double t68 = t67*t67;

		s[0*ng+i] = 0.0;
		s[1*ng+i] = 0.0;
		s[2*ng+i] = 0.0;
		s[3*ng+i] = 0.0;
		s[4*ng+i] = (deltaU*deltaU)*ft1*param3*u1+cb2*t9*t36*u1*((t2*t2)*t3+t3*(t5*t5))-cb1*mu[i]*t15*t52*(ft2-1.0)-mu[i]*t4*t9*t36*(t2*t4*u6+t4*t5*u11)*(c+t14+1.0)+(mu[i]*mu[i])*t4*t9*t16*t27*(cb1*ft2*t26+cw1*t66*pow((t39*t40+1.0)/(t39*t40+t67*t68),1.0/6.0));
	}
    
//     for (int i = 0; i <ng*5; i++) {
//         s[i] /= param3;
//     }

	if (computeJacobian == 1) {

		for (int i = 0; i <ng; i++) {
			double x1 = pg[0*ng+i];
			double x2 = pg[1*ng+i];

			double wallDistance = max(pg[3*ng+i], 1.0e-8);

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

			double t4 = 1.0/u1;
			double t8 = t4*u5*u6;
			double t2 = -t8+u10;
			double t3 = 1.0/(u1*u1);
			double t11 = t4*u5*u11;
			double t5 = -t11+u15;
			double t6 = 1.0/param3;
			double t7 = 1.0/sigma;
			double t9 = t2*t2;
			double t10 = 1.0/(u1*u1*u1);
			double t12 = t5*t5;
			double t13 = 1.0/(u1*u1*u1*u1);
			double t14 = 1.0/mu[i];
			double t15 = 1.0/3.141592653589793;
			double t16 = b*t14*u5;
			double t17 = atan(t16);
			double t18 = t15*t17;
			double t19 = t18+1.0/2.0;
			double t20 = t14*t19*u5;
			double t21 = c+t20;
			double t22 = t21*t21;
			double t23 = cv1*cv1;
			double t25 = t4*u3*u6;
			double t26 = t25-u8;
			double t27 = t4*t26;
			double t28 = t4*u2*u11;
			double t33 = t28-u12;
			double t29 = t4*t33;
			double t24 = -t27+t29;
			double t30 = t24*t24;
			double t31 = t30+1.0E-20;
			double t32 = 1.0/sqrt(t31);
			double t34 = 1.0/(kappa*kappa);
			double t35 = 1.0/(wallDistance*wallDistance);
			double t36 = t22*t22;
			double t37 = cv1*t23+t21*t22;
			double t38 = 1.0/t37;
			double t39 = t36*t38;
			double t40 = t39+1.0;
			double t41 = 1.0/t40;
			double t42 = t21*t41;
			double t43 = t42-1.0;
			double t44 = t4*(t28-u12);
			double t45 = -t27+t44;
			double t46 = t3*t26;
			double t47 = t10*u3*u6;
			double t48 = t45*t45;
			double t49 = t48+1.0E-20;
			double t50 = 1.0/sqrt(t49);
			double t59 = mu[i]*t4*t6*t21*t34*t35*t43*t50;
			double t51 = t59-9.0/1.0E1;
			double t53 = t3*t33;
			double t54 = t10*u2*u11;
			double t52 = t46+t47-t53-t54;
			double t55 = c+t20+1.0;
			double t56 = cw3*cw3;
			double t57 = t56*t56;
			double t58 = sqrt(t49);
			double t60 = t58*(9.0/1.0E1);
			double t67 = mu[i]*t4*t6*t21*t34*t35*t43;
			double t61 = t60-t67;
			double t62 = t58*(1.0/1.0E1);
			double t63 = b*t51;
			double t64 = atan(t63);
			double t65 = t15*t64;
			double t66 = t65-1.0/2.0;
			double t68 = c*t58;
			double t75 = t61*t66;
			double t69 = t62+t68-t75;
			double t70 = 1.0/t69;
			double t76 = mu[i]*t4*t6*t21*t34*t35*t70;
			double t71 = rMax-t76;
			double t77 = b*t71;
			double t78 = atan(t77);
			double t79 = t15*t78;
			double t80 = t79+1.0/2.0;
			double t81 = t71*t80;
			double t72 = c-rMax+t81;
			double t73 = t72*t72;
			double t74 = t73*t73;
			double t117 = c-rMax+t81+t73*t74;
			double t118 = cw2*t117;
			double t82 = c-rMax+t81-t118;
			double t83 = t82*t82;
			double t84 = t83*t83;
			double t85 = mu[i]*t3*t6*t21*t34*t35*t43;
			double t86 = mu[i]*t3*t6*t21*t34*t35*t43*t50;
			double t87 = 1.0/pow(t49,3.0/2.0);
			double t88 = b*b;
			double t89 = t51*t51;
			double t90 = t88*t89;
			double t91 = t90+1.0;
			double t92 = 1.0/t91;
			double t93 = mu[i]*t3*t6*t21*t34*t35*t70;
			double t94 = 1.0/(t69*t69);
			double t95 = t45*t50*t52*(9.0/1.0E1);
			double t96 = t85+t95;
			double t98 = t3*(t28-u12);
			double t99 = t46+t47-t54-t98;
			double t97 = t45*t50*t99*(1.0/1.0E1);
			double t100 = c*t45*t50*(t46+t47-t54-t98);
			double t101 = mu[i]*t4*t6*t21*t34*t35*t43*t45*t87*(t46+t47-t54-t98);
			double t102 = t86+t101;
			double t103 = b*t15*t61*t92*t102;
			double t104 = t45*t50*t99*(9.0/1.0E1);
			double t105 = t85+t104;
			double t106 = t45*t50*(t46+t47-t54-t98)*(1.0/1.0E1);
			double t110 = t66*t105;
			double t107 = t100+t103+t106-t110;
			double t108 = mu[i]*t4*t6*t21*t34*t35*t94*t107;
			double t109 = t93+t108;
			double t111 = t71*t71;
			double t112 = t88*t111;
			double t113 = t112+1.0;
			double t114 = 1.0/t113;
			double t115 = b*t15*t71*t109*t114;
			double t116 = t56*t57+1.0;
			double t119 = t56*t57+t83*t84;
			double t120 = 1.0/t119;
			double t121 = t116*t120;
			double t122 = t80*t109;
			double t123 = t115+t122;
			double t124 = t72*t74*t123*6.0;
			double t125 = mu[i]*mu[i];
			double t126 = pow(t121,1.0/6.0);
			double t127 = ft2-1.0;
			double t128 = t3*t45*t50*u11*(1.0/1.0E1);
			double t129 = c*t3*t45*t50*u11;
			double t130 = b*mu[i]*t6*t10*t15*t21*t34*t35*t43*t45*t61*t87*t92*u11;
			double t132 = t3*t45*t50*t66*u11*(9.0/1.0E1);
			double t131 = t128+t129+t130-t132;
			double t133 = mu[i]*t4*t6*t21*t34*t35*t80*t94*t131;
			double t134 = b*mu[i]*t4*t6*t15*t21*t34*t35*t71*t94*t114*t131;
			double t135 = 1.0/pow(t121,5.0/6.0);
			double t136 = 1.0/(t119*t119);
			double t137 = t133+t134;
			double t138 = t72*t74*t137*6.0;
			double t139 = t133+t134+t138;
			double t140 = t133+t134-cw2*t139;
			double t141 = t3*t45*t50*u6*(1.0/1.0E1);
			double t142 = c*t3*t45*t50*u6;
			double t143 = b*mu[i]*t6*t10*t15*t21*t34*t35*t43*t45*t61*t87*t92*u6;
			double t145 = t3*t45*t50*t66*u6*(9.0/1.0E1);
			double t144 = t141+t142+t143-t145;
			double t146 = mu[i]*t4*t6*t21*t34*t35*t80*t94*t144;
			double t147 = b*mu[i]*t4*t6*t15*t21*t34*t35*t71*t94*t114*t144;
			double t148 = t146+t147;
			double t149 = t72*t74*t148*6.0;
			double t150 = t146+t147+t149;
			double t151 = t146+t147-cw2*t150;
			double t152 = 1.0/(mu[i]*mu[i]);
			double t153 = t14*t19;
			double t154 = u5*u5;
			double t155 = t88*t152*t154;
			double t156 = t155+1.0;
			double t157 = 1.0/t156;
			double t158 = b*t15*t152*t157*u5;
			double t159 = t153+t158;
			double t160 = t41*t159;
			double t161 = 1.0/(t40*t40);
			double t162 = t21*t22*t38*t159*4.0;
			double t163 = 1.0/(t37*t37);
			double t167 = t22*t36*t159*t163*3.0;
			double t164 = t162-t167;
			double t168 = t21*t161*t164;
			double t165 = t160-t168;
			double t166 = mu[i]*t4*t6*t34*t35*t43*t159;
			double t169 = mu[i]*t4*t6*t21*t34*t35*t165;
			double t170 = t166+t169;
			double t171 = t66*t170;
			double t172 = mu[i]*t4*t6*t34*t35*t43*t50*t159;
			double t173 = mu[i]*t4*t6*t21*t34*t35*t50*t165;
			double t174 = t172+t173;
			double t177 = b*t15*t61*t92*t174;
			double t175 = t171-t177;
			double t176 = mu[i]*t4*t6*t34*t35*t70*t159;
			double t180 = mu[i]*t4*t6*t21*t34*t35*t94*t175;
			double t178 = t176-t180;
			double t179 = t80*t178;
			double t181 = b*t15*t71*t114*t178;
			double t182 = t179+t181;
			double t183 = t72*t74*t182*6.0;
			double t184 = t179+t181+t183;
			double t185 = t179+t181-cw2*t184;
			double t186 = t2*t4*u6;
			double t187 = t4*t5*u11;
			double t188 = t186+t187;
			double t189 = u6*u6;
			double t190 = u11*u11;
			double t191 = cw1*t82*t126;
			double t192 = cb1*ft2*t34;
			double t193 = t191+t192;
			double t194 = t3*t45*t50*u3*(1.0/1.0E1);
			double t195 = c*t3*t45*t50*u3;
			double t196 = b*mu[i]*t6*t10*t15*t21*t34*t35*t43*t45*t61*t87*t92*u3;
			double t198 = t3*t45*t50*t66*u3*(9.0/1.0E1);
			double t197 = t194+t195+t196-t198;
			double t199 = mu[i]*t4*t6*t21*t34*t35*t80*t94*t197;
			double t200 = b*mu[i]*t4*t6*t15*t21*t34*t35*t71*t94*t114*t197;
			double t201 = t199+t200;
			double t202 = t72*t74*t201*6.0;
			double t203 = t199+t200+t202;
			double t204 = t199+t200-cw2*t203;
			double t205 = t4*t45*t50*(1.0/1.0E1);
			double t206 = c*t4*t45*t50;
			double t207 = b*mu[i]*t3*t6*t15*t21*t34*t35*t43*t45*t61*t87*t92;
			double t209 = t4*t45*t50*t66*(9.0/1.0E1);
			double t208 = t205+t206+t207-t209;
			double t210 = mu[i]*t4*t6*t21*t34*t35*t80*t94*t208;
			double t211 = b*mu[i]*t4*t6*t15*t21*t34*t35*t71*t94*t114*t208;
			double t212 = t210+t211;
			double t213 = t72*t74*t212*6.0;
			double t214 = t210+t211+t213;
			double t227 = cw2*t214;
			double t215 = t210+t211-t227;
			double t216 = t3*t45*t50*u2*(1.0/1.0E1);
			double t217 = c*t3*t45*t50*u2;
			double t218 = b*mu[i]*t6*t10*t15*t21*t34*t35*t43*t45*t61*t87*t92*u2;
			double t220 = t3*t45*t50*t66*u2*(9.0/1.0E1);
			double t219 = t216+t217+t218-t220;
			double t221 = mu[i]*t4*t6*t21*t34*t35*t80*t94*t219;
			double t222 = b*mu[i]*t4*t6*t15*t21*t34*t35*t71*t94*t114*t219;
			double t223 = t221+t222;
			double t224 = t72*t74*t223*6.0;
			double t225 = t221+t222+t224;
			double t226 = t221+t222-cw2*t225;
			double t228 = cw1*t126*t215;
			double t229 = t228-cw1*t83*t84*t116*t135*t136*t215;
			double t230 = t4*t6*t22*t35*t125*t229;

			s_udg[0*ng+i] = 0.0;
			s_udg[1*ng+i] = 0.0;
			s_udg[2*ng+i] = 0.0;
			s_udg[3*ng+i] = 0.0;
			s_udg[4*ng+i] = (deltaU*deltaU)*ft1*param3+cb2*t6*t7*(t3*t9+t3*t12)-cb2*t6*t7*u1*(t9*t10*2.0+t10*t12*2.0-t2*t13*u5*u6*2.0-t5*t13*u5*u11*2.0)-cb1*mu[i]*t21*t127*(-(t85+t24*t32*t52*(9.0/1.0E1))*(t15*atan(b*(mu[i]*t4*t6*t21*t32*t34*t35*t43-9.0/1.0E1))-1.0/2.0)+t45*t50*t52*(1.0/1.0E1)+c*t45*t50*t52+b*t15*t61*t92*(t86+mu[i]*t4*t6*t21*t34*t35*t43*t45*t52*t87))+mu[i]*t3*t6*t7*t55*t188-t3*t6*t22*t35*t125*t193+t4*t6*t22*t35*t125*(cw1*t126*(t115-cw2*(t115+t124+t80*(t93+mu[i]*t4*t6*t21*t34*t35*t94*(t97+t100+t103-t66*t96)))+t80*(t93+mu[i]*t4*t6*t21*t34*t35*t94*(t97-t66*t96+c*t45*t50*(t46+t47-t54-t3*(t28-u12))+b*t15*t61*t92*(t86+mu[i]*t4*t6*t21*t34*t35*t43*t45*t87*(t46+t47-t54-t3*(t28-u12))))))-cw1*t83*t84*t116*t135*t136*(t115+t122-cw2*(t115+t122+t124)))+mu[i]*t4*t6*t7*t55*(t2*t3*u6+t3*t5*u11-t10*t189*u5-t10*t190*u5);
			s_udg[5*ng+i] = 0.0;
			s_udg[6*ng+i] = 0.0;
			s_udg[7*ng+i] = 0.0;
			s_udg[8*ng+i] = 0.0;
			s_udg[9*ng+i] = -cb1*mu[i]*t21*t127*t131+t4*t6*t22*t35*t125*(cw1*t126*t140-cw1*t83*t84*t116*t135*t136*t140);
			s_udg[10*ng+i] = 0.0;
			s_udg[11*ng+i] = 0.0;
			s_udg[12*ng+i] = 0.0;
			s_udg[13*ng+i] = 0.0;
			s_udg[14*ng+i] = cb1*mu[i]*t21*t127*t144-t4*t6*t22*t35*t125*(cw1*t126*t151-cw1*t83*t84*t116*t135*t136*t151);
			s_udg[15*ng+i] = 0.0;
			s_udg[16*ng+i] = 0.0;
			s_udg[17*ng+i] = 0.0;
			s_udg[18*ng+i] = 0.0;
			s_udg[19*ng+i] = 0.0;
			s_udg[20*ng+i] = 0.0;
			s_udg[21*ng+i] = 0.0;
			s_udg[22*ng+i] = 0.0;
			s_udg[23*ng+i] = 0.0;
			s_udg[24*ng+i] = -cb1*mu[i]*t21*t127*t175-cb1*mu[i]*t69*t127*t159-cb2*t6*t7*u1*(t2*t10*u6*2.0+t5*t10*u11*2.0)-mu[i]*t4*t6*t7*t159*t188-t4*t6*t22*t35*t125*(cw1*t126*t185-cw1*t83*t84*t116*t135*t136*t185)+mu[i]*t4*t6*t7*t55*(t3*t189+t3*t190)+t4*t6*t21*t35*t125*t159*t193*2.0;
			s_udg[25*ng+i] = 0.0;
			s_udg[26*ng+i] = 0.0;
			s_udg[27*ng+i] = 0.0;
			s_udg[28*ng+i] = 0.0;
			s_udg[29*ng+i] = cb1*mu[i]*t21*t127*t197-mu[i]*t4*t6*t7*t55*(t2*t4-t3*u5*u6)-cb2*t2*t3*t6*t7*u5*2.0-t4*t6*t22*t35*t125*(cw1*t126*t204-cw1*t83*t84*t116*t135*t136*t204);
			s_udg[30*ng+i] = 0.0;
			s_udg[31*ng+i] = 0.0;
			s_udg[32*ng+i] = 0.0;
			s_udg[33*ng+i] = 0.0;
			s_udg[34*ng+i] = 0.0;
			s_udg[35*ng+i] = 0.0;
			s_udg[36*ng+i] = 0.0;
			s_udg[37*ng+i] = 0.0;
			s_udg[38*ng+i] = 0.0;
			s_udg[39*ng+i] = t230-cb1*mu[i]*t21*t127*t208;
			s_udg[40*ng+i] = 0.0;
			s_udg[41*ng+i] = 0.0;
			s_udg[42*ng+i] = 0.0;
			s_udg[43*ng+i] = 0.0;
			s_udg[44*ng+i] = 0.0;
			s_udg[45*ng+i] = 0.0;
			s_udg[46*ng+i] = 0.0;
			s_udg[47*ng+i] = 0.0;
			s_udg[48*ng+i] = 0.0;
			s_udg[49*ng+i] = cb2*t4*t6*t7*(u10*2.0-t4*u5*u6*2.0)-mu[i]*t3*t6*t7*t55*u6;
			s_udg[50*ng+i] = 0.0;
			s_udg[51*ng+i] = 0.0;
			s_udg[52*ng+i] = 0.0;
			s_udg[53*ng+i] = 0.0;
			s_udg[54*ng+i] = -cb1*mu[i]*t21*t127*t219-mu[i]*t4*t6*t7*t55*(t4*t5-t3*u5*u11)-cb2*t3*t5*t6*t7*u5*2.0+t4*t6*t22*t35*t125*(cw1*t126*t226-cw1*t83*t84*t116*t135*t136*t226);
			s_udg[55*ng+i] = 0.0;
			s_udg[56*ng+i] = 0.0;
			s_udg[57*ng+i] = 0.0;
			s_udg[58*ng+i] = 0.0;
			s_udg[59*ng+i] = -t230+cb1*mu[i]*t21*t127*t208;
			s_udg[60*ng+i] = 0.0;
			s_udg[61*ng+i] = 0.0;
			s_udg[62*ng+i] = 0.0;
			s_udg[63*ng+i] = 0.0;
			s_udg[64*ng+i] = 0.0;
			s_udg[65*ng+i] = 0.0;
			s_udg[66*ng+i] = 0.0;
			s_udg[67*ng+i] = 0.0;
			s_udg[68*ng+i] = 0.0;
			s_udg[69*ng+i] = 0.0;
			s_udg[70*ng+i] = 0.0;
			s_udg[71*ng+i] = 0.0;
			s_udg[72*ng+i] = 0.0;
			s_udg[73*ng+i] = 0.0;
			s_udg[74*ng+i] = cb2*t4*t6*t7*(u15*2.0-t4*u5*u11*2.0)-mu[i]*t3*t6*t7*t55*u11;
		}

		if (app.viscosityModel != 0) {

			for (int i = 0; i <ng; i++) {
				double x1 = pg[0*ng+i];
				double x2 = pg[1*ng+i];

				double wallDistance = max(pg[3*ng+i], 1.0e-8);

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

				double t2 = 1.0/mu[i];
				double t3 = 1.0/u1;
				double t14 = t3*u3*u6;
				double t15 = t14-u8;
				double t16 = t3*t15;
				double t17 = t3*u2*u11;
				double t18 = t3*(t17-u12);
				double t4 = -t16+t18;
				double t5 = 1.0/3.141592653589793;
				double t6 = b*t2*u5;
				double t7 = atan(t6);
				double t8 = t5*t7;
				double t9 = t8+1.0/2.0;
				double t10 = t2*t9*u5;
				double t11 = c+t10;
				double t12 = t11*t11;
				double t13 = cv1*cv1;
				double t19 = t4*t4;
				double t20 = t19+1.0E-20;
				double t21 = sqrt(t20);
				double t22 = 1.0/(kappa*kappa);
				double t23 = 1.0/param3;
				double t24 = 1.0/(wallDistance*wallDistance);
				double t25 = t12*t12;
				double t26 = cv1*t13+t11*t12;
				double t27 = 1.0/t26;
				double t28 = t25*t27;
				double t29 = t28+1.0;
				double t30 = 1.0/t29;
				double t31 = t11*t30;
				double t32 = t31-1.0;
				double t33 = 1.0/sqrt(t20);
				double t34 = mu[i]*t3*t11*t22*t23*t24*t32*t33;
				double t35 = t34-9.0/1.0E1;
				double t36 = b*t35;
				double t37 = atan(t36);
				double t38 = t5*t37;
				double t39 = t38-1.0/2.0;
				double t40 = 1.0/(mu[i]*mu[i]);
				double t41 = u5*u5;
				double t42 = t9*t40*u5;
				double t43 = 1.0/(mu[i]*mu[i]*mu[i]);
				double t44 = b*b;
				double t45 = t40*t41*t44;
				double t46 = t45+1.0;
				double t47 = 1.0/t46;
				double t48 = b*t5*t41*t43*t47;
				double t49 = t42+t48;
				double t50 = t21*(9.0/1.0E1);
				double t60 = mu[i]*t3*t11*t22*t23*t24*t32;
				double t51 = t50-t60;
				double t52 = t30*t49;
				double t53 = t11*t12*t27*t49*4.0;
				double t54 = 1.0/(t26*t26);
				double t94 = t12*t25*t49*t54*3.0;
				double t55 = t53-t94;
				double t56 = 1.0/(t29*t29);
				double t95 = t11*t55*t56;
				double t57 = t52-t95;
				double t58 = ft2-1.0;
				double t59 = t21*(1.0/1.0E1);
				double t61 = c*t21;
				double t71 = t39*t51;
				double t62 = t59+t61-t71;
				double t63 = 1.0/sigma;
				double t64 = u10-t3*u5*u6;
				double t65 = t3*t64*u6;
				double t66 = u15-t3*u5*u11;
				double t67 = t3*t66*u11;
				double t68 = t65+t67;
				double t69 = cw3*cw3;
				double t70 = t69*t69;
				double t72 = 1.0/t62;
				double t77 = mu[i]*t3*t11*t22*t23*t24*t72;
				double t73 = rMax-t77;
				double t78 = b*t73;
				double t79 = atan(t78);
				double t80 = t5*t79;
				double t81 = t80+1.0/2.0;
				double t82 = t73*t81;
				double t74 = c-rMax+t82;
				double t75 = t74*t74;
				double t76 = t75*t75;
				double t86 = c-rMax+t82+t75*t76;
				double t87 = cw2*t86;
				double t83 = c-rMax+t82-t87;
				double t84 = t83*t83;
				double t85 = t84*t84;
				double t88 = t69*t70+1.0;
				double t89 = t69*t70+t84*t85;
				double t90 = 1.0/t89;
				double t91 = t88*t90;
				double t92 = pow(t91,1.0/6.0);
				double t93 = mu[i]*t3*t22*t23*t24*t32*t49;
				double t96 = mu[i]*t3*t11*t22*t23*t24*t57;
				double t108 = t3*t11*t22*t23*t24*t32;
				double t97 = t93+t96-t108;
				double t98 = t39*t97;
				double t99 = t35*t35;
				double t100 = t44*t99;
				double t101 = t100+1.0;
				double t102 = 1.0/t101;
				double t103 = mu[i]*t3*t11*t22*t23*t24*t33*t57;
				double t104 = mu[i]*t3*t22*t23*t24*t32*t33*t49;
				double t109 = t3*t11*t22*t23*t24*t32*t33;
				double t105 = t103+t104-t109;
				double t110 = b*t5*t51*t102*t105;
				double t106 = t98-t110;
				double t107 = t3*t11*t22*t23*t24*t72;
				double t111 = 1.0/(t62*t62);
				double t112 = mu[i]*t3*t11*t22*t23*t24*t106*t111;
				double t115 = mu[i]*t3*t22*t23*t24*t49*t72;
				double t113 = t107+t112-t115;
				double t114 = t81*t113;
				double t116 = t73*t73;
				double t117 = t44*t116;
				double t118 = t117+1.0;
				double t119 = 1.0/t118;
				double t120 = b*t5*t73*t113*t119;
				double t121 = t114+t120;
				double t122 = t74*t76*t121*6.0;
				double t123 = t114+t120+t122;
				double t124 = t114+t120-cw2*t123;
				double t125 = mu[i]*mu[i];
				double t126 = cw1*t83*t92;
				double t127 = cb1*ft2*t22;
				double t128 = t126+t127;

				s_mu[0*ng+i] = 0.0;
				s_mu[1*ng+i] = 0.0;
				s_mu[2*ng+i] = 0.0;
				s_mu[3*ng+i] = 0.0;
				s_mu[4*ng+i] = -cb1*t11*t58*t62-t3*t23*t63*t68*(c+t10+1.0)+cb1*mu[i]*t49*t58*t62+cb1*mu[i]*t11*t58*t106+mu[i]*t3*t12*t23*t24*t128*2.0+mu[i]*t3*t23*t49*t63*t68-t3*t12*t23*t24*t125*(cw1*t92*t124-cw1*t84*t85*t88*1.0/(t89*t89)*1.0/pow(t91,5.0/6.0)*t124)-t3*t11*t23*t24*t49*t125*t128*2.0;
			}

			chainRule_s_udg(s_udg, s_mu, mu_udg, ng, nc, ncu, nd);
		}
        
//         for (int i = 0; i <ng*5*15; i++) {
//             s_udg[i] /= param3;
//         }
	}
    
	/* Deallocate dynamic memory */
	delete[] mu;
	if (computeJacobian == 1 && app.viscosityModel != 0) {
		delete[] mu_udg;
		delete[] s_mu;
	}
}
