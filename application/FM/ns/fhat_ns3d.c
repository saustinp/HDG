void fhat_ns3d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = 1.0/uh1;
		double t3 = uh2*uh2;
		double t4 = 1.0/param3;
		double t5 = 1.0/(uh1*uh1);
		double t6 = uh3*uh3;
		double t7 = u19*uh1*2.0;
		double t8 = param1-1.0;
		double t9 = uh4*uh4;
		double t18 = uh1*uh5*2.0;
		double t10 = t3+t6+t9-t18;
		double t11 = t2*t8*t10*(1.0/2.0);
		double t12 = u8*uh1;
		double t13 = u12*uh1;
		double t14 = param3*uh1*uh2*uh3;
		double t15 = t12+t13+t14-u6*uh3-u11*uh2;
		double t16 = u7*uh1*2.0;
		double t17 = u13*uh1*2.0;
		double t19 = u9*uh1;
		double t20 = u17*uh1;
		double t21 = param3*uh1*uh2*uh4;
		double t22 = t19+t20+t21-u6*uh4-u16*uh2;
		double t23 = u14*uh1;
		double t24 = u18*uh1;
		double t25 = param3*uh1*uh3*uh4;
		double t26 = t23+t24+t25-u11*uh4-u16*uh3;
		double t27 = uh1*uh1;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+nl3*uh4+param6*(u1-uh1);
		fh[1*ng+i] = param6*(u2-uh2)-nl1*(t11-t2*t3+t4*t5*(t7+t17+u6*uh2*4.0-u7*uh1*4.0-u11*uh3*2.0-u16*uh4*2.0)*(1.0/3.0))+nl2*t4*t5*t15+nl3*t4*t5*t22;
		fh[2*ng+i] = param6*(u3-uh3)-nl2*(t11-t2*t6+t4*t5*(t7+t16-u6*uh2*2.0+u11*uh3*4.0-u13*uh1*4.0-u16*uh4*2.0)*(1.0/3.0))+nl1*t4*t5*t15+nl3*t4*t5*t26;
		fh[3*ng+i] = param6*(u4-uh4)-nl3*(t11-t2*t9+t4*t5*(t16+t17-u6*uh2*2.0-u11*uh3*2.0+u16*uh4*4.0-u19*uh1*4.0)*(1.0/3.0))+nl1*t4*t5*t22+nl2*t4*t5*t26;
		fh[4*ng+i] = (t4*1.0/(uh1*uh1*uh1)*(nl1*param1*t3*u6*-6.0-nl1*param1*t6*u6*6.0+nl1*param4*t3*u6*8.0-nl1*param1*t9*u6*6.0+nl1*param4*t6*u6*6.0-nl2*param1*t3*u11*6.0+nl1*param4*t9*u6*6.0-nl2*param1*t6*u11*6.0+nl2*param4*t3*u11*6.0-nl2*param1*t9*u11*6.0+nl2*param4*t6*u11*8.0-nl3*param1*t3*u16*6.0+nl2*param4*t9*u11*6.0-nl3*param1*t6*u16*6.0+nl3*param4*t3*u16*6.0-nl3*param1*t9*u16*6.0+nl3*param4*t6*u16*6.0+nl3*param4*t9*u16*8.0-nl1*param1*t27*u10*6.0-nl2*param1*t27*u15*6.0-nl3*param1*t27*u20*6.0+nl1*param1*u7*uh1*uh2*6.0+nl1*param1*u6*uh1*uh5*6.0+nl1*param1*u8*uh1*uh3*6.0-nl1*param4*u7*uh1*uh2*8.0+nl1*param1*u9*uh1*uh4*6.0-nl1*param4*u8*uh1*uh3*6.0+nl2*param4*u6*uh2*uh3*2.0+nl2*param4*u7*uh1*uh3*4.0-nl2*param4*u8*uh1*uh2*6.0+nl2*param1*u12*uh1*uh2*6.0-nl1*param4*u9*uh1*uh4*6.0+nl3*param4*u6*uh2*uh4*2.0+nl3*param4*u7*uh1*uh4*4.0-nl3*param4*u9*uh1*uh2*6.0+nl2*param1*u11*uh1*uh5*6.0+nl2*param1*u13*uh1*uh3*6.0+nl1*param4*u11*uh2*uh3*2.0-nl1*param4*u12*uh1*uh3*6.0+nl1*param4*u13*uh1*uh2*4.0-nl2*param4*u12*uh1*uh2*6.0+nl2*param1*u14*uh1*uh4*6.0-nl2*param4*u13*uh1*uh3*8.0+nl3*param1*u17*uh1*uh2*6.0-nl2*param4*u14*uh1*uh4*6.0+nl3*param4*u11*uh3*uh4*2.0+nl3*param4*u13*uh1*uh4*4.0-nl3*param4*u14*uh1*uh3*6.0+nl3*param1*u16*uh1*uh5*6.0+nl3*param1*u18*uh1*uh3*6.0+nl1*param4*u16*uh2*uh4*2.0-nl1*param4*u17*uh1*uh4*6.0+nl1*param4*u19*uh1*uh2*4.0-nl3*param4*u17*uh1*uh2*6.0+nl3*param1*u19*uh1*uh4*6.0+nl2*param4*u16*uh3*uh4*2.0-nl2*param4*u18*uh1*uh4*6.0+nl2*param4*u19*uh1*uh3*4.0-nl3*param4*u18*uh1*uh3*6.0-nl3*param4*u19*uh1*uh4*8.0-nl1*param3*param4*t3*uh1*uh2*3.0-nl2*param3*param4*t3*uh1*uh3*3.0-nl1*param3*param4*t6*uh1*uh2*3.0-nl3*param3*param4*t3*uh1*uh4*3.0-nl2*param3*param4*t6*uh1*uh3*3.0-nl1*param3*param4*t9*uh1*uh2*3.0-nl3*param3*param4*t6*uh1*uh4*3.0-nl2*param3*param4*t9*uh1*uh3*3.0-nl3*param3*param4*t9*uh1*uh4*3.0-param3*param4*param6*t27*u5*uh1*6.0+param3*param4*param6*t27*uh1*uh5*6.0+nl1*param1*param3*param4*t3*uh1*uh2*3.0+nl2*param1*param3*param4*t3*uh1*uh3*3.0+nl1*param1*param3*param4*t6*uh1*uh2*3.0+nl3*param1*param3*param4*t3*uh1*uh4*3.0+nl2*param1*param3*param4*t6*uh1*uh3*3.0+nl1*param1*param3*param4*t9*uh1*uh2*3.0+nl3*param1*param3*param4*t6*uh1*uh4*3.0+nl2*param1*param3*param4*t9*uh1*uh3*3.0+nl3*param1*param3*param4*t9*uh1*uh4*3.0-nl1*param1*param3*param4*t27*uh2*uh5*6.0-nl2*param1*param3*param4*t27*uh3*uh5*6.0-nl3*param1*param3*param4*t27*uh4*uh5*6.0)*(-1.0/6.0))/param4;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = 1.0/param3;
		double t3 = 1.0/(uh1*uh1);
		double t4 = uh2*uh2;
		double t5 = uh3*uh3;
		double t6 = uh4*uh4;
		double t7 = 1.0/uh1;
		double t8 = 1.0/param4;
		double t9 = 1.0/(uh1*uh1*uh1);
		double t10 = nl1*t2*t7;
		double t11 = nl2*t2*t7;
		double t12 = nl1*param4*uh1*uh3*6.0;
		double t13 = nl2*param4*uh1*uh2*6.0;
		double t14 = nl3*param4*uh1*uh4*4.0;
		double t15 = nl3*t2*t7;
		double t16 = nl1*param4*uh1*uh4*6.0;
		double t17 = nl3*param4*uh1*uh2*6.0;
		double t18 = nl2*param4*uh1*uh4*6.0;
		double t19 = nl3*param4*uh1*uh3*6.0;
		double t20 = nl1*param4*uh1*uh2*4.0;
		double t21 = nl2*param4*uh1*uh3*4.0;
		fh_udg[0*ng+i] = param6;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = 0.0;
		fh_udg[5*ng+i] = 0.0;
		fh_udg[6*ng+i] = param6;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = 0.0;
		fh_udg[9*ng+i] = 0.0;
		fh_udg[10*ng+i] = 0.0;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = param6;
		fh_udg[13*ng+i] = 0.0;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = 0.0;
		fh_udg[16*ng+i] = 0.0;
		fh_udg[17*ng+i] = 0.0;
		fh_udg[18*ng+i] = param6;
		fh_udg[19*ng+i] = 0.0;
		fh_udg[20*ng+i] = 0.0;
		fh_udg[21*ng+i] = 0.0;
		fh_udg[22*ng+i] = 0.0;
		fh_udg[23*ng+i] = 0.0;
		fh_udg[24*ng+i] = param6;
		fh_udg[25*ng+i] = 0.0;
		fh_udg[26*ng+i] = nl1*t2*t3*uh2*(-4.0/3.0)-nl2*t2*t3*uh3-nl3*t2*t3*uh4;
		fh_udg[27*ng+i] = -nl1*t2*t3*uh3+nl2*t2*t3*uh2*(2.0/3.0);
		fh_udg[28*ng+i] = -nl1*t2*t3*uh4+nl3*t2*t3*uh2*(2.0/3.0);
		fh_udg[29*ng+i] = t2*t8*t9*(nl1*param1*t4*-6.0-nl1*param1*t5*6.0-nl1*param1*t6*6.0+nl1*param4*t4*8.0+nl1*param4*t5*6.0+nl1*param4*t6*6.0+nl1*param1*uh1*uh5*6.0+nl2*param4*uh2*uh3*2.0+nl3*param4*uh2*uh4*2.0)*(-1.0/6.0);
		fh_udg[30*ng+i] = 0.0;
		fh_udg[31*ng+i] = nl1*t2*t7*(4.0/3.0);
		fh_udg[32*ng+i] = nl2*t2*t7*(-2.0/3.0);
		fh_udg[33*ng+i] = nl3*t2*t7*(-2.0/3.0);
		fh_udg[34*ng+i] = t2*t8*t9*(t14+t21+nl1*param1*uh1*uh2*6.0-nl1*param4*uh1*uh2*8.0)*(-1.0/6.0);
		fh_udg[35*ng+i] = 0.0;
		fh_udg[36*ng+i] = t11;
		fh_udg[37*ng+i] = t10;
		fh_udg[38*ng+i] = 0.0;
		fh_udg[39*ng+i] = t2*t8*t9*(t12+t13-nl1*param1*uh1*uh3*6.0)*(1.0/6.0);
		fh_udg[40*ng+i] = 0.0;
		fh_udg[41*ng+i] = t15;
		fh_udg[42*ng+i] = 0.0;
		fh_udg[43*ng+i] = t10;
		fh_udg[44*ng+i] = t2*t8*t9*(t16+t17-nl1*param1*uh1*uh4*6.0)*(1.0/6.0);
		fh_udg[45*ng+i] = 0.0;
		fh_udg[46*ng+i] = 0.0;
		fh_udg[47*ng+i] = 0.0;
		fh_udg[48*ng+i] = 0.0;
		fh_udg[49*ng+i] = nl1*param1*t2*t7*t8;
		fh_udg[50*ng+i] = 0.0;
		fh_udg[51*ng+i] = nl1*t2*t3*uh3*(2.0/3.0)-nl2*t2*t3*uh2;
		fh_udg[52*ng+i] = -nl1*t2*t3*uh2-nl2*t2*t3*uh3*(4.0/3.0)-nl3*t2*t3*uh4;
		fh_udg[53*ng+i] = -nl2*t2*t3*uh4+nl3*t2*t3*uh3*(2.0/3.0);
		fh_udg[54*ng+i] = t2*t8*t9*(nl2*param1*t4*-6.0-nl2*param1*t5*6.0-nl2*param1*t6*6.0+nl2*param4*t4*6.0+nl2*param4*t5*8.0+nl2*param4*t6*6.0+nl2*param1*uh1*uh5*6.0+nl1*param4*uh2*uh3*2.0+nl3*param4*uh3*uh4*2.0)*(-1.0/6.0);
		fh_udg[55*ng+i] = 0.0;
		fh_udg[56*ng+i] = t11;
		fh_udg[57*ng+i] = t10;
		fh_udg[58*ng+i] = 0.0;
		fh_udg[59*ng+i] = t2*t8*t9*(t12+t13-nl2*param1*uh1*uh2*6.0)*(1.0/6.0);
		fh_udg[60*ng+i] = 0.0;
		fh_udg[61*ng+i] = nl1*t2*t7*(-2.0/3.0);
		fh_udg[62*ng+i] = nl2*t2*t7*(4.0/3.0);
		fh_udg[63*ng+i] = nl3*t2*t7*(-2.0/3.0);
		fh_udg[64*ng+i] = t2*t8*t9*(t14+t20+nl2*param1*uh1*uh3*6.0-nl2*param4*uh1*uh3*8.0)*(-1.0/6.0);
		fh_udg[65*ng+i] = 0.0;
		fh_udg[66*ng+i] = 0.0;
		fh_udg[67*ng+i] = t15;
		fh_udg[68*ng+i] = t11;
		fh_udg[69*ng+i] = t2*t8*t9*(t18+t19-nl2*param1*uh1*uh4*6.0)*(1.0/6.0);
		fh_udg[70*ng+i] = 0.0;
		fh_udg[71*ng+i] = 0.0;
		fh_udg[72*ng+i] = 0.0;
		fh_udg[73*ng+i] = 0.0;
		fh_udg[74*ng+i] = nl2*param1*t2*t7*t8;
		fh_udg[75*ng+i] = 0.0;
		fh_udg[76*ng+i] = nl1*t2*t3*uh4*(2.0/3.0)-nl3*t2*t3*uh2;
		fh_udg[77*ng+i] = nl2*t2*t3*uh4*(2.0/3.0)-nl3*t2*t3*uh3;
		fh_udg[78*ng+i] = -nl1*t2*t3*uh2-nl2*t2*t3*uh3-nl3*t2*t3*uh4*(4.0/3.0);
		fh_udg[79*ng+i] = t2*t8*t9*(nl3*param1*t4*-6.0-nl3*param1*t5*6.0-nl3*param1*t6*6.0+nl3*param4*t4*6.0+nl3*param4*t5*6.0+nl3*param4*t6*8.0+nl3*param1*uh1*uh5*6.0+nl1*param4*uh2*uh4*2.0+nl2*param4*uh3*uh4*2.0)*(-1.0/6.0);
		fh_udg[80*ng+i] = 0.0;
		fh_udg[81*ng+i] = t15;
		fh_udg[82*ng+i] = 0.0;
		fh_udg[83*ng+i] = t10;
		fh_udg[84*ng+i] = t2*t8*t9*(t16+t17-nl3*param1*uh1*uh2*6.0)*(1.0/6.0);
		fh_udg[85*ng+i] = 0.0;
		fh_udg[86*ng+i] = 0.0;
		fh_udg[87*ng+i] = t15;
		fh_udg[88*ng+i] = t11;
		fh_udg[89*ng+i] = t2*t8*t9*(t18+t19-nl3*param1*uh1*uh3*6.0)*(1.0/6.0);
		fh_udg[90*ng+i] = 0.0;
		fh_udg[91*ng+i] = nl1*t2*t7*(-2.0/3.0);
		fh_udg[92*ng+i] = nl2*t2*t7*(-2.0/3.0);
		fh_udg[93*ng+i] = nl3*t2*t7*(4.0/3.0);
		fh_udg[94*ng+i] = t2*t8*t9*(t20+t21+nl3*param1*uh1*uh4*6.0-nl3*param4*uh1*uh4*8.0)*(-1.0/6.0);
		fh_udg[95*ng+i] = 0.0;
		fh_udg[96*ng+i] = 0.0;
		fh_udg[97*ng+i] = 0.0;
		fh_udg[98*ng+i] = 0.0;
		fh_udg[99*ng+i] = nl3*param1*t2*t7*t8;

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = 1.0/(uh1*uh1);
		double t3 = uh2*uh2;
		double t4 = param1-1.0;
		double t5 = 1.0/param3;
		double t6 = 1.0/(uh1*uh1*uh1);
		double t7 = uh3*uh3;
		double t8 = u19*uh1*2.0;
		double t9 = uh4*uh4;
		double t23 = uh1*uh5*2.0;
		double t10 = t3+t7+t9-t23;
		double t11 = t2*t4*t10*(1.0/2.0);
		double t12 = 1.0/uh1;
		double t13 = t4*t12*uh5;
		double t14 = u19*2.0;
		double t15 = u8*uh1;
		double t16 = u12*uh1;
		double t17 = param3*uh1*uh2*uh3;
		double t18 = t15+t16+t17-u6*uh3-u11*uh2;
		double t19 = param3*uh2*uh3;
		double t20 = t19+u8+u12;
		double t21 = u7*uh1*2.0;
		double t22 = u13*uh1*2.0;
		double t24 = u7*2.0;
		double t25 = u13*2.0;
		double t26 = u9*uh1;
		double t27 = u17*uh1;
		double t28 = param3*uh1*uh2*uh4;
		double t29 = t26+t27+t28-u6*uh4-u16*uh2;
		double t30 = u14*uh1;
		double t31 = u18*uh1;
		double t32 = param3*uh1*uh3*uh4;
		double t33 = t30+t31+t32-u11*uh4-u16*uh3;
		double t34 = param3*uh2*uh4;
		double t35 = t34+u9+u17;
		double t36 = param3*uh3*uh4;
		double t37 = t36+u14+u18;
		double t38 = uh1*uh1;
		double t39 = 1.0/param4;
		double t40 = t4*t12*uh2;
		double t51 = param3*uh1*uh3;
		double t41 = -t51+u11;
		double t42 = t40-t2*t5*u6*(2.0/3.0);
		double t45 = param3*uh1*uh4;
		double t43 = -t45+u16;
		double t48 = param3*uh1*uh2;
		double t44 = -t48+u6;
		double t46 = t2*t5*u11*(2.0/3.0);
		double t47 = t4*t12*uh3;
		double t49 = t2*t5*u16*(2.0/3.0);
		double t52 = t4*t12*uh4;
		double t50 = t49-t52;
		fh_uh[0*ng+i] = -param6;
		fh_uh[1*ng+i] = nl1*(t11+t13-t2*t3-t2*t5*(t14+t25-u7*4.0)*(1.0/3.0)+t5*t6*(t8+t22+u6*uh2*4.0-u7*uh1*4.0-u11*uh3*2.0-u16*uh4*2.0)*(2.0/3.0))+nl2*t2*t5*t20-nl2*t5*t6*t18*2.0-nl3*t5*t6*t29*2.0+nl3*t2*t5*t35;
		fh_uh[2*ng+i] = nl2*(t11+t13-t2*t7-t2*t5*(t14+t24-u13*4.0)*(1.0/3.0)+t5*t6*(t8+t21-u6*uh2*2.0+u11*uh3*4.0-u13*uh1*4.0-u16*uh4*2.0)*(2.0/3.0))+nl1*t2*t5*t20-nl1*t5*t6*t18*2.0+nl3*t2*t5*t37-nl3*t5*t6*t33*2.0;
		fh_uh[3*ng+i] = nl3*(t11+t13-t2*t9-t2*t5*(t24+t25-u19*4.0)*(1.0/3.0)+t5*t6*(t21+t22-u6*uh2*2.0-u11*uh3*2.0+u16*uh4*4.0-u19*uh1*4.0)*(2.0/3.0))-nl1*t5*t6*t29*2.0+nl1*t2*t5*t35+nl2*t2*t5*t37-nl2*t5*t6*t33*2.0;
		fh_uh[4*ng+i] = t5*t6*t39*(nl1*param1*u7*uh2*-6.0-nl1*param1*u6*uh5*6.0-nl1*param1*u8*uh3*6.0+nl1*param1*u10*uh1*1.2E1+nl1*param4*u7*uh2*8.0-nl1*param1*u9*uh4*6.0+nl1*param4*u8*uh3*6.0-nl2*param4*u7*uh3*4.0+nl2*param4*u8*uh2*6.0-nl2*param1*u12*uh2*6.0+nl1*param4*u9*uh4*6.0-nl3*param4*u7*uh4*4.0+nl3*param4*u9*uh2*6.0-nl2*param1*u11*uh5*6.0-nl2*param1*u13*uh3*6.0+nl2*param1*u15*uh1*1.2E1+nl1*param4*u12*uh3*6.0-nl1*param4*u13*uh2*4.0+nl2*param4*u12*uh2*6.0-nl2*param1*u14*uh4*6.0+nl2*param4*u13*uh3*8.0-nl3*param1*u17*uh2*6.0+nl2*param4*u14*uh4*6.0-nl3*param4*u13*uh4*4.0+nl3*param4*u14*uh3*6.0-nl3*param1*u16*uh5*6.0-nl3*param1*u18*uh3*6.0+nl3*param1*u20*uh1*1.2E1+nl1*param4*u17*uh4*6.0-nl1*param4*u19*uh2*4.0+nl3*param4*u17*uh2*6.0-nl3*param1*u19*uh4*6.0+nl2*param4*u18*uh4*6.0-nl2*param4*u19*uh3*4.0+nl3*param4*u18*uh3*6.0+nl3*param4*u19*uh4*8.0+nl1*param3*param4*t3*uh2*3.0+nl2*param3*param4*t3*uh3*3.0+nl1*param3*param4*t7*uh2*3.0+nl3*param3*param4*t3*uh4*3.0+nl1*param3*param4*t9*uh2*3.0+nl2*param3*param4*t7*uh3*3.0+nl2*param3*param4*t9*uh3*3.0+nl3*param3*param4*t7*uh4*3.0+nl3*param3*param4*t9*uh4*3.0+param3*param4*param6*t38*u5*1.8E1-param3*param4*param6*t38*uh5*1.8E1-nl1*param1*param3*param4*t3*uh2*3.0-nl2*param1*param3*param4*t3*uh3*3.0-nl1*param1*param3*param4*t7*uh2*3.0-nl3*param1*param3*param4*t3*uh4*3.0-nl1*param1*param3*param4*t9*uh2*3.0-nl2*param1*param3*param4*t7*uh3*3.0-nl2*param1*param3*param4*t9*uh3*3.0-nl3*param1*param3*param4*t7*uh4*3.0-nl3*param1*param3*param4*t9*uh4*3.0+nl1*param1*param3*param4*uh1*uh2*uh5*1.2E1+nl2*param1*param3*param4*uh1*uh3*uh5*1.2E1+nl3*param1*param3*param4*uh1*uh4*uh5*1.2E1)*(1.0/6.0)+t5*t39*1.0/(uh1*uh1*uh1*uh1)*(nl1*param1*t3*u6*-6.0+nl1*param4*t3*u6*8.0-nl1*param1*t7*u6*6.0-nl1*param1*t9*u6*6.0-nl2*param1*t3*u11*6.0+nl1*param4*t7*u6*6.0+nl1*param4*t9*u6*6.0+nl2*param4*t3*u11*6.0-nl2*param1*t7*u11*6.0-nl2*param1*t9*u11*6.0-nl3*param1*t3*u16*6.0+nl2*param4*t7*u11*8.0+nl2*param4*t9*u11*6.0+nl3*param4*t3*u16*6.0-nl3*param1*t7*u16*6.0-nl3*param1*t9*u16*6.0+nl3*param4*t7*u16*6.0+nl3*param4*t9*u16*8.0-nl1*param1*t38*u10*6.0-nl2*param1*t38*u15*6.0-nl3*param1*t38*u20*6.0+nl1*param1*u7*uh1*uh2*6.0+nl1*param1*u6*uh1*uh5*6.0+nl1*param1*u8*uh1*uh3*6.0-nl1*param4*u7*uh1*uh2*8.0+nl1*param1*u9*uh1*uh4*6.0-nl1*param4*u8*uh1*uh3*6.0+nl2*param4*u6*uh2*uh3*2.0+nl2*param4*u7*uh1*uh3*4.0-nl2*param4*u8*uh1*uh2*6.0+nl2*param1*u12*uh1*uh2*6.0-nl1*param4*u9*uh1*uh4*6.0+nl3*param4*u6*uh2*uh4*2.0+nl3*param4*u7*uh1*uh4*4.0-nl3*param4*u9*uh1*uh2*6.0+nl2*param1*u11*uh1*uh5*6.0+nl2*param1*u13*uh1*uh3*6.0+nl1*param4*u11*uh2*uh3*2.0-nl1*param4*u12*uh1*uh3*6.0+nl1*param4*u13*uh1*uh2*4.0-nl2*param4*u12*uh1*uh2*6.0+nl2*param1*u14*uh1*uh4*6.0-nl2*param4*u13*uh1*uh3*8.0+nl3*param1*u17*uh1*uh2*6.0-nl2*param4*u14*uh1*uh4*6.0+nl3*param4*u11*uh3*uh4*2.0+nl3*param4*u13*uh1*uh4*4.0-nl3*param4*u14*uh1*uh3*6.0+nl3*param1*u16*uh1*uh5*6.0+nl3*param1*u18*uh1*uh3*6.0+nl1*param4*u16*uh2*uh4*2.0-nl1*param4*u17*uh1*uh4*6.0+nl1*param4*u19*uh1*uh2*4.0-nl3*param4*u17*uh1*uh2*6.0+nl3*param1*u19*uh1*uh4*6.0+nl2*param4*u16*uh3*uh4*2.0-nl2*param4*u18*uh1*uh4*6.0+nl2*param4*u19*uh1*uh3*4.0-nl3*param4*u18*uh1*uh3*6.0-nl3*param4*u19*uh1*uh4*8.0-nl1*param3*param4*t3*uh1*uh2*3.0-nl2*param3*param4*t3*uh1*uh3*3.0-nl1*param3*param4*t7*uh1*uh2*3.0-nl3*param3*param4*t3*uh1*uh4*3.0-nl1*param3*param4*t9*uh1*uh2*3.0-nl2*param3*param4*t7*uh1*uh3*3.0-nl2*param3*param4*t9*uh1*uh3*3.0-nl3*param3*param4*t7*uh1*uh4*3.0-nl3*param3*param4*t9*uh1*uh4*3.0-param3*param4*param6*t38*u5*uh1*6.0+param3*param4*param6*t38*uh1*uh5*6.0+nl1*param1*param3*param4*t3*uh1*uh2*3.0+nl2*param1*param3*param4*t3*uh1*uh3*3.0+nl1*param1*param3*param4*t7*uh1*uh2*3.0+nl3*param1*param3*param4*t3*uh1*uh4*3.0+nl1*param1*param3*param4*t9*uh1*uh2*3.0+nl2*param1*param3*param4*t7*uh1*uh3*3.0+nl2*param1*param3*param4*t9*uh1*uh3*3.0+nl3*param1*param3*param4*t7*uh1*uh4*3.0+nl3*param1*param3*param4*t9*uh1*uh4*3.0-nl1*param1*param3*param4*t38*uh2*uh5*6.0-nl2*param1*param3*param4*t38*uh3*uh5*6.0-nl3*param1*param3*param4*t38*uh4*uh5*6.0)*(1.0/2.0);
		fh_uh[5*ng+i] = nl1;
		fh_uh[6*ng+i] = -param6-nl1*(t40-t12*uh2*2.0+t2*t5*u6*(4.0/3.0))-nl2*t2*t5*t41-nl3*t2*t5*t43;
		fh_uh[7*ng+i] = -nl2*t42-nl1*t2*t5*t41;
		fh_uh[8*ng+i] = -nl3*t42-nl1*t2*t5*t43;
		fh_uh[9*ng+i] = t5*t6*t39*(nl1*param1*u6*uh2*-1.2E1+nl1*param1*u7*uh1*6.0+nl1*param4*u6*uh2*1.6E1-nl1*param4*u7*uh1*8.0+nl2*param4*u6*uh3*2.0-nl2*param4*u8*uh1*6.0-nl2*param1*u11*uh2*1.2E1+nl2*param1*u12*uh1*6.0+nl3*param4*u6*uh4*2.0-nl3*param4*u9*uh1*6.0+nl1*param4*u11*uh3*2.0+nl1*param4*u13*uh1*4.0+nl2*param4*u11*uh2*1.2E1-nl2*param4*u12*uh1*6.0-nl3*param1*u16*uh2*1.2E1+nl3*param1*u17*uh1*6.0+nl1*param4*u16*uh4*2.0+nl1*param4*u19*uh1*4.0+nl3*param4*u16*uh2*1.2E1-nl3*param4*u17*uh1*6.0-nl1*param3*param4*t3*uh1*9.0-nl1*param3*param4*t7*uh1*3.0-nl1*param3*param4*t9*uh1*3.0+nl1*param1*param3*param4*t3*uh1*9.0+nl1*param1*param3*param4*t7*uh1*3.0+nl1*param1*param3*param4*t9*uh1*3.0-nl1*param1*param3*param4*t38*uh5*6.0-nl2*param3*param4*uh1*uh2*uh3*6.0-nl3*param3*param4*uh1*uh2*uh4*6.0+nl2*param1*param3*param4*uh1*uh2*uh3*6.0+nl3*param1*param3*param4*uh1*uh2*uh4*6.0)*(-1.0/6.0);
		fh_uh[10*ng+i] = nl2;
		fh_uh[11*ng+i] = nl1*(t46-t4*t12*uh3)-nl2*t2*t5*t44;
		fh_uh[12*ng+i] = -param6-nl2*(t47-t12*uh3*2.0+t2*t5*u11*(4.0/3.0))-nl1*t2*t5*t44-nl3*t2*t5*t43;
		fh_uh[13*ng+i] = nl3*(t46-t47)-nl2*t2*t5*t43;
		fh_uh[14*ng+i] = t5*t6*t39*(nl1*param1*u6*uh3*-1.2E1+nl1*param1*u8*uh1*6.0+nl1*param4*u6*uh3*1.2E1-nl1*param4*u8*uh1*6.0+nl2*param4*u6*uh2*2.0+nl2*param4*u7*uh1*4.0-nl2*param1*u11*uh3*1.2E1+nl2*param1*u13*uh1*6.0+nl1*param4*u11*uh2*2.0-nl1*param4*u12*uh1*6.0+nl2*param4*u11*uh3*1.6E1-nl2*param4*u13*uh1*8.0+nl3*param4*u11*uh4*2.0-nl3*param4*u14*uh1*6.0-nl3*param1*u16*uh3*1.2E1+nl3*param1*u18*uh1*6.0+nl2*param4*u16*uh4*2.0+nl2*param4*u19*uh1*4.0+nl3*param4*u16*uh3*1.2E1-nl3*param4*u18*uh1*6.0-nl2*param3*param4*t3*uh1*3.0-nl2*param3*param4*t7*uh1*9.0-nl2*param3*param4*t9*uh1*3.0+nl2*param1*param3*param4*t3*uh1*3.0+nl2*param1*param3*param4*t7*uh1*9.0+nl2*param1*param3*param4*t9*uh1*3.0-nl2*param1*param3*param4*t38*uh5*6.0-nl1*param3*param4*uh1*uh2*uh3*6.0-nl3*param3*param4*uh1*uh3*uh4*6.0+nl1*param1*param3*param4*uh1*uh2*uh3*6.0+nl3*param1*param3*param4*uh1*uh3*uh4*6.0)*(-1.0/6.0);
		fh_uh[15*ng+i] = nl3;
		fh_uh[16*ng+i] = nl1*t50-nl3*t2*t5*t44;
		fh_uh[17*ng+i] = nl2*t50-nl3*t2*t5*t41;
		fh_uh[18*ng+i] = -param6-nl3*(t52-t12*uh4*2.0+t2*t5*u16*(4.0/3.0))-nl2*t2*t5*t41-nl1*t2*t5*t44;
		fh_uh[19*ng+i] = t5*t6*t39*(nl1*param1*u6*uh4*-1.2E1+nl1*param1*u9*uh1*6.0+nl1*param4*u6*uh4*1.2E1-nl1*param4*u9*uh1*6.0+nl3*param4*u6*uh2*2.0+nl3*param4*u7*uh1*4.0-nl2*param1*u11*uh4*1.2E1+nl2*param1*u14*uh1*6.0+nl2*param4*u11*uh4*1.2E1-nl2*param4*u14*uh1*6.0+nl3*param4*u11*uh3*2.0+nl3*param4*u13*uh1*4.0+nl1*param4*u16*uh2*2.0-nl1*param4*u17*uh1*6.0-nl3*param1*u16*uh4*1.2E1+nl3*param1*u19*uh1*6.0+nl2*param4*u16*uh3*2.0-nl2*param4*u18*uh1*6.0+nl3*param4*u16*uh4*1.6E1-nl3*param4*u19*uh1*8.0-nl3*param3*param4*t3*uh1*3.0-nl3*param3*param4*t7*uh1*3.0-nl3*param3*param4*t9*uh1*9.0+nl3*param1*param3*param4*t3*uh1*3.0+nl3*param1*param3*param4*t7*uh1*3.0+nl3*param1*param3*param4*t9*uh1*9.0-nl3*param1*param3*param4*t38*uh5*6.0-nl1*param3*param4*uh1*uh2*uh4*6.0-nl2*param3*param4*uh1*uh3*uh4*6.0+nl1*param1*param3*param4*uh1*uh2*uh4*6.0+nl2*param1*param3*param4*uh1*uh3*uh4*6.0)*(-1.0/6.0);
		fh_uh[20*ng+i] = 0.0;
		fh_uh[21*ng+i] = nl1*t4;
		fh_uh[22*ng+i] = nl2*t4;
		fh_uh[23*ng+i] = nl3*t4;
		fh_uh[24*ng+i] = t5*t6*t39*(nl1*param1*u6*uh1*6.0+nl2*param1*u11*uh1*6.0+nl3*param1*u16*uh1*6.0+param3*param4*param6*t38*uh1*6.0-nl1*param1*param3*param4*t38*uh2*6.0-nl2*param1*param3*param4*t38*uh3*6.0-nl3*param1*param3*param4*t38*uh4*6.0)*(-1.0/6.0);

	}
}


void fhatonly_ns3d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double uh4 = uh[3*ng+i];
		double uh5 = uh[4*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];

		double t2 = 1.0/uh1;
		double t3 = uh2*uh2;
		double t4 = 1.0/param3;
		double t5 = 1.0/(uh1*uh1);
		double t6 = uh3*uh3;
		double t7 = u19*uh1*2.0;
		double t8 = param1-1.0;
		double t9 = uh4*uh4;
		double t18 = uh1*uh5*2.0;
		double t10 = t3+t6+t9-t18;
		double t11 = t2*t8*t10*(1.0/2.0);
		double t12 = u8*uh1;
		double t13 = u12*uh1;
		double t14 = param3*uh1*uh2*uh3;
		double t15 = t12+t13+t14-u6*uh3-u11*uh2;
		double t16 = u7*uh1*2.0;
		double t17 = u13*uh1*2.0;
		double t19 = u9*uh1;
		double t20 = u17*uh1;
		double t21 = param3*uh1*uh2*uh4;
		double t22 = t19+t20+t21-u6*uh4-u16*uh2;
		double t23 = u14*uh1;
		double t24 = u18*uh1;
		double t25 = param3*uh1*uh3*uh4;
		double t26 = t23+t24+t25-u11*uh4-u16*uh3;
		double t27 = uh1*uh1;
		fh[0*ng+i] = nl1*uh2+nl2*uh3+nl3*uh4+param6*(u1-uh1);
		fh[1*ng+i] = param6*(u2-uh2)-nl1*(t11-t2*t3+t4*t5*(t7+t17+u6*uh2*4.0-u7*uh1*4.0-u11*uh3*2.0-u16*uh4*2.0)*(1.0/3.0))+nl2*t4*t5*t15+nl3*t4*t5*t22;
		fh[2*ng+i] = param6*(u3-uh3)-nl2*(t11-t2*t6+t4*t5*(t7+t16-u6*uh2*2.0+u11*uh3*4.0-u13*uh1*4.0-u16*uh4*2.0)*(1.0/3.0))+nl1*t4*t5*t15+nl3*t4*t5*t26;
		fh[3*ng+i] = param6*(u4-uh4)-nl3*(t11-t2*t9+t4*t5*(t16+t17-u6*uh2*2.0-u11*uh3*2.0+u16*uh4*4.0-u19*uh1*4.0)*(1.0/3.0))+nl1*t4*t5*t22+nl2*t4*t5*t26;
		fh[4*ng+i] = (t4*1.0/(uh1*uh1*uh1)*(nl1*param1*t3*u6*-6.0-nl1*param1*t6*u6*6.0+nl1*param4*t3*u6*8.0-nl1*param1*t9*u6*6.0+nl1*param4*t6*u6*6.0-nl2*param1*t3*u11*6.0+nl1*param4*t9*u6*6.0-nl2*param1*t6*u11*6.0+nl2*param4*t3*u11*6.0-nl2*param1*t9*u11*6.0+nl2*param4*t6*u11*8.0-nl3*param1*t3*u16*6.0+nl2*param4*t9*u11*6.0-nl3*param1*t6*u16*6.0+nl3*param4*t3*u16*6.0-nl3*param1*t9*u16*6.0+nl3*param4*t6*u16*6.0+nl3*param4*t9*u16*8.0-nl1*param1*t27*u10*6.0-nl2*param1*t27*u15*6.0-nl3*param1*t27*u20*6.0+nl1*param1*u7*uh1*uh2*6.0+nl1*param1*u6*uh1*uh5*6.0+nl1*param1*u8*uh1*uh3*6.0-nl1*param4*u7*uh1*uh2*8.0+nl1*param1*u9*uh1*uh4*6.0-nl1*param4*u8*uh1*uh3*6.0+nl2*param4*u6*uh2*uh3*2.0+nl2*param4*u7*uh1*uh3*4.0-nl2*param4*u8*uh1*uh2*6.0+nl2*param1*u12*uh1*uh2*6.0-nl1*param4*u9*uh1*uh4*6.0+nl3*param4*u6*uh2*uh4*2.0+nl3*param4*u7*uh1*uh4*4.0-nl3*param4*u9*uh1*uh2*6.0+nl2*param1*u11*uh1*uh5*6.0+nl2*param1*u13*uh1*uh3*6.0+nl1*param4*u11*uh2*uh3*2.0-nl1*param4*u12*uh1*uh3*6.0+nl1*param4*u13*uh1*uh2*4.0-nl2*param4*u12*uh1*uh2*6.0+nl2*param1*u14*uh1*uh4*6.0-nl2*param4*u13*uh1*uh3*8.0+nl3*param1*u17*uh1*uh2*6.0-nl2*param4*u14*uh1*uh4*6.0+nl3*param4*u11*uh3*uh4*2.0+nl3*param4*u13*uh1*uh4*4.0-nl3*param4*u14*uh1*uh3*6.0+nl3*param1*u16*uh1*uh5*6.0+nl3*param1*u18*uh1*uh3*6.0+nl1*param4*u16*uh2*uh4*2.0-nl1*param4*u17*uh1*uh4*6.0+nl1*param4*u19*uh1*uh2*4.0-nl3*param4*u17*uh1*uh2*6.0+nl3*param1*u19*uh1*uh4*6.0+nl2*param4*u16*uh3*uh4*2.0-nl2*param4*u18*uh1*uh4*6.0+nl2*param4*u19*uh1*uh3*4.0-nl3*param4*u18*uh1*uh3*6.0-nl3*param4*u19*uh1*uh4*8.0-nl1*param3*param4*t3*uh1*uh2*3.0-nl2*param3*param4*t3*uh1*uh3*3.0-nl1*param3*param4*t6*uh1*uh2*3.0-nl3*param3*param4*t3*uh1*uh4*3.0-nl2*param3*param4*t6*uh1*uh3*3.0-nl1*param3*param4*t9*uh1*uh2*3.0-nl3*param3*param4*t6*uh1*uh4*3.0-nl2*param3*param4*t9*uh1*uh3*3.0-nl3*param3*param4*t9*uh1*uh4*3.0-param3*param4*param6*t27*u5*uh1*6.0+param3*param4*param6*t27*uh1*uh5*6.0+nl1*param1*param3*param4*t3*uh1*uh2*3.0+nl2*param1*param3*param4*t3*uh1*uh3*3.0+nl1*param1*param3*param4*t6*uh1*uh2*3.0+nl3*param1*param3*param4*t3*uh1*uh4*3.0+nl2*param1*param3*param4*t6*uh1*uh3*3.0+nl1*param1*param3*param4*t9*uh1*uh2*3.0+nl3*param1*param3*param4*t6*uh1*uh4*3.0+nl2*param1*param3*param4*t9*uh1*uh3*3.0+nl3*param1*param3*param4*t9*uh1*uh4*3.0-nl1*param1*param3*param4*t27*uh2*uh5*6.0-nl2*param1*param3*param4*t27*uh3*uh5*6.0-nl3*param1*param3*param4*t27*uh4*uh5*6.0)*(-1.0/6.0))/param4;

	}
}

