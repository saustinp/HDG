#include <math.h>

void source_streamer2d(double *s, double *s_udg, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
{
	double param1 = param[0];
	double param2 = param[1];
	double param3 = param[2];
	double param4 = param[3];
	double param5 = param[4];
	double param6 = param[5];
	double param7 = param[6];
	double param8 = param[7];
	double param9 = param[8];
	double param10 = param[9];

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/param2;
		double t5 = 1.0/param3;
		double t8 = param1*3.4075E+2;
		double t6 = t5*t5*t5;
		double t7 = t2+t3;
		double t9 = sqrt(t7);
		double t10 = 1.0/t9;
		double t12 = param3*t9;
		double t11 = t10*t10*t10;
		double t13 = 1.0/pow(t12,1.3E+1/5.0E+1);
		double t14 = t5*t10*2.73E+7;
		double t15 = -t14;
		double t17 = t6*t11*4.3666E+26;
		double t16 = exp(t15);
		double t18 = t17+1.1944E+6;
		double t19 = param1*t16*t18;
		double t20 = -t19;
		double t21 = t8+t20;
		double t22 = t4*t9*t13*t21*u1*x1*2.3987;
		double t23 = -t22;
		s[0*ng+i] = t23;
		s[1*ng+i] = t23;
		s[2*ng+i] = -1.0/(param1*param1)*param4*t5*x1*(u1-u2);

	}

	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double u4 = udg[3*ng+i];
		double u5 = udg[4*ng+i];
		double u6 = udg[5*ng+i];
		double u7 = udg[6*ng+i];
		double u8 = udg[7*ng+i];
		double u9 = udg[8*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/(param1*param1);
		double t5 = 1.0/param2;
		double t6 = 1.0/param3;
		double t9 = param1*3.4075E+2;
		double t7 = t6*t6*t6;
		double t8 = t2+t3;
		double t11 = param4*t4*t6*x1;
		double t10 = sqrt(t8);
		double t12 = 1.0/t10;
		double t15 = param3*t10;
		double t13 = t12*t12*t12;
		double t14 = t12*t12*t12*t12*t12;
		double t16 = 1.0/pow(t15,1.3E+1/5.0E+1);
		double t17 = 1.0/pow(t15,6.3E+1/5.0E+1);
		double t18 = t6*t12*2.73E+7;
		double t19 = -t18;
		double t21 = t7*t13*4.3666E+26;
		double t20 = exp(t19);
		double t22 = t21+1.1944E+6;
		double t23 = param1*t7*t14*t20*u6*1.30998E+27;
		double t24 = param1*t7*t14*t20*u9*1.30998E+27;
		double t25 = param1*t20*t22;
		double t26 = -t25;
		double t28 = t6*t13*t25*u6*2.73E+7;
		double t29 = t6*t13*t25*u9*2.73E+7;
		double t27 = t9+t26;
		double t30 = -t28;
		double t31 = -t29;
		double t32 = param3*t5*t17*t27*u1*u6*x1*6.23662E-1;
		double t33 = param3*t5*t17*t27*u1*u9*x1*6.23662E-1;
		double t34 = t5*t10*t16*t27*x1*2.3987;
		double t36 = t5*t12*t16*t27*u1*u6*x1*2.3987;
		double t37 = t5*t12*t16*t27*u1*u9*x1*2.3987;
		double t40 = t23+t30;
		double t41 = t24+t31;
		double t35 = -t34;
		double t38 = -t36;
		double t39 = -t37;
		double t42 = t5*t10*t16*t40*u1*x1*2.3987;
		double t43 = t5*t10*t16*t41*u1*x1*2.3987;
		double t44 = -t42;
		double t45 = -t43;
		double t46 = t32+t38+t44;
		double t47 = t33+t39+t45;
		s_udg[0*ng+i] = t35;
		s_udg[1*ng+i] = t35;
		s_udg[2*ng+i] = -t11;
		s_udg[3*ng+i] = 0.0;
		s_udg[4*ng+i] = 0.0;
		s_udg[5*ng+i] = t11;
		s_udg[6*ng+i] = 0.0;
		s_udg[7*ng+i] = 0.0;
		s_udg[8*ng+i] = 0.0;
		s_udg[9*ng+i] = 0.0;
		s_udg[10*ng+i] = 0.0;
		s_udg[11*ng+i] = 0.0;
		s_udg[12*ng+i] = 0.0;
		s_udg[13*ng+i] = 0.0;
		s_udg[14*ng+i] = 0.0;
		s_udg[15*ng+i] = t46;
		s_udg[16*ng+i] = t46;
		s_udg[17*ng+i] = 0.0;
		s_udg[18*ng+i] = 0.0;
		s_udg[19*ng+i] = 0.0;
		s_udg[20*ng+i] = 0.0;
		s_udg[21*ng+i] = 0.0;
		s_udg[22*ng+i] = 0.0;
		s_udg[23*ng+i] = 0.0;
		s_udg[24*ng+i] = t47;
		s_udg[25*ng+i] = t47;
		s_udg[26*ng+i] = 0.0;

	}
}

// void sourceonly_streamer2d(double *s, double *pg, double *udg, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
// {
// 	double param1 = param[0];
// 	double param2 = param[1];
// 	double param3 = param[2];
// 	double param4 = param[3];
// 	double param5 = param[4];
// 	double param6 = param[5];
// 	double param7 = param[6];
// 	double param8 = param[7];
// 	double param9 = param[8];
// 	double param10 = param[9];

// 	for (int i = 0; i <ng; i++) {
// 		double x1 = pg[0*ng+i];
// 		double x2 = pg[1*ng+i];
// 		double u1 = udg[0*ng+i];
// 		double u2 = udg[1*ng+i];
// 		double u3 = udg[2*ng+i];
// 		double u4 = udg[3*ng+i];
// 		double u5 = udg[4*ng+i];
// 		double u6 = udg[5*ng+i];
// 		double u7 = udg[6*ng+i];
// 		double u8 = udg[7*ng+i];
// 		double u9 = udg[8*ng+i];

// 		double t2 = u6*u6;
// 		double t3 = u9*u9;
// 		double t4 = 1.0/param2;
// 		double t5 = 1.0/param3;
// 		double t8 = param1*3.4075E+2;
// 		double t6 = t5*t5*t5;
// 		double t7 = t2+t3;
// 		double t9 = sqrt(t7);
// 		double t10 = 1.0/t9;
// 		double t12 = param3*t9;
// 		double t11 = t10*t10*t10;
// 		double t13 = 1.0/pow(t12,1.3E+1/5.0E+1);
// 		double t14 = t5*t10*2.73E+7;
// 		double t15 = -t14;
// 		double t17 = t6*t11*4.3666E+26;
// 		double t16 = exp(t15);
// 		double t18 = t17+1.1944E+6;
// 		double t19 = param1*t16*t18;
// 		double t20 = -t19;
// 		double t21 = t8+t20;
// 		double t22 = t4*t9*t13*t21*u1*x1*2.3987;
// 		double t23 = -t22;
// 		s[0*ng+i] = t23;
// 		s[1*ng+i] = t23;
// 		s[2*ng+i] = -1.0/(param1*param1)*param4*t5*x1*(u1-u2);

// 	}
// }

