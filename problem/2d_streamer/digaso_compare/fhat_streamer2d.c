#include <math.h>

void fhat_streamer2d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/param1;
		double t5 = 1.0/param2;
		double t6 = 1.0/param3;
		double t7 = t2+t3;
		double t8 = sqrt(t7);
		double t9 = param3*t8;
		double t10 = pow(t9,1.1E+1/5.0E+1);
		double t11 = 1.0/pow(t9,1.3E+1/5.0E+1);
		fh[0*ng+i] = -nl1*x1*(t5*t11*u6*uh1*2.3987-t4*t5*t6*t10*u4*4.3628E-3)-nl2*x1*(t5*t11*u9*uh1*2.3987-t4*t5*t6*t10*u7*4.3628E-3)+param10*x1*(u1-uh1);
		fh[1*ng+i] = param10*x1*(u2-uh2);
		fh[2*ng+i] = nl1*u6*x1+nl2*u9*x1+param10*x1*(u3-uh3);

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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = param10*x1;
		double t3 = u6*u6;
		double t4 = u9*u9;
		double t5 = 1.0/param1;
		double t6 = 1.0/param2;
		double t7 = 1.0/param3;
		double t8 = t3+t4;
		double t9 = sqrt(t8);
		double t10 = 1.0/t9;
		double t11 = param3*t9;
		double t12 = pow(t11,1.1E+1/5.0E+1);
		double t13 = 1.0/pow(t11,1.3E+1/5.0E+1);
		double t15 = 1.0/pow(t11,6.3E+1/5.0E+1);
		double t14 = t13*t13*t13;
		double t16 = t6*t13*uh1*2.3987;
		double t18 = param3*t6*t10*t15*u6*u9*uh1*6.23662E-1;
		double t17 = -t16;
		fh_udg[0*ng+i] = t2;
		fh_udg[1*ng+i] = 0.0;
		fh_udg[2*ng+i] = 0.0;
		fh_udg[3*ng+i] = 0.0;
		fh_udg[4*ng+i] = t2;
		fh_udg[5*ng+i] = 0.0;
		fh_udg[6*ng+i] = 0.0;
		fh_udg[7*ng+i] = 0.0;
		fh_udg[8*ng+i] = t2;
		fh_udg[9*ng+i] = nl1*t5*t6*t7*t12*x1*4.3628E-3;
		fh_udg[10*ng+i] = 0.0;
		fh_udg[11*ng+i] = 0.0;
		fh_udg[12*ng+i] = 0.0;
		fh_udg[13*ng+i] = 0.0;
		fh_udg[14*ng+i] = 0.0;
		fh_udg[15*ng+i] = nl2*x1*(t18+t5*t6*t10*t14*u6*u7*9.59816E-4)+nl1*x1*(t17+param3*t3*t6*t10*t15*uh1*6.23662E-1+t5*t6*t10*t14*u4*u6*9.59816E-4);
		fh_udg[16*ng+i] = 0.0;
		fh_udg[17*ng+i] = nl1*x1;
		fh_udg[18*ng+i] = nl2*t5*t6*t7*t12*x1*4.3628E-3;
		fh_udg[19*ng+i] = 0.0;
		fh_udg[20*ng+i] = 0.0;
		fh_udg[21*ng+i] = 0.0;
		fh_udg[22*ng+i] = 0.0;
		fh_udg[23*ng+i] = 0.0;
		fh_udg[24*ng+i] = nl1*x1*(t18+t5*t6*t10*t14*u4*u9*9.59816E-4)+nl2*x1*(t17+param3*t4*t6*t10*t15*uh1*6.23662E-1+t5*t6*t10*t14*u7*u9*9.59816E-4);
		fh_udg[25*ng+i] = 0.0;
		fh_udg[26*ng+i] = nl2*x1;

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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = param10*x1;
		double t3 = u6*u6;
		double t4 = u9*u9;
		double t5 = 1.0/param2;
		double t6 = -t2;
		double t7 = t3+t4;
		double t8 = sqrt(t7);
		double t9 = param3*t8;
		double t10 = 1.0/pow(t9,1.3E+1/5.0E+1);
		fh_uh[0*ng+i] = t6-nl1*t5*t10*u6*x1*2.3987-nl2*t5*t10*u9*x1*2.3987;
		fh_uh[1*ng+i] = 0.0;
		fh_uh[2*ng+i] = 0.0;
		fh_uh[3*ng+i] = 0.0;
		fh_uh[4*ng+i] = t6;
		fh_uh[5*ng+i] = 0.0;
		fh_uh[6*ng+i] = 0.0;
		fh_uh[7*ng+i] = 0.0;
		fh_uh[8*ng+i] = t6;

	}
}

void fhatonly_streamer2d(double *fh, double *pg, double *udg, double *uh, double *nl, double *param, double time, int ng, int nc, int ncu, int nd, int ncd)
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
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];

		double t2 = u6*u6;
		double t3 = u9*u9;
		double t4 = 1.0/param1;
		double t5 = 1.0/param2;
		double t6 = 1.0/param3;
		double t7 = t2+t3;
		double t8 = sqrt(t7);
		double t9 = param3*t8;
		double t10 = pow(t9,1.1E+1/5.0E+1);
		double t11 = 1.0/pow(t9,1.3E+1/5.0E+1);
		fh[0*ng+i] = -nl1*x1*(t5*t11*u6*uh1*2.3987-t4*t5*t6*t10*u4*4.3628E-3)-nl2*x1*(t5*t11*u9*uh1*2.3987-t4*t5*t6*t10*u7*4.3628E-3)+param10*x1*(u1-uh1);
		fh[1*ng+i] = param10*x1*(u2-uh2);
		fh[2*ng+i] = nl1*u6*x1+nl2*u9*x1+param10*x1*(u3-uh3);

	}
}

