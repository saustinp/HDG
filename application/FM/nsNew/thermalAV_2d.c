		double t2 = 1.0/(r*r);
		double t3 = ru*ru;
		double t4 = t2*t3*(1.0/2.0);
		double t5 = rv*rv;
		double t6 = t2*t5*(1.0/2.0);
		double t7 = t4+t6;
		double t8 = 1.0/r;
		double t9 = 1.0/gam1;
		double t10 = rE-r*t7;
		f[0*ng+i] += 0.0;
		f[1*ng+i] += 0.0;
		f[2*ng+i] += 0.0;
		f[3*ng+i] += -av*t8*t9*(gam1*r*(-rEx+rx*t7+r*(ru*t2*(rux-ru*rx*t8)+rv*t2*(rvx-rv*rx*t8)))+gam1*rx*t10);
		f[4*ng+i] += 0.0;
		f[5*ng+i] += 0.0;
		f[6*ng+i] += 0.0;
		f[7*ng+i] += -av*t8*t9*(gam1*r*(-rEy+ry*t7+r*(ru*t2*(ruy-ru*ry*t8)+rv*t2*(rvy-rv*ry*t8)))+gam1*ry*t10);

			double t2 = 1.0/(r*r);
			double t3 = 1.0/r;
			double t4 = ru*ru;
			double t5 = 1.0/(r*r*r*r);
			double t6 = rv*rv;
			double t10 = ru*rx*t3;
			double t7 = rux-t10;
			double t8 = 1.0/(r*r*r);
			double t12 = rv*rx*t3;
			double t9 = rvx-t12;
			double t11 = ru*t2*t7;
			double t13 = rv*t2*t9;
			double t14 = t2*t4*(1.0/2.0);
			double t15 = t2*t6*(1.0/2.0);
			double t16 = t4*t8;
			double t17 = t6*t8;
			double t18 = t16+t17;
			double t19 = 1.0/gam1;
			double t20 = t14+t15;
			double t21 = rx*t20;
			double t22 = t11+t13;
			double t23 = r*t22;
			double t24 = -rEx+t21+t23;
			double t27 = ru*ry*t3;
			double t25 = ruy-t27;
			double t29 = rv*ry*t3;
			double t26 = rvy-t29;
			double t28 = ru*t2*t25;
			double t30 = rv*t2*t26;
			double t38 = r*t18;
			double t31 = t14+t15-t38;
			double t37 = r*t20;
			double t32 = rE-t37;
			double t33 = ry*t20;
			double t34 = t28+t30;
			double t35 = r*t34;
			double t36 = -rEy+t33+t35;
			double t39 = gam1*t32;
			double t40 = gam1*r*t31;
			double t41 = t39+t40;
			f_udg[0*ng+i] += 0.0;
			f_udg[1*ng+i] += 0.0;
			f_udg[2*ng+i] += 0.0;
			f_udg[3*ng+i] += av*t2*t19*(gam1*r*t24+gam1*rx*t32)-av*t3*t19*(gam1*t24-gam1*rx*t31+gam1*r*(t11+t13-r*(ru*t7*t8*2.0+rv*t8*t9*2.0-rx*t4*t5-rx*t5*t6)-rx*t18));
			f_udg[4*ng+i] += 0.0;
			f_udg[5*ng+i] += 0.0;
			f_udg[6*ng+i] += 0.0;
			f_udg[7*ng+i] += av*t2*t19*(gam1*r*t36+gam1*ry*t32)-av*t3*t19*(gam1*t36-gam1*ry*t31+gam1*r*(t28+t30-r*(ru*t8*t25*2.0+rv*t8*t26*2.0-ry*t4*t5-ry*t5*t6)-ry*t18));
			f_udg[8*ng+i] += 0.0;
			f_udg[9*ng+i] += 0.0;
			f_udg[10*ng+i] += 0.0;
			f_udg[11*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t7-ru*rx*t8)+ru*rx*t2)-gam1*ru*rx*t3);
			f_udg[12*ng+i] += 0.0;
			f_udg[13*ng+i] += 0.0;
			f_udg[14*ng+i] += 0.0;
			f_udg[15*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t25-ru*ry*t8)+ru*ry*t2)-gam1*ru*ry*t3);
			f_udg[16*ng+i] += 0.0;
			f_udg[17*ng+i] += 0.0;
			f_udg[18*ng+i] += 0.0;
			f_udg[19*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t9-rv*rx*t8)+rv*rx*t2)-gam1*rv*rx*t3);
			f_udg[20*ng+i] += 0.0;
			f_udg[21*ng+i] += 0.0;
			f_udg[22*ng+i] += 0.0;
			f_udg[23*ng+i] += -av*t3*t19*(gam1*r*(r*(t2*t26-rv*ry*t8)+rv*ry*t2)-gam1*rv*ry*t3);
			f_udg[24*ng+i] += 0.0;
			f_udg[25*ng+i] += 0.0;
			f_udg[26*ng+i] += 0.0;
			f_udg[27*ng+i] += -av*rx*t3;
			f_udg[28*ng+i] += 0.0;
			f_udg[29*ng+i] += 0.0;
			f_udg[30*ng+i] += 0.0;
			f_udg[31*ng+i] += -av*ry*t3;
			f_udg[32*ng+i] += 0.0;
			f_udg[33*ng+i] += 0.0;
			f_udg[34*ng+i] += 0.0;
			f_udg[35*ng+i] += -av*t3*t19*t41;
			f_udg[36*ng+i] += 0.0;
			f_udg[37*ng+i] += 0.0;
			f_udg[38*ng+i] += 0.0;
			f_udg[39*ng+i] += 0.0;
			f_udg[40*ng+i] += 0.0;
			f_udg[41*ng+i] += 0.0;
			f_udg[42*ng+i] += 0.0;
			f_udg[43*ng+i] += -av*ru*t3;
			f_udg[44*ng+i] += 0.0;
			f_udg[45*ng+i] += 0.0;
			f_udg[46*ng+i] += 0.0;
			f_udg[47*ng+i] += 0.0;
			f_udg[48*ng+i] += 0.0;
			f_udg[49*ng+i] += 0.0;
			f_udg[50*ng+i] += 0.0;
			f_udg[51*ng+i] += -av*rv*t3;
			f_udg[52*ng+i] += 0.0;
			f_udg[53*ng+i] += 0.0;
			f_udg[54*ng+i] += 0.0;
			f_udg[55*ng+i] += 0.0;
			f_udg[56*ng+i] += 0.0;
			f_udg[57*ng+i] += 0.0;
			f_udg[58*ng+i] += 0.0;
			f_udg[59*ng+i] += av;
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
			f_udg[71*ng+i] += -av*t3*t19*t41;
			f_udg[72*ng+i] += 0.0;
			f_udg[73*ng+i] += 0.0;
			f_udg[74*ng+i] += 0.0;
			f_udg[75*ng+i] += 0.0;
			f_udg[76*ng+i] += 0.0;
			f_udg[77*ng+i] += 0.0;
			f_udg[78*ng+i] += 0.0;
			f_udg[79*ng+i] += -av*ru*t3;
			f_udg[80*ng+i] += 0.0;
			f_udg[81*ng+i] += 0.0;
			f_udg[82*ng+i] += 0.0;
			f_udg[83*ng+i] += 0.0;
			f_udg[84*ng+i] += 0.0;
			f_udg[85*ng+i] += 0.0;
			f_udg[86*ng+i] += 0.0;
			f_udg[87*ng+i] += -av*rv*t3;
			f_udg[88*ng+i] += 0.0;
			f_udg[89*ng+i] += 0.0;
			f_udg[90*ng+i] += 0.0;
			f_udg[91*ng+i] += 0.0;
			f_udg[92*ng+i] += 0.0;
			f_udg[93*ng+i] += 0.0;
			f_udg[94*ng+i] += 0.0;
			f_udg[95*ng+i] += av;
