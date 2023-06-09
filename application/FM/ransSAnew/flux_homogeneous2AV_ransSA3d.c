
// Written by: C. Nguyen & P. Fernandez

void flux_homogeneous2AV_ransSA3d(double *f, double *f_udg, double *pg, double *udg, appstruct &app, double *param,
                                  double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)
{
    // This function assumes the non-dimensionalization is such that v_inf = 1

    double r, ru, rv, rw, rE, rx, ry, rz, rux, rvy, rwz, u, v, w, ux, vy, wz, p;
    double p_reg, D_p_reg_D_p, r_reg, D_r_reg_D_r, trash;
    double c, g, x, l, wallDistance, ff, D_ff_D_x;

    double param1 = param[0];
    double param2 = param[1];
    double param3 = param[2];
    double param4 = param[3];
    double param5 = param[4];
    double param6 = param[5];
    double c_inf = 1/param5;

    double alpha = 1.0e4;
    double beta = 1.0e-2;
    double alpha_ff = 10.0;
    double beta_ff = 0.5;
    double eps_v = 1.0e-8;
    double rRegMin = 1.0e-2;
    double pRegMin = 1.0e-3;
    double rampFactor = app.rampFactor;
    double epsilon0 = 1.0;
    double h0 = 1.0e-3;

    for (int i = 0; i <ng; i++) {
        double x1 = pg[0*ng+i];
        double x2 = pg[1*ng+i];
        double x3 = pg[2*ng+i];

        wallDistance = pg[3*ng+i];
        error("Need to determine the initial index of wallDistance\n");

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
        double u21 = udg[20*ng+i];
        double u22 = udg[21*ng+i];
        double u23 = udg[22*ng+i];
        double u24 = udg[23*ng+i];

        r = u1;
        ru = u2;
        rv = u3;
        rw = u4;
        rE = u5;
        rx = u7;
        rux = u8;
        ry = u13;
        rvy = u15;
        rz = u19;
        rwz = u22;

        u = ru/r;
        v = rv/r;
        w = rw/r;
        p = (param1-1)*(rE-r*(u*u+v*v+w*w+eps_v*eps_v)/2);
        surrogateMaximum(&p_reg, &trash, &D_p_reg_D_p, pRegMin, p, alpha);
        surrogateMaximum(&r_reg, &trash, &D_r_reg_D_r, rRegMin, r, alpha);

        ux = (rux - u*rx)/r;
        vy = (rvy - v*ry)/r;
        wz = (rwz - w*rz)/r;

        c = sqrt((param1*p_reg)/r_reg);
        l = min(h0,10*wallDistance);
        g = c_inf*sqrt(1.0 + 0.5*log(1.0 + exp(2*((c/c_inf)*(c/c_inf)-1.0))));
        x = (l*(ux+vy+wz))/g;
        ff = log(1 + exp(alpha_ff*(x-beta_ff)))/alpha_ff;

        double t2 = 1.0/(u1*u1);
        double t3 = 1.0/u1;
        double t4 = param1*(1.0/2.0);
        double t5 = t4-1.0/2.0;
        double t6 = u2*u2;
        double t7 = t2*t6;
        double t8 = u3*u3;
        double t9 = t2*t8;
        double t10 = u4*u4;
        double t11 = t2*t10;
        double t12 = eps_v*eps_v;
        double t13 = t7+t9+t11+t12;

        f[0*ng+i] += epsilon0*ff*rampFactor*u7;
        f[1*ng+i] += epsilon0*ff*rampFactor*u8;
        f[2*ng+i] += epsilon0*ff*rampFactor*u9;
        f[3*ng+i] += epsilon0*ff*rampFactor*u10;
        f[4*ng+i] += -epsilon0*ff*rampFactor*(t5*(t13*u7+t3*u2*(u8-t3*u2*u7)*2.0+t3*u3*(u9-t3*u3*u7)*2.0+t3*u4*(u10-t3*u4*u7)*2.0)-param1*u11);
        f[5*ng+i] += epsilon0*ff*rampFactor*u12;
        f[6*ng+i] += epsilon0*ff*rampFactor*u13;
        f[7*ng+i] += epsilon0*ff*rampFactor*u14;
        f[8*ng+i] += epsilon0*ff*rampFactor*u15;
        f[9*ng+i] += epsilon0*ff*rampFactor*u16;
        f[10*ng+i] += -epsilon0*ff*rampFactor*(t5*(t13*u13+t3*u2*(u14-t3*u2*u13)*2.0+t3*u3*(u15-t3*u3*u13)*2.0+t3*u4*(u16-t3*u4*u13)*2.0)-param1*u17);
        f[11*ng+i] += epsilon0*ff*rampFactor*u18;
        f[12*ng+i] += epsilon0*ff*rampFactor*u19;
        f[13*ng+i] += epsilon0*ff*rampFactor*u20;
        f[14*ng+i] += epsilon0*ff*rampFactor*u21;
        f[15*ng+i] += epsilon0*ff*rampFactor*u22;
        f[16*ng+i] += -epsilon0*ff*rampFactor*(t5*(t13*u19+t3*u2*(u20-t3*u2*u19)*2.0+t3*u3*(u21-t3*u3*u19)*2.0+t3*u4*(u22-t3*u4*u19)*2.0)-param1*u23);
        f[17*ng+i] += epsilon0*ff*rampFactor*u24;
    }

    if (computeJacobian == 1) {

        for (int i = 0; i <ng; i++) {
            double x1 = pg[0*ng+i];
            double x2 = pg[1*ng+i];
            double x3 = pg[2*ng+i];

            wallDistance = pg[3*ng+i];
            error("Need to determine the initial index of wallDistance\n");

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
            double u21 = udg[20*ng+i];
            double u22 = udg[21*ng+i];
            double u23 = udg[22*ng+i];
            double u24 = udg[23*ng+i];

            r = u1;
            ru = u2;
            rv = u3;
            rw = u4;
            rE = u5;
            rx = u7;
            rux = u8;
            ry = u13;
            rvy = u15;
            rz = u19;
            rwz = u22;

            u = ru/r;
            v = rv/r;
            w = rw/r;
            p = (param1-1)*(rE-r*(u*u+v*v+w*w+eps_v*eps_v)/2);
            surrogateMaximum(&p_reg, &trash, &D_p_reg_D_p, pRegMin, p, alpha);
            surrogateMaximum(&r_reg, &trash, &D_r_reg_D_r, rRegMin, r, alpha);

            ux = (rux - u*rx)/r;
            vy = (rvy - v*ry)/r;
            wz = (rwz - w*rz)/r;

            c = sqrt((param1*p_reg)/r_reg);
            l = min(h0,10*wallDistance);
            g = c_inf*sqrt(1.0 + 0.5*log(1.0 + exp(2*((c/c_inf)*(c/c_inf)-1.0))));
            x = (l*(ux+vy+wz))/g;
            ff = log(1 + exp(alpha_ff*(x-beta_ff)))/alpha_ff;
            D_ff_D_x = exp(alpha_ff*(x-beta_ff))/(1 + exp(alpha_ff*(x-beta_ff)));

            double t2 = 1.0/(u1*u1);
            double t3 = 1.0/u1;
            double t4 = 1.0/(u1*u1*u1);
            double t5 = param5*param5;
            double t6 = 1.0/r_reg;
            double t7 = p_reg*param1*t5*t6*2.0;
            double t8 = t7-2.0;
            double t9 = exp(t8);
            double t10 = t9+1.0;
            double t11 = log(t10);
            double t12 = t11*(1.0/2.0);
            double t13 = t12+1.0;
            double t19 = t3*u2*u7;
            double t14 = -t19+u8;
            double t21 = t3*u3*u13;
            double t15 = -t21+u15;
            double t23 = t3*u4*u19;
            double t16 = -t23+u22;
            double t17 = 1.0/t10;
            double t18 = 1.0/pow(t13,3.0/2.0);
            double t20 = t3*t14;
            double t22 = t3*t15;
            double t24 = t3*t16;
            double t25 = t20+t22+t24;
            double t26 = u2*u2;
            double t27 = u3*u3;
            double t28 = u4*u4;
            double t29 = 1.0/sqrt(t13);
            double t30 = t2*t14;
            double t31 = t2*t15;
            double t32 = t2*t16;
            double t49 = t4*u2*u7;
            double t50 = t4*u3*u13;
            double t51 = t4*u4*u19;
            double t33 = t30+t31+t32-t49-t50-t51;
            double t34 = 1.0/(r_reg*r_reg);
            double t35 = D_r_reg_D_r*l*p_reg*param1*param5*t5*t9*t17*t18*t25*t34*(1.0/2.0);
            double t36 = param1-1.0;
            double t37 = t2*t26*(1.0/2.0);
            double t38 = t2*t27*(1.0/2.0);
            double t39 = t2*t28*(1.0/2.0);
            double t40 = t4*t26*2.0;
            double t41 = t4*t27*2.0;
            double t42 = t4*t28*2.0;
            double t43 = t40+t41+t42;
            double t44 = eps_v*eps_v;
            double t45 = t44*(1.0/2.0);
            double t53 = t43*u1*(1.0/2.0);
            double t46 = t37+t38+t39+t45-t53;
            double t47 = D_p_reg_D_p*l*param1*param5*t5*t6*t9*t17*t18*t25*t36*t46*(1.0/2.0);
            double t52 = l*param5*t29*t33;
            double t48 = t35+t47-t52;
            double t54 = param1*(1.0/2.0);
            double t55 = t54-1.0/2.0;
            double t72 = t3*u3*u7;
            double t56 = -t72+u9;
            double t74 = t3*u4*u7;
            double t57 = -t74+u10;
            double t58 = t2*t26;
            double t59 = t2*t27;
            double t60 = t2*t28;
            double t61 = t44+t58+t59+t60;
            double t80 = t3*u2*u13;
            double t62 = -t80+u14;
            double t83 = t3*u4*u13;
            double t63 = -t83+u16;
            double t89 = t3*u2*u19;
            double t64 = -t89+u20;
            double t91 = t3*u3*u19;
            double t65 = -t91+u21;
            double t66 = l*param5*t2*t29*u7;
            double t68 = D_p_reg_D_p*l*param1*param5*t3*t5*t6*t9*t17*t18*t25*t36*u2*(1.0/2.0);
            double t67 = t66-t68;
            double t69 = param1*u11;
            double t70 = t61*u7;
            double t71 = t3*t14*u2*2.0;
            double t73 = t3*t56*u3*2.0;
            double t75 = t3*t57*u4*2.0;
            double t76 = t70+t71+t73+t75;
            double t99 = t55*t76;
            double t77 = t69-t99;
            double t78 = param1*u17;
            double t79 = t61*u13;
            double t81 = t3*t62*u2*2.0;
            double t82 = t3*t15*u3*2.0;
            double t84 = t3*t63*u4*2.0;
            double t85 = t79+t81+t82+t84;
            double t100 = t55*t85;
            double t86 = t78-t100;
            double t87 = param1*u23;
            double t88 = t61*u19;
            double t90 = t3*t64*u2*2.0;
            double t92 = t3*t65*u3*2.0;
            double t93 = t3*t16*u4*2.0;
            double t94 = t88+t90+t92+t93;
            double t101 = t55*t94;
            double t95 = t87-t101;
            double t96 = l*param5*t2*t29*u13;
            double t98 = D_p_reg_D_p*l*param1*param5*t3*t5*t6*t9*t17*t18*t25*t36*u3*(1.0/2.0);
            double t97 = t96-t98;
            double t102 = l*param5*t2*t29*u19;
            double t104 = D_p_reg_D_p*l*param1*param5*t3*t5*t6*t9*t17*t18*t25*t36*u4*(1.0/2.0);
            double t103 = t102-t104;
            double t105 = epsilon0*ff*rampFactor;
            double t106 = -t44+t58+t59+t60;
            double t107 = epsilon0*ff*rampFactor*t55*t106;
            double t108 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u7;
            double t109 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u8;
            double t110 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u9;
            double t111 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u10;
            double t112 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*t77;
            double t113 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u12;
            double t114 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u13;
            double t115 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u14;
            double t116 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u15;
            double t117 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u16;
            double t118 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*t86;
            double t119 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u18;
            double t120 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u19;
            double t121 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u20;
            double t122 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u21;
            double t123 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u22;
            double t124 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*t95;
            double t125 = D_ff_D_x*epsilon0*l*param5*rampFactor*t3*t29*u24;
            double t126 = epsilon0*ff*param1*rampFactor;

            f_udg[0*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u7;
            f_udg[1*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u8;
            f_udg[2*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u9;
            f_udg[3*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u10;
            f_udg[4*ng+i] += epsilon0*ff*rampFactor*t55*(t43*u7+t2*t14*u2*2.0-t4*t26*u7*2.0-t4*t27*u7*2.0-t4*t28*u7*2.0+t2*t56*u3*2.0+t2*t57*u4*2.0)+D_ff_D_x*epsilon0*rampFactor*t48*t77;
            f_udg[5*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u12;
            f_udg[6*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u13;
            f_udg[7*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u14;
            f_udg[8*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u15;
            f_udg[9*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u16;
            f_udg[10*ng+i] += epsilon0*ff*rampFactor*t55*(t43*u13+t2*t15*u3*2.0-t4*t26*u13*2.0-t4*t27*u13*2.0-t4*t28*u13*2.0+t2*t62*u2*2.0+t2*t63*u4*2.0)+D_ff_D_x*epsilon0*rampFactor*t48*t86;
            f_udg[11*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u18;
            f_udg[12*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u19;
            f_udg[13*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u20;
            f_udg[14*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u21;
            f_udg[15*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u22;
            f_udg[16*ng+i] += epsilon0*ff*rampFactor*t55*(t43*u19+t2*t16*u4*2.0-t4*t26*u19*2.0-t4*t27*u19*2.0-t4*t28*u19*2.0+t2*t64*u2*2.0+t2*t65*u3*2.0)+D_ff_D_x*epsilon0*rampFactor*t48*t95;
            f_udg[17*ng+i] += D_ff_D_x*epsilon0*rampFactor*t48*u24;
            f_udg[18*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u7;
            f_udg[19*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u8;
            f_udg[20*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u9;
            f_udg[21*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u10;
            f_udg[22*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*t77-epsilon0*ff*rampFactor*t3*t14*t55*2.0;
            f_udg[23*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u12;
            f_udg[24*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u13;
            f_udg[25*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u14;
            f_udg[26*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u15;
            f_udg[27*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u16;
            f_udg[28*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*t86-epsilon0*ff*rampFactor*t3*t55*t62*2.0;
            f_udg[29*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u18;
            f_udg[30*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u19;
            f_udg[31*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u20;
            f_udg[32*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u21;
            f_udg[33*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u22;
            f_udg[34*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*t95-epsilon0*ff*rampFactor*t3*t55*t64*2.0;
            f_udg[35*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t67*u24;
            f_udg[36*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u7;
            f_udg[37*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u8;
            f_udg[38*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u9;
            f_udg[39*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u10;
            f_udg[40*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t77*t97-epsilon0*ff*rampFactor*t3*t55*t56*2.0;
            f_udg[41*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u12;
            f_udg[42*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u13;
            f_udg[43*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u14;
            f_udg[44*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u15;
            f_udg[45*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u16;
            f_udg[46*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t86*t97-epsilon0*ff*rampFactor*t3*t15*t55*2.0;
            f_udg[47*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u18;
            f_udg[48*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u19;
            f_udg[49*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u20;
            f_udg[50*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u21;
            f_udg[51*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u22;
            f_udg[52*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t95*t97-epsilon0*ff*rampFactor*t3*t55*t65*2.0;
            f_udg[53*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t97*u24;
            f_udg[54*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u7;
            f_udg[55*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u8;
            f_udg[56*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u9;
            f_udg[57*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u10;
            f_udg[58*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t77*t103-epsilon0*ff*rampFactor*t3*t55*t57*2.0;
            f_udg[59*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u12;
            f_udg[60*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u13;
            f_udg[61*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u14;
            f_udg[62*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u15;
            f_udg[63*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u16;
            f_udg[64*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t86*t103-epsilon0*ff*rampFactor*t3*t55*t63*2.0;
            f_udg[65*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u18;
            f_udg[66*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u19;
            f_udg[67*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u20;
            f_udg[68*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u21;
            f_udg[69*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u22;
            f_udg[70*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t95*t103-epsilon0*ff*rampFactor*t3*t16*t55*2.0;
            f_udg[71*ng+i] += -D_ff_D_x*epsilon0*rampFactor*t103*u24;
            f_udg[72*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u7*(-1.0/2.0);
            f_udg[73*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u8*(-1.0/2.0);
            f_udg[74*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u9*(-1.0/2.0);
            f_udg[75*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u10*(-1.0/2.0);
            f_udg[76*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*t77*(-1.0/2.0);
            f_udg[77*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u12*(-1.0/2.0);
            f_udg[78*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u13*(-1.0/2.0);
            f_udg[79*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u14*(-1.0/2.0);
            f_udg[80*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u15*(-1.0/2.0);
            f_udg[81*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u16*(-1.0/2.0);
            f_udg[82*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*t86*(-1.0/2.0);
            f_udg[83*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u18*(-1.0/2.0);
            f_udg[84*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u19*(-1.0/2.0);
            f_udg[85*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u20*(-1.0/2.0);
            f_udg[86*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u21*(-1.0/2.0);
            f_udg[87*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u22*(-1.0/2.0);
            f_udg[88*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*t95*(-1.0/2.0);
            f_udg[89*ng+i] += D_ff_D_x*D_p_reg_D_p*epsilon0*l*param1*param5*rampFactor*t5*t6*t9*t17*t18*t25*t36*u24*(-1.0/2.0);
            f_udg[90*ng+i] += 0.0;
            f_udg[91*ng+i] += 0.0;
            f_udg[92*ng+i] += 0.0;
            f_udg[93*ng+i] += 0.0;
            f_udg[94*ng+i] += 0.0;
            f_udg[95*ng+i] += 0.0;
            f_udg[96*ng+i] += 0.0;
            f_udg[97*ng+i] += 0.0;
            f_udg[98*ng+i] += 0.0;
            f_udg[99*ng+i] += 0.0;
            f_udg[100*ng+i] += 0.0;
            f_udg[101*ng+i] += 0.0;
            f_udg[102*ng+i] += 0.0;
            f_udg[103*ng+i] += 0.0;
            f_udg[104*ng+i] += 0.0;
            f_udg[105*ng+i] += 0.0;
            f_udg[106*ng+i] += 0.0;
            f_udg[107*ng+i] += 0.0;
            f_udg[108*ng+i] += t105-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u7;
            f_udg[109*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u8;
            f_udg[110*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u9;
            f_udg[111*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u10;
            f_udg[112*ng+i] += t107-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t77*u2;
            f_udg[113*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u12;
            f_udg[114*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u13;
            f_udg[115*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u14;
            f_udg[116*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u15;
            f_udg[117*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u16;
            f_udg[118*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t86*u2;
            f_udg[119*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u18;
            f_udg[120*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u19;
            f_udg[121*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u20;
            f_udg[122*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u21;
            f_udg[123*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u22;
            f_udg[124*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t95*u2;
            f_udg[125*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u2*u24;
            f_udg[126*ng+i] += t108;
            f_udg[127*ng+i] += t105+t109;
            f_udg[128*ng+i] += t110;
            f_udg[129*ng+i] += t111;
            f_udg[130*ng+i] += t112-epsilon0*ff*rampFactor*t3*t55*u2*2.0;
            f_udg[131*ng+i] += t113;
            f_udg[132*ng+i] += t114;
            f_udg[133*ng+i] += t115;
            f_udg[134*ng+i] += t116;
            f_udg[135*ng+i] += t117;
            f_udg[136*ng+i] += t118;
            f_udg[137*ng+i] += t119;
            f_udg[138*ng+i] += t120;
            f_udg[139*ng+i] += t121;
            f_udg[140*ng+i] += t122;
            f_udg[141*ng+i] += t123;
            f_udg[142*ng+i] += t124;
            f_udg[143*ng+i] += t125;
            f_udg[144*ng+i] += 0.0;
            f_udg[145*ng+i] += 0.0;
            f_udg[146*ng+i] += t105;
            f_udg[147*ng+i] += 0.0;
            f_udg[148*ng+i] += epsilon0*ff*rampFactor*t3*t55*u3*-2.0;
            f_udg[149*ng+i] += 0.0;
            f_udg[150*ng+i] += 0.0;
            f_udg[151*ng+i] += 0.0;
            f_udg[152*ng+i] += 0.0;
            f_udg[153*ng+i] += 0.0;
            f_udg[154*ng+i] += 0.0;
            f_udg[155*ng+i] += 0.0;
            f_udg[156*ng+i] += 0.0;
            f_udg[157*ng+i] += 0.0;
            f_udg[158*ng+i] += 0.0;
            f_udg[159*ng+i] += 0.0;
            f_udg[160*ng+i] += 0.0;
            f_udg[161*ng+i] += 0.0;
            f_udg[162*ng+i] += 0.0;
            f_udg[163*ng+i] += 0.0;
            f_udg[164*ng+i] += 0.0;
            f_udg[165*ng+i] += t105;
            f_udg[166*ng+i] += epsilon0*ff*rampFactor*t3*t55*u4*-2.0;
            f_udg[167*ng+i] += 0.0;
            f_udg[168*ng+i] += 0.0;
            f_udg[169*ng+i] += 0.0;
            f_udg[170*ng+i] += 0.0;
            f_udg[171*ng+i] += 0.0;
            f_udg[172*ng+i] += 0.0;
            f_udg[173*ng+i] += 0.0;
            f_udg[174*ng+i] += 0.0;
            f_udg[175*ng+i] += 0.0;
            f_udg[176*ng+i] += 0.0;
            f_udg[177*ng+i] += 0.0;
            f_udg[178*ng+i] += 0.0;
            f_udg[179*ng+i] += 0.0;
            f_udg[180*ng+i] += 0.0;
            f_udg[181*ng+i] += 0.0;
            f_udg[182*ng+i] += 0.0;
            f_udg[183*ng+i] += 0.0;
            f_udg[184*ng+i] += t126;
            f_udg[185*ng+i] += 0.0;
            f_udg[186*ng+i] += 0.0;
            f_udg[187*ng+i] += 0.0;
            f_udg[188*ng+i] += 0.0;
            f_udg[189*ng+i] += 0.0;
            f_udg[190*ng+i] += 0.0;
            f_udg[191*ng+i] += 0.0;
            f_udg[192*ng+i] += 0.0;
            f_udg[193*ng+i] += 0.0;
            f_udg[194*ng+i] += 0.0;
            f_udg[195*ng+i] += 0.0;
            f_udg[196*ng+i] += 0.0;
            f_udg[197*ng+i] += 0.0;
            f_udg[198*ng+i] += 0.0;
            f_udg[199*ng+i] += 0.0;
            f_udg[200*ng+i] += 0.0;
            f_udg[201*ng+i] += 0.0;
            f_udg[202*ng+i] += 0.0;
            f_udg[203*ng+i] += t105;
            f_udg[204*ng+i] += 0.0;
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
            f_udg[216*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u7;
            f_udg[217*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u8;
            f_udg[218*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u9;
            f_udg[219*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u10;
            f_udg[220*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t77*u3;
            f_udg[221*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u12;
            f_udg[222*ng+i] += t105-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u13;
            f_udg[223*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u14;
            f_udg[224*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u15;
            f_udg[225*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u16;
            f_udg[226*ng+i] += t107-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t86*u3;
            f_udg[227*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u18;
            f_udg[228*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u19;
            f_udg[229*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u20;
            f_udg[230*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u21;
            f_udg[231*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u22;
            f_udg[232*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t95*u3;
            f_udg[233*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u3*u24;
            f_udg[234*ng+i] += 0.0;
            f_udg[235*ng+i] += 0.0;
            f_udg[236*ng+i] += 0.0;
            f_udg[237*ng+i] += 0.0;
            f_udg[238*ng+i] += 0.0;
            f_udg[239*ng+i] += 0.0;
            f_udg[240*ng+i] += 0.0;
            f_udg[241*ng+i] += t105;
            f_udg[242*ng+i] += 0.0;
            f_udg[243*ng+i] += 0.0;
            f_udg[244*ng+i] += epsilon0*ff*rampFactor*t3*t55*u2*-2.0;
            f_udg[245*ng+i] += 0.0;
            f_udg[246*ng+i] += 0.0;
            f_udg[247*ng+i] += 0.0;
            f_udg[248*ng+i] += 0.0;
            f_udg[249*ng+i] += 0.0;
            f_udg[250*ng+i] += 0.0;
            f_udg[251*ng+i] += 0.0;
            f_udg[252*ng+i] += t108;
            f_udg[253*ng+i] += t109;
            f_udg[254*ng+i] += t110;
            f_udg[255*ng+i] += t111;
            f_udg[256*ng+i] += t112;
            f_udg[257*ng+i] += t113;
            f_udg[258*ng+i] += t114;
            f_udg[259*ng+i] += t115;
            f_udg[260*ng+i] += t105+t116;
            f_udg[261*ng+i] += t117;
            f_udg[262*ng+i] += t118-epsilon0*ff*rampFactor*t3*t55*u3*2.0;
            f_udg[263*ng+i] += t119;
            f_udg[264*ng+i] += t120;
            f_udg[265*ng+i] += t121;
            f_udg[266*ng+i] += t122;
            f_udg[267*ng+i] += t123;
            f_udg[268*ng+i] += t124;
            f_udg[269*ng+i] += t125;
            f_udg[270*ng+i] += 0.0;
            f_udg[271*ng+i] += 0.0;
            f_udg[272*ng+i] += 0.0;
            f_udg[273*ng+i] += 0.0;
            f_udg[274*ng+i] += 0.0;
            f_udg[275*ng+i] += 0.0;
            f_udg[276*ng+i] += 0.0;
            f_udg[277*ng+i] += 0.0;
            f_udg[278*ng+i] += 0.0;
            f_udg[279*ng+i] += t105;
            f_udg[280*ng+i] += epsilon0*ff*rampFactor*t3*t55*u4*-2.0;
            f_udg[281*ng+i] += 0.0;
            f_udg[282*ng+i] += 0.0;
            f_udg[283*ng+i] += 0.0;
            f_udg[284*ng+i] += 0.0;
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
            f_udg[298*ng+i] += t126;
            f_udg[299*ng+i] += 0.0;
            f_udg[300*ng+i] += 0.0;
            f_udg[301*ng+i] += 0.0;
            f_udg[302*ng+i] += 0.0;
            f_udg[303*ng+i] += 0.0;
            f_udg[304*ng+i] += 0.0;
            f_udg[305*ng+i] += 0.0;
            f_udg[306*ng+i] += 0.0;
            f_udg[307*ng+i] += 0.0;
            f_udg[308*ng+i] += 0.0;
            f_udg[309*ng+i] += 0.0;
            f_udg[310*ng+i] += 0.0;
            f_udg[311*ng+i] += 0.0;
            f_udg[312*ng+i] += 0.0;
            f_udg[313*ng+i] += 0.0;
            f_udg[314*ng+i] += 0.0;
            f_udg[315*ng+i] += 0.0;
            f_udg[316*ng+i] += 0.0;
            f_udg[317*ng+i] += t105;
            f_udg[318*ng+i] += 0.0;
            f_udg[319*ng+i] += 0.0;
            f_udg[320*ng+i] += 0.0;
            f_udg[321*ng+i] += 0.0;
            f_udg[322*ng+i] += 0.0;
            f_udg[323*ng+i] += 0.0;
            f_udg[324*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u7;
            f_udg[325*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u8;
            f_udg[326*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u9;
            f_udg[327*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u10;
            f_udg[328*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t77*u4;
            f_udg[329*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u12;
            f_udg[330*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u13;
            f_udg[331*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u14;
            f_udg[332*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u15;
            f_udg[333*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u16;
            f_udg[334*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t86*u4;
            f_udg[335*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u18;
            f_udg[336*ng+i] += t105-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u19;
            f_udg[337*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u20;
            f_udg[338*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u21;
            f_udg[339*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u22;
            f_udg[340*ng+i] += t107-D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*t95*u4;
            f_udg[341*ng+i] += -D_ff_D_x*epsilon0*l*param5*rampFactor*t2*t29*u4*u24;
            f_udg[342*ng+i] += 0.0;
            f_udg[343*ng+i] += 0.0;
            f_udg[344*ng+i] += 0.0;
            f_udg[345*ng+i] += 0.0;
            f_udg[346*ng+i] += 0.0;
            f_udg[347*ng+i] += 0.0;
            f_udg[348*ng+i] += 0.0;
            f_udg[349*ng+i] += 0.0;
            f_udg[350*ng+i] += 0.0;
            f_udg[351*ng+i] += 0.0;
            f_udg[352*ng+i] += 0.0;
            f_udg[353*ng+i] += 0.0;
            f_udg[354*ng+i] += 0.0;
            f_udg[355*ng+i] += t105;
            f_udg[356*ng+i] += 0.0;
            f_udg[357*ng+i] += 0.0;
            f_udg[358*ng+i] += epsilon0*ff*rampFactor*t3*t55*u2*-2.0;
            f_udg[359*ng+i] += 0.0;
            f_udg[360*ng+i] += 0.0;
            f_udg[361*ng+i] += 0.0;
            f_udg[362*ng+i] += 0.0;
            f_udg[363*ng+i] += 0.0;
            f_udg[364*ng+i] += 0.0;
            f_udg[365*ng+i] += 0.0;
            f_udg[366*ng+i] += 0.0;
            f_udg[367*ng+i] += 0.0;
            f_udg[368*ng+i] += 0.0;
            f_udg[369*ng+i] += 0.0;
            f_udg[370*ng+i] += 0.0;
            f_udg[371*ng+i] += 0.0;
            f_udg[372*ng+i] += 0.0;
            f_udg[373*ng+i] += 0.0;
            f_udg[374*ng+i] += t105;
            f_udg[375*ng+i] += 0.0;
            f_udg[376*ng+i] += epsilon0*ff*rampFactor*t3*t55*u3*-2.0;
            f_udg[377*ng+i] += 0.0;
            f_udg[378*ng+i] += t108;
            f_udg[379*ng+i] += t109;
            f_udg[380*ng+i] += t110;
            f_udg[381*ng+i] += t111;
            f_udg[382*ng+i] += t112;
            f_udg[383*ng+i] += t113;
            f_udg[384*ng+i] += t114;
            f_udg[385*ng+i] += t115;
            f_udg[386*ng+i] += t116;
            f_udg[387*ng+i] += t117;
            f_udg[388*ng+i] += t118;
            f_udg[389*ng+i] += t119;
            f_udg[390*ng+i] += t120;
            f_udg[391*ng+i] += t121;
            f_udg[392*ng+i] += t122;
            f_udg[393*ng+i] += t105+t123;
            f_udg[394*ng+i] += t124-epsilon0*ff*rampFactor*t3*t55*u4*2.0;
            f_udg[395*ng+i] += t125;
            f_udg[396*ng+i] += 0.0;
            f_udg[397*ng+i] += 0.0;
            f_udg[398*ng+i] += 0.0;
            f_udg[399*ng+i] += 0.0;
            f_udg[400*ng+i] += 0.0;
            f_udg[401*ng+i] += 0.0;
            f_udg[402*ng+i] += 0.0;
            f_udg[403*ng+i] += 0.0;
            f_udg[404*ng+i] += 0.0;
            f_udg[405*ng+i] += 0.0;
            f_udg[406*ng+i] += 0.0;
            f_udg[407*ng+i] += 0.0;
            f_udg[408*ng+i] += 0.0;
            f_udg[409*ng+i] += 0.0;
            f_udg[410*ng+i] += 0.0;
            f_udg[411*ng+i] += 0.0;
            f_udg[412*ng+i] += t126;
            f_udg[413*ng+i] += 0.0;
            f_udg[414*ng+i] += 0.0;
            f_udg[415*ng+i] += 0.0;
            f_udg[416*ng+i] += 0.0;
            f_udg[417*ng+i] += 0.0;
            f_udg[418*ng+i] += 0.0;
            f_udg[419*ng+i] += 0.0;
            f_udg[420*ng+i] += 0.0;
            f_udg[421*ng+i] += 0.0;
            f_udg[422*ng+i] += 0.0;
            f_udg[423*ng+i] += 0.0;
            f_udg[424*ng+i] += 0.0;
            f_udg[425*ng+i] += 0.0;
            f_udg[426*ng+i] += 0.0;
            f_udg[427*ng+i] += 0.0;
            f_udg[428*ng+i] += 0.0;
            f_udg[429*ng+i] += 0.0;
            f_udg[430*ng+i] += 0.0;
            f_udg[431*ng+i] += t105;
        }
    }
}
