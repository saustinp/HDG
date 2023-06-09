#include <math.h>

void saflux2d (
  double U[15],
  double Fix[5],
  double Fiy[5],
  double Fvx[5],
  double Fvy[5],
  double Jix[75],
  double Jiy[75],
  double Jvx[75],
  double Jvy[75])
{
  double t1;
  double t105;
  double t107;
  double t109;
  double t110;
  double t114;
  double t117;
  double t119;
  double t120;
  double t125;
  double t127;
  double t128;
  double t13;
  double t145;
  double t146;
  double t147;
  double t149;
  double t15;
  double t151;
  double t153;
  double t154;
  double t155;
  double t156;
  double t157;
  double t158;
  double t16;
  double t167;
  double t17;
  double t170;
  double t171;
  double t172;
  double t174;
  double t177;
  double t18;
  double t186;
  double t187;
  double t189;
  double t190;
  double t196;
  double t2;
  double t20;
  double t204;
  double t207;
  double t21;
  double t211;
  double t212;
  double t215;
  double t217;
  double t218;
  double t220;
  double t223;
  double t227;
  double t229;
  double t23;
  double t231;
  double t237;
  double t238;
  double t24;
  double t242;
  double t26;
  double t27;
  double t3;
  double t30;
  double t31;
  double t33;
  double t34;
  double t36;
  double t37;
  double t38;
  double t39;
  double t40;
  double t43;
  double t44;
  double t46;
  double t47;
  double t5;
  double t50;
  double t52;
  double t54;
  double t57;
  double t59;
  double t6;
  double t61;
  double t62;
  double t64;
  double t67;
  double t69;
  double t7;
  double t71;
  double t73;
  double t74;
  double t75;
  double t77;
  double t8;
  double t81;
  double t82;
  double t83;
  double t85;
  double t87;
  double t88;
  double t9;
  double t91;
  double t92;
  double t94;
  double t97;
  double t99;
  Fix[0] = U[1];
  t1 = pow(Fix[0], 0.2e1);
  t2 = U[0];
  t3 = 0.1e1 / t2;
  t5 = gam - 0.1e1;
  t6 = U[3];
  t7 = U[2];
  t8 = t7 * t7;
  t9 = t1 + t8;
  t13 = t5 * (t6 - t9 * t3 / 0.2e1);
  Fix[1] = t1 * t3 + t13;
  Fix[2] = Fix[0] * t7 * t3;
  t15 = t6 + t13;
  t16 = Fix[0] * t15;
  Fix[3] = t16 * t3;
  t17 = U[4];
  t18 = t17 * Fix[0];
  Fix[4] = t18 * t3;
  Fiy[0] = t7;
  Fiy[1] = Fix[2];
  Fiy[2] = t8 * t3 + t13;
  t20 = Fiy[0] * t15;
  Fiy[3] = t20 * t3;
  t21 = t17 * Fiy[0];
  Fiy[4] = t21 * t3;
  Fvx[0] = 0.0e0;
  t23 = U[5];
  t24 = t23 * Fix[0];
  t26 = U[6] - t24 * t3;
  t27 = t26 * t3;
  t30 = U[10];
  t31 = t30 * Fiy[0];
  t33 = U[12] - t31 * t3;
  t34 = t33 * t3;
  t36 = 0.4e1 / 0.3e1 * t27 - 0.2e1 / 0.3e1 * t34;
  t37 = t17 * t17;
  t38 = t37 * t37;
  t39 = pow(Re, 0.2e1);
  t40 = t39 * Re;
  t43 = t37 * t17 * t40;
  t44 = cv1 * cv1;
  t46 = t43 + t44 * cv1;
  t47 = 0.1e1 / t46;
  t50 = t38 * t40 * t47 + 0.1e1 / Re;
  Fvx[1] = t36 * t50;
  t52 = t23 * Fiy[0];
  t54 = U[7] - t52 * t3;
  t57 = t30 * Fix[0];
  t59 = U[11] - t57 * t3;
  t61 = t54 * t3 + t59 * t3;
  Fvx[2] = t61 * t50;
  t62 = t36 * Fix[0];
  t64 = t61 * Fiy[0];
  t67 = gam / Pr;
  t69 = t23 * t6;
  t71 = U[8] - t69 * t3;
  t73 = t2 * t2;
  t74 = 0.1e1 / t73;
  t75 = Fix[0] * t74;
  t77 = Fiy[0] * t74;
  t81 = t62 * t3 + t64 * t3 + t67 * (t71 * t3 - t75 * t26 - t77 * t54);
  Fvx[3] = t81 * t50;
  t82 = 0.1e1 / sigm;
  t83 = t50 * t82;
  t85 = t23 * t17;
  t87 = U[9] - t85 * t3;
  t88 = t87 * t3;
  Fvx[4] = t83 * t88;
  Fvy[0] = 0.0e0;
  Fvy[1] = Fvx[2];
  t91 = 0.4e1 / 0.3e1 * t34 - 0.2e1 / 0.3e1 * t27;
  Fvy[2] = t91 * t50;
  t92 = t61 * Fix[0];
  t94 = t91 * Fiy[0];
  t97 = t30 * t6;
  t99 = U[13] - t97 * t3;
  t105 = t92 * t3 + t94 * t3 + t67 * (t99 * t3 - t75 * t59 - t77 * t33);
  Fvy[3] = t105 * t50;
  t107 = t30 * t17;
  t109 = U[14] - t107 * t3;
  t110 = t109 * t3;
  Fvy[4] = t83 * t110;
  Jix[0] = 0.0e0;
  t114 = t5 * t9 * t74 / 0.2e1;
  Jix[1] = -t1 * t74 + t114;
  Jix[2] = -Fix[0] * Fiy[0] * t74;
  t117 = Fix[0] * t5;
  t119 = 0.1e1 / t73 / t2;
  t120 = t9 * t119;
  Jix[3] = t117 * t120 / 0.2e1 - t16 * t74;
  Jix[4] = -t18 * t74;
  Jix[5] = 0.1e1;
  t125 = Fix[0] * t3;
  t127 = t117 * t3;
  Jix[6] = 0.2e1 * t125 - t127;
  Jix[7] = Fiy[0] * t3;
  t128 = t15 * t3;
  Jix[8] = t128 - t1 * t5 * t74;
  Jix[9] = t17 * t3;
  Jix[10] = 0.0e0;
  Jix[11] = -t5 * Fiy[0] * t3;
  Jix[12] = t125;
  Jix[13] = -t117 * t77;
  Jix[14] = 0.0e0;
  Jix[15] = 0.0e0;
  Jix[16] = t5;
  Jix[17] = 0.0e0;
  Jix[18] = Fix[0] * gam * t3;
  Jix[19] = 0.0e0;
  Jix[20] = 0.0e0;
  Jix[21] = 0.0e0;
  Jix[22] = 0.0e0;
  Jix[23] = 0.0e0;
  Jix[24] = Jix[12];
  Jix[25] = 0.0e0;
  Jix[26] = 0.0e0;
  Jix[27] = 0.0e0;
  Jix[28] = 0.0e0;
  Jix[29] = 0.0e0;
  Jix[30] = 0.0e0;
  Jix[31] = 0.0e0;
  Jix[32] = 0.0e0;
  Jix[33] = 0.0e0;
  Jix[34] = 0.0e0;
  Jix[35] = 0.0e0;
  Jix[36] = 0.0e0;
  Jix[37] = 0.0e0;
  Jix[38] = 0.0e0;
  Jix[39] = 0.0e0;
  Jix[40] = 0.0e0;
  Jix[41] = 0.0e0;
  Jix[42] = 0.0e0;
  Jix[43] = 0.0e0;
  Jix[44] = 0.0e0;
  Jix[45] = 0.0e0;
  Jix[46] = 0.0e0;
  Jix[47] = 0.0e0;
  Jix[48] = 0.0e0;
  Jix[49] = 0.0e0;
  Jix[50] = 0.0e0;
  Jix[51] = 0.0e0;
  Jix[52] = 0.0e0;
  Jix[53] = 0.0e0;
  Jix[54] = 0.0e0;
  Jix[55] = 0.0e0;
  Jix[56] = 0.0e0;
  Jix[57] = 0.0e0;
  Jix[58] = 0.0e0;
  Jix[59] = 0.0e0;
  Jix[60] = 0.0e0;
  Jix[61] = 0.0e0;
  Jix[62] = 0.0e0;
  Jix[63] = 0.0e0;
  Jix[64] = 0.0e0;
  Jix[65] = 0.0e0;
  Jix[66] = 0.0e0;
  Jix[67] = 0.0e0;
  Jix[68] = 0.0e0;
  Jix[69] = 0.0e0;
  Jix[70] = 0.0e0;
  Jix[71] = 0.0e0;
  Jix[72] = 0.0e0;
  Jix[73] = 0.0e0;
  Jix[74] = 0.0e0;
  Jiy[0] = 0.0e0;
  Jiy[1] = Jix[2];
  Jiy[2] = -t8 * t74 + t114;
  Jiy[3] = Fiy[0] * Jix[16] * t120 / 0.2e1 - t20 * t74;
  Jiy[4] = -t21 * t74;
  Jiy[5] = 0.0e0;
  Jiy[6] = Jix[7];
  Jiy[7] = -t127;
  Jiy[8] = Jix[13];
  Jiy[9] = 0.0e0;
  Jiy[10] = 0.1e1;
  Jiy[11] = Jix[24];
  Jiy[12] = 0.2e1 * Jiy[6] + Jix[11];
  Jiy[13] = t128 - t8 * Jix[16] * t74;
  Jiy[14] = Jix[9];
  Jiy[15] = 0.0e0;
  Jiy[16] = 0.0e0;
  Jiy[17] = Jix[16];
  Jiy[18] = Fiy[0] * gam * t3;
  Jiy[19] = 0.0e0;
  Jiy[20] = 0.0e0;
  Jiy[21] = 0.0e0;
  Jiy[22] = 0.0e0;
  Jiy[23] = 0.0e0;
  Jiy[24] = Jiy[6];
  Jiy[25] = 0.0e0;
  Jiy[26] = 0.0e0;
  Jiy[27] = 0.0e0;
  Jiy[28] = 0.0e0;
  Jiy[29] = 0.0e0;
  Jiy[30] = 0.0e0;
  Jiy[31] = 0.0e0;
  Jiy[32] = 0.0e0;
  Jiy[33] = 0.0e0;
  Jiy[34] = 0.0e0;
  Jiy[35] = 0.0e0;
  Jiy[36] = 0.0e0;
  Jiy[37] = 0.0e0;
  Jiy[38] = 0.0e0;
  Jiy[39] = 0.0e0;
  Jiy[40] = 0.0e0;
  Jiy[41] = 0.0e0;
  Jiy[42] = 0.0e0;
  Jiy[43] = 0.0e0;
  Jiy[44] = 0.0e0;
  Jiy[45] = 0.0e0;
  Jiy[46] = 0.0e0;
  Jiy[47] = 0.0e0;
  Jiy[48] = 0.0e0;
  Jiy[49] = 0.0e0;
  Jiy[50] = 0.0e0;
  Jiy[51] = 0.0e0;
  Jiy[52] = 0.0e0;
  Jiy[53] = 0.0e0;
  Jiy[54] = 0.0e0;
  Jiy[55] = 0.0e0;
  Jiy[56] = 0.0e0;
  Jiy[57] = 0.0e0;
  Jiy[58] = 0.0e0;
  Jiy[59] = 0.0e0;
  Jiy[60] = 0.0e0;
  Jiy[61] = 0.0e0;
  Jiy[62] = 0.0e0;
  Jiy[63] = 0.0e0;
  Jiy[64] = 0.0e0;
  Jiy[65] = 0.0e0;
  Jiy[66] = 0.0e0;
  Jiy[67] = 0.0e0;
  Jiy[68] = 0.0e0;
  Jiy[69] = 0.0e0;
  Jiy[70] = 0.0e0;
  Jiy[71] = 0.0e0;
  Jiy[72] = 0.0e0;
  Jiy[73] = 0.0e0;
  Jiy[74] = 0.0e0;
  Jvx[0] = 0.0e0;
  t145 = t24 * t119;
  t146 = 0.4e1 / 0.3e1 * t145;
  t147 = t26 * t74;
  t149 = t31 * t119;
  t151 = t33 * t74;
  t153 = t146 - 0.4e1 / 0.3e1 * t147 - 0.2e1 / 0.3e1 * t149 + 0.2e1 / 0.3e1 * t151;
  Jvx[1] = t153 * t50;
  t154 = t52 * t119;
  t155 = t54 * t74;
  t156 = t57 * t119;
  t157 = t59 * t74;
  t158 = t154 - t155 + t156 - t157;
  Jvx[2] = t158 * t50;
  t167 = Fix[0] * t119;
  t170 = t73 * t73;
  t171 = 0.1e1 / t170;
  t172 = t1 * t171;
  t174 = Fiy[0] * t119;
  t177 = t8 * t171;
  Jvx[3] = (t153 * Fix[0] * t3 - t62 * t74 + t158 * Fiy[0] * t3 - t64 * t74 + t67 * (t69 * t119 - t71 * t74 + 0.2e1 * t167 * t26 - t172 * t23 + 0.2e1 * t174 * t54 - t177 * t23)) * t50;
  Jvx[4] = t83 * t85 * t119 - t83 * t87 * t74;
  Jvx[5] = 0.0e0;
  t186 = t23 * t74;
  t187 = t186 * t50;
  Jvx[6] = -0.4e1 / 0.3e1 * t187;
  t189 = t30 * t74;
  t190 = t189 * t50;
  Jvx[7] = -t190;
  Jvx[8] = (-t146 + t36 * t3 - t149 + t67 * (t145 - t147)) * t50;
  Jvx[9] = 0.0e0;
  Jvx[10] = 0.0e0;
  Jvx[11] = 0.2e1 / 0.3e1 * t190;
  Jvx[12] = -t187;
  t196 = t61 * t3;
  Jvx[13] = (0.2e1 / 0.3e1 * t156 - t154 + t196 + t67 * (t154 - t155)) * t50;
  Jvx[14] = 0.0e0;
  Jvx[15] = 0.0e0;
  Jvx[16] = 0.0e0;
  Jvx[17] = 0.0e0;
  Jvx[18] = -t67 * t187;
  Jvx[19] = 0.0e0;
  Jvx[20] = 0.0e0;
  t204 = t39 * t39;
  t207 = t46 * t46;
  t211 = 0.4e1 * t43 * t47 - 0.3e1 * t38 * t37 * t204 * t39 / t207;
  Jvx[21] = t36 * t211;
  Jvx[22] = t61 * t211;
  Jvx[23] = t81 * t211;
  t212 = t211 * t82;
  Jvx[24] = t212 * t88 - t83 * t186;
  Jvx[25] = 0.0e0;
  t215 = t75 * t50;
  Jvx[26] = -0.4e1 / 0.3e1 * t215;
  t217 = t77 * t50;
  Jvx[27] = -t217;
  t218 = t1 * t119;
  t220 = t8 * t119;
  t223 = t67 * (-t6 * t74 + t218 + t220);
  Jvx[28] = (-0.4e1 / 0.3e1 * t218 - t220 + t223) * t50;
  Jvx[29] = -t83 * t17 * t74;
  Jvx[30] = 0.0e0;
  t227 = t3 * t50;
  Jvx[31] = 0.4e1 / 0.3e1 * t227;
  Jvx[32] = 0.0e0;
  t229 = t67 * t75;
  Jvx[33] = (0.4e1 / 0.3e1 * t75 - t229) * t50;
  Jvx[34] = 0.0e0;
  Jvx[35] = 0.0e0;
  Jvx[36] = 0.0e0;
  Jvx[37] = t227;
  t231 = t67 * t77;
  Jvx[38] = (t77 - t231) * t50;
  Jvx[39] = 0.0e0;
  Jvx[40] = 0.0e0;
  Jvx[41] = 0.0e0;
  Jvx[42] = 0.0e0;
  Jvx[43] = t67 * Jvx[37];
  Jvx[44] = 0.0e0;
  Jvx[45] = 0.0e0;
  Jvx[46] = 0.0e0;
  Jvx[47] = 0.0e0;
  Jvx[48] = 0.0e0;
  Jvx[49] = t83 * t3;
  Jvx[50] = 0.0e0;
  Jvx[51] = 0.2e1 / 0.3e1 * t217;
  Jvx[52] = -t215;
  Jvx[53] = -t174 * Fix[0] * t50 / 0.3e1;
  Jvx[54] = 0.0e0;
  Jvx[55] = 0.0e0;
  Jvx[56] = 0.0e0;
  Jvx[57] = Jvx[37];
  Jvx[58] = t217;
  Jvx[59] = 0.0e0;
  Jvx[60] = 0.0e0;
  Jvx[61] = -0.2e1 / 0.3e1 * Jvx[57];
  Jvx[62] = 0.0e0;
  t237 = 0.2e1 / 0.3e1 * t215;
  Jvx[63] = -t237;
  Jvx[64] = 0.0e0;
  Jvx[65] = 0.0e0;
  Jvx[66] = 0.0e0;
  Jvx[67] = 0.0e0;
  Jvx[68] = 0.0e0;
  Jvx[69] = 0.0e0;
  Jvx[70] = 0.0e0;
  Jvx[71] = 0.0e0;
  Jvx[72] = 0.0e0;
  Jvx[73] = 0.0e0;
  Jvx[74] = 0.0e0;
  Jvy[0] = 0.0e0;
  Jvy[1] = Jvx[2];
  t238 = 0.4e1 / 0.3e1 * t149;
  t242 = t238 - 0.4e1 / 0.3e1 * t151 - 0.2e1 / 0.3e1 * t145 + 0.2e1 / 0.3e1 * t147;
  Jvy[2] = t242 * t50;
  Jvy[3] = (t158 * Fix[0] * t3 - t92 * t74 + t242 * Fiy[0] * t3 - t94 * t74 + t67 * (t97 * t119 - t99 * t74 + 0.2e1 * t167 * t59 - t172 * t30 + 0.2e1 * t174 * t33 - t177 * t30)) * t50;
  Jvy[4] = t83 * t107 * t119 - t83 * t109 * t74;
  Jvy[5] = 0.0e0;
  Jvy[6] = Jvx[7];
  Jvy[7] = 0.2e1 / 0.3e1 * t187;
  Jvy[8] = (-t156 + t196 + 0.2e1 / 0.3e1 * t154 + t67 * (t156 - t157)) * t50;
  Jvy[9] = 0.0e0;
  Jvy[10] = 0.0e0;
  Jvy[11] = Jvx[12];
  Jvy[12] = -0.4e1 / 0.3e1 * t190;
  Jvy[13] = (-t145 - t238 + t91 * t3 + t67 * (t149 - t151)) * t50;
  Jvy[14] = 0.0e0;
  Jvy[15] = 0.0e0;
  Jvy[16] = 0.0e0;
  Jvy[17] = 0.0e0;
  Jvy[18] = -t67 * t190;
  Jvy[19] = 0.0e0;
  Jvy[20] = 0.0e0;
  Jvy[21] = Jvx[22];
  Jvy[22] = t91 * t211;
  Jvy[23] = t105 * t211;
  Jvy[24] = t212 * t110 - t83 * t189;
  Jvy[25] = 0.0e0;
  Jvy[26] = Jvx[27];
  Jvy[27] = t237;
  Jvy[28] = Jvx[53];
  Jvy[29] = 0.0e0;
  Jvy[30] = 0.0e0;
  Jvy[31] = 0.0e0;
  Jvy[32] = Jvx[61];
  Jvy[33] = -Jvx[51];
  Jvy[34] = 0.0e0;
  Jvy[35] = 0.0e0;
  Jvy[36] = Jvx[57];
  Jvy[37] = 0.0e0;
  Jvy[38] = t215;
  Jvy[39] = 0.0e0;
  Jvy[40] = 0.0e0;
  Jvy[41] = 0.0e0;
  Jvy[42] = 0.0e0;
  Jvy[43] = 0.0e0;
  Jvy[44] = 0.0e0;
  Jvy[45] = 0.0e0;
  Jvy[46] = 0.0e0;
  Jvy[47] = 0.0e0;
  Jvy[48] = 0.0e0;
  Jvy[49] = 0.0e0;
  Jvy[50] = 0.0e0;
  Jvy[51] = Jvx[52];
  Jvy[52] = -0.4e1 / 0.3e1 * Jvx[58];
  Jvy[53] = (-t218 - 0.4e1 / 0.3e1 * t220 + t223) * t50;
  Jvy[54] = Jvx[29];
  Jvy[55] = 0.0e0;
  Jvy[56] = Jvy[36];
  Jvy[57] = 0.0e0;
  Jvy[58] = (t75 - t229) * t50;
  Jvy[59] = 0.0e0;
  Jvy[60] = 0.0e0;
  Jvy[61] = 0.0e0;
  Jvy[62] = Jvx[31];
  Jvy[63] = (0.4e1 / 0.3e1 * t77 - t231) * t50;
  Jvy[64] = 0.0e0;
  Jvy[65] = 0.0e0;
  Jvy[66] = 0.0e0;
  Jvy[67] = 0.0e0;
  Jvy[68] = Jvx[43];
  Jvy[69] = 0.0e0;
  Jvy[70] = 0.0e0;
  Jvy[71] = 0.0e0;
  Jvy[72] = 0.0e0;
  Jvy[73] = 0.0e0;
  Jvy[74] = Jvx[49];
  return;
}
