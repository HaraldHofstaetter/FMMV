/*
! +-----------------------------------------------+
! |        Generated by SPL Compiler 3.30         |
! +-----------------------------------------------+
! Command-line options: -l -B 32 -C -xprecision=double -xposw -xsubname=FFT_M2X_48 
! 
*/
/*
! The SPL Program: (compose (tensor (compose (..)(..))(I 12))(compose (T 48 12)(compose (..)(..))))
! node size: 48 X 48
*/


	void FFT_M2X_48(double *y, double *x) {
	int i0;
	double f0;
	double f1;
	double f2;
	double f3;
	double f4;
	double f5;
	double f6;
	double f7;
	double f8;
	double f9;
	double f10;
	double f11;
	double f12;
	double f13;
	double f14;
	double f15;
	double f16;
	double f17;
	double f18;
	double f19;
	double f20;
	double f21;
	double f22;
	double f23;
	double f24;
	double f25;
	double f26;
	double f27;
	double f28;
	double f29;
	double f30;
	double f31;
	double f32;
	double f33;
	double f34;
	double f35;
	double f36;
	double f37;
	double f38;
	double f39;
	double f40;
	double f41;
	double f42;
	double f43;
	double f44;
	double f45;
	double f46;
	double f47;
	double f48;
	double f49;
	double f50;
	double f51;
	double f54;
	double f55;
	double f56;
	double f57;
	double f58;
	double f59;
	double f64;
	double f65;
	double f66;
	double f67;
	double f70;
	double f71;
	double f72;
	double f73;
	double f74;
	double f75;
	double f80;
	double f81;
	double f82;
	double f83;
	double f86;
	double f87;
	double f88;
	double f89;
	double f90;
	double f91;
	double f96;
	double f97;
	double f98;
	double f99;
	double f102;
	double f103;
	double f104;
	double f105;
	double f106;
	double f107;
	double f112;
	double f113;
	double f114;
	double f115;
	double f118;
	double f119;
	double f120;
	double f121;
	double f122;
	double f123;
	double f124;
	double f125;
	  static double t31[96];
	  static double t32[96];
	  static double t33[96];
	  extern double T_48_12 [96];

	  /* deftemp t31 ( 96 ) */
	  /* deftemp t32 ( 96 ) */
	  /* deftemp t33 ( 96 ) */
	  t33[0] = x[0];
	  t33[1] = x[1];
	  t33[2] = x[8];
	  t33[3] = x[9];
	  t33[4] = x[16];
	  t33[5] = x[17];
	  t33[6] = x[24];
	  t33[7] = x[25];
	  t33[8] = x[32];
	  t33[9] = x[33];
	  t33[10] = x[40];
	  t33[11] = x[41];
	  t33[12] = x[48];
	  t33[13] = x[49];
	  t33[14] = x[56];
	  t33[15] = x[57];
	  t33[16] = x[64];
	  t33[17] = x[65];
	  t33[18] = x[72];
	  t33[19] = x[73];
	  t33[20] = x[80];
	  t33[21] = x[81];
	  t33[22] = x[88];
	  t33[23] = x[89];
	  t33[24] = x[2];
	  t33[25] = x[3];
	  t33[26] = x[10];
	  t33[27] = x[11];
	  t33[28] = x[18];
	  t33[29] = x[19];
	  t33[30] = x[26];
	  t33[31] = x[27];
	  t33[32] = x[34];
	  t33[33] = x[35];
	  t33[34] = x[42];
	  t33[35] = x[43];
	  t33[36] = x[50];
	  t33[37] = x[51];
	  t33[38] = x[58];
	  t33[39] = x[59];
	  t33[40] = x[66];
	  t33[41] = x[67];
	  t33[42] = x[74];
	  t33[43] = x[75];
	  t33[44] = x[82];
	  t33[45] = x[83];
	  t33[46] = x[90];
	  t33[47] = x[91];
	  t33[48] = x[4];
	  t33[49] = x[5];
	  t33[50] = x[12];
	  t33[51] = x[13];
	  t33[52] = x[20];
	  t33[53] = x[21];
	  t33[54] = x[28];
	  t33[55] = x[29];
	  t33[56] = x[36];
	  t33[57] = x[37];
	  t33[58] = x[44];
	  t33[59] = x[45];
	  t33[60] = x[52];
	  t33[61] = x[53];
	  t33[62] = x[60];
	  t33[63] = x[61];
	  t33[64] = x[68];
	  t33[65] = x[69];
	  t33[66] = x[76];
	  t33[67] = x[77];
	  t33[68] = x[84];
	  t33[69] = x[85];
	  t33[70] = x[92];
	  t33[71] = x[93];
	  t33[72] = x[6];
	  t33[73] = x[7];
	  t33[74] = x[14];
	  t33[75] = x[15];
	  t33[76] = x[22];
	  t33[77] = x[23];
	  t33[78] = x[30];
	  t33[79] = x[31];
	  t33[80] = x[38];
	  t33[81] = x[39];
	  t33[82] = x[46];
	  t33[83] = x[47];
	  t33[84] = x[54];
	  t33[85] = x[55];
	  t33[86] = x[62];
	  t33[87] = x[63];
	  t33[88] = x[70];
	  t33[89] = x[71];
	  t33[90] = x[78];
	  t33[91] = x[79];
	  t33[92] = x[86];
	  t33[93] = x[87];
	  t33[94] = x[94];
	  t33[95] = x[95];
	  for (i0 = 0; i0 < 4; i0++) {  
	    f0 = t33[24*i0+18] - t33[24*i0+6];
	    f1 = t33[24*i0+19] - t33[24*i0+7];
	    f2 = t33[24*i0+18] + t33[24*i0+6];
	    f3 = t33[24*i0+19] + t33[24*i0+7];
	    f4 = t33[24*i0] - t33[24*i0+12];
	    f5 = t33[24*i0+1] - t33[24*i0+13];
	    f6 = t33[24*i0] + t33[24*i0+12];
	    f7 = t33[24*i0+1] + t33[24*i0+13];
	    f8 = f6 - f2;
	    f9 = f7 - f3;
	    f10 = f6 + f2;
	    f11 = f7 + f3;
	    f12 = f4 + f1;
	    f13 = f5 - f0;
	    f14 = f4 - f1;
	    f15 = f5 + f0;
	    f16 = t33[24*i0+2] - t33[24*i0+14];
	    f17 = t33[24*i0+3] - t33[24*i0+15];
	    f18 = t33[24*i0+2] + t33[24*i0+14];
	    f19 = t33[24*i0+3] + t33[24*i0+15];
	    f20 = t33[24*i0+8] - t33[24*i0+20];
	    f21 = t33[24*i0+9] - t33[24*i0+21];
	    f22 = t33[24*i0+8] + t33[24*i0+20];
	    f23 = t33[24*i0+9] + t33[24*i0+21];
	    f24 = f22 - f18;
	    f25 = f23 - f19;
	    f26 = f22 + f18;
	    f27 = f23 + f19;
	    f28 = f20 + f17;
	    f29 = f21 - f16;
	    f30 = f20 - f17;
	    f31 = f21 + f16;
	    f32 = t33[24*i0+10] - t33[24*i0+22];
	    f33 = t33[24*i0+11] - t33[24*i0+23];
	    f34 = t33[24*i0+10] + t33[24*i0+22];
	    f35 = t33[24*i0+11] + t33[24*i0+23];
	    f36 = t33[24*i0+16] - t33[24*i0+4];
	    f37 = t33[24*i0+17] - t33[24*i0+5];
	    f38 = t33[24*i0+16] + t33[24*i0+4];
	    f39 = t33[24*i0+17] + t33[24*i0+5];
	    f40 = f38 - f34;
	    f41 = f39 - f35;
	    f42 = f38 + f34;
	    f43 = f39 + f35;
	    f44 = f36 + f33;
	    f45 = f37 - f32;
	    f46 = f36 - f33;
	    f47 = f37 + f32;
	    f48 = f26 - f42;
	    f49 = f27 - f43;
	    f50 = f26 + f42;
	    f51 = f27 + f43;
	    t32[24*i0] = f10 + f50;
	    t32[24*i0+1] = f11 + f51;
	    f54 = 5.0000000000000000e-01 * f50;
	    f55 = 5.0000000000000000e-01 * f51;
	    f56 = f10 - f54;
	    f57 = f11 - f55;
	    f58 = 8.6602540378443860e-01 * f49;
	    f59 = 8.6602540378443860e-01 * f48;
	    t32[24*i0+16] = f56 + f58;
	    t32[24*i0+17] = f57 - f59;
	    t32[24*i0+8] = f56 - f58;
	    t32[24*i0+9] = f57 + f59;
	    f64 = f30 - f46;
	    f65 = f31 - f47;
	    f66 = f30 + f46;
	    f67 = f31 + f47;
	    t32[24*i0+6] = f14 + f66;
	    t32[24*i0+7] = f15 + f67;
	    f70 = 5.0000000000000000e-01 * f66;
	    f71 = 5.0000000000000000e-01 * f67;
	    f72 = f14 - f70;
	    f73 = f15 - f71;
	    f74 = 8.6602540378443860e-01 * f65;
	    f75 = 8.6602540378443860e-01 * f64;
	    t32[24*i0+22] = f72 + f74;
	    t32[24*i0+23] = f73 - f75;
	    t32[24*i0+14] = f72 - f74;
	    t32[24*i0+15] = f73 + f75;
	    f80 = f24 - f40;
	    f81 = f25 - f41;
	    f82 = f24 + f40;
	    f83 = f25 + f41;
	    t32[24*i0+12] = f8 + f82;
	    t32[24*i0+13] = f9 + f83;
	    f86 = 5.0000000000000000e-01 * f82;
	    f87 = 5.0000000000000000e-01 * f83;
	    f88 = f8 - f86;
	    f89 = f9 - f87;
	    f90 = 8.6602540378443860e-01 * f81;
	    f91 = 8.6602540378443860e-01 * f80;
	    t32[24*i0+4] = f88 + f90;
	    t32[24*i0+5] = f89 - f91;
	    t32[24*i0+20] = f88 - f90;
	    t32[24*i0+21] = f89 + f91;
	    f96 = f28 - f44;
	    f97 = f29 - f45;
	    f98 = f28 + f44;
	    f99 = f29 + f45;
	    t32[24*i0+18] = f12 + f98;
	    t32[24*i0+19] = f13 + f99;
	    f102 = 5.0000000000000000e-01 * f98;
	    f103 = 5.0000000000000000e-01 * f99;
	    f104 = f12 - f102;
	    f105 = f13 - f103;
	    f106 = 8.6602540378443860e-01 * f97;
	    f107 = 8.6602540378443860e-01 * f96;
	    t32[24*i0+10] = f104 + f106;
	    t32[24*i0+11] = f105 - f107;
	    t32[24*i0+2] = f104 - f106;
	    t32[24*i0+3] = f105 + f107;
	    }
	  /* undeftemp t33 ( 96 ) */
	  for (i0 = 0; i0 < 48; i0++) {  
	    f112 = T_48_12[2*i0] * t32[2*i0];
	    f113 = T_48_12[2*i0+1] * t32[2*i0+1];
	    f114 = T_48_12[2*i0] * t32[2*i0+1];
	    f115 = T_48_12[2*i0+1] * t32[2*i0];
	    t31[2*i0] = f112 - f113;
	    t31[2*i0+1] = f114 + f115;
	    }
	  /* undeftemp t32 ( 96 ) */
	  for (i0 = 0; i0 < 12; i0++) {  
	    f118 = t31[2*i0+24] - t31[2*i0+72];
	    f119 = t31[2*i0+25] - t31[2*i0+73];
	    f120 = t31[2*i0+24] + t31[2*i0+72];
	    f121 = t31[2*i0+25] + t31[2*i0+73];
	    f122 = t31[2*i0] - t31[2*i0+48];
	    f123 = t31[2*i0+1] - t31[2*i0+49];
	    f124 = t31[2*i0] + t31[2*i0+48];
	    f125 = t31[2*i0+1] + t31[2*i0+49];
	    y[2*i0+48] = f124 - f120;
	    y[2*i0+49] = f125 - f121;
	    y[2*i0] = f124 + f120;
	    y[2*i0+1] = f125 + f121;
	    y[2*i0+72] = f122 + f119;
	    y[2*i0+73] = f123 - f118;
	    y[2*i0+24] = f122 - f119;
	    y[2*i0+25] = f123 + f118;
	    }
	  /* undeftemp t31 ( 96 ) */

	  /* undecl double t33[96]; */
	  /* undecl double t32[96]; */
	  /* undecl double t31[96]; */

	}
#include <math.h>
#ifndef NO_GLOBAL
	double T_48_12 [96];
#endif

	void init_FFT_M2X_48( void )
	{
	int i;
	int j;
	  extern double T_48_12 [96];
	  for (i=0; i<=3; i++) {
	    for (j=0; j<=11; j++) {
	      T_48_12[24*i+2*j] = (cos(6.2831853071795862e+00*i*j/48));
	      T_48_12[24*i+2*j+1] = (sin(6.2831853071795862e+00*i*j/48));
	    }
	  }

	}

