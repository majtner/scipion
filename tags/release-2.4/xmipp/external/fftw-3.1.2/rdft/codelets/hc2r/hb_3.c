/*
 * Copyright (c) 2003, 2006 Matteo Frigo
 * Copyright (c) 2003, 2006 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Sun Jul  2 16:30:59 EDT 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 3 -dif -name hb_3 -include hb.h */

/*
 * This function contains 16 FP additions, 14 FP multiplications,
 * (or, 6 additions, 4 multiplications, 10 fused multiply/add),
 * 27 stack variables, and 12 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.9 2006-02-12 23:34:12 athena Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.16 2006-02-12 23:34:12 athena Exp $
 */

#include "hb.h"

static const R *hb_3(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 4, MAKE_VOLATILE_STRIDE(ios)) {
	  E Tm, Tc, Tn, Th, Td, Ti, Tl, Te;
	  {
	       E T1, T5, T6, T7, T4, Tg, T2, T3;
	       T1 = rio[0];
	       T2 = rio[WS(ios, 1)];
	       T3 = iio[-WS(ios, 2)];
	       T5 = iio[0];
	       T6 = iio[-WS(ios, 1)];
	       T7 = rio[WS(ios, 2)];
	       T4 = T2 + T3;
	       Tg = T2 - T3;
	       {
		    E Tj, Tf, Tb, T8, Ta, Tk, To, T9;
		    Tj = W[2];
		    Tb = T7 + T6;
		    T8 = T6 - T7;
		    rio[0] = T1 + T4;
		    Ta = FNMS(KP500000000, T4, T1);
		    Tm = W[3];
		    iio[-WS(ios, 2)] = T5 + T8;
		    Tf = FNMS(KP500000000, T8, T5);
		    Tc = FNMS(KP866025403, Tb, Ta);
		    Tn = FMA(KP866025403, Tb, Ta);
		    T9 = W[0];
		    Th = FMA(KP866025403, Tg, Tf);
		    Tk = FNMS(KP866025403, Tg, Tf);
		    To = Tj * Tn;
		    Td = T9 * Tc;
		    Ti = T9 * Th;
		    Tl = Tj * Tk;
		    rio[WS(ios, 2)] = FNMS(Tm, Tk, To);
		    Te = W[1];
	       }
	  }
	  iio[0] = FMA(Tm, Tn, Tl);
	  iio[-WS(ios, 1)] = FMA(Te, Tc, Ti);
	  rio[WS(ios, 1)] = FNMS(Te, Th, Td);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 3},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 3, "hb_3", twinstr, &GENUS, {6, 4, 10, 0}, 0, 0, 0 };

void X(codelet_hb_3) (planner *p) {
     X(khc2hc_register) (p, hb_3, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -sign 1 -n 3 -dif -name hb_3 -include hb.h */

/*
 * This function contains 16 FP additions, 12 FP multiplications,
 * (or, 10 additions, 6 multiplications, 6 fused multiply/add),
 * 15 stack variables, and 12 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.9 2006-02-12 23:34:12 athena Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.16 2006-02-12 23:34:12 athena Exp $
 */

#include "hb.h"

static const R *hb_3(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 4, MAKE_VOLATILE_STRIDE(ios)) {
	  E T1, T4, Ta, Te, T5, T8, Tb, Tf;
	  {
	       E T2, T3, T6, T7;
	       T1 = rio[0];
	       T2 = rio[WS(ios, 1)];
	       T3 = iio[-WS(ios, 2)];
	       T4 = T2 + T3;
	       Ta = FNMS(KP500000000, T4, T1);
	       Te = KP866025403 * (T2 - T3);
	       T5 = iio[0];
	       T6 = rio[WS(ios, 2)];
	       T7 = iio[-WS(ios, 1)];
	       T8 = T6 - T7;
	       Tb = KP866025403 * (T6 + T7);
	       Tf = FMA(KP500000000, T8, T5);
	  }
	  rio[0] = T1 + T4;
	  iio[-WS(ios, 2)] = T5 - T8;
	  {
	       E Ti, Tk, Th, Tj;
	       Ti = Tf - Te;
	       Tk = Ta + Tb;
	       Th = W[2];
	       Tj = W[3];
	       iio[0] = FMA(Th, Ti, Tj * Tk);
	       rio[WS(ios, 2)] = FNMS(Tj, Ti, Th * Tk);
	  }
	  {
	       E Tc, Tg, T9, Td;
	       Tc = Ta - Tb;
	       Tg = Te + Tf;
	       T9 = W[0];
	       Td = W[1];
	       rio[WS(ios, 1)] = FNMS(Td, Tg, T9 * Tc);
	       iio[-WS(ios, 1)] = FMA(T9, Tg, Td * Tc);
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 3},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 3, "hb_3", twinstr, &GENUS, {10, 6, 6, 0}, 0, 0, 0 };

void X(codelet_hb_3) (planner *p) {
     X(khc2hc_register) (p, hb_3, &desc);
}
#endif				/* HAVE_FMA */
