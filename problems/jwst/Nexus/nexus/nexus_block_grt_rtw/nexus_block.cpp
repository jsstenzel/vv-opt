/*
 * nexus_block.cpp
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "nexus_block".
 *
 * Model version              : 1.259
 * Simulink Coder version : 24.2 (R2024b) 21-Jun-2024
 * C++ source code generated on : Sun Mar  2 14:58:31 2025
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Emulation hardware selection:
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "nexus_block.h"
#include "rtwtypes.h"
#include <cmath>
#include "nexus_block_private.h"
#include "cmath"

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
void nexus_block::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t { rtsiGetT(si) };

  time_T tnew { rtsiGetSolverStopTime(si) };

  time_T h { rtsiGetStepSize(si) };

  real_T *x { rtsiGetContStates(si) };

  ODE4_IntgData *id { static_cast<ODE4_IntgData *>(rtsiGetSolverData(si)) };

  real_T *y { id->y };

  real_T *f0 { id->f[0] };

  real_T *f1 { id->f[1] };

  real_T *f2 { id->f[2] };

  real_T *f3 { id->f[3] };

  real_T temp;
  int_T i;
  int_T nXc { 240 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  nexus_block_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  this->step();
  nexus_block_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  this->step();
  nexus_block_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  this->step();
  nexus_block_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T hi;
  uint32_T lo;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return static_cast<real_T>(*u) * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T si;
  real_T sr;
  real_T y;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = std::sqrt(-2.0 * std::log(si) / si) * sr;
  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (std::isinf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

/* Model step function */
void nexus_block::step()
{
  real_T rtb_NEXUSPlantDynamics[36];
  real_T rtb_ACSController[3];
  real_T rtb_Angles[3];
  real_T rtb_CryoNoiseFilters[3];
  real_T rtb_u[2];
  real_T GSNoise_CSTATE;
  real_T GSNoise_CSTATE_0;
  real_T GSNoise_CSTATE_1;
  real_T STNoise_CSTATE;
  real_T rtb_WhiteNoise;
  real_T u;
  int_T ci;
  int_T iy;
  int_T rtb_RW1_tmp;
  boolean_T tmp;
  if (rtmIsMajorTimeStep((&nexus_block_M))) {
    /* set solver stop time */
    if (!((&nexus_block_M)->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&(&nexus_block_M)->solverInfo, (((&nexus_block_M)
        ->Timing.clockTickH0 + 1) * (&nexus_block_M)->Timing.stepSize0 *
        4294967296.0));
    } else {
      rtsiSetSolverStopTime(&(&nexus_block_M)->solverInfo, (((&nexus_block_M)
        ->Timing.clockTick0 + 1) * (&nexus_block_M)->Timing.stepSize0 +
        (&nexus_block_M)->Timing.clockTickH0 * (&nexus_block_M)
        ->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep((&nexus_block_M))) {
    (&nexus_block_M)->Timing.t[0] = rtsiGetT(&(&nexus_block_M)->solverInfo);
  }

  /* StateSpace: '<Root>/NEXUS Plant Dynamics' */
  for (iy = 0; iy < 36; iy++) {
    rtb_WhiteNoise = 0.0;
    for (ci = 0; ci < 158; ci++) {
      rtb_WhiteNoise += nexus_block_P.Cp[ci * 36 + iy] *
        nexus_block_X.NEXUSPlantDynamics_CSTATE[ci];
    }

    rtb_NEXUSPlantDynamics[iy] = rtb_WhiteNoise;
  }

  /* End of StateSpace: '<Root>/NEXUS Plant Dynamics' */
  for (iy = 0; iy < 134; iy++) {
    /* Gain: '<Root>/WFE Sensitivity' */
    rtb_WhiteNoise = 0.0;
    for (ci = 0; ci < 30; ci++) {
      rtb_WhiteNoise += nexus_block_P.dwdu[134 * ci + iy] *
        rtb_NEXUSPlantDynamics[ci];
    }

    nexus_block_B.WFESensitivity[iy] = rtb_WhiteNoise;

    /* End of Gain: '<Root>/WFE Sensitivity' */
  }

  tmp = rtmIsMajorTimeStep((&nexus_block_M));
  if (tmp) {
    /* S-Function (sdspstatfcns): '<S4>/Mean' */
    nexus_block_DW.Mean_AccVal = nexus_block_B.WFESensitivity[0];
    iy = 1;

    /* S-Function (sdspstatfcns): '<S4>/Variance' */
    nexus_block_DW.Variance_AccVal = nexus_block_B.WFESensitivity[0];
    nexus_block_DW.Variance_SqData = nexus_block_B.WFESensitivity[0] *
      nexus_block_B.WFESensitivity[0];
    for (ci = 132; ci >= 0; ci--) {
      /* S-Function (sdspstatfcns): '<S4>/Mean' */
      nexus_block_DW.Mean_AccVal += nexus_block_B.WFESensitivity[133 - ci];

      /* S-Function (sdspstatfcns): '<S4>/Variance' */
      nexus_block_DW.Variance_AccVal += nexus_block_B.WFESensitivity[iy];
      nexus_block_DW.Variance_SqData += nexus_block_B.WFESensitivity[iy] *
        nexus_block_B.WFESensitivity[iy];
      iy++;
    }

    /* S-Function (sdspstatfcns): '<S4>/Mean' */
    rtb_WhiteNoise = nexus_block_DW.Mean_AccVal / 134.0;

    /* Signum: '<S4>/Sign' */
    u = rtb_WhiteNoise;

    /* Sum: '<S4>/Sum' incorporates:
     *  Fcn: '<S4>/Fcn2'
     *  S-Function (sdspstatfcns): '<S4>/Variance'
     */
    rtb_WhiteNoise = (nexus_block_DW.Variance_SqData -
                      nexus_block_DW.Variance_AccVal *
                      nexus_block_DW.Variance_AccVal / 134.0) / 133.0 +
      rt_powd_snf(rtb_WhiteNoise, 2.0);

    /* Signum: '<S4>/Sign' */
    if (std::isnan(u)) {
      STNoise_CSTATE = (rtNaN);
    } else if (u < 0.0) {
      STNoise_CSTATE = -1.0;
    } else {
      STNoise_CSTATE = (u > 0.0);
    }

    /* Math: '<S4>/sqrt'
     *
     * About '<S4>/sqrt':
     *  Operator: sqrt
     */
    if (rtb_WhiteNoise < 0.0) {
      rtb_WhiteNoise = -std::sqrt(std::abs(rtb_WhiteNoise));
    } else {
      rtb_WhiteNoise = std::sqrt(rtb_WhiteNoise);
    }

    /* Gain: '<S4>/Gain' incorporates:
     *  Math: '<S4>/sqrt'
     *  Product: '<S4>/Product'
     *  Signum: '<S4>/Sign'
     *
     * About '<S4>/sqrt':
     *  Operator: sqrt
     */
    rtb_WhiteNoise = STNoise_CSTATE * rtb_WhiteNoise * nexus_block_P.m2nm;
  }

  /* StateSpace: '<Root>/ACS Controller' */
  for (iy = 0; iy < 3; iy++) {
    rtb_WhiteNoise = 0.0;
    for (ci = 0; ci < 9; ci++) {
      rtb_WhiteNoise += nexus_block_P.Cca[ci * 3 + iy] *
        nexus_block_X.ACSController_CSTATE[ci];
    }

    rtb_ACSController[iy] = rtb_WhiteNoise;
  }

  /* End of StateSpace: '<Root>/ACS Controller' */
  if (tmp) {
    /* Gain: '<S6>/Output' incorporates:
     *  RandomNumber: '<S6>/White Noise'
     */
    nexus_block_B.Output = std::sqrt(nexus_block_P.BandLimitedWhiteNoise6_Cov) /
      0.031622776601683791 * nexus_block_DW.NextOutput;

    /* Gain: '<S7>/Output' incorporates:
     *  RandomNumber: '<S7>/White Noise'
     */
    nexus_block_B.Output_l = std::sqrt(nexus_block_P.BandLimitedWhiteNoise7_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_j;

    /* Gain: '<S8>/Output' incorporates:
     *  RandomNumber: '<S8>/White Noise'
     */
    nexus_block_B.Output_k = std::sqrt(nexus_block_P.BandLimitedWhiteNoise8_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_l;
  }

  /* StateSpace: '<S1>/ST Noise' */
  rtb_WhiteNoise = nexus_block_X.STNoise_CSTATE[1];
  u = nexus_block_X.STNoise_CSTATE[0];
  STNoise_CSTATE = nexus_block_X.STNoise_CSTATE[2];
  for (iy = 0; iy < 3; iy++) {
    rtb_Angles[iy] = (nexus_block_P.Cds[iy + 3] * rtb_WhiteNoise +
                      nexus_block_P.Cds[iy] * u) + nexus_block_P.Cds[iy + 6] *
      STNoise_CSTATE;
  }

  /* End of StateSpace: '<S1>/ST Noise' */

  /* SignalConversion generated from: '<S1>/ST Noise' */
  nexus_block_B.TmpSignalConversionAtSTNoiseInp[0] = nexus_block_B.Output;
  nexus_block_B.TmpSignalConversionAtSTNoiseInp[1] = nexus_block_B.Output_l;
  nexus_block_B.TmpSignalConversionAtSTNoiseInp[2] = nexus_block_B.Output_k;

  /* StateSpace: '<Root>/FSM Controller' */
  for (iy = 0; iy < 2; iy++) {
    rtb_WhiteNoise = 0.0;
    for (ci = 0; ci < 6; ci++) {
      rtb_WhiteNoise += -nexus_block_P.Ccf[(ci << 1) + iy] *
        nexus_block_X.FSMController_CSTATE[ci];
    }

    rtb_u[iy] = rtb_WhiteNoise;
  }

  /* End of StateSpace: '<Root>/FSM Controller' */

  /* Gain: '<Root>/FSM Plant' */
  rtb_WhiteNoise = nexus_block_P.FSMPlant_Gain[0] * rtb_u[0] + rtb_u[1] *
    nexus_block_P.FSMPlant_Gain[2];
  u = rtb_u[0] * nexus_block_P.FSMPlant_Gain[1] + rtb_u[1] *
    nexus_block_P.FSMPlant_Gain[3];
  for (iy = 0; iy < 2; iy++) {
    /* Sum: '<Root>/Centroid' incorporates:
     *  Gain: '<Root>/Centroid Sensitivity'
     *  Gain: '<Root>/FSM Coupling'
     *  Gain: '<Root>/FSM Plant'
     */
    STNoise_CSTATE = 0.0;
    for (ci = 0; ci < 30; ci++) {
      STNoise_CSTATE += nexus_block_P.dcdu[(ci << 1) + iy] *
        rtb_NEXUSPlantDynamics[ci];
    }

    rtb_u[iy] = STNoise_CSTATE - (nexus_block_P.Kfsm[iy + 2] * u +
      nexus_block_P.Kfsm[iy] * rtb_WhiteNoise);

    /* End of Sum: '<Root>/Centroid' */
  }

  if (tmp) {
    /* Gain: '<S9>/Output' incorporates:
     *  RandomNumber: '<S9>/White Noise'
     */
    nexus_block_B.Output_g = std::sqrt(nexus_block_P.BandLimitedWhiteNoise5_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_m;
  }

  /* StateSpace: '<S2>/Cryo Noise Filters' */
  for (iy = 0; iy < 3; iy++) {
    STNoise_CSTATE = 0.0;
    for (ci = 0; ci < 12; ci++) {
      STNoise_CSTATE += nexus_block_P.Cdc[ci * 3 + iy] *
        nexus_block_X.CryoNoiseFilters_CSTATE[ci];
    }

    rtb_CryoNoiseFilters[iy] = STNoise_CSTATE;
  }

  /* End of StateSpace: '<S2>/Cryo Noise Filters' */
  if (tmp) {
    /* Gain: '<S10>/Output' incorporates:
     *  RandomNumber: '<S10>/White Noise'
     */
    nexus_block_B.Output_j = std::sqrt(nexus_block_P.BandLimitedWhiteNoise10_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_i;

    /* Gain: '<S11>/Output' incorporates:
     *  RandomNumber: '<S11>/White Noise'
     */
    nexus_block_B.Output_d = std::sqrt(nexus_block_P.BandLimitedWhiteNoise9_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_h;
  }

  /* SignalConversion generated from: '<S3>/GS  Noise' */
  nexus_block_B.TmpSignalConversionAtGSNoiseInp[0] = nexus_block_B.Output_d;
  nexus_block_B.TmpSignalConversionAtGSNoiseInp[1] = nexus_block_B.Output_j;

  /* StateSpace: '<S3>/GS  Noise' */
  STNoise_CSTATE = nexus_block_X.GSNoise_CSTATE[1];
  GSNoise_CSTATE = nexus_block_X.GSNoise_CSTATE[0];
  GSNoise_CSTATE_0 = nexus_block_X.GSNoise_CSTATE[2];
  GSNoise_CSTATE_1 = nexus_block_X.GSNoise_CSTATE[3];
  for (iy = 0; iy < 2; iy++) {
    /* Sum: '<Root>/Measured Centroid' incorporates:
     *  Gain: '<S3>/Gain'
     */
    nexus_block_B.MeasuredCentroid[iy] = (((nexus_block_P.Cdg[iy + 2] *
      nexus_block_P.psc * STNoise_CSTATE + nexus_block_P.psc *
      nexus_block_P.Cdg[iy] * GSNoise_CSTATE) + nexus_block_P.Cdg[iy + 4] *
      nexus_block_P.psc * GSNoise_CSTATE_0) + nexus_block_P.Cdg[iy + 6] *
      nexus_block_P.psc * GSNoise_CSTATE_1) * nexus_block_P.Gain_Gain_n +
      rtb_u[iy];
  }

  /* End of StateSpace: '<S3>/GS  Noise' */
  if (tmp) {
    /* Gain: '<S12>/Output' incorporates:
     *  RandomNumber: '<S12>/White Noise'
     */
    nexus_block_B.Output_f = std::sqrt(nexus_block_P.BandLimitedWhiteNoise1_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_k;

    /* Gain: '<S13>/Output' incorporates:
     *  RandomNumber: '<S13>/White Noise'
     */
    nexus_block_B.Output_d5 = std::sqrt(nexus_block_P.BandLimitedWhiteNoise2_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_g;

    /* Gain: '<S14>/Output' incorporates:
     *  RandomNumber: '<S14>/White Noise'
     */
    nexus_block_B.Output_g0 = std::sqrt(nexus_block_P.BandLimitedWhiteNoise3_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_gq;

    /* Gain: '<S15>/Output' incorporates:
     *  RandomNumber: '<S15>/White Noise'
     */
    nexus_block_B.Output_jf = std::sqrt(nexus_block_P.BandLimitedWhiteNoise4_Cov)
      / 0.031622776601683791 * nexus_block_DW.NextOutput_ig;
  }

  /* SignalConversion generated from: '<Root>/ACS Controller' incorporates:
   *  Gain: '<Root>/FSM Plant'
   *  Gain: '<S1>/Gain'
   *  Sum: '<Root>/Angles'
   */
  nexus_block_B.TmpSignalConversionAtACSControl[0] = rtb_NEXUSPlantDynamics[33];
  nexus_block_B.TmpSignalConversionAtACSControl[3] = nexus_block_P.Gain_Gain *
    rtb_Angles[0] + rtb_NEXUSPlantDynamics[30];
  nexus_block_B.TmpSignalConversionAtACSControl[1] = rtb_NEXUSPlantDynamics[34];
  nexus_block_B.TmpSignalConversionAtACSControl[4] = nexus_block_P.Gain_Gain *
    rtb_Angles[1] + rtb_NEXUSPlantDynamics[31];
  nexus_block_B.TmpSignalConversionAtACSControl[2] = rtb_NEXUSPlantDynamics[35];
  nexus_block_B.TmpSignalConversionAtACSControl[5] = nexus_block_P.Gain_Gain *
    rtb_Angles[2] + rtb_NEXUSPlantDynamics[32];
  nexus_block_B.TmpSignalConversionAtACSControl[6] = rtb_WhiteNoise;
  nexus_block_B.TmpSignalConversionAtACSControl[7] = u;
  for (iy = 0; iy < 6; iy++) {
    /* StateSpace: '<S5>/RW1' */
    rtb_WhiteNoise = 0.0;

    /* StateSpace: '<S5>/RW2' */
    u = 0.0;

    /* StateSpace: '<S5>/RW3' */
    STNoise_CSTATE = 0.0;

    /* StateSpace: '<S5>/RW4' */
    GSNoise_CSTATE = 0.0;
    for (ci = 0; ci < 12; ci++) {
      /* StateSpace: '<S5>/RW1' incorporates:
       *  StateSpace: '<S5>/RW2'
       *  StateSpace: '<S5>/RW3'
       *  StateSpace: '<S5>/RW4'
       */
      rtb_RW1_tmp = ci * 6 + iy;
      rtb_WhiteNoise += nexus_block_P.RW1_C[rtb_RW1_tmp] *
        nexus_block_X.RW1_CSTATE[ci];

      /* StateSpace: '<S5>/RW2' */
      u += nexus_block_P.RW2_C[rtb_RW1_tmp] * nexus_block_X.RW2_CSTATE[ci];

      /* StateSpace: '<S5>/RW3' */
      STNoise_CSTATE += nexus_block_P.RW3_C[rtb_RW1_tmp] *
        nexus_block_X.RW3_CSTATE[ci];

      /* StateSpace: '<S5>/RW4' */
      GSNoise_CSTATE += nexus_block_P.RW4_C[rtb_RW1_tmp] *
        nexus_block_X.RW4_CSTATE[ci];
    }

    /* SignalConversion generated from: '<Root>/NEXUS Plant Dynamics' incorporates:
     *  Gain: '<S5>/Gain'
     *  StateSpace: '<S5>/RW1'
     *  StateSpace: '<S5>/RW2'
     *  StateSpace: '<S5>/RW3'
     *  StateSpace: '<S5>/RW4'
     */
    nexus_block_B.TmpSignalConversionAtNEXUSPlant[iy] =
      nexus_block_P.Gain_Gain_o * rtb_WhiteNoise;
    nexus_block_B.TmpSignalConversionAtNEXUSPlant[iy + 6] =
      nexus_block_P.Gain_Gain_o * u;
    nexus_block_B.TmpSignalConversionAtNEXUSPlant[iy + 12] =
      nexus_block_P.Gain_Gain_o * STNoise_CSTATE;
    nexus_block_B.TmpSignalConversionAtNEXUSPlant[iy + 18] =
      nexus_block_P.Gain_Gain_o * GSNoise_CSTATE;
  }

  /* SignalConversion generated from: '<Root>/NEXUS Plant Dynamics' incorporates:
   *  Gain: '<S2>/Gain'
   */
  nexus_block_B.TmpSignalConversionAtNEXUSPlant[24] = nexus_block_P.Gain_Gain_i *
    rtb_CryoNoiseFilters[0];
  nexus_block_B.TmpSignalConversionAtNEXUSPlant[27] = rtb_ACSController[0];
  nexus_block_B.TmpSignalConversionAtNEXUSPlant[25] = nexus_block_P.Gain_Gain_i *
    rtb_CryoNoiseFilters[1];
  nexus_block_B.TmpSignalConversionAtNEXUSPlant[28] = rtb_ACSController[1];
  nexus_block_B.TmpSignalConversionAtNEXUSPlant[26] = nexus_block_P.Gain_Gain_i *
    rtb_CryoNoiseFilters[2];
  nexus_block_B.TmpSignalConversionAtNEXUSPlant[29] = rtb_ACSController[2];
  if (rtmIsMajorTimeStep((&nexus_block_M))) {
    if (rtmIsMajorTimeStep((&nexus_block_M))) {
      /* Update for RandomNumber: '<S6>/White Noise' */
      nexus_block_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed) * nexus_block_P.WhiteNoise_StdDev +
        nexus_block_P.WhiteNoise_Mean;

      /* Update for RandomNumber: '<S7>/White Noise' */
      nexus_block_DW.NextOutput_j = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_a) * nexus_block_P.WhiteNoise_StdDev_k +
        nexus_block_P.WhiteNoise_Mean_p;

      /* Update for RandomNumber: '<S8>/White Noise' */
      nexus_block_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_e) * nexus_block_P.WhiteNoise_StdDev_b +
        nexus_block_P.WhiteNoise_Mean_o;

      /* Update for RandomNumber: '<S9>/White Noise' */
      nexus_block_DW.NextOutput_m = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_o) * nexus_block_P.WhiteNoise_StdDev_e +
        nexus_block_P.WhiteNoise_Mean_g;

      /* Update for RandomNumber: '<S10>/White Noise' */
      nexus_block_DW.NextOutput_i = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_c) * nexus_block_P.WhiteNoise_StdDev_f +
        nexus_block_P.WhiteNoise_Mean_c;

      /* Update for RandomNumber: '<S11>/White Noise' */
      nexus_block_DW.NextOutput_h = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_d) * nexus_block_P.WhiteNoise_StdDev_n +
        nexus_block_P.WhiteNoise_Mean_p1;

      /* Update for RandomNumber: '<S12>/White Noise' */
      nexus_block_DW.NextOutput_k = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_f) * nexus_block_P.WhiteNoise_StdDev_eg +
        nexus_block_P.WhiteNoise_Mean_k;

      /* Update for RandomNumber: '<S13>/White Noise' */
      nexus_block_DW.NextOutput_g = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_eu) * nexus_block_P.WhiteNoise_StdDev_g +
        nexus_block_P.WhiteNoise_Mean_l;

      /* Update for RandomNumber: '<S14>/White Noise' */
      nexus_block_DW.NextOutput_gq = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_k) * nexus_block_P.WhiteNoise_StdDev_bq +
        nexus_block_P.WhiteNoise_Mean_h;

      /* Update for RandomNumber: '<S15>/White Noise' */
      nexus_block_DW.NextOutput_ig = rt_nrand_Upu32_Yd_f_pw_snf
        (&nexus_block_DW.RandSeed_p) * nexus_block_P.WhiteNoise_StdDev_p +
        nexus_block_P.WhiteNoise_Mean_l3;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep((&nexus_block_M))) {
    rt_ertODEUpdateContinuousStates(&(&nexus_block_M)->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++(&nexus_block_M)->Timing.clockTick0)) {
      ++(&nexus_block_M)->Timing.clockTickH0;
    }

    (&nexus_block_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&nexus_block_M)
      ->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      (&nexus_block_M)->Timing.clockTick1++;
      if (!(&nexus_block_M)->Timing.clockTick1) {
        (&nexus_block_M)->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void nexus_block::nexus_block_derivatives()
{
  XDot_nexus_block_T *_rtXdot;
  real_T NEXUSPlantDynamics_CSTATE;
  real_T STNoise_CSTATE;
  real_T STNoise_CSTATE_0;
  real_T TmpSignalConversionAtSTNoiseI_0;
  real_T TmpSignalConversionAtSTNoiseI_1;
  real_T TmpSignalConversionAtSTNoiseInp;
  int_T ci;
  int_T is;
  _rtXdot = ((XDot_nexus_block_T *) (&nexus_block_M)->derivs);

  /* Derivatives for StateSpace: '<Root>/NEXUS Plant Dynamics' */
  for (is = 0; is < 158; is++) {
    NEXUSPlantDynamics_CSTATE = 0.0;
    for (ci = 0; ci < 158; ci++) {
      NEXUSPlantDynamics_CSTATE += nexus_block_P.Ap[ci * 158 + is] *
        nexus_block_X.NEXUSPlantDynamics_CSTATE[ci];
    }

    for (ci = 0; ci < 30; ci++) {
      NEXUSPlantDynamics_CSTATE += nexus_block_P.Bp[ci * 158 + is] *
        nexus_block_B.TmpSignalConversionAtNEXUSPlant[ci];
    }

    _rtXdot->NEXUSPlantDynamics_CSTATE[is] = NEXUSPlantDynamics_CSTATE;
  }

  /* End of Derivatives for StateSpace: '<Root>/NEXUS Plant Dynamics' */

  /* Derivatives for StateSpace: '<Root>/ACS Controller' */
  for (is = 0; is < 9; is++) {
    NEXUSPlantDynamics_CSTATE = 0.0;
    for (ci = 0; ci < 9; ci++) {
      NEXUSPlantDynamics_CSTATE += nexus_block_P.Aca[ci * 9 + is] *
        nexus_block_X.ACSController_CSTATE[ci];
    }

    for (ci = 0; ci < 8; ci++) {
      NEXUSPlantDynamics_CSTATE += nexus_block_P.Bca[ci * 9 + is] *
        nexus_block_B.TmpSignalConversionAtACSControl[ci];
    }

    _rtXdot->ACSController_CSTATE[is] = NEXUSPlantDynamics_CSTATE;
  }

  /* End of Derivatives for StateSpace: '<Root>/ACS Controller' */

  /* Derivatives for StateSpace: '<S1>/ST Noise' */
  NEXUSPlantDynamics_CSTATE = nexus_block_X.STNoise_CSTATE[1];
  STNoise_CSTATE = nexus_block_X.STNoise_CSTATE[0];
  STNoise_CSTATE_0 = nexus_block_X.STNoise_CSTATE[2];
  TmpSignalConversionAtSTNoiseInp =
    nexus_block_B.TmpSignalConversionAtSTNoiseInp[0];
  TmpSignalConversionAtSTNoiseI_0 =
    nexus_block_B.TmpSignalConversionAtSTNoiseInp[1];
  TmpSignalConversionAtSTNoiseI_1 =
    nexus_block_B.TmpSignalConversionAtSTNoiseInp[2];
  for (is = 0; is < 3; is++) {
    _rtXdot->STNoise_CSTATE[is] = ((((nexus_block_P.Ads[is + 3] *
      NEXUSPlantDynamics_CSTATE + nexus_block_P.Ads[is] * STNoise_CSTATE) +
      nexus_block_P.Ads[is + 6] * STNoise_CSTATE_0) + nexus_block_P.Bds[is] *
      TmpSignalConversionAtSTNoiseInp) + nexus_block_P.Bds[is + 3] *
      TmpSignalConversionAtSTNoiseI_0) + nexus_block_P.Bds[is + 6] *
      TmpSignalConversionAtSTNoiseI_1;
  }

  /* End of Derivatives for StateSpace: '<S1>/ST Noise' */

  /* Derivatives for StateSpace: '<Root>/FSM Controller' */
  for (is = 0; is < 6; is++) {
    NEXUSPlantDynamics_CSTATE = 0.0;
    for (ci = 0; ci < 6; ci++) {
      NEXUSPlantDynamics_CSTATE += nexus_block_P.Acf[ci * 6 + is] *
        nexus_block_X.FSMController_CSTATE[ci];
    }

    _rtXdot->FSMController_CSTATE[is] = (nexus_block_P.Bcf[is] *
      nexus_block_B.MeasuredCentroid[0] + NEXUSPlantDynamics_CSTATE) +
      nexus_block_P.Bcf[is + 6] * nexus_block_B.MeasuredCentroid[1];
  }

  /* End of Derivatives for StateSpace: '<Root>/FSM Controller' */

  /* Derivatives for StateSpace: '<S2>/Cryo Noise Filters' */
  for (is = 0; is < 12; is++) {
    NEXUSPlantDynamics_CSTATE = 0.0;
    for (ci = 0; ci < 12; ci++) {
      NEXUSPlantDynamics_CSTATE += nexus_block_P.Adc[ci * 12 + is] *
        nexus_block_X.CryoNoiseFilters_CSTATE[ci];
    }

    _rtXdot->CryoNoiseFilters_CSTATE[is] = nexus_block_P.Bdc[is] *
      nexus_block_B.Output_g + NEXUSPlantDynamics_CSTATE;
  }

  /* End of Derivatives for StateSpace: '<S2>/Cryo Noise Filters' */

  /* Derivatives for StateSpace: '<S3>/GS  Noise' */
  NEXUSPlantDynamics_CSTATE = nexus_block_X.GSNoise_CSTATE[1];
  STNoise_CSTATE = nexus_block_X.GSNoise_CSTATE[0];
  STNoise_CSTATE_0 = nexus_block_X.GSNoise_CSTATE[2];
  TmpSignalConversionAtSTNoiseInp = nexus_block_X.GSNoise_CSTATE[3];
  TmpSignalConversionAtSTNoiseI_0 =
    nexus_block_B.TmpSignalConversionAtGSNoiseInp[0];
  TmpSignalConversionAtSTNoiseI_1 =
    nexus_block_B.TmpSignalConversionAtGSNoiseInp[1];
  for (is = 0; is < 4; is++) {
    _rtXdot->GSNoise_CSTATE[is] = ((((nexus_block_P.Adg[is + 4] *
      NEXUSPlantDynamics_CSTATE + nexus_block_P.Adg[is] * STNoise_CSTATE) +
      nexus_block_P.Adg[is + 8] * STNoise_CSTATE_0) + nexus_block_P.Adg[is + 12]
      * TmpSignalConversionAtSTNoiseInp) + nexus_block_P.Bdg[is] *
      TmpSignalConversionAtSTNoiseI_0) + nexus_block_P.Bdg[is + 4] *
      TmpSignalConversionAtSTNoiseI_1;
  }

  /* End of Derivatives for StateSpace: '<S3>/GS  Noise' */
  for (is = 0; is < 12; is++) {
    /* Derivatives for StateSpace: '<S5>/RW1' */
    NEXUSPlantDynamics_CSTATE = 0.0;

    /* Derivatives for StateSpace: '<S5>/RW2' */
    STNoise_CSTATE = 0.0;

    /* Derivatives for StateSpace: '<S5>/RW3' */
    STNoise_CSTATE_0 = 0.0;

    /* Derivatives for StateSpace: '<S5>/RW4' */
    TmpSignalConversionAtSTNoiseInp = 0.0;
    for (ci = 0; ci < 12; ci++) {
      /* Derivatives for StateSpace: '<S5>/RW1' */
      TmpSignalConversionAtSTNoiseI_0 = nexus_block_P.Adw[ci * 12 + is];
      NEXUSPlantDynamics_CSTATE += TmpSignalConversionAtSTNoiseI_0 *
        nexus_block_X.RW1_CSTATE[ci];

      /* Derivatives for StateSpace: '<S5>/RW2' */
      STNoise_CSTATE += TmpSignalConversionAtSTNoiseI_0 *
        nexus_block_X.RW2_CSTATE[ci];

      /* Derivatives for StateSpace: '<S5>/RW3' */
      STNoise_CSTATE_0 += TmpSignalConversionAtSTNoiseI_0 *
        nexus_block_X.RW3_CSTATE[ci];

      /* Derivatives for StateSpace: '<S5>/RW4' */
      TmpSignalConversionAtSTNoiseInp += TmpSignalConversionAtSTNoiseI_0 *
        nexus_block_X.RW4_CSTATE[ci];
    }

    /* Derivatives for StateSpace: '<S5>/RW1' */
    TmpSignalConversionAtSTNoiseI_0 = nexus_block_P.Bdw[is];
    _rtXdot->RW1_CSTATE[is] = TmpSignalConversionAtSTNoiseI_0 *
      nexus_block_B.Output_f + NEXUSPlantDynamics_CSTATE;

    /* Derivatives for StateSpace: '<S5>/RW2' */
    _rtXdot->RW2_CSTATE[is] = TmpSignalConversionAtSTNoiseI_0 *
      nexus_block_B.Output_d5 + STNoise_CSTATE;

    /* Derivatives for StateSpace: '<S5>/RW3' */
    _rtXdot->RW3_CSTATE[is] = TmpSignalConversionAtSTNoiseI_0 *
      nexus_block_B.Output_g0 + STNoise_CSTATE_0;

    /* Derivatives for StateSpace: '<S5>/RW4' */
    _rtXdot->RW4_CSTATE[is] = TmpSignalConversionAtSTNoiseI_0 *
      nexus_block_B.Output_jf + TmpSignalConversionAtSTNoiseInp;
  }
}

/* Model initialize function */
void nexus_block::initialize()
{
  /* Registration code */
  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&(&nexus_block_M)->solverInfo, &(&nexus_block_M)
                          ->Timing.simTimeStep);
    rtsiSetTPtr(&(&nexus_block_M)->solverInfo, &rtmGetTPtr((&nexus_block_M)));
    rtsiSetStepSizePtr(&(&nexus_block_M)->solverInfo, &(&nexus_block_M)
                       ->Timing.stepSize0);
    rtsiSetdXPtr(&(&nexus_block_M)->solverInfo, &(&nexus_block_M)->derivs);
    rtsiSetContStatesPtr(&(&nexus_block_M)->solverInfo, (real_T **)
                         &(&nexus_block_M)->contStates);
    rtsiSetNumContStatesPtr(&(&nexus_block_M)->solverInfo, &(&nexus_block_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&nexus_block_M)->solverInfo,
      &(&nexus_block_M)->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&nexus_block_M)->solverInfo,
      &(&nexus_block_M)->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&nexus_block_M)->solverInfo,
      &(&nexus_block_M)->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&nexus_block_M)->solverInfo, (boolean_T**)
      &(&nexus_block_M)->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&nexus_block_M)->solverInfo, (&rtmGetErrorStatus
      ((&nexus_block_M))));
    rtsiSetRTModelPtr(&(&nexus_block_M)->solverInfo, (&nexus_block_M));
  }

  rtsiSetSimTimeStep(&(&nexus_block_M)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&nexus_block_M)->solverInfo, false);
  rtsiSetIsContModeFrozen(&(&nexus_block_M)->solverInfo, false);
  (&nexus_block_M)->intgData.y = (&nexus_block_M)->odeY;
  (&nexus_block_M)->intgData.f[0] = (&nexus_block_M)->odeF[0];
  (&nexus_block_M)->intgData.f[1] = (&nexus_block_M)->odeF[1];
  (&nexus_block_M)->intgData.f[2] = (&nexus_block_M)->odeF[2];
  (&nexus_block_M)->intgData.f[3] = (&nexus_block_M)->odeF[3];
  (&nexus_block_M)->contStates = ((X_nexus_block_T *) &nexus_block_X);
  (&nexus_block_M)->contStateDisabled = ((XDis_nexus_block_T *)
    &nexus_block_XDis);
  (&nexus_block_M)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&nexus_block_M)->solverInfo, static_cast<void *>
                    (&(&nexus_block_M)->intgData));
  rtsiSetSolverName(&(&nexus_block_M)->solverInfo,"ode4");
  rtmSetTPtr((&nexus_block_M), &(&nexus_block_M)->Timing.tArray[0]);
  (&nexus_block_M)->Timing.stepSize0 = 0.001;

  {
    real_T tmp;
    int32_T t;
    int_T is;
    uint32_T tseed;

    /* InitializeConditions for StateSpace: '<Root>/NEXUS Plant Dynamics' */
    for (is = 0; is < 158; is++) {
      nexus_block_X.NEXUSPlantDynamics_CSTATE[is] =
        nexus_block_P.NEXUSPlantDynamics_InitialCondi;
    }

    /* End of InitializeConditions for StateSpace: '<Root>/NEXUS Plant Dynamics' */

    /* InitializeConditions for StateSpace: '<Root>/ACS Controller' */
    for (is = 0; is < 9; is++) {
      nexus_block_X.ACSController_CSTATE[is] =
        nexus_block_P.ACSController_InitialCondition;
    }

    /* End of InitializeConditions for StateSpace: '<Root>/ACS Controller' */

    /* InitializeConditions for RandomNumber: '<S6>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[5]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed = tseed;
    nexus_block_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed) * nexus_block_P.WhiteNoise_StdDev +
      nexus_block_P.WhiteNoise_Mean;

    /* End of InitializeConditions for RandomNumber: '<S6>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S7>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[6]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_a = tseed;
    nexus_block_DW.NextOutput_j = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_a) * nexus_block_P.WhiteNoise_StdDev_k +
      nexus_block_P.WhiteNoise_Mean_p;

    /* End of InitializeConditions for RandomNumber: '<S7>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S8>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[7]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_e = tseed;
    nexus_block_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_e) * nexus_block_P.WhiteNoise_StdDev_b +
      nexus_block_P.WhiteNoise_Mean_o;

    /* End of InitializeConditions for RandomNumber: '<S8>/White Noise' */

    /* InitializeConditions for StateSpace: '<S1>/ST Noise' */
    nexus_block_X.STNoise_CSTATE[0] = nexus_block_P.STNoise_InitialCondition;
    nexus_block_X.STNoise_CSTATE[1] = nexus_block_P.STNoise_InitialCondition;
    nexus_block_X.STNoise_CSTATE[2] = nexus_block_P.STNoise_InitialCondition;

    /* InitializeConditions for StateSpace: '<Root>/FSM Controller' */
    for (is = 0; is < 6; is++) {
      nexus_block_X.FSMController_CSTATE[is] =
        nexus_block_P.FSMController_InitialCondition;
    }

    /* End of InitializeConditions for StateSpace: '<Root>/FSM Controller' */

    /* InitializeConditions for RandomNumber: '<S9>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[4]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_o = tseed;
    nexus_block_DW.NextOutput_m = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_o) * nexus_block_P.WhiteNoise_StdDev_e +
      nexus_block_P.WhiteNoise_Mean_g;

    /* End of InitializeConditions for RandomNumber: '<S9>/White Noise' */

    /* InitializeConditions for StateSpace: '<S2>/Cryo Noise Filters' */
    for (is = 0; is < 12; is++) {
      nexus_block_X.CryoNoiseFilters_CSTATE[is] =
        nexus_block_P.CryoNoiseFilters_InitialConditi;
    }

    /* End of InitializeConditions for StateSpace: '<S2>/Cryo Noise Filters' */

    /* InitializeConditions for RandomNumber: '<S10>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[9]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>(static_cast<
      uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_c = tseed;
    nexus_block_DW.NextOutput_i = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_c) * nexus_block_P.WhiteNoise_StdDev_f +
      nexus_block_P.WhiteNoise_Mean_c;

    /* End of InitializeConditions for RandomNumber: '<S10>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S11>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[8]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_d = tseed;
    nexus_block_DW.NextOutput_h = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_d) * nexus_block_P.WhiteNoise_StdDev_n +
      nexus_block_P.WhiteNoise_Mean_p1;

    /* End of InitializeConditions for RandomNumber: '<S11>/White Noise' */

    /* InitializeConditions for StateSpace: '<S3>/GS  Noise' */
    nexus_block_X.GSNoise_CSTATE[0] = nexus_block_P.GSNoise_InitialCondition;
    nexus_block_X.GSNoise_CSTATE[1] = nexus_block_P.GSNoise_InitialCondition;
    nexus_block_X.GSNoise_CSTATE[2] = nexus_block_P.GSNoise_InitialCondition;
    nexus_block_X.GSNoise_CSTATE[3] = nexus_block_P.GSNoise_InitialCondition;

    /* InitializeConditions for RandomNumber: '<S12>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[0]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_f = tseed;
    nexus_block_DW.NextOutput_k = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_f) * nexus_block_P.WhiteNoise_StdDev_eg +
      nexus_block_P.WhiteNoise_Mean_k;

    /* End of InitializeConditions for RandomNumber: '<S12>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S13>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[1]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_eu = tseed;
    nexus_block_DW.NextOutput_g = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_eu) * nexus_block_P.WhiteNoise_StdDev_g +
      nexus_block_P.WhiteNoise_Mean_l;

    /* End of InitializeConditions for RandomNumber: '<S13>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S14>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[2]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_k = tseed;
    nexus_block_DW.NextOutput_gq = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_k) * nexus_block_P.WhiteNoise_StdDev_bq +
      nexus_block_P.WhiteNoise_Mean_h;

    /* End of InitializeConditions for RandomNumber: '<S14>/White Noise' */

    /* InitializeConditions for RandomNumber: '<S15>/White Noise' */
    tmp = std::floor(nexus_block_P.Seed[3]);
    if (std::isnan(tmp) || std::isinf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = std::fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? static_cast<uint32_T>(-static_cast<int32_T>
      (static_cast<uint32_T>(-tmp))) : static_cast<uint32_T>(tmp);
    is = static_cast<int32_T>(tseed >> 16U);
    t = static_cast<int32_T>(tseed & 32768U);
    tseed = ((((tseed - (static_cast<uint32_T>(is) << 16U)) +
               static_cast<uint32_T>(t)) << 16U) + static_cast<uint32_T>(t)) +
      static_cast<uint32_T>(is);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    }

    nexus_block_DW.RandSeed_p = tseed;
    nexus_block_DW.NextOutput_ig = rt_nrand_Upu32_Yd_f_pw_snf
      (&nexus_block_DW.RandSeed_p) * nexus_block_P.WhiteNoise_StdDev_p +
      nexus_block_P.WhiteNoise_Mean_l3;

    /* End of InitializeConditions for RandomNumber: '<S15>/White Noise' */
    for (is = 0; is < 12; is++) {
      /* InitializeConditions for StateSpace: '<S5>/RW1' */
      nexus_block_X.RW1_CSTATE[is] = nexus_block_P.RW1_InitialCondition;

      /* InitializeConditions for StateSpace: '<S5>/RW2' */
      nexus_block_X.RW2_CSTATE[is] = nexus_block_P.RW2_InitialCondition;

      /* InitializeConditions for StateSpace: '<S5>/RW3' */
      nexus_block_X.RW3_CSTATE[is] = nexus_block_P.RW3_InitialCondition;

      /* InitializeConditions for StateSpace: '<S5>/RW4' */
      nexus_block_X.RW4_CSTATE[is] = nexus_block_P.RW4_InitialCondition;
    }
  }
}

/* Model terminate function */
void nexus_block::terminate()
{
  /* (no terminate code required) */
}

/* Constructor */
nexus_block::nexus_block() :
  nexus_block_B(),
  nexus_block_DW(),
  nexus_block_X(),
  nexus_block_XDis(),
  nexus_block_M()
{
  /* Currently there is no constructor body generated.*/
}

/* Destructor */
/* Currently there is no destructor body generated.*/
nexus_block::~nexus_block() = default;

/* Real-Time Model get method */
RT_MODEL_nexus_block_T * nexus_block::getRTM()
{
  return (&nexus_block_M);
}
