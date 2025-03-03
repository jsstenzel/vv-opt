/*
 * nexus_block.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "nexus_block".
 *
 * Model version              : 1.265
 * Simulink Coder version : 24.2 (R2024b) 21-Jun-2024
 * C++ source code generated on : Sun Mar  2 17:48:36 2025
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Emulation hardware selection:
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef nexus_block_h_
#define nexus_block_h_
#include <stdlib.h>
#include <cmath>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "rt_nonfinite.h"
#include "nexus_block_types.h"

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rtGetNaN.h"

}

#include <cfloat>
#include <cstring>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

/* Block signals (default storage) */
struct B_nexus_block_T {
  real_T WFESensitivity[134];          /* '<Root>/WFE Sensitivity' */
  real_T m2mic[2];                     /* '<Root>/m2mic' */
  real_T Output;                       /* '<S6>/Output' */
  real_T Output_l;                     /* '<S7>/Output' */
  real_T Output_k;                     /* '<S8>/Output' */
  real_T TmpSignalConversionAtSTNoiseInp[3];
  real_T Output_g;                     /* '<S9>/Output' */
  real_T Output_j;                     /* '<S10>/Output' */
  real_T Output_d;                     /* '<S11>/Output' */
  real_T TmpSignalConversionAtGSNoiseInp[2];
  real_T MeasuredCentroid[2];          /* '<Root>/Measured Centroid' */
  real_T Output_f;                     /* '<S12>/Output' */
  real_T Output_d5;                    /* '<S13>/Output' */
  real_T Output_g0;                    /* '<S14>/Output' */
  real_T Output_jf;                    /* '<S15>/Output' */
  real_T TmpSignalConversionAtACSControl[8];
  real_T TmpSignalConversionAtNEXUSPlant[30];
};

/* Block states (default storage) for system '<Root>' */
struct DW_nexus_block_T {
  real_T Mean_AccVal;                  /* '<S4>/Mean' */
  real_T Variance_AccVal;              /* '<S4>/Variance' */
  real_T Variance_SqData;              /* '<S4>/Variance' */
  real_T NextOutput;                   /* '<S6>/White Noise' */
  real_T NextOutput_j;                 /* '<S7>/White Noise' */
  real_T NextOutput_l;                 /* '<S8>/White Noise' */
  real_T NextOutput_m;                 /* '<S9>/White Noise' */
  real_T NextOutput_i;                 /* '<S10>/White Noise' */
  real_T NextOutput_h;                 /* '<S11>/White Noise' */
  real_T NextOutput_k;                 /* '<S12>/White Noise' */
  real_T NextOutput_g;                 /* '<S13>/White Noise' */
  real_T NextOutput_gq;                /* '<S14>/White Noise' */
  real_T NextOutput_ig;                /* '<S15>/White Noise' */
  struct {
    void *LoggedData;
  } Performance1_PWORK;                /* '<Root>/Performance 1' */

  struct {
    void *LoggedData;
  } Performance2_PWORK;                /* '<Root>/Performance 2' */

  uint32_T RandSeed;                   /* '<S6>/White Noise' */
  uint32_T RandSeed_a;                 /* '<S7>/White Noise' */
  uint32_T RandSeed_e;                 /* '<S8>/White Noise' */
  uint32_T RandSeed_o;                 /* '<S9>/White Noise' */
  uint32_T RandSeed_c;                 /* '<S10>/White Noise' */
  uint32_T RandSeed_d;                 /* '<S11>/White Noise' */
  uint32_T RandSeed_f;                 /* '<S12>/White Noise' */
  uint32_T RandSeed_eu;                /* '<S13>/White Noise' */
  uint32_T RandSeed_k;                 /* '<S14>/White Noise' */
  uint32_T RandSeed_p;                 /* '<S15>/White Noise' */
};

/* Continuous states (default storage) */
struct X_nexus_block_T {
  real_T NEXUSPlantDynamics_CSTATE[158];/* '<Root>/NEXUS Plant Dynamics' */
  real_T FSMController_CSTATE[6];      /* '<Root>/FSM Controller' */
  real_T ACSController_CSTATE[9];      /* '<Root>/ACS Controller' */
  real_T STNoise_CSTATE[3];            /* '<S1>/ST Noise' */
  real_T CryoNoiseFilters_CSTATE[12];  /* '<S2>/Cryo Noise Filters' */
  real_T GSNoise_CSTATE[4];            /* '<S3>/GS  Noise' */
  real_T RW1_CSTATE[12];               /* '<S5>/RW1' */
  real_T RW2_CSTATE[12];               /* '<S5>/RW2' */
  real_T RW3_CSTATE[12];               /* '<S5>/RW3' */
  real_T RW4_CSTATE[12];               /* '<S5>/RW4' */
};

/* State derivatives (default storage) */
struct XDot_nexus_block_T {
  real_T NEXUSPlantDynamics_CSTATE[158];/* '<Root>/NEXUS Plant Dynamics' */
  real_T FSMController_CSTATE[6];      /* '<Root>/FSM Controller' */
  real_T ACSController_CSTATE[9];      /* '<Root>/ACS Controller' */
  real_T STNoise_CSTATE[3];            /* '<S1>/ST Noise' */
  real_T CryoNoiseFilters_CSTATE[12];  /* '<S2>/Cryo Noise Filters' */
  real_T GSNoise_CSTATE[4];            /* '<S3>/GS  Noise' */
  real_T RW1_CSTATE[12];               /* '<S5>/RW1' */
  real_T RW2_CSTATE[12];               /* '<S5>/RW2' */
  real_T RW3_CSTATE[12];               /* '<S5>/RW3' */
  real_T RW4_CSTATE[12];               /* '<S5>/RW4' */
};

/* State disabled  */
struct XDis_nexus_block_T {
  boolean_T NEXUSPlantDynamics_CSTATE[158];/* '<Root>/NEXUS Plant Dynamics' */
  boolean_T FSMController_CSTATE[6];   /* '<Root>/FSM Controller' */
  boolean_T ACSController_CSTATE[9];   /* '<Root>/ACS Controller' */
  boolean_T STNoise_CSTATE[3];         /* '<S1>/ST Noise' */
  boolean_T CryoNoiseFilters_CSTATE[12];/* '<S2>/Cryo Noise Filters' */
  boolean_T GSNoise_CSTATE[4];         /* '<S3>/GS  Noise' */
  boolean_T RW1_CSTATE[12];            /* '<S5>/RW1' */
  boolean_T RW2_CSTATE[12];            /* '<S5>/RW2' */
  boolean_T RW3_CSTATE[12];            /* '<S5>/RW3' */
  boolean_T RW4_CSTATE[12];            /* '<S5>/RW4' */
};

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
struct ODE4_IntgData {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
};

#endif

/* Parameters (default storage) */
struct P_nexus_block_T_ {
  real_T Aca[81];                      /* Variable: Aca
                                        * Referenced by: '<Root>/ACS Controller'
                                        */
  real_T Acf[36];                      /* Variable: Acf
                                        * Referenced by: '<Root>/FSM Controller'
                                        */
  real_T Adc[144];                     /* Variable: Adc
                                        * Referenced by: '<S2>/Cryo Noise Filters'
                                        */
  real_T Adg[16];                      /* Variable: Adg
                                        * Referenced by: '<S3>/GS  Noise'
                                        */
  real_T Ads[9];                       /* Variable: Ads
                                        * Referenced by: '<S1>/ST Noise'
                                        */
  real_T Adw[144];                     /* Variable: Adw
                                        * Referenced by:
                                        *   '<S5>/RW1'
                                        *   '<S5>/RW2'
                                        *   '<S5>/RW3'
                                        *   '<S5>/RW4'
                                        */
  real_T Ap[24964];                    /* Variable: Ap
                                        * Referenced by: '<Root>/NEXUS Plant Dynamics'
                                        */
  real_T Bca[72];                      /* Variable: Bca
                                        * Referenced by: '<Root>/ACS Controller'
                                        */
  real_T Bcf[12];                      /* Variable: Bcf
                                        * Referenced by: '<Root>/FSM Controller'
                                        */
  real_T Bdc[12];                      /* Variable: Bdc
                                        * Referenced by: '<S2>/Cryo Noise Filters'
                                        */
  real_T Bdg[8];                       /* Variable: Bdg
                                        * Referenced by: '<S3>/GS  Noise'
                                        */
  real_T Bds[9];                       /* Variable: Bds
                                        * Referenced by: '<S1>/ST Noise'
                                        */
  real_T Bdw[12];                      /* Variable: Bdw
                                        * Referenced by:
                                        *   '<S5>/RW1'
                                        *   '<S5>/RW2'
                                        *   '<S5>/RW3'
                                        *   '<S5>/RW4'
                                        */
  real_T Bp[4740];                     /* Variable: Bp
                                        * Referenced by: '<Root>/NEXUS Plant Dynamics'
                                        */
  real_T Cca[27];                      /* Variable: Cca
                                        * Referenced by: '<Root>/ACS Controller'
                                        */
  real_T Ccf[12];                      /* Variable: Ccf
                                        * Referenced by: '<Root>/FSM Controller'
                                        */
  real_T Cdc[36];                      /* Variable: Cdc
                                        * Referenced by: '<S2>/Cryo Noise Filters'
                                        */
  real_T Cdg[8];                       /* Variable: Cdg
                                        * Referenced by: '<S3>/GS  Noise'
                                        */
  real_T Cds[9];                       /* Variable: Cds
                                        * Referenced by: '<S1>/ST Noise'
                                        */
  real_T Cp[5688];                     /* Variable: Cp
                                        * Referenced by: '<Root>/NEXUS Plant Dynamics'
                                        */
  real_T Kfsm[4];                      /* Variable: Kfsm
                                        * Referenced by: '<Root>/FSM Coupling'
                                        */
  real_T Seed[10];                     /* Variable: Seed
                                        * Referenced by:
                                        *   '<S6>/White Noise'
                                        *   '<S7>/White Noise'
                                        *   '<S8>/White Noise'
                                        *   '<S9>/White Noise'
                                        *   '<S10>/White Noise'
                                        *   '<S11>/White Noise'
                                        *   '<S12>/White Noise'
                                        *   '<S13>/White Noise'
                                        *   '<S14>/White Noise'
                                        *   '<S15>/White Noise'
                                        */
  real_T dcdu[60];                     /* Variable: dcdu
                                        * Referenced by: '<Root>/Centroid Sensitivity'
                                        */
  real_T dwdu[4020];                   /* Variable: dwdu
                                        * Referenced by: '<Root>/WFE Sensitivity'
                                        */
  real_T m2micron;                     /* Variable: m2micron
                                        * Referenced by: '<Root>/m2mic'
                                        */
  real_T m2nm;                         /* Variable: m2nm
                                        * Referenced by: '<S4>/Gain'
                                        */
  real_T psc;                          /* Variable: psc
                                        * Referenced by: '<S3>/GS  Noise'
                                        */
  real_T BandLimitedWhiteNoise6_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise6_Cov
                                    * Referenced by: '<S6>/Output'
                                    */
  real_T BandLimitedWhiteNoise7_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise7_Cov
                                    * Referenced by: '<S7>/Output'
                                    */
  real_T BandLimitedWhiteNoise8_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise8_Cov
                                    * Referenced by: '<S8>/Output'
                                    */
  real_T BandLimitedWhiteNoise5_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise5_Cov
                                    * Referenced by: '<S9>/Output'
                                    */
  real_T BandLimitedWhiteNoise10_Cov;
                                  /* Mask Parameter: BandLimitedWhiteNoise10_Cov
                                   * Referenced by: '<S10>/Output'
                                   */
  real_T BandLimitedWhiteNoise9_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise9_Cov
                                    * Referenced by: '<S11>/Output'
                                    */
  real_T BandLimitedWhiteNoise1_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise1_Cov
                                    * Referenced by: '<S12>/Output'
                                    */
  real_T BandLimitedWhiteNoise2_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise2_Cov
                                    * Referenced by: '<S13>/Output'
                                    */
  real_T BandLimitedWhiteNoise3_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise3_Cov
                                    * Referenced by: '<S14>/Output'
                                    */
  real_T BandLimitedWhiteNoise4_Cov;
                                   /* Mask Parameter: BandLimitedWhiteNoise4_Cov
                                    * Referenced by: '<S15>/Output'
                                    */
  real_T NEXUSPlantDynamics_InitialCondi;/* Expression: 0
                                          * Referenced by: '<Root>/NEXUS Plant Dynamics'
                                          */
  real_T FSMController_InitialCondition;/* Expression: 0
                                         * Referenced by: '<Root>/FSM Controller'
                                         */
  real_T FSMPlant_Gain[4];             /* Expression: eye(2,2)
                                        * Referenced by: '<Root>/FSM Plant'
                                        */
  real_T ACSController_InitialCondition;/* Expression: 0
                                         * Referenced by: '<Root>/ACS Controller'
                                         */
  real_T WhiteNoise_Mean;              /* Expression: 0
                                        * Referenced by: '<S6>/White Noise'
                                        */
  real_T WhiteNoise_StdDev;            /* Computed Parameter: WhiteNoise_StdDev
                                        * Referenced by: '<S6>/White Noise'
                                        */
  real_T WhiteNoise_Mean_p;            /* Expression: 0
                                        * Referenced by: '<S7>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_k;         /* Computed Parameter: WhiteNoise_StdDev_k
                                       * Referenced by: '<S7>/White Noise'
                                       */
  real_T WhiteNoise_Mean_o;            /* Expression: 0
                                        * Referenced by: '<S8>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_b;         /* Computed Parameter: WhiteNoise_StdDev_b
                                       * Referenced by: '<S8>/White Noise'
                                       */
  real_T STNoise_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S1>/ST Noise'
                                        */
  real_T Gain_Gain;                    /* Expression: 1
                                        * Referenced by: '<S1>/Gain'
                                        */
  real_T WhiteNoise_Mean_g;            /* Expression: 0
                                        * Referenced by: '<S9>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_e;         /* Computed Parameter: WhiteNoise_StdDev_e
                                       * Referenced by: '<S9>/White Noise'
                                       */
  real_T CryoNoiseFilters_InitialConditi;/* Expression: 0
                                          * Referenced by: '<S2>/Cryo Noise Filters'
                                          */
  real_T Gain_Gain_i;                  /* Expression: 1
                                        * Referenced by: '<S2>/Gain'
                                        */
  real_T WhiteNoise_Mean_c;            /* Expression: 0
                                        * Referenced by: '<S10>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_f;         /* Computed Parameter: WhiteNoise_StdDev_f
                                       * Referenced by: '<S10>/White Noise'
                                       */
  real_T WhiteNoise_Mean_p1;           /* Expression: 0
                                        * Referenced by: '<S11>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_n;         /* Computed Parameter: WhiteNoise_StdDev_n
                                       * Referenced by: '<S11>/White Noise'
                                       */
  real_T GSNoise_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S3>/GS  Noise'
                                        */
  real_T Gain_Gain_n;                  /* Expression: 1
                                        * Referenced by: '<S3>/Gain'
                                        */
  real_T WhiteNoise_Mean_k;            /* Expression: 0
                                        * Referenced by: '<S12>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_eg;       /* Computed Parameter: WhiteNoise_StdDev_eg
                                      * Referenced by: '<S12>/White Noise'
                                      */
  real_T WhiteNoise_Mean_l;            /* Expression: 0
                                        * Referenced by: '<S13>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_g;         /* Computed Parameter: WhiteNoise_StdDev_g
                                       * Referenced by: '<S13>/White Noise'
                                       */
  real_T WhiteNoise_Mean_h;            /* Expression: 0
                                        * Referenced by: '<S14>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_bq;       /* Computed Parameter: WhiteNoise_StdDev_bq
                                      * Referenced by: '<S14>/White Noise'
                                      */
  real_T WhiteNoise_Mean_l3;           /* Expression: 0
                                        * Referenced by: '<S15>/White Noise'
                                        */
  real_T WhiteNoise_StdDev_p;         /* Computed Parameter: WhiteNoise_StdDev_p
                                       * Referenced by: '<S15>/White Noise'
                                       */
  real_T RW1_C[72];                    /* Expression: R1*Cdw
                                        * Referenced by: '<S5>/RW1'
                                        */
  real_T RW1_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S5>/RW1'
                                        */
  real_T RW2_C[72];                    /* Expression: R2*Cdw
                                        * Referenced by: '<S5>/RW2'
                                        */
  real_T RW2_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S5>/RW2'
                                        */
  real_T RW3_C[72];                    /* Expression: R3*Cdw
                                        * Referenced by: '<S5>/RW3'
                                        */
  real_T RW3_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S5>/RW3'
                                        */
  real_T RW4_C[72];                    /* Expression: R4*Cdw
                                        * Referenced by: '<S5>/RW4'
                                        */
  real_T RW4_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S5>/RW4'
                                        */
  real_T Gain_Gain_o;                  /* Expression: 1
                                        * Referenced by: '<S5>/Gain'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_nexus_block_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_nexus_block_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_nexus_block_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[240];
  real_T odeF[4][240];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tStart;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Class declaration for model nexus_block */
class nexus_block final
{
  /* public data and function members */
 public:
  /* Copy Constructor */
  nexus_block(nexus_block const&) = delete;

  /* Assignment Operator */
  nexus_block& operator= (nexus_block const&) & = delete;

  /* Move Constructor */
  nexus_block(nexus_block &&) = delete;

  /* Move Assignment Operator */
  nexus_block& operator= (nexus_block &&) = delete;

  /* Real-Time Model get method */
  RT_MODEL_nexus_block_T * getRTM();

  /* model start function */
  void start();

  /* Initial conditions function */
  void initialize();

  /* model step function */
  void step();

  /* model terminate function */
  static void terminate();

  /* Constructor */
  nexus_block();

  /* Destructor */
  ~nexus_block();

  /* private data and function members */
 private:
  /* Block signals */
  B_nexus_block_T nexus_block_B;

  /* Block states */
  DW_nexus_block_T nexus_block_DW;

  /* Tunable parameters */
  static P_nexus_block_T nexus_block_P;

  /* Block continuous states */
  X_nexus_block_T nexus_block_X;

  /* Block Continuous state disabled vector */
  XDis_nexus_block_T nexus_block_XDis;

  /* Global mass matrix */

  /* Continuous states update member function*/
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  /* Derivatives member function */
  void nexus_block_derivatives();

  /* Real-Time Model */
  RT_MODEL_nexus_block_T nexus_block_M;
};

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'nexus_block'
 * '<S1>'   : 'nexus_block/ACS Noise'
 * '<S2>'   : 'nexus_block/Cryo Noise'
 * '<S3>'   : 'nexus_block/GS Noise'
 * '<S4>'   : 'nexus_block/RMMS'
 * '<S5>'   : 'nexus_block/RWA Noise'
 * '<S6>'   : 'nexus_block/ACS Noise/Band-Limited White Noise6'
 * '<S7>'   : 'nexus_block/ACS Noise/Band-Limited White Noise7'
 * '<S8>'   : 'nexus_block/ACS Noise/Band-Limited White Noise8'
 * '<S9>'   : 'nexus_block/Cryo Noise/Band-Limited White Noise5'
 * '<S10>'  : 'nexus_block/GS Noise/Band-Limited White Noise10'
 * '<S11>'  : 'nexus_block/GS Noise/Band-Limited White Noise9'
 * '<S12>'  : 'nexus_block/RWA Noise/Band-Limited White Noise1'
 * '<S13>'  : 'nexus_block/RWA Noise/Band-Limited White Noise2'
 * '<S14>'  : 'nexus_block/RWA Noise/Band-Limited White Noise3'
 * '<S15>'  : 'nexus_block/RWA Noise/Band-Limited White Noise4'
 */
#endif                                 /* nexus_block_h_ */
