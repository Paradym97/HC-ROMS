/*
** svn $Id: np9km.h 585 2012-01-03 18:44:28Z arango $
*******************************************************************************
** Copyright (c) 2002-2012 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for North Pacific 9km exp.
**
** Application flag:   np9km
** Input script:       np9km.in
*/

#define NP_9km_CFSR_exp_new03

/*Options associated with momentum equations*/
#define UV_ADV                  /* active advection terms*/
#define UV_COR                  /* active Coriolis terms*/
#define UV_QDRAG                /* active quadratic bottom friction*/
#define UV_VIS4                 /* active biharmonic horizontal mixing*/
#undef  UV_SMAGORINSKY

/*Options associated with tracers equations*/
#define TS_MPDATA               /*recursive MPDATA 3D advection*/
#define TS_DIF4                 /*turn ON or OFF harmonic horizontal mixing*/
#define SALINITY                /*use if having salinity*/
#define NONLIN_EOS              /*using nonlinear equation of state*/

#undef  QCORRECTION            /*use if net heat flux correction*/
#define SCORRECTION            /*use if freshwater flux correction*/
#undef  SRELAXATION

/*Options for pressure gradient algorithm*/
#define DJ_GRADPS               /*use if splines density Jacobian*/

/*Options for horizontal mixing of momentum*/
#define VISC_GRID              /*use to scale viscosity coefficient by grid size*/
#define MIX_S_UV               /*use if mixing along constant S-surfaces*/

/*Options for horizontal mixing of tracers*/
#define DIFF_GRID              /*use to scale diffusion coefficients by grid size*/
#define MIX_GEO_TS             /*use if mixing along geopotential(constant depth)surfaces*/

/*Options for model configuration*/
#define SOLVE3D                /*solving 3D primitive equations*/
/*activate parabolic splines vertical derivatives*/
!#define SPLINES
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define RI_SPLINES

#define MASKING                /*use if land/sea masking*/
#define AVERAGES               /*use if writing out time-averaged data*/

#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV
#define AVERAGES

/*
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV
#undef AVERAGES_AKV
#undef AVERAGES_AKT
#undef AVERAGES_AKS
*/


/*surface forcing*/
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE_OUT
# define EMINUSP
# define SOLAR_SOURCE
# define COOL_SKIN
# undef DIURNAL_SRFLUX
#endif

#define WIND_MINUS_CURRENT 

#  undef ANA_SSFLUX
#  undef ANA_SMFLUX
#  undef ANA_STFLUX

/*Options for vertical mixing of momentum and tracers*/
/*
** BVF_MIXING
** GLS_MIXING
** MY25_MIXING
** LMD_MIXING
*/

/*turbulence mixing*/
#undef MY25_MIXING
# ifdef MY25_MIXING
#  define N2S2_HORAVG
#  define KANTHA_CLAYSON
# endif

#define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_BKPP
#  define LMD_DDMIX
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  define RI_SPLINES
# endif

# undef  GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  undef  CANUTO_A
#  define N2S2_HORAVG
# endif


/*Options for reading and processing of climatological fields*/
#define TCLIMATOLOGY            /*use if processing tracers climatology*/

/*Options to nudge climatology data*/
#define TCLM_NUDGING            /*use if nudging tracers climatology*/


#define ANA_BSFLUX
#define ANA_BTFLUX

/*define River flows*/
/*
!#define UV_PSOURCE
!#define TS_PSOURCE
*/


/*define NetCDF input/output OPTIONS*/
#define RST_SINGLE 
#define INLINE_2DIO
!#define PIO_LIB
!#define PNETCDF
#undef  PERFECT_RESTART         /* use to include perfect restart variables  */

!#define HCGCM /* by yman*/
#ifdef HCGCM
#define HCGCM_T
#define HCGCM_Q
#define HCGCM_U
#endif
!#define DISTRIBUTE



