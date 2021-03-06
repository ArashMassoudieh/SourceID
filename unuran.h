/* file automatically generated by unuran/scripts/merge_h.pl                 */

/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran.h                                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#  define __BEGIN_DECLS extern "C" {
#  define __END_DECLS }
#else
#  define __BEGIN_DECLS /* empty */
#  define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


#ifndef UNURAN_H_IN_SEEN
#define UNURAN_H_IN_SEEN

#include <limits.h>


#include <stdio.h>


#include <stdlib.h>


#include "unuran_config.h"


/*-----*/
/* <1> `unur_typedefs.h' */
#ifndef UNUR_TYPEDEFS_H_SEEN
#define UNUR_TYPEDEFS_H_SEEN
struct unur_distr;                       
typedef struct unur_distr UNUR_DISTR;
struct unur_par;                         
typedef struct unur_par   UNUR_PAR;
struct unur_gen;                         
typedef struct unur_gen   UNUR_GEN;
typedef double UNUR_FUNCT_CONT(double x, const struct unur_distr *distr);
typedef double UNUR_FUNCT_DISCR(int x, const struct unur_distr *distr);
typedef double UNUR_FUNCT_CVEC(const double *x, const struct unur_distr *distr);
typedef int UNUR_VFUNCT_CVEC(double *result, const double *x, const struct unur_distr *distr);
struct unur_slist;         
#endif  
/* end of `unur_typedefs.h' */
/*-----*/
/*-----*/
/* <1> `unur_uniform.h' */
#ifndef UNUR_UNIFORM_H_SEEN
#define UNUR_UNIFORM_H_SEEN
double unur_urng_MRG31k3p (void);
int unur_urng_MRG31k3p_seed (long seed);
int unur_urng_MRG31k3p_reset (void);
double unur_urng_fish (void);
int unur_urng_fish_seed (long seed);
int unur_urng_fish_reset (void);
double unur_urng_mstd (void);
int unur_urng_mstd_seed (long seed);
int unur_urng_mstd_reset (void);
#endif  
/* end of `unur_uniform.h' */
/*-----*/
/*-----*/
/* <1> `x_urng.h' */
#ifndef X_URNG_H_SEEN
#define X_URNG_H_SEEN
UNUR_URNG *unur_get_default_urng( void );
UNUR_URNG *unur_set_default_urng( UNUR_URNG *urng_new );
UNUR_URNG *unur_set_default_urng_aux( UNUR_URNG *urng_new );
UNUR_URNG *unur_get_default_urng_aux( void );
int unur_set_urng( UNUR_PAR *parameters, UNUR_URNG *urng );
UNUR_URNG *unur_chg_urng( UNUR_GEN *generator, UNUR_URNG *urng );
UNUR_URNG *unur_get_urng( UNUR_GEN *generator );
int unur_set_urng_aux( UNUR_PAR *parameters, UNUR_URNG *urng_aux );
int unur_use_urng_aux_default( UNUR_PAR *parameters );
int unur_chgto_urng_aux_default( UNUR_GEN *generator );
UNUR_URNG *unur_chg_urng_aux( UNUR_GEN *generator, UNUR_URNG *urng_aux );
UNUR_URNG *unur_get_urng_aux( UNUR_GEN *generator );
#endif  
/* end of `x_urng.h' */
/*-----*/
/*-----*/
/* <1> `distr.h' */
enum {
  UNUR_DISTR_CONT  = 0x010u,      
  UNUR_DISTR_CEMP  = 0x011u,      
  UNUR_DISTR_CVEC  = 0x110u,      
  UNUR_DISTR_CVEMP = 0x111u,      
  UNUR_DISTR_DISCR = 0x020u       
};
void unur_distr_free( UNUR_DISTR *distribution );
int unur_distr_set_name( UNUR_DISTR *distribution, const char *name );
const char *unur_distr_get_name( const UNUR_DISTR *distribution );
int unur_distr_get_dim( const UNUR_DISTR *distribution );
unsigned int unur_distr_get_type( const UNUR_DISTR *distribution );
int unur_distr_is_cont( const UNUR_DISTR *distribution );
int unur_distr_is_cvec( const UNUR_DISTR *distribution );
int unur_distr_is_cemp( const UNUR_DISTR *distribution );
int unur_distr_is_cvemp( const UNUR_DISTR *distribution );
int unur_distr_is_discr( const UNUR_DISTR *distribution );
/* end of `distr.h' */
/*-----*/
/*-----*/
/* <1> `cemp.h' */
UNUR_DISTR *unur_distr_cemp_new( void );
int unur_distr_cemp_set_data( UNUR_DISTR *distribution, const double *sample, int n_sample );
int unur_distr_cemp_read_data( UNUR_DISTR *distribution, const char *filename );
int unur_distr_cemp_get_data( const UNUR_DISTR *distribution, const double **sample );
/* end of `cemp.h' */
/*-----*/
/*-----*/
/* <1> `cont.h' */
UNUR_DISTR *unur_distr_cont_new( void );
int unur_distr_cont_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *pdf );
int unur_distr_cont_set_dpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *dpdf );
int unur_distr_cont_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *cdf );
UNUR_FUNCT_CONT *unur_distr_cont_get_pdf( const UNUR_DISTR *distribution );
UNUR_FUNCT_CONT *unur_distr_cont_get_dpdf( const UNUR_DISTR *distribution );
UNUR_FUNCT_CONT *unur_distr_cont_get_cdf( const UNUR_DISTR *distribution );
double unur_distr_cont_eval_pdf( double x, const UNUR_DISTR *distribution );
double unur_distr_cont_eval_dpdf( double x, const UNUR_DISTR *distribution );
double unur_distr_cont_eval_cdf( double x, const UNUR_DISTR *distribution );
int unur_distr_cont_set_pdfstr( UNUR_DISTR *distribution, const char *pdfstr );
int unur_distr_cont_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
char *unur_distr_cont_get_pdfstr( const UNUR_DISTR *distribution );
char *unur_distr_cont_get_dpdfstr( const UNUR_DISTR *distribution );
char *unur_distr_cont_get_cdfstr( const UNUR_DISTR *distribution );
int unur_distr_cont_set_pdfparams( UNUR_DISTR *distribution, const double *params, int n_params );
int unur_distr_cont_get_pdfparams( const UNUR_DISTR *distribution, const double **params );
int unur_distr_cont_set_domain( UNUR_DISTR *distribution, double left, double right );
int unur_distr_cont_get_domain( const UNUR_DISTR *distribution, double *left, double *right );
int unur_distr_cont_get_truncated( const UNUR_DISTR *distribution, double *left, double *right );
int unur_distr_cont_set_hr( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *hazard );
UNUR_FUNCT_CONT *unur_distr_cont_get_hr( const UNUR_DISTR *distribution );
double unur_distr_cont_eval_hr( double x, const UNUR_DISTR *distribution );
int unur_distr_cont_set_hrstr( UNUR_DISTR *distribution, const char *hrstr );
char *unur_distr_cont_get_hrstr( const UNUR_DISTR *distribution );
int unur_distr_cont_set_mode( UNUR_DISTR *distribution, double mode );
int unur_distr_cont_upd_mode( UNUR_DISTR *distribution );
double unur_distr_cont_get_mode( UNUR_DISTR *distribution );
int unur_distr_cont_set_pdfarea( UNUR_DISTR *distribution, double area );
int unur_distr_cont_upd_pdfarea( UNUR_DISTR *distribution );
double unur_distr_cont_get_pdfarea( UNUR_DISTR *distribution );
/* end of `cont.h' */
/*-----*/
/*-----*/
/* <1> `corder.h' */
UNUR_DISTR *unur_distr_corder_new( const UNUR_DISTR *distribution, int n, int k );
const UNUR_DISTR *unur_distr_corder_get_distribution( const UNUR_DISTR *distribution );
int unur_distr_corder_set_rank( UNUR_DISTR *distribution, int n, int k );
int unur_distr_corder_get_rank( const UNUR_DISTR *distribution, int *n, int *k );
#define unur_distr_corder_get_pdf(distr)   unur_distr_cont_get_pdf((distr))
#define unur_distr_corder_get_dpdf(distr)  unur_distr_cont_get_dpdf((distr))
#define unur_distr_corder_get_cdf(distr)   unur_distr_cont_get_cdf((distr))
#define unur_distr_corder_eval_pdf(x,distr)  unur_distr_cont_eval_pdf((x),(distr))
#define unur_distr_corder_eval_dpdf(x,distr) unur_distr_cont_eval_dpdf((x),(distr))
#define unur_distr_corder_eval_cdf(x,distr)  unur_distr_cont_eval_cdf((x),(distr))
#define unur_distr_corder_set_pdfparams(distr,params,n)  unur_distr_cont_set_pdfparams((distr),(params),(n))
#define unur_distr_corder_get_pdfparams(distr,params)  unur_distr_cont_get_pdfparams((distr),(params))
#define unur_distr_corder_set_domain(distr,left,right)  unur_distr_cont_set_domain((distr),(left),(right))
#define unur_distr_corder_get_domain(distr,left,right)  unur_distr_cont_get_domain((distr),(left),(right))
#define unur_distr_corder_get_truncated(distr,left,right)  unur_distr_cont_get_truncated((distr),(left),(right))
#define unur_distr_corder_set_mode(distr,mode)   unur_distr_cont_set_mode((distr),(mode))
#define unur_distr_corder_upd_mode(distr)   unur_distr_cont_upd_mode((distr))
#define unur_distr_corder_get_mode(distr)   unur_distr_cont_get_mode((distr))
#define unur_distr_corder_set_pdfarea(distr,area)   unur_distr_cont_set_pdfarea((distr),(area))
#define unur_distr_corder_upd_pdfarea(distr)   unur_distr_cont_upd_pdfarea((distr))
#define unur_distr_corder_get_pdfarea(distr)   unur_distr_cont_get_pdfarea((distr))
/* end of `corder.h' */
/*-----*/
/*-----*/
/* <1> `cvec.h' */
UNUR_DISTR *unur_distr_cvec_new( int dim );
int unur_distr_cvec_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *pdf );
int unur_distr_cvec_set_dpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dpdf );
UNUR_FUNCT_CVEC *unur_distr_cvec_get_pdf( const UNUR_DISTR *distribution );
UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dpdf( const UNUR_DISTR *distribution );
double unur_distr_cvec_eval_pdf( const double *x, const UNUR_DISTR *distribution );
int unur_distr_cvec_eval_dpdf( double *result, const double *x, const UNUR_DISTR *distribution );
int unur_distr_cvec_set_mean( UNUR_DISTR *distribution, const double *mean );
const double *unur_distr_cvec_get_mean( const UNUR_DISTR *distribution );
int unur_distr_cvec_set_covar( UNUR_DISTR *distribution, const double *covar );
const double *unur_distr_cvec_get_covar( const UNUR_DISTR *distribution );
int unur_distr_cvec_set_pdfparams( UNUR_DISTR *distribution, int par, const double *params, int n_params );
int unur_distr_cvec_get_pdfparams( const UNUR_DISTR *distribution, int par, const double **params );
int unur_distr_cvec_set_mode( UNUR_DISTR *distribution, const double *mode );
const double *unur_distr_cvec_get_mode( const UNUR_DISTR *distribution );
int unur_distr_cvec_set_pdfvol( UNUR_DISTR *distribution, double volume );
double unur_distr_cvec_get_pdfvol( const UNUR_DISTR *distribution );
/* end of `cvec.h' */
/*-----*/
/*-----*/
/* <1> `cvemp.h' */
UNUR_DISTR *unur_distr_cvemp_new( int dim ); 
int unur_distr_cvemp_set_data( UNUR_DISTR *distribution, const double *sample, int n_sample );
int unur_distr_cvemp_read_data( UNUR_DISTR *distribution, const char *filename );
int unur_distr_cvemp_get_data( const UNUR_DISTR *distribution, const double **sample );
/* end of `cvemp.h' */
/*-----*/
/*-----*/
/* <1> `discr.h' */
UNUR_DISTR *unur_distr_discr_new( void );
int unur_distr_discr_set_pv( UNUR_DISTR *distribution, const double *pv, int n_pv );
int unur_distr_discr_make_pv( UNUR_DISTR *distribution );
int unur_distr_discr_get_pv( const UNUR_DISTR *distribution, const double **pv );
int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( const UNUR_DISTR *distribution );
UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( const UNUR_DISTR *distribution );
double unur_distr_discr_eval_pv(int k, const UNUR_DISTR *distribution );
double unur_distr_discr_eval_pmf( int k, const UNUR_DISTR *distribution );
double unur_distr_discr_eval_cdf( int k, const UNUR_DISTR *distribution );
int unur_distr_discr_set_pmfstr( UNUR_DISTR *distribution, const char *pmfstr );
int unur_distr_discr_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
char *unur_distr_discr_get_pmfstr( const UNUR_DISTR *distribution );
char *unur_distr_discr_get_cdfstr( const UNUR_DISTR *distribution );
int unur_distr_discr_set_pmfparams( UNUR_DISTR *distribution, const double *params, int n_params );
int unur_distr_discr_get_pmfparams( const UNUR_DISTR *distribution, const double **params );
int unur_distr_discr_set_domain( UNUR_DISTR *distribution, int left, int right );
int unur_distr_discr_get_domain( const UNUR_DISTR *distribution, int *left, int *right );
int unur_distr_discr_set_mode( UNUR_DISTR *distribution, int mode );
int unur_distr_discr_upd_mode( UNUR_DISTR *distribution );
int unur_distr_discr_get_mode( UNUR_DISTR *distribution );
int unur_distr_discr_set_pmfsum( UNUR_DISTR *distribution, double sum );
int unur_distr_discr_upd_pmfsum( UNUR_DISTR *distribution );
double unur_distr_discr_get_pmfsum( UNUR_DISTR *distribution );
/* end of `discr.h' */
/*-----*/
/*-----*/
/* <1> `auto.h' */
UNUR_PAR *unur_auto_new( const UNUR_DISTR *distribution );
int unur_auto_set_logss( UNUR_PAR *parameters, int logss );
/* end of `auto.h' */
/*-----*/
/*-----*/
/* <1> `dari.h' */
UNUR_PAR *unur_dari_new( const UNUR_DISTR *distribution );
int unur_dari_reinit( UNUR_GEN *generator );
int unur_dari_set_squeeze( UNUR_PAR *parameters, int squeeze );
int unur_dari_set_tablesize( UNUR_PAR *parameters, int size );
int unur_dari_set_cpfactor( UNUR_PAR *parameters, double cp_factor );
int unur_dari_set_verify( UNUR_PAR *parameters, int verify );
int unur_dari_chg_verify( UNUR_GEN *generator, int verify );
int unur_dari_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_dari_chg_domain( UNUR_GEN *generator, int left, int right );
int unur_dari_chg_mode( UNUR_GEN *generator, int mode );
int unur_dari_upd_mode( UNUR_GEN *generator );
int unur_dari_chg_pmfsum( UNUR_GEN *generator, double sum );
int unur_dari_upd_pmfsum( UNUR_GEN *generator );
/* end of `dari.h' */
/*-----*/
/*-----*/
/* <1> `dau.h' */
UNUR_PAR *unur_dau_new( const UNUR_DISTR *distribution );
int unur_dau_set_urnfactor( UNUR_PAR *parameters, double factor );
/* end of `dau.h' */
/*-----*/
/*-----*/
/* <1> `dgt.h' */
UNUR_PAR *unur_dgt_new( const UNUR_DISTR *distribution );
int unur_dgt_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_dgt_set_variant( UNUR_PAR *parameters, unsigned variant );
/* end of `dgt.h' */
/*-----*/
/*-----*/
/* <1> `dsrou.h' */
UNUR_PAR *unur_dsrou_new( const UNUR_DISTR *distribution );
int unur_dsrou_reinit( UNUR_GEN *generator );
int unur_dsrou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
int unur_dsrou_set_verify( UNUR_PAR *parameters, int verify );
int unur_dsrou_chg_verify( UNUR_GEN *generator, int verify );
int unur_dsrou_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_dsrou_chg_domain( UNUR_GEN *generator, int left, int right );
int unur_dsrou_chg_mode( UNUR_GEN *generator, int mode );
int unur_dsrou_upd_mode( UNUR_GEN *generator );
int unur_dsrou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
int unur_dsrou_chg_pmfsum( UNUR_GEN *generator, double sum );
int unur_dsrou_upd_pmfsum( UNUR_GEN *generator );
/* end of `dsrou.h' */
/*-----*/
/*-----*/
/* <1> `arou.h' */
UNUR_PAR *unur_arou_new( const UNUR_DISTR *distribution );
int unur_arou_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
double unur_arou_get_sqhratio( const UNUR_GEN *generator );
double unur_arou_get_hatarea( const UNUR_GEN *generator );
double unur_arou_get_squeezearea( const UNUR_GEN *generator );
int unur_arou_set_max_segments( UNUR_PAR *parameters, int max_segs );
int unur_arou_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
int unur_arou_set_center( UNUR_PAR *parameters, double center );
int unur_arou_set_usecenter( UNUR_PAR *parameters, int usecenter );
int unur_arou_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_arou_set_verify( UNUR_PAR *parameters, int verify );
int unur_arou_chg_verify( UNUR_GEN *generator, int verify );
int unur_arou_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* end of `arou.h' */
/*-----*/
/*-----*/
/* <1> `hinv.h' */
UNUR_PAR *unur_hinv_new( const UNUR_DISTR *distribution );
int unur_hinv_set_order( UNUR_PAR *parameters, int order);
int unur_hinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
int unur_hinv_set_cpoints( UNUR_PAR *parameters, const double *stp, int n_stp );
int unur_hinv_set_boundary( UNUR_PAR *parameters, double left, double right );
int unur_hinv_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_hinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_hinv_get_n_intervals( const UNUR_GEN *generator );
double unur_hinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
int unur_hinv_chg_truncated( UNUR_GEN *generator, double left, double right );
int unur_hinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
/* end of `hinv.h' */
/*-----*/
/*-----*/
/* <1> `hrb.h' */
UNUR_PAR *unur_hrb_new( const UNUR_DISTR *distribution );
int unur_hrb_set_upperbound( UNUR_PAR *parameters, double upperbound );
int unur_hrb_set_verify( UNUR_PAR *parameters, int verify );
int unur_hrb_chg_verify( UNUR_GEN *generator, int verify );
/* end of `hrb.h' */
/*-----*/
/*-----*/
/* <1> `hrd.h' */
UNUR_PAR *unur_hrd_new( const UNUR_DISTR *distribution );
int unur_hrd_set_verify( UNUR_PAR *parameters, int verify );
int unur_hrd_chg_verify( UNUR_GEN *generator, int verify );
/* end of `hrd.h' */
/*-----*/
/*-----*/
/* <1> `hri.h' */
UNUR_PAR *unur_hri_new( const UNUR_DISTR *distribution );
int unur_hri_set_p0( UNUR_PAR *parameters, double p0 );
int unur_hri_set_verify( UNUR_PAR *parameters, int verify );
int unur_hri_chg_verify( UNUR_GEN *generator, int verify );
/* end of `hri.h' */
/*-----*/
/*-----*/
/* <1> `ninv.h' */
UNUR_PAR *unur_ninv_new( const UNUR_DISTR *distribution );
int unur_ninv_set_useregula( UNUR_PAR *parameters );
int unur_ninv_set_usenewton( UNUR_PAR *parameters );
int unur_ninv_set_max_iter( UNUR_PAR *parameters, int max_iter );
int unur_ninv_set_x_resolution( UNUR_PAR *parameters, double x_resolution);
int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
int unur_ninv_set_table(UNUR_PAR *parameters, int no_of_points);
int unur_ninv_chg_max_iter(UNUR_GEN *generator, int max_iter);
int unur_ninv_chg_x_resolution(UNUR_GEN *generator, double x_resolution);
int unur_ninv_chg_start(UNUR_GEN *gen, double left, double right);
int unur_ninv_chg_table(UNUR_GEN *gen, int no_of_points);
int unur_ninv_chg_truncated(UNUR_GEN *gen, double left, double right);
int unur_ninv_chg_pdfparams(UNUR_GEN *generator, double *params, int n_params);
/* end of `ninv.h' */
/*-----*/
/*-----*/
/* <1> `srou.h' */
UNUR_PAR *unur_srou_new( const UNUR_DISTR *distribution );
int unur_srou_reinit( UNUR_GEN *generator );
int unur_srou_set_r( UNUR_PAR *parameters, double r );
int unur_srou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
int unur_srou_set_pdfatmode( UNUR_PAR *parameters, double fmode );
int unur_srou_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
int unur_srou_set_usemirror( UNUR_PAR *parameters, int usemirror );
int unur_srou_set_verify( UNUR_PAR *parameters, int verify );
int unur_srou_chg_verify( UNUR_GEN *generator, int verify );
int unur_srou_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
int unur_srou_upd_mode( UNUR_GEN *generator );
int unur_srou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
int unur_srou_chg_pdfatmode( UNUR_GEN *generator, double fmode );
int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
int unur_srou_upd_pdfarea( UNUR_GEN *generator );
/* end of `srou.h' */
/*-----*/
/*-----*/
/* <1> `ssr.h' */
UNUR_PAR *unur_ssr_new( const UNUR_DISTR *distribution );
int unur_ssr_reinit( UNUR_GEN *generator );
int unur_ssr_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
int unur_ssr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
int unur_ssr_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
int unur_ssr_set_verify( UNUR_PAR *parameters, int verify );
int unur_ssr_chg_verify( UNUR_GEN *generator, int verify );
int unur_ssr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_ssr_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_ssr_chg_mode( UNUR_GEN *generator, double mode );
int unur_ssr_upd_mode( UNUR_GEN *generator );
int unur_ssr_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
int unur_ssr_chg_pdfatmode( UNUR_GEN *generator, double fmode );
int unur_ssr_chg_pdfarea( UNUR_GEN *generator, double area );
int unur_ssr_upd_pdfarea( UNUR_GEN *generator );
/* end of `ssr.h' */
/*-----*/
/*-----*/
/* <1> `tabl.h' */
UNUR_PAR *unur_tabl_new( const UNUR_DISTR* distribution );
int unur_tabl_set_usedars( UNUR_PAR *parameters, int usedars );
int unur_tabl_set_variant_splitmode( UNUR_PAR *parameters, unsigned splitmode );
int unur_tabl_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
double unur_tabl_get_sqhratio( const UNUR_GEN *generator );
double unur_tabl_get_hatarea( const UNUR_GEN *generator );
double unur_tabl_get_squeezearea( const UNUR_GEN *generator );
int unur_tabl_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_tabl_set_areafraction( UNUR_PAR *parameters, double fraction );
int unur_tabl_set_nstp( UNUR_PAR *parameters, int n_stp );
int unur_tabl_set_slopes( UNUR_PAR *parameters, const double *slopes, int n_slopes );
int unur_tabl_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_tabl_set_boundary( UNUR_PAR *parameters, double left, double right );
int unur_tabl_set_verify( UNUR_PAR *parameters, int verify );
int unur_tabl_chg_verify( UNUR_GEN *generator, int verify );
/* end of `tabl.h' */
/*-----*/
/*-----*/
/* <1> `tdr.h' */
UNUR_PAR *unur_tdr_new( const UNUR_DISTR* distribution );
int unur_tdr_set_c( UNUR_PAR *parameters, double c );
int unur_tdr_set_variant_gw( UNUR_PAR *parameters );
int unur_tdr_set_variant_ps( UNUR_PAR *parameters );
int unur_tdr_set_variant_ia( UNUR_PAR *parameters );
int unur_tdr_set_usedars( UNUR_PAR *parameters, int usedars );
int unur_tdr_set_darsfactor( UNUR_PAR *parameters, double factor );
int unur_tdr_chg_truncated(UNUR_GEN *gen, double left, double right);
int unur_tdr_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
double unur_tdr_get_sqhratio( const UNUR_GEN *generator );
double unur_tdr_get_hatarea( const UNUR_GEN *generator );
double unur_tdr_get_squeezearea( const UNUR_GEN *generator );
int _unur_tdr_is_ARS_running( const UNUR_GEN *generator );
int unur_tdr_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
int unur_tdr_set_center( UNUR_PAR *parameters, double center );
int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
int unur_tdr_set_verify( UNUR_PAR *parameters, int verify );
int unur_tdr_chg_verify( UNUR_GEN *generator, int verify );
int unur_tdr_set_pedantic( UNUR_PAR *parameters, int pedantic );
double unur_tdr_eval_invcdfhat( const UNUR_GEN *generator, double u, 
double *hx, double *fx, double *sqx );
/* end of `tdr.h' */
/*-----*/
/*-----*/
/* <1> `utdr.h' */
UNUR_PAR *unur_utdr_new( const UNUR_DISTR *distribution );
int unur_utdr_reinit( UNUR_GEN *generator );
int unur_utdr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
int unur_utdr_set_cpfactor( UNUR_PAR *parameters, double cp_factor );
int unur_utdr_set_deltafactor( UNUR_PAR *parameters, double delta );
int unur_utdr_set_verify( UNUR_PAR *parameters, int verify );
int unur_utdr_chg_verify( UNUR_GEN *generator, int verify );
int unur_utdr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_utdr_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_utdr_chg_mode( UNUR_GEN *generator, double mode );
int unur_utdr_upd_mode( UNUR_GEN *generator );
int unur_utdr_chg_pdfatmode( UNUR_GEN *generator, double fmode );
int unur_utdr_chg_pdfarea( UNUR_GEN *generator, double area );
int unur_utdr_upd_pdfarea( UNUR_GEN *generator );
/* end of `utdr.h' */
/*-----*/
/*-----*/
/* <1> `empk.h' */
UNUR_PAR *unur_empk_new( const UNUR_DISTR *distribution );
int unur_empk_set_kernel( UNUR_PAR *parameters, unsigned kernel);
int unur_empk_set_kernelgen( UNUR_PAR *parameters, const UNUR_GEN *kernelgen, double alpha, double kernelvar );
int unur_empk_set_beta( UNUR_PAR *parameters, double beta );
int unur_empk_set_smoothing( UNUR_PAR *parameters, double smoothing );
int unur_empk_chg_smoothing( UNUR_GEN *generator, double smoothing );
int unur_empk_set_varcor( UNUR_PAR *parameters, int varcor );
int unur_empk_chg_varcor( UNUR_GEN *generator, int varcor );
int unur_empk_set_positive( UNUR_PAR *parameters, int positive );
/* end of `empk.h' */
/*-----*/
/*-----*/
/* <1> `empl.h' */
UNUR_PAR *unur_empl_new( const UNUR_DISTR *distribution );
/* end of `empl.h' */
/*-----*/
/*-----*/
/* <1> `vmt.h' */
UNUR_PAR *unur_vmt_new( const UNUR_DISTR *distribution );
int unur_vmt_set_marginalgen( UNUR_PAR *parameters, const UNUR_GEN *uvgen );
/* end of `vmt.h' */
/*-----*/
/*-----*/
/* <1> `vempk.h' */
UNUR_PAR *unur_vempk_new( const UNUR_DISTR *distribution );
int unur_vempk_set_smoothing( UNUR_PAR *parameters, double smoothing );
int unur_vempk_chg_smoothing( UNUR_GEN *generator, double smoothing );
int unur_vempk_set_varcor( UNUR_PAR *parameters, int varcor );
int unur_vempk_chg_varcor( UNUR_GEN *generator, int varcor );
/* end of `vempk.h' */
/*-----*/
/*-----*/
/* <1> `cstd.h' */
#define UNUR_STDGEN_DEFAULT   0        
#define UNUR_STDGEN_INVERSION (~0u)    
#define UNUR_STDGEN_FAST      (0)      
UNUR_PAR *unur_cstd_new( const UNUR_DISTR *distribution );
int unur_cstd_set_variant( UNUR_PAR *parameters, unsigned variant );
int unur_cstd_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_cstd_chg_truncated( UNUR_GEN *generator, double left, double right );
/* end of `cstd.h' */
/*-----*/
/*-----*/
/* <1> `dstd.h' */
UNUR_PAR *unur_dstd_new( const UNUR_DISTR *distribution );
int unur_dstd_set_variant( UNUR_PAR *parameters, unsigned variant );
int unur_dstd_chg_pmfparams( UNUR_GEN *gen, double *params, int n_params );
/* end of `dstd.h' */
/*-----*/
/*-----*/
/* <1> `unif.h' */
UNUR_PAR *unur_unif_new( const UNUR_DISTR *dummy );
/* end of `unif.h' */
/*-----*/
/*-----*/
/* <1> `parser.h' */
UNUR_GEN *unur_str2gen( const char *string );
UNUR_DISTR *unur_str2distr( const char *string );
UNUR_PAR *_unur_str2par( const UNUR_DISTR *distribution, const char *string, struct unur_slist **mlist );
/* end of `parser.h' */
/*-----*/
/*-----*/
/* <1> `x_gen.h' */
UNUR_GEN *unur_init( UNUR_PAR *parameters );
int    unur_sample_discr(UNUR_GEN *generator);
double unur_sample_cont(UNUR_GEN *generator);
void   unur_sample_vec(UNUR_GEN *generator, double *vector);
void  unur_free( UNUR_GEN *generator );
int unur_get_dimension( const UNUR_GEN *generator );
const char *unur_get_genid( const UNUR_GEN *generator );
const UNUR_DISTR *unur_get_distr( const UNUR_GEN *generator );
UNUR_GEN *unur_gen_clone( const UNUR_GEN *gen );
/* end of `x_gen.h' */
/*-----*/
/*-----*/
/* <1> `unur_distributions.h' */
#ifndef UNURAN_DISTRIBUTIONS_H_SEEN
#define UNURAN_DISTRIBUTIONS_H_SEEN
/*-----*/
/* <2> `unur_stddistr.h' */
#ifndef UNUR_STDDISTR_H_SEEN
#define UNUR_STDDISTR_H_SEEN
enum {
  UNUR_DISTR_GENERIC  = 0x0u,
  UNUR_DISTR_CORDER,              
  UNUR_DISTR_BETA,                
  UNUR_DISTR_CAUCHY,              
  UNUR_DISTR_CHI,                 
  UNUR_DISTR_CHISQUARE,           
  UNUR_DISTR_EPANECHNIKOV,        
  UNUR_DISTR_EXPONENTIAL,         
  UNUR_DISTR_EXTREME_I,           
  UNUR_DISTR_EXTREME_II,          
  UNUR_DISTR_GAMMA,               
  UNUR_DISTR_GIG,                 
  UNUR_DISTR_LAPLACE,             
  UNUR_DISTR_LOGISTIC,            
  UNUR_DISTR_LOGNORMAL,           
  UNUR_DISTR_LOMAX,               
  UNUR_DISTR_NORMAL   = 0x0100u,  
   UNUR_DISTR_GAUSSIAN = 0x0100u, 
  UNUR_DISTR_PARETO,              
  UNUR_DISTR_POWEREXPONENTIAL,    
  UNUR_DISTR_RAYLEIGH,            
  UNUR_DISTR_SLASH,               
  UNUR_DISTR_STUDENT,             
  UNUR_DISTR_TRIANGULAR,          
  UNUR_DISTR_UNIFORM = 0x0200u,   
   UNUR_DISTR_BOXCAR = 0x0200u,   
  UNUR_DISTR_WEIBULL,             
  UNUR_DISTR_BURR_I,              
  UNUR_DISTR_BURR_II,             
  UNUR_DISTR_BURR_III,            
  UNUR_DISTR_BURR_IV,             
  UNUR_DISTR_BURR_V,              
  UNUR_DISTR_BURR_VI,             
  UNUR_DISTR_BURR_VII,            
  UNUR_DISTR_BURR_VIII,           
  UNUR_DISTR_BURR_IX,             
  UNUR_DISTR_BURR_X,              
  UNUR_DISTR_BURR_XI,             
  UNUR_DISTR_BURR_XII,            
  UNUR_DISTR_BINOMIAL,            
  UNUR_DISTR_GEOMETRIC,           
  UNUR_DISTR_HYPERGEOMETRIC,      
  UNUR_DISTR_LOGARITHMIC,         
  UNUR_DISTR_NEGATIVEBINOMIAL,    
  UNUR_DISTR_POISSON,             
  UNUR_DISTR_ZIPF,                
  UNUR_DISTR_MNORMAL              
};
#endif  
/* end of `unur_stddistr.h' */
/*-----*/
UNUR_DISTR *unur_distr_beta(const double *params, int n_params);
UNUR_DISTR *unur_distr_burr(const double *params, int n_params);
UNUR_DISTR *unur_distr_cauchy(const double *params, int n_params);
UNUR_DISTR *unur_distr_chi(const double *params, int n_params);
UNUR_DISTR *unur_distr_chisquare(const double *params, int n_params);
UNUR_DISTR *unur_distr_exponential(const double *params, int n_params);
UNUR_DISTR *unur_distr_extremeI(const double *params, int n_params);
UNUR_DISTR *unur_distr_extremeII(const double *params, int n_params);
UNUR_DISTR *unur_distr_gamma(const double *params, int n_params);
UNUR_DISTR *unur_distr_gig(const double *params, int n_params);
UNUR_DISTR *unur_distr_laplace(const double *params, int n_params);
UNUR_DISTR *unur_distr_logistic(const double *params, int n_params);
UNUR_DISTR *unur_distr_lognormal(const double *params, int n_params);
UNUR_DISTR *unur_distr_lomax(const double *params, int n_params);
UNUR_DISTR *unur_distr_normal( const double *params, int n_params );
UNUR_DISTR *unur_distr_pareto( const double *params, int n_params );
UNUR_DISTR *unur_distr_powerexponential(const double *params, int n_params);
UNUR_DISTR *unur_distr_rayleigh(const double *params, int n_params);
UNUR_DISTR *unur_distr_slash(const double *params, int n_params);
UNUR_DISTR *unur_distr_student(const double *params, int n_params);
UNUR_DISTR *unur_distr_triangular(const double *params, int n_params);
UNUR_DISTR *unur_distr_uniform(const double *params, int n_params);
UNUR_DISTR *unur_distr_weibull(const double *params, int n_params);
UNUR_DISTR *unur_distr_multinormal(int dim, const double *mean, const double *covar);
UNUR_DISTR *unur_distr_binomial(const double *params, int n_params);
UNUR_DISTR *unur_distr_geometric(const double *params, int n_params);
UNUR_DISTR *unur_distr_hypergeometric(const double *params, int n_params);
UNUR_DISTR *unur_distr_logarithmic(const double *params, int n_params);
UNUR_DISTR *unur_distr_negativebinomial(const double *params, int n_params);
UNUR_DISTR *unur_distr_poisson(const double *params, int n_params);
UNUR_DISTR *unur_distr_zipf(const double *params, int n_params);
#endif  
/* end of `unur_distributions.h' */
/*-----*/
/*-----*/
/* <1> `unur_errno.h' */
#ifndef UNUR_ERRNO_H_SEEN
#define UNUR_ERRNO_H_SEEN
enum { 
  UNUR_ERR_DISTR_SET      = 0x11u,    
  UNUR_ERR_DISTR_GET      = 0x12u,    
  UNUR_ERR_DISTR_NPARAMS  = 0x13u,    
  UNUR_ERR_DISTR_DOMAIN   = 0x14u,    
  UNUR_ERR_DISTR_GEN      = 0x15u,    
  UNUR_ERR_DISTR_REQUIRED = 0x16u,    
  UNUR_ERR_DISTR_UNKNOWN  = 0x17u,    
  UNUR_ERR_DISTR_INVALID  = 0x18u,    
  UNUR_ERR_DISTR_DATA     = 0x19u,    
  UNUR_ERR_PAR_SET        = 0x21u,    
  UNUR_ERR_PAR_VARIANT    = 0x22u,    
  UNUR_ERR_PAR_INVALID    = 0x23u,    
  UNUR_ERR_GEN            = 0x31u,    
  UNUR_ERR_GEN_DATA       = 0x32u,    
  UNUR_ERR_GEN_CONDITION  = 0x33u,    
  UNUR_ERR_GEN_INVALID    = 0x34u,    
  UNUR_ERR_GEN_SAMPLING   = 0x35u,    
  UNUR_ERR_STR            = 0x41u,    
  UNUR_ERR_STR_UNKNOWN    = 0x42u,    
  UNUR_ERR_STR_SYNTAX     = 0x43u,    
  UNUR_ERR_STR_INVALID    = 0x44u,    
  UNUR_ERR_FSTR_SYNTAX    = 0x45u,    
  UNUR_ERR_FSTR_DERIV     = 0x46u,    
  UNUR_ERR_DOMAIN         = 0x01u,    
  UNUR_ERR_ROUNDOFF       = 0x02u,    
  UNUR_ERR_MALLOC         = 0x03u,    
  UNUR_ERR_NULL           = 0x04u,     
  UNUR_ERR_COOKIE         = 0x05u,    
  UNUR_ERR_GENERIC        = 0x06u,    
  UNUR_ERR_COMPILE        = 0x0eu,    
  UNUR_ERR_SHOULD_NOT_HAPPEN = 0x0fu  
};
extern unsigned unur_errno;
const char *unur_get_strerror ( const int unur_errno );
FILE *unur_set_stream( FILE *new_stream );
FILE *unur_get_stream( void );
#endif  
/* end of `unur_errno.h' */
/*-----*/
/*-----*/
/* <1> `debug.h' */
#define UNUR_DEBUG_OFF     (0u)           
#define UNUR_DEBUG_ALL     (~0u)      
#define UNUR_DEBUG_INIT    0x00000001u    
#define UNUR_DEBUG_SETUP   0x00000fffu    
#define UNUR_DEBUG_ADAPT   0x00fff000u    
#define UNUR_DEBUG_SAMPLE  0xff000000u    
int unur_set_debug( UNUR_PAR *parameters, unsigned debug );
int unur_chg_debug( UNUR_GEN *generator, unsigned debug );
int unur_set_default_debug( unsigned debug );
/* end of `debug.h' */
/*-----*/
/*-----*/
/* <1> `umath.h' */
#ifndef MATH_H_SEEN
#define MATH_H_SEEN

#include <math.h>

#define UNUR_INFINITY      HUGE_VAL     
#ifndef TRUE
#define TRUE   (1)
#endif
#ifndef FALSE
#define FALSE  (0)
#endif
#endif  
/* end of `umath.h' */
/*-----*/
/*-----*/
/* <1> `slist.h' */
#ifndef SLIST_H_SEEN
#define SLIST_H_SEEN
struct unur_slist *_unur_slist_new( void );
void _unur_slist_append( struct unur_slist *slist, void *element );
int _unur_slist_length( const struct unur_slist *slist );
void *_unur_slist_get( const struct unur_slist *slist, int n );
void _unur_slist_free( struct unur_slist *slist );
#endif  
/* end of `slist.h' */
/*-----*/
#ifndef TRUE
#define TRUE   (1)
#endif
#ifndef FALSE
#define FALSE  (0)
#endif
#endif  
__END_DECLS
