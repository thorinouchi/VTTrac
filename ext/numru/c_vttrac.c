/**
 * @file  c_vttrac.c
 * @brief  The core part of VTTrac, Velocimetry by Template Tracking
 * @author Takeshi Horinouchi
 * @date   2021-05-26
 
Copyright 2021 (C) Takeshi Horinouchi. All rights reserved.

LICENSE: see ../../LICENSE.txt

This file provides the core part of the tracking for a time sequence
of monochromatic image-like data. It is meant to fit to wrapping for
object-oriented languages. The entire data and parameters are stored
(or referenced) in the "object" in the type vttrac_obj_t, so multiple
instances are used for tracking simultanously.

The axes of the image arrays z[tid] (tid is time index) are named as x
and y, where the x axis is the fastest varying dimension consecutive
on memory.

The supported image-like data are single-precision floats, as
double-precision data are rarely required in practive. All the other
numeric values, including score-related ones, are double-precisioned
to avoid possible severe run-off.

To use this C library, one needs to set the entire parameters of the
vttrac_obj_t instances properly. The Ruby wrapper in ../../lib/numru
provides a more easy-to-use user interface.

**/

#include<stdio.h>
#include<stdlib.h>
#include<alloca.h>
#include<math.h>
#include "c_vttrac.h"

#define XCOR 1
#define NCOV 2

/**
 * @fn sets data for tracking (you need to set parameters separately)
 *
 * @param [out] o  the object to be initialized.
 * @param [in] nx  size of the first (fastest varying) dimension of the "images"
 * @param [in] ny  size of the second dimension of the "images"
 * @param [in] nt  number of the times
 * @param [in] z   array of image-like data. z[i] contains i-th image data.
 * @param [in] chk_zmiss  if z may have data missing.
 * @param [in] zmiss  missing value used if chk_zmiss.
 * @param [in] t  times at which the images are for (size must be nt).
 * @param [in] tunit  unit of t (such as "s" or ""..  Unused, just kept.)
 * @param [in] imiss  missing value to be set for integers
 * @param [in] rmiss  missing value to be set for doubles
**/
void
c_vttrac_initialize(vttrac_obj_t *o, int nx, int ny, int nt,
                    float  **z,  int chk_zmiss, float zmiss,
                    double *t, char *tunit, int imiss, double rmiss)
{
    o->nx = nx;
    o->ny = ny;
    o->nt = nt;
    o->z = z;
    o->chk_zmiss = chk_zmiss;
    o->zmiss = zmiss;
    o->t = t;
    o->tunit = tunit;
    o->imiss = imiss;
    o->rmiss = rmiss;
    o->dtmean = (t[nt-1] - t[0])/(nt-1.0);
}

/**
 * @fn sets the tracking parameters i[xy]hw
 *
 * @param [in,out] o  The object.
 * @param [in] ixhw  The range over which next x is searched around initial guess
 * @param [in] iyhw  The range over which next y is searched around initial guess
 **/
void
c_vttrac_set_ixyhw_directly(vttrac_obj_t *o, int ixhw, int iyhw)
{
    o->ixhw = ixhw;
    o->iyhw = iyhw;
    o->vxhw = ixhw/o->dtmean - 1; // max displacement
    o->vyhw = iyhw/o->dtmean - 1; // -1 is from margin to find peak
}

/**
 * @fn sets the tracking parameters i[xy]hw from velocities (v[xy]hh)
 *
 * @param [in,out] o  The object.
 * @param [in] vxhw  The range over which vx is searched around initial guess
 * @param [in] vyhw  The range over which vy is searched around initial guess
**/
void
c_vttrac_set_ixyhw_from_v(vttrac_obj_t *o, double vxhw, double vyhw)
{
    o->vxhw = vxhw;
    o->vyhw = vyhw;
    o->ixhw = ceil( fabs(vxhw * o->dtmean) ) + 1; // max displacement
    o->iyhw = ceil( fabs(vyhw * o->dtmean) ) + 1; // +1 is margin to find peak
}

/**
 * @fn sets basic (mandatory) tracking parameters. Also sets default vals for optional params.
 *
 * @param [in,out] o  The object.
 * @param [in] nsx  The template subimage size (x)
 * @param [in] nsy  The template subimage size (y)
 * @param [in] itstep  Index-based time step for tracking (can be negative)
 * @param [in] ntrac  Number of times for each initial template is tracked
**/
void
c_vttrac_set_basic(vttrac_obj_t *o,
                   int nsx, int nsy,
                   int itstep, int ntrac)
{
    o->nsx = nsx;
    o->nsy = nsy;
    o->itstep = itstep;
    o->ntrac = ntrac;

    // optional parameters (default exists)
    o->subgrid = 1;
    o->subgrid_gaus = 0;
    o->score_method = XCOR;
    o->score_th0 = 0.8;
    o->score_th1 = 0.7;
    o->peak_inside_th = 0.03;   // unused if <0
    o->min_contrast = -999.0;   // unused if <0
    o->vxch = -999.0;   // unused if <0
    o->vych = -999.0;   // unused if <0
    o->use_init_temp = 0; // false
}

/**
 * @fn sets optional tracking parameters.
 *
 * @param [in,out] o  The object.
 * @param [in] subgrid  Whether to conduct subgrid tracking
 * @param [in] subgrid_gaus  Whether subgrid peak finding is by gaussian
 * @param [in] score_method  Scoring method (such as XCOR for cross-correlation)
 * @param [in] score_th0  (Result screening parameter) Minimum score required for the first-time tracking
 * @param [in] score_th1  (Result screening parameter) Minimum score required for the subsequent tracking
 * @param [in] peak_inside_th  (Template screening parameter) If positive, an initial template is used only when it is peaked (max or min) inside, exceeding the max or min along the sides by the ratio specified by its value.
 * @param [in] min_contrast  (Template screening parameter) If positive, an initial template is used only when it has a difference in max and min greater than its value.
 * @param [in] vxch  (Result screening parameter) If positive, tracking result is rejected if the vx chnages along trajecty greather than this value (thus used only when ntrac>=2). As a special case, if the result of the second tracking is rejected, the first one is also rejected, since there is no consecutive consistent result in this case.
 * @param [in] vych  (Result screening parameter) As vxch but for the y-component.
 **/
void
c_vttrac_set_optional(vttrac_obj_t *o,
                      int subgrid, int subgrid_gauss, int score_method,
                      double score_th0, double score_th1,
                      float peak_inside_th, float min_contrast,
                      double vxch, double vych, int use_init_temp)
{
    o->subgrid = subgrid;
    o->score_method = score_method;
    o->score_th0 = score_th0;
    o->score_th1 = score_th1;
    o->peak_inside_th = peak_inside_th;
    o->min_contrast = min_contrast;
    o->vxch = vxch;   // unused if <0
    o->vych = vxch;
    o->use_init_temp = use_init_temp;
}

/* (private function) to check whether a time index is valid

 Returns 0 if valid, 1 if not.
 */
static int
inspect_t_index(vttrac_obj_t *o, int tid)
{
    int stat;
    stat = !( tid >= 0 && tid < o->nt);
    return stat;
}


/* (private function) print a 2D double array */
static void
print_dary2d(char *ttl, char *fmt, double *a, int nx, int ny)
{
    int i,j;
    printf("%s\n", ttl);
    for(j=0 ; j<ny ; j++){
        for(i=0 ; i<nx ; i++){
            printf(fmt,a[i+nx*j]);
        }
        printf("\n");
    }
}

/* (private function) multiply a constant to a subimage  */
static void
zsub_fact(vttrac_obj_t *o, float *zs, float fact)
{
    int i;
    for (i = 0; i < o->nsx*o->nsy; i++){
        zs[i] *= fact;
    }
}

/* (private function) add a subimage onto another */
static void
zsub_add(vttrac_obj_t *o, float *zs, float *zs2)
{
    int i;
    for ( i = 0; i < o->nsx*o->nsy; i++){
        zs[i] += zs2[i];
    }
}

/* (private function) copy a subimage */
static void
zsub_copy(vttrac_obj_t *o, float *zsto, float *zsfrom)
{
    int i;
    for ( i = 0; i < o->nsx*o->nsy; i++){
        zsto[i] += zsfrom[i];
    }
}

/* (private function) derive subimage's mean and the deviation from it */
static void
zsub_mean_dev(vttrac_obj_t *o, float *zs, double *mean, float *dev)
{
    int i, nsxy;
    *mean = 0.0;
    nsxy = o->nsx * o->nsy;
    for (i = 0; i < nsxy; i++){
        *mean += zs[i];
    }
    *mean /= nsxy;
    for (i = 0; i < nsxy; i++){
        dev[i] = zs[i] - *mean;
    }
}

/* (private function) 
   compute subimage's sigma (not unbiased) from the deviation from mean
 */
static double
zsub_sig_from_dev(vttrac_obj_t *o, float *zs)
{
    int i, nsxy;
    double var = 0.0;
    nsxy = o->nsx * o->nsy;
    for (i = 0; i < nsxy; i++){
       var += zs[i] * zs[i];
    }
    var /= nsxy;
    return( sqrt(var) );
}

/* (private function) derive max and min in subimage */
static void
zsub_min_max(vttrac_obj_t *o, float *zs, float *min, float *max)
{
    int i, j, nsx;
    nsx = o->nsx;
    *max = *min = zs[0];   // initialization
    for(j = 0; j < o->nsy; j++){
        for(i = 0; i < nsx; i++){
            if (zs[i + nsx*j] > *max) {*max = zs[i + nsx*j];}
            if (zs[i + nsx*j] < *min) {*min = zs[i + nsx*j];}
        }
    }
}

/* (private function)  read out a template subimage from the image at tid.

   The sub-image positons are specified at its center. (If the sub-image
   size is even, with one more pix on the "left" / "bottom", by starting
   from the index xi-nsx/2, yi-nsy/2).

   Returns 0 if successful (specified region is valid and, if chk_zmiss, 
   no data missing), 1 if not.
 */
static int
get_zsub(vttrac_obj_t *o, int tid, int xi, int yi, float *zs)
{
    float **z;
    int i, j, xi0, yi0, nsx2, nsy2;
    int stat=0;
    z = o->z;
    nsx2 = o->nsx/2;
    nsy2 = o->nsy/2;
    xi0 = xi - nsx2;
    yi0 = yi - nsy2;
    if (xi0 < 0 || xi0 + o->nsx-1 >= o->nx ||
        yi0 < 0 || yi0 + o->nsy-1 >= o->ny) {
        stat = 1;    // sub-image is not within the original image
        return(stat);
    }
    for (j = 0; j < o->nsy; j++){
        for (i = 0; i < o->nsx; i++){
            zs[i + o->nsx*j] = z[tid][(xi0+i) + o->nx*(yi0+j)];
        }
    }
    if (o->chk_zmiss){
        for (i = 0; i < o->nsx*o->nsy; i++){
            if (zs[i] == o->zmiss){stat=1;}
        }
        if(stat){return(stat);}
    }
    return(stat);
}

/* (private function)
   read out a template submimage from the image at tid.
   Possibly at subgrid: Linearly interpolated, if
   x or y has deviation from integer (bilinear if x and y).
   Efficient: no unnecessary read-out is made.

   Returns 0 if successful (specified region is valid and, if chk_zmiss, 
   no data missing), 1 if not.
*/
static int
get_zsub_subgrid(vttrac_obj_t *o,
                 int tid, double x, double y,
                 float *zs)
{
    float **z;
    int stat, isx=0, isy=0;
    float *zsw1;
    int xi, yi;
    double dx, dy;

    xi = round(x);
    yi = round(y);
    dx = x - xi;
    dy = y - yi;
    z = o->z;
    zsw1 =  (float *) alloca( sizeof(float)* o->nsx * o->nsy );
    if (dx==0.0 && dy==0.0) {
        stat = get_zsub(o, tid, xi, yi, zs);
    } else {
        if (dx>0) {isx=1;} else if(dx<0) {isx=-1; dx=-dx;}
        if (dy>0) {isy=1;} else if(dy<0) {isy=-1; dy=-dy;}
        stat = get_zsub(o, tid, xi, yi, zs);
        if (stat) {return(stat);}
        zsub_fact(o, zs, (1.0-dx)*(1.0-dy));
        if (isx){
            stat = get_zsub(o, tid, xi+isx, yi, zsw1);
            if (stat) {return(stat);}
            zsub_fact(o, zsw1, dx*(1.0-dy));
            zsub_add(o, zs, zsw1);
            if (isy){
                stat = get_zsub(o, tid, xi+isx, yi+isy, zsw1);
                if (stat) {return(stat);}
                zsub_fact(o, zsw1, dx*dy);
                zsub_add(o, zs, zsw1);
            }
        }
        if (isy){
            stat = get_zsub(o, tid, xi, yi+isy, zsw1);
            if (stat) {return(stat);}
            zsub_fact(o, zsw1, (1.0-dx)*dy);
            zsub_add(o, zs, zsw1);
        }
            
    }
    return(stat);
}

/* (private function)
   check if there is data missing in the specified region at tid.

   Returns 0 if there is no data missing, 1 if not.
 */
static int
chk_zmiss_region(vttrac_obj_t *o, int tid, int k0, int k1, int l0, int l1)
{
    int k,l, stat=0, nx;
    nx = o->nx;
    for ( l=l0; l<=l1; l++ ) {
        for ( k=k0; k<=k1; k++ ) {
            stat = (o->z[tid][k + nx*l] == o->zmiss);
            if (stat) {return(stat);}
        }
    }
    return(stat);
}


/* (private function)
   Check whether the template subimage is peaked (maximized or minimized)
   inside, and the peak is conspicuous enough, having a difference from the 
   max or min on the sides greater than 
   peak_inside_th*(inside_max - inside_min).

  CAUTION:
  If o->peak_inside_th < 0, no checking is conducted.

  RETURN VALUE
  * 0 if passed the check, 1 if not.
 */
static int
chk_zsub_peak_inside(vttrac_obj_t *o, float *zs)
{
    float side_max, side_min, inner_max, inner_min;
    int i, j, nsx, nsy;
    int stat=1;

    if ( o->peak_inside_th < 0){
        stat = 0;   // do not check --> no problem
        return(stat);
    }
    nsx = o->nsx;
    nsy = o->nsy;
        
    // find max and min along sides
    side_max = side_min = zs[0];   // initialization
    for(j = 0; j < nsy; j += nsy-1){
        for( i = 0; i < nsx; i++ ){
            if (zs[i + nsx*j] > side_max) {side_max = zs[i + nsx*j];}
            if (zs[i + nsx*j] < side_min) {side_min = zs[i + nsx*j];}
        }
    }
    for(j = 1; j < nsy-1; j++){
        for(i = 0; i < nsx; i += nsx-1){
            if (zs[i + nsx*j] > side_max) {side_max = zs[i + nsx*j];}
            if (zs[i + nsx*j] < side_min) {side_min = zs[i + nsx*j];}
        }
    }

    // find max and min inside the sides
    inner_max = inner_min = zs[1 + nsx];   // initialization
    for(j = 1; j < nsy-1; j++){
        for(i = 1; i < nsx-1; i++){
            if (zs[i + nsx*j] > inner_max) {inner_max = zs[i + nsx*j];}
            if (zs[i + nsx*j] < inner_min) {inner_min = zs[i + nsx*j];}
        }
    }
    if ( inner_max > side_max + o->peak_inside_th*(inner_max-inner_min) ||
         inner_min < side_min - o->peak_inside_th*(inner_max-inner_min) ) {
        stat = 0; // OK, because the max or min is inside and the difference
                  // from the max or min on the sides is not too tiny
    }
    return(stat);
}


/* (private function)
   sliding cross-correlation between the sugimage and image at tid.

   RETURN VALUE
   * 0 if all the relevant data and regions are valid, so all the scores
     (x-cor) are defined at all tested center locations; 1 if not.
 */
static int
sliding_xcor(vttrac_obj_t *o, double sigx, float *xd,
             int tid, int k0, int k1, int l0, int l1,
             double *xcor)
{
    int i, j, k, l, nx, ny, nk, nl, stat;
    int nsx, nsy, nsxy, nsx2, nsy2;
    double xysum, vxy, vyy, ysum, yysum, ymean, yd;
    float y;

    nx = o->nx;
    ny = o->ny;
    nsx = o->nsx;
    nsy = o->nsy;
    nsx2 = nsx/2;
    nsy2 = nsy/2;
    nk = k1-k0+1;
    nl = l1-l0+1;
    nsxy = nsx * nsy;
    k0 = k0 - nsx2;
    l0 = l0 - nsy2;
    stat = ( k0 < 0 || k1+nsx-1 > nx-1 || l0 < 0 || l1+nsy-1 > ny-1 );
    if (stat) {return(stat);}
    if (o->chk_zmiss){
        stat = chk_zmiss_region(o, tid, k0, k1+nsx-1, l0, l1+nsy-1);
        if (stat) {return(stat);}
    }
    for ( l=0; l<nl; l++ ) {
        k=0;
        ysum = 0.0;
        for ( j=0; j<nsy; j++ ) {
            for ( i=0; i<nsx; i++ ) {
                ysum += o->z[tid][(k0+k+i) + nx*(l0+l+j)];
            }
        }
        ymean = ysum/nsxy;
        yysum = 0.0;
        xysum = 0.0;
        for ( j=0; j<nsy; j++ ) {
            for ( i=0; i<nsx; i++ ) {
                yd = o->z[tid][(k0+k+i) + nx*(l0+l+j)] - ymean;
                yysum += yd*yd;
                xysum += xd[i + nsx*j]*yd;
            }
        }
        vyy = yysum/nsxy;
        vxy = xysum/nsxy;
        xcor[k + nk*l] = vxy/sqrt(vyy)/sigx;  // cross-correlataion coef
        for ( k=1; k<nk; k++ ) {
            vyy = vyy + ymean*ymean;   // mean(y^2) for previous k
            yysum = ysum = 0.0;
            for ( j=0; j<nsy; j++ ) {
                i=-1;
                y = o->z[tid][(k0+k+i) + nx*(l0+l+j)];
                ysum -= y;   // to subtract the first in the previous summation
                yysum -= y*y;
                i=nsx-1;
                y = o->z[tid][(k0+k+i) + nx*(l0+l+j)];
                ysum += y;
                yysum += y*y;
            }
            ymean = ymean + ysum/nsxy;     // ymean is renewed.
            vyy = vyy + yysum/nsxy - ymean*ymean; //new mean(y^2) - new ymean^2
            xysum = 0.0;
            for ( j=0; j<nsy; j++ ) {
                for ( i=0; i<nsx; i++ ) {
                    yd = o->z[tid][(k0+k+i) + nx*(l0+l+j)] - ymean;
                    xysum += xd[i + nsx*j]*yd;
                }
            }
            vxy = xysum/nsxy;
            xcor[k + nk*l] = vxy/sqrt(vyy)/sigx;  // cross-correlataion coef
        }
    }
    return(stat);
}

/* (private function)
   sliding normalized covariance between the sugimage and image at tid.

   normalization is done by the sigma of the fist image : cov(x',y')/sigx^2
   (in contrast to cov(x',y')/sigx/sigy in the correlation coefficient).

   RETURN VALUE
   * 0 if all the relevant data and regions are valid, so all the scores
     (x-cor) are defined at all tested center locations; 1 if not.
 */
static int
sliding_ncov(vttrac_obj_t *o, double sigx, float *xd,
             int tid, int k0, int k1, int l0, int l1,
             double *ncov)
{
    int i, j, k, l, nx, ny, nk, nl, stat;
    int nsx, nsy, nsxy, nsx2, nsy2;
    double xysum, vxy, ysum, ymean, yd;
    float y;

    nx = o->nx;
    ny = o->ny;
    nsx = o->nsx;
    nsy = o->nsy;
    nsx2 = nsx/2;
    nsy2 = nsy/2;
    nk = k1-k0+1;
    nl = l1-l0+1;
    nsxy = nsx * nsy;
    k0 = k0 - nsx2;
    l0 = l0 - nsy2;
    stat = ( k0 < 0 || k1+nsx-1 > nx-1 || l0 < 0 || l1+nsy-1 > ny-1 );
    if (stat) {return(stat);}
    if (o->chk_zmiss){
        stat = chk_zmiss_region(o, tid, k0, k1+nsx-1, l0, l1+nsy-1);
        if (stat) {return(stat);}
    }
    for ( l=0; l<nl; l++ ) {
        k=0;
        ysum = 0.0;
        for ( j=0; j<nsy; j++ ) {
            for ( i=0; i<nsx; i++ ) {
                ysum += o->z[tid][(k0+k+i) + nx*(l0+l+j)];
            }
        }
        ymean = ysum/nsxy;
        xysum = 0.0;
        for ( j=0; j<nsy; j++ ) {
            for ( i=0; i<nsx; i++ ) {
                yd = o->z[tid][(k0+k+i) + nx*(l0+l+j)] - ymean;
                xysum += xd[i + nsx*j]*yd;
            }
        }
        vxy = xysum/nsxy;
        ncov[k + nk*l] = vxy/(sigx*sigx);
        for ( k=1; k<nk; k++ ) {
            ysum = 0.0;
            for ( j=0; j<nsy; j++ ) {
                i=-1;
                y = o->z[tid][(k0+k+i) + nx*(l0+l+j)];
                ysum -= y;   // to subtract the first in the previous summation
                i=nsx-1;
                y = o->z[tid][(k0+k+i) + nx*(l0+l+j)];
                ysum += y;
            }
            ymean = ymean + ysum/nsxy;     // ymean is renewed.
            xysum = 0.0;
            for ( j=0; j<nsy; j++ ) {
                for ( i=0; i<nsx; i++ ) {
                    yd = o->z[tid][(k0+k+i) + nx*(l0+l+j)] - ymean;
                    xysum += xd[i + nsx*j]*yd;
                }
            }
            vxy = xysum/nsxy;
            ncov[k + nk*l] = vxy/(sigx*sigx);
        }
    }
    return(stat);
}

/* (private function)
   conduct template matching, scoring by cross-correlation

   RETURN VALUE
   * 0 if passed the check, 1 if not.
 */
static int
get_score_xcor(vttrac_obj_t *o, float *x, int tid,
          int k0, int k1, int l0, int l1,
          double *scr)
{
    int stat = 0;
    int nsxy;
    double xmean, sigx, sigy, xy;
    float *xd;
    nsxy = o->nsx * o->nsy;
    xd =  (float *) alloca( sizeof(float)* nsxy );
    zsub_mean_dev(o, x, &xmean, xd);
    sigx = zsub_sig_from_dev(o, xd);
    stat = sliding_xcor(o, sigx, xd, tid, k0, k1, l0, l1, scr);
    
    return(stat);
}

/* (private function)
   conduct template matching, scoring by cross-correlation

   RETURN VALUE
   * 0 if passed the check, 1 if not.
 */
static int
get_score_ncov(vttrac_obj_t *o, float *x, int tid,
          int k0, int k1, int l0, int l1,
          double *scr)
{
    int stat = 0;
    int nsxy;
    double xmean, sigx, sigy, xy;
    float *xd;
    nsxy = o->nsx * o->nsy;
    xd =  (float *) alloca( sizeof(float)* nsxy );
    zsub_mean_dev(o, x, &xmean, xd);
    sigx = zsub_sig_from_dev(o, xd);
    stat = sliding_ncov(o, sigx, xd, tid, k0, k1, l0, l1, scr);
    
    return(stat);
}

/* (private function)
 conduct template matching driver
 */
static int
get_score(vttrac_obj_t *o, float *zs0, int tid,
          int k0, int k1, int l0, int l1,
          double *scr)
{
    int stat;
    switch(o->score_method){
      case XCOR:
          stat = get_score_xcor(o, zs0, tid,
                                k0, k1, l0, l1,
                                scr);
          break;
      case NCOV:
          stat = get_score_ncov(o, zs0, tid,
                                k0, k1, l0, l1,
                                scr);
          break;
      default:
          stat = 1;
    }
    return(stat);
}

/* (private function)
 Find subgrid peak from 5 points with elliptic paraboloid
 
 input c(enter), l(eft), r(ight), b(ottom), t(op) at
 (0,0), (-1,0), (0,1), (0,-1), (0,1), respectively.
 c must be greater than any of l,r,b,t.

   equation: z = -p(x-x0)^2 + -q(y-y0)^2 + r
               = -p x^2 + 2p x0 x - q y^2 + 2q y0 y + c

 where u = p x0, v = q y0, c = -p x0^2 - q y0^2 + r
 */
static int
find_subgrid_peak_5pt_epara(double c, double l, double r, double b, double t,
                            double *x0, double *y0, double *max)
{
    double p, q, stat;
    l = l-c;
    r = r-c;
    b = b-c;
    t = t-c;
    stat = !( l<=0.0 && r<=0.0 && b<=0.0 && t<=0.0 &&
              ( l<0.0 || r<0.0 ) && ( b<0.0 || t<=0.0) );
    if (stat) return(stat);
    p = -(l+r)/2.0;
    q = -(b+t)/2.0;
    *x0 = (r-l)/(4.0*p);     // --> |x0| < 0.5, if c >= [l,r,b,t] > 0
    *y0 = (t-b)/(4.0*q);     // --> |y0| < 0.5, if c >= [l,r,b,t] > 0
    *max = c +  p * *x0 * *x0 + q * *y0 * *y0; // r= c + p x0^2 + q y0^2
    return(stat);
}

/* (private function)
 Find subgrid peak from 5 points by interpolating with a 2D gaussian (for positive scores).

 It is simply a log-version of the elliptic-paraboloid method
 (find_subgrid_peak_5pt_epara). It appears that this method is
 preferred in many PIVs over the elliptic-paraboloid method
 (find_subgrid_peak_5pt_epara). In a test, there were not much
 difference between them, though.

 input c(enter), l(eft), r(ight), b(ottom), t(op) at
 (0,0), (-1,0), (0,1), (0,-1), (0,1), respectively.
 all of them must be positive
 c must be greater than any of l,r,b,t.
*/
static int
find_subgrid_peak_5pt_gaus(double c, double l, double r, double b, double t,
                           double *x0, double *y0, double *max)
{
    double p, q, stat;
    stat = !( c>0.0 && l>0.0 && r>0.0 && b>0.0 && t>0.0 ); // all must be >0
    if (stat) return(stat);
    c = log(c);
    l = log(l) - c;
    r = log(r) - c;
    b = log(b) - c;
    t = log(t) - c;
    stat = !( l<=0.0 && r<=0.0 && b<=0.0 && t<=0.0 &&
              ( l<0.0 || r<0.0 ) && ( b<0.0 || t<=0.0) );
    if (stat) return(stat);
    p = -(l+r)/2.0;
    q = -(b+t)/2.0;
    *x0 = (r-l)/(4.0*p);     // --> |x0| < 0.5, if c >= [l,r,b,t] > 0
    *y0 = (t-b)/(4.0*q);     // --> |y0| < 0.5, if c >= [l,r,b,t] > 0
    *max = exp(c +  p * *x0 * *x0 + q * *y0 * *y0);
    return(stat);
}

/* (private function)
 Find the score peak and its location

 RETURN VALUE
 * stat : 0 if the peak is inside; 1 if not
 */
static int
find_score_peak(vttrac_obj_t *o, double *scr, int kw, int lw,
                double *kp, double *lp, double *scrp)
{
    int k, l, stat;
    int kpi, lpi;
    double sp;

    // find the max and its index
    kpi = lpi = 0;   // initialization
    sp = scr[0];     // initialization
    for (l=0; l<lw; l++){
        for (k=0; k<kw; k++){
            if (scr[k + kw*l] > sp){
                kpi = k;
                lpi = l;
                sp = scr[k + kw*l];
            }
        }
    }
    *kp = (double) kpi;
    *lp = (double) lpi;
    *scrp = sp;

    // whether on the sides or not
    stat = ( kpi==0 || kpi==kw-1 || lpi==0 || lpi == lw -1);
    if (stat) return(stat);
    
    // subgrid determination
    if (o->subgrid) {
        if (o->subgrid_gaus) {
            stat = find_subgrid_peak_5pt_gaus( scr[kpi+kw*lpi],
                                      scr[kpi-1+kw*lpi], scr[kpi+1+kw*lpi],
                                      scr[kpi+kw*(lpi-1)], scr[kpi+kw*(lpi+1)],
                                      kp, lp, scrp); //[kl]p: relative to [kl]pi
        } else {
            stat = find_subgrid_peak_5pt_epara( scr[kpi+kw*lpi],
                                      scr[kpi-1+kw*lpi], scr[kpi+1+kw*lpi],
                                      scr[kpi+kw*(lpi-1)], scr[kpi+kw*(lpi+1)],
                                      kp, lp, scrp); //[kl]p: relative to [kl]pi
        }
        if (stat) return(stat);
        *kp += kpi;
        *lp += lpi;
    }

    return(stat);
}
    


/**
 * @fn  Conduct tracking
 * 
 * @param [in] o  The traking oject
 * @param [in] len  # of initial time/loc.. (lengths of input arrays from tid to vxg)
 * @param [in] tid0  Tracking initial time indices
 * @param [in] x0  Tracking initial template-center x location (index-based; non-integer for subgrid)
 * @param [in] y0  Tracking initial template-center y location (index-based; non-integer for subgrid)
 * @param [in] vx0g  First guess of vx (to search around it). Can be 0.
 * @param [in] vy0g  First guess of vy (to search around it). Can be 0.
 * @param [out] count (len: len)  The number of successful tracking for each initial template.
 * @param [out] tid  (len: (ntrac+1)*len)  time index of the trajectories (tid0 and subsequent ones)
 * @param [out] x  (len: (ntrac+1)*len)  x locations of the trajectories (x0 and derived ones)
 * @param [out] y  (len: (ntrac+1)*len)  y locations of trajectories (x0 and derived ones)
 * @param [out] vx  (len: ntrac*len)  Derived x-velocity
 * @param [out] vy  (len: ntrac*len)  Derived y-velocity
 * @param [out] score  (len: ntrac*len)  Scores along the trajectory (max values, possibly at subgrid)
 * @param [out]  zss (optional, if non-NULL)  (Diagnosis output if wanted) The subimages along the track (1D pointer for 4D array; nsx * nsy * (ntrac+1) * len
 * @param [out]  score_arry (optional, if non-NULL) (Diagnosis output if wanted) the entire scores (1D pointer for 4D array; (x-sliding size) * (y-sliding size) * (ntrac+1) * len
 **/
void
c_vttrac_trac(vttrac_obj_t *o,
              int len, int *tid0, double *x0, double *y0,   // IN
              double *vx0g, double *vy0g,                   // IN
              int *count, int *tid,                         // OUT
              double *x, double *y, double *vx,  double *vy, // OUT
              double *score,   // OUT
              float *zss, double *score_ary)  // (Optional) OUT
{
    int j, m, stat, out_zss, out_score_ary, ntr0, ntr1, itstep;
    int *alive, nsxy, tidf, tidl, ixhw, iyhw, kc, lc;
    int kw, lw, chk_vchange;
    float *zs0, zmin, zmax, *zsw;
    double dt, *t, *scr, xcur, ycur, xp, yp, sp, vxg, vyg;
    double xw, yw, vxw, vyw;  // for temporary storage for x,y,vx,vy
    float **z;

    out_zss = (zss != NULL);
    out_score_ary = (score_ary != NULL);
    ntr0 = o->ntrac;
    ntr1 = o->ntrac+1;
    nsxy = o->nsx*o->nsy;
    ixhw = o->ixhw;
    iyhw = o->iyhw;
    kw = 2*ixhw + 1;
    lw = 2*iyhw + 1;
    t = o->t;
    z = o->z;
    itstep = o->itstep;
    chk_vchange = ( o->vxch > 0.0 && o->vych > 0.0);

    zs0 =  (float *) alloca( sizeof(float)* o->nsx * o->nsy );
    if (o->use_init_temp)
        zsw =  (float *) alloca( sizeof(float)* o->nsx * o->nsy );
    alive =  (int *) alloca( sizeof(int)* len );
    if (!out_score_ary)
        scr =  (double *) alloca( sizeof(double)* kw*lw );

    j=0;
    for (m=0; m<len; m++){

        // record initial data
        count[m] = 0;
        alive[m] = 1;
        if (o->subgrid) { // initial position
            x[j + m*ntr1] = x0[m];
            y[j + m*ntr1] = y0[m];
        } else {
            x[j + m*ntr1] = round(x0[m]);
            y[j + m*ntr1] = round(y0[m]);
        }
        tid[j + m*ntr1] = tid0[m];
    }

    for (m=0; m<len; m++){
        for (j=0; j<ntr0; j++){
            if (alive[m]) {
                xcur = x[j + m*ntr1];
                ycur = y[j + m*ntr1];

                // record initial data
                tidf = tid0[m] + j*itstep;  // index of the tracking start time
                tidl = tidf + itstep;      // index of the tracking end time
                stat = inspect_t_index(o, tidf);
                if (stat) {alive[m]=0;  continue;}
                if (j==0 || !o->use_init_temp) {
                    stat = get_zsub_subgrid(o, tidf, xcur, ycur, zs0);
                }
                if (stat) {alive[m]=0;  continue;}
                if (o->min_contrast > 0.0) {
                    zsub_min_max(o, zs0, &zmin, &zmax);
                    if ( zmax-zmin < o->min_contrast) {alive[m]=0;  continue;}
                }
                if (o->peak_inside_th > 0.0) {
                    stat = chk_zsub_peak_inside(o, zs0);
                    if (stat) {alive[m]=0;  continue;}
                }
                if (out_zss) {
                    if (j==0 || !o->use_init_temp) {
                        zsub_copy(o, zss + (j + m*ntr1)*nsxy, zs0);
                    } else {
                        stat = get_zsub_subgrid(o, tidf, xcur, ycur, zsw);
                        zsub_copy(o, zss + (j + m*ntr1)*nsxy, zsw);
                    }
                }

                // inspect the tracking end time
                stat = inspect_t_index(o, tidl);
                if (stat) {alive[m]=0;  continue;}
                dt = t[tidl] - t[tidf];    // time diff. can be negative

                if (j==0) {
                    vxg = vx0g[m];
                    vyg = vy0g[m];
                } else {
                    vxg = vx[j-1 + m*ntr0];   // previous result
                    vyg = vy[j-1 + m*ntr0];   // previous result
                }
                kc = round(xcur + vxg*dt);
                lc = round(ycur + vyg*dt);
                if (out_score_ary) {
                    scr = score_ary + (j + m*ntr0)*kw*lw;
                }
                stat = get_score(o, zs0, tidl,
                                 kc-ixhw, kc+ixhw, lc-iyhw, lc+iyhw,
                                 scr);
                alive[m] = !stat;
                //print_dary2d("**score", "%7.3f", scr, kw, lw );
                stat = find_score_peak(o, scr, kw, lw, &xp, &yp, &sp);
                if (stat) {alive[m]=0;  continue;}
                if ( (j==0 && sp < o->score_th0) ||
                     (j!=0 && sp < o->score_th1) ){
                    alive[m]=0;  continue;
                }
                xw = xp+kc-ixhw;    // next position (x)
                yw = yp+lc-iyhw;    // next position (y)
                score[j + m*ntr0] = sp;
                vxw = (xw - x[j + m*ntr1])/dt;    // velocity (x)
                vyw = (yw - y[j + m*ntr1])/dt;    // velocity (y)
                if (chk_vchange && j>0){
                    stat = ( fabs( vxw - vx[j-1 + m*ntr0] ) > o->vxch ||
                             fabs( vyw - vy[j-1 + m*ntr0] ) > o->vych);
                    if (stat) {
                        if (j==1) {   // invalidate the j==0 results too
                            count[m] = 0;
                            x[1 + m*ntr1] = y[1 + m*ntr1] = o->rmiss;
                            vx[m*ntr0] = vy[m*ntr0] = o->rmiss;
                        }
                        alive[m]=0;
                        continue;
                    }
                }
                count[m] = j+1;   // # of valid tracking
                tid[j+1 + m*ntr1] = tidl;
                x[j+1 + m*ntr1] = xw;
                y[j+1 + m*ntr1] = yw;
                vx[j + m*ntr0] = vxw;
                vy[j + m*ntr0] = vyw;
                if (out_zss && j==ntr0-1) { // last sub image
                    xcur = x[j+1 + m*ntr1];
                    ycur = y[j+1 + m*ntr1];
                    stat = get_zsub_subgrid(o, tidf, xcur, ycur, zs0);
                    zsub_copy(o, zss + (j+1 + m*ntr1)*nsxy, zs0);
                }
            }
        }
    }
}
