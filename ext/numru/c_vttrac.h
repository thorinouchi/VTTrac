typedef struct vttrac_obj {
    // data on which tracking is made
    int nx, ny;   // image size
    int nt;    // time length
    float  **z;
    double *t;
    char *tunit;
    double dtmean;

    // tracking parameters
    int imiss; // = -999        // missing value to be set in int vars
    double rmiss; // = -999.0  // missing value to be set in real vars
    int chk_zmiss; // = 0 // if true (!=0), check missing values in z (image)
    float zmiss;        // missing value in z
    int nsx, nsy; // sub-image size
    double vxhw, vyhw;  // velocities corresponding to ixhw, iyhw throug dtmean
    int ixhw, iyhw;     // max displacement for template matching
    double vxch, vych;
    int itstep, ntrac;

    int subgrid;
    int subgrid_gaus; //1: subgrid peak finding is by gaussian; 0: e-paraboloid
    int score_method;
    double score_th0, score_th1;
    float peak_inside_th; //threshold for the peak-inside screening(unused if<0)
    float min_contrast; // minimum contrast in the template (unused if=<0)
    int use_init_temp; // if true, always use initial template submimages
} vttrac_obj_t;

extern void
c_vttrac_initialize(vttrac_obj_t *o, int nx, int ny, int nt,
                    float  **z,  int chk_zmiss, float zmiss,
                    double *t, char *tunit, int imiss, double rmiss);

extern void
c_vttrac_set_ixyhw_directly(vttrac_obj_t *o, int ixhw, int iyhw);

extern void
c_vttrac_set_ixyhw_from_v(vttrac_obj_t *o, double vxhw, double vyhw);

extern void
c_vttrac_set_basic(vttrac_obj_t *o,
                   int nsx, int nsy,
                   int itstep, int ntrac);

extern void
c_vttrac_set_optional(vttrac_obj_t *o,
                      int subgrid, int subgrid_gaus, int score_method,
                      double score_th0, double score_th1,
                      float peak_inside_th, float min_contrast,
                      double vxch, double vych, int use_init_temp);

extern void
c_vttrac_trac(vttrac_obj_t *o,
              int len, int *tid0, double *x0, double *y0,   // IN
              double *vx0g, double *vy0g,                   // IN
              int *count, int *tid,                         // OUT
              double *x, double *y, double *vx,  double *vy, // OUT
              double *score,   // OUT
              float *zss, double *score_ary  // (Optional) OUT
              );
