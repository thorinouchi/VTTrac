/**
 * @file  rb_vttrac_core.c
 * @brief  Ruby interface of VTTrac, definifng a pure C-based class, VTTracCore
 * @author Takeshi Horinouchi
 * @date   2021-05-26

Copyright 2021 (C) Takeshi Horinouchi. All rights reserved.

LICENSE: see ../../LICENSE.txt

A more user-friendly Ruby class, VTTrac, which inherits VTTracCore, is
defined in ../../lib/numru/vttrac.rb

**/

#include<stdio.h>
#include "ruby.h"
#include "narray.h"
#include "c_vttrac.h"

/* for compatibility for NArray and NArray with big memory patch */
#ifndef NARRAY_BIGMEM
typedef int    na_shape_t;
#endif

/* default missing values */
#define IMISS -999
#define RMISS -999.0

struct vttrac_cdata {
    vttrac_obj_t o;
    int nzna;    // number of z-representing NArrays to be marked
    VALUE *zna;  // z-representing NArrays to prevent unexpected GC
    VALUE tna;   // time-representing NArray
};

static void
vttrac_cdata_mark(void *ptr)   // for GC
{
    int i;
    struct vttrac_cdata *op = ptr;
    for (i=0; i<op->nzna; i++) {
        rb_gc_mark(op->zna[i]);
    }
    rb_gc_mark(op->tna);
}

static void
vttrac_cdata_free(void *ptr)
{
    struct vttrac_cdata *op = ptr;

    xfree(op->o.z);
    xfree(op->zna);
    xfree(op);
}

static size_t
vttrac_cdata_memsize(const void *ptr)
{
    return sizeof(struct vttrac_cdata);
}

static const rb_data_type_t vttrac_cdata_type = {
    "vttrac_cdata",
    {vttrac_cdata_mark, vttrac_cdata_free, vttrac_cdata_memsize,},
    0, 0, RUBY_TYPED_FREE_IMMEDIATELY
};

static VALUE
rb_vttrac_core_alloc(VALUE klass)
{
    struct vttrac_cdata *op;
    VALUE obj = TypedData_Make_Struct(klass, struct vttrac_cdata,
                                      &vttrac_cdata_type, op);
    op->o.z = NULL;
    op->o.t = NULL;
    op->o.tunit = NULL;
    return obj;
}

static VALUE
rb_vttrac_core_initialize(VALUE obj,
                          VALUE rb_z, VALUE rb_chk_zmiss, VALUE rb_zmiss,
                          VALUE rb_t,  VALUE rb_tunit
                         )
{
    struct vttrac_cdata *op;
    double *t;
    float  **z;
    int nt, it, nx, ny;   // won't be very long -> int is enough
    int chk_zmiss;
    VALUE zs;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    if ( !rb_obj_is_kind_of(rb_t, cNArray) ) {
	rb_raise(rb_eArgError, "t must be an NArray");
    }
    rb_t = na_cast_object(rb_t, NA_DFLOAT);
    t = NA_PTR_TYPE(rb_t, double *);
    nt = NA_TOTAL(rb_t);
    op->tna = rb_t;
    
    z = (float **) xmalloc( sizeof(float *)*nt );
    if (rb_obj_is_kind_of(rb_z, rb_cArray)) {
        if (RARRAY_LEN(rb_z) != nt) {
            rb_raise(rb_eArgError, "len of z (Array) != t.length");
        }
        op->nzna = nt;
        op->zna = (VALUE *) xmalloc( sizeof(VALUE)*nt );
        for (it=0; it<nt; it++){
            zs = na_cast_object(RARRAY_AREF(rb_z, it), NA_SFLOAT);
            nx = NA_SHAPE0(zs);
            ny = NA_SHAPE1(zs);
            z[it] = NA_PTR_TYPE(zs, float *);  // do not copy, but share ptr
            op->zna[it] = zs;                   // --> need to pretect from GC
        }
    } else if (rb_obj_is_kind_of(rb_z, cNArray)) {
        zs = na_cast_object(rb_z, NA_SFLOAT);
        op->nzna = 1;
        op->zna = (VALUE *) xmalloc( sizeof(VALUE) );
        op->zna[0] = zs;                   // --> need to pretect from GC
        nx = NA_SHAPE0(zs);
        ny = NA_SHAPE1(zs);
        for (it=0; it<nt; it++){
            z[it] = NA_PTR_TYPE(zs, float *) + (size_t) it*nx*ny;
        }
    } else {
        rb_raise(rb_eArgError,"z must be an Array or NArray");
    }

    chk_zmiss = ( rb_chk_zmiss == Qtrue );
    c_vttrac_initialize(&op->o, nx, ny, nt,
                        z,  chk_zmiss, (float) NUM2DBL(rb_zmiss),
                        t, StringValuePtr(rb_tunit), IMISS, RMISS );
    return obj;
}

static VALUE
rb_vttrac_core_set_ixyhw_from_v(VALUE obj,
                                VALUE rb_vxhw, VALUE rb_vyhw)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    c_vttrac_set_ixyhw_from_v(&op->o, NUM2DBL(rb_vxhw), NUM2DBL(rb_vyhw));
}

static VALUE
rb_vttrac_core_set_ixyhw_directly(VALUE obj,
                                  VALUE rb_vxhw, VALUE rb_vyhw)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    c_vttrac_set_ixyhw_directly(&op->o, NUM2INT(rb_vxhw), NUM2INT(rb_vyhw));
}

static VALUE
rb_vttrac_core_set_basic(VALUE obj,
                         VALUE rb_nsx, VALUE rb_nsy,
                         VALUE rb_itstep, VALUE rb_ntrac)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    c_vttrac_set_basic(&op->o,
                       NUM2INT(rb_nsx), NUM2INT(rb_nsy),
                       NUM2INT(rb_itstep), NUM2INT(rb_ntrac)
                       );

    return(Qnil);
}

static VALUE
rb_vttrac_core_set_optional(VALUE obj,
                            VALUE rb_subgrid, VALUE rb_subgrid_gaus,
                            VALUE rb_score_method,
                            VALUE rb_score_th0, VALUE rb_score_th1,
                            VALUE rb_peak_inside_th, VALUE rb_min_contrast,
                            VALUE rb_vxch, VALUE rb_vych,VALUE rb_use_init_temp)
{
    struct vttrac_cdata *op;
    int subgrid, subgrid_gaus, use_init_temp;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    subgrid = ( rb_subgrid == Qtrue );
    subgrid_gaus = ( rb_subgrid_gaus == Qtrue );
    use_init_temp = ( rb_use_init_temp == Qtrue );
    c_vttrac_set_optional(&op->o,
                          subgrid, subgrid_gaus, NUM2INT(rb_score_method),
                          NUM2DBL(rb_score_th0), NUM2DBL(rb_score_th1),
                          NUM2DBL(rb_peak_inside_th), NUM2DBL(rb_min_contrast),
                          NUM2DBL(rb_vxch), NUM2DBL(rb_vych), use_init_temp);

    return(Qnil);
}

static VALUE
rb_vttrac_core_t(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(op->tna);
}

static VALUE
rb_vttrac_core_z(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    if (op->nzna == 1){
        return(op->zna[0]);
    } else {
        return( rb_ary_new_from_values( op->nzna, op->zna ) );
    }
}

static VALUE
rb_vttrac_core_trac(VALUE obj,
                    VALUE rb_tid0, VALUE rb_x0, VALUE rb_y0,
                    VALUE rb_vx0g, VALUE rb_vy0g,
                    VALUE rb_out_subimage,  VALUE rb_out_score_ary)
{
    VALUE chk, zs, rb_count, rb_tid, rb_x, rb_y, rb_vx, rb_vy, rb_score, na;
    VALUE rb_rmiss, rb_zss, rb_score_ary;
    double *t;
    float  **z;
    int nt, it, len, nx, ny, i;   // won't be very long -> int is enough
    int *tid0;         // won't be very long -> int is enough
    double *x0, *y0, *vx0g, *vy0g;
    na_shape_t shape0[2], shape1[2], shape2[4];
    int *count, *tid;
    double *x, *y, *vx, *vy, *score, *score_ary;
    float *zss;

    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);

    if (!rb_obj_is_kind_of(rb_tid0, cNArray) ||
        na_get_typecode(rb_tid0) != NA_LINT ) {
	rb_raise(rb_eArgError, "tid0 must be an Integer NArray");
    }
    tid0 = NA_PTR_TYPE(rb_tid0, int *);
    len = NA_TOTAL(rb_tid0);

    if (!rb_obj_is_kind_of(rb_x0, cNArray) )
        rb_raise(rb_eArgError, "x0 != NArray");
    if (NA_TOTAL(rb_x0) != len) rb_raise(rb_eArgError, "invalid x0 length");
    rb_x0 = na_cast_object(rb_x0, NA_DFLOAT);
    x0 = NA_PTR_TYPE(rb_x0, double *);

    if (!rb_obj_is_kind_of(rb_y0, cNArray) )
	rb_raise(rb_eArgError, "y0 != NArray");
    if (NA_TOTAL(rb_y0) != len) rb_raise(rb_eArgError, "invalid y0 length");
    rb_y0 = na_cast_object(rb_y0, NA_DFLOAT);
    y0 = NA_PTR_TYPE(rb_y0, double *);

    if (!rb_obj_is_kind_of(rb_vx0g, cNArray) )
	rb_raise(rb_eArgError, "vx0g != NArray");
    if (NA_TOTAL(rb_vx0g) != len) rb_raise(rb_eArgError,"invalid vx0g length");
    rb_vx0g = na_cast_object(rb_vx0g, NA_DFLOAT);
    vx0g = NA_PTR_TYPE(rb_vx0g, double *);

    if (!rb_obj_is_kind_of(rb_vy0g, cNArray) ) {
	rb_raise(rb_eArgError, "vy0g != NArray");
    }
    if (NA_TOTAL(rb_vy0g) != len) rb_raise(rb_eArgError,"invalid vy0g length");
    rb_vy0g = na_cast_object(rb_vy0g, NA_DFLOAT);
    vy0g = NA_PTR_TYPE(rb_vy0g, double *);

    shape0[0] = op->o.ntrac;
    shape0[1] = len;
    shape1[0] = op->o.ntrac+1;
    shape1[1] = len;

    rb_rmiss = DBL2NUM(op->o.rmiss);

    rb_count = na_make_object(NA_LINT, 1, &len, cNArray);
    count = NA_PTR_TYPE(rb_count, int *);

    rb_tid = na_make_object(NA_LINT, 2, shape1, cNArray);
    na_fill(rb_tid, rb_rmiss);
    tid = NA_PTR_TYPE(rb_tid, int *);
    rb_x = na_make_object(NA_DFLOAT, 2, shape1, cNArray);
    na_fill(rb_x, rb_rmiss);
    x = NA_PTR_TYPE(rb_x, double *);
    rb_y = na_make_object(NA_DFLOAT, 2, shape1, cNArray);
    na_fill(rb_y, rb_rmiss);
    y = NA_PTR_TYPE(rb_y, double *);
    rb_vx = na_make_object(NA_DFLOAT, 2, shape0, cNArray);
    na_fill(rb_vx, rb_rmiss);
    vx = NA_PTR_TYPE(rb_vx, double *);
    rb_vy = na_make_object(NA_DFLOAT, 2, shape0, cNArray);
    na_fill(rb_vy, rb_rmiss);
    vy = NA_PTR_TYPE(rb_vy, double *);
    rb_score = na_make_object(NA_DFLOAT, 2, shape0, cNArray);
    na_fill(rb_score, rb_rmiss);
    score = NA_PTR_TYPE(rb_score, double *);

    if (rb_out_subimage == Qtrue){
        shape2[0] = op->o.nsx;
        shape2[1] = op->o.nsy;
        shape2[2] = shape1[0];
        shape2[3] = shape1[1];
        rb_zss = na_make_object(NA_SFLOAT, 4, shape2, cNArray);
        zss = NA_PTR_TYPE(rb_zss, float *);
    } else {
        rb_zss = Qnil;
        zss = NULL;
    }

    if (rb_out_score_ary == Qtrue){
        shape2[0] = 2 * op->o.ixhw + 1;
        shape2[1] = 2 * op->o.iyhw + 1;
        shape2[2] = shape0[0];
        shape2[3] = shape0[1];
        rb_score_ary = na_make_object(NA_DFLOAT, 4, shape2, cNArray);
        score_ary = NA_PTR_TYPE(rb_score_ary, double *);
    } else {
        rb_score_ary = Qnil;
        score_ary = NULL;
    }

    c_vttrac_trac(&op->o,
                  len, tid0, x0, y0, vx0g, vy0g,
                  count, tid, x, y, vx, vy, score, zss, score_ary);

    return( rb_ary_new_from_args(9,
                                 rb_count, rb_tid, rb_x, rb_y, rb_vx, rb_vy,
                                 rb_score, rb_zss, rb_score_ary) );
}

static VALUE
rb_vttrac_core_get_rmiss(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM(op->o.rmiss));
}

static VALUE
rb_vttrac_core_set_rmiss(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.rmiss = NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_imiss(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(INT2NUM(op->o.imiss));
}

static VALUE
rb_vttrac_core_set_imiss(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.imiss = NUM2INT(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_nsx(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(INT2NUM(op->o.nsx));
}

static VALUE
rb_vttrac_core_set_nsx(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.nsx = NUM2INT(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_nsy(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(INT2NUM(op->o.nsy));
}

static VALUE
rb_vttrac_core_set_nsy(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.nsy = NUM2INT(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_vxch(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM(op->o.vxch));
}

static VALUE
rb_vttrac_core_set_vxch(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.vxch = NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_vych(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM(op->o.vych));
}

static VALUE
rb_vttrac_core_set_vych(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.vych = NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_itstep(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(INT2NUM(op->o.itstep));
}

static VALUE
rb_vttrac_core_set_itstep(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.itstep = NUM2INT(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_ntrac(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(INT2NUM(op->o.ntrac));
}

static VALUE
rb_vttrac_core_set_ntrac(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.ntrac = NUM2INT(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_subgrid(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(op->o.subgrid ? Qtrue : Qfalse);
}

static VALUE
rb_vttrac_core_set_subgrid(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.subgrid = ( val == Qtrue );
    return(val);
}

static VALUE
rb_vttrac_core_get_subgrid_gaus(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(op->o.subgrid_gaus ? Qtrue : Qfalse);
}

static VALUE
rb_vttrac_core_set_subgrid_gaus(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.subgrid_gaus = ( val == Qtrue );
    return(val);
}

static VALUE
rb_vttrac_core_get_score_method(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(INT2NUM(op->o.score_method));
}

static VALUE
rb_vttrac_core_set_score_method(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.score_method = NUM2INT(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_score_th0(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM(op->o.score_th0));
}

static VALUE
rb_vttrac_core_set_score_th0(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.score_th0 = NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_score_th1(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM(op->o.score_th1));
}

static VALUE
rb_vttrac_core_set_score_th1(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.score_th1 = NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_peak_inside_th(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM( (double) op->o.peak_inside_th));
}

static VALUE
rb_vttrac_core_set_peak_inside_th(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.peak_inside_th = (float) NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_min_contrast(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(DBL2NUM( (double) op->o.min_contrast));
}

static VALUE
rb_vttrac_core_set_min_contrast(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.min_contrast = (float) NUM2DBL(val);
    return(val);
}

static VALUE
rb_vttrac_core_get_use_init_temp(VALUE obj)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    return(op->o.use_init_temp ? Qtrue : Qfalse);
}

static VALUE
rb_vttrac_core_set_use_init_temp(VALUE obj, VALUE val)
{
    struct vttrac_cdata *op;
    TypedData_Get_Struct(obj, struct vttrac_cdata, &vttrac_cdata_type, op);
    op->o.use_init_temp = ( val == Qtrue );
    return(val);
}

void
Init_vttrac_core()
{
    static VALUE mNumRu, cVTTracCore;

    mNumRu = rb_define_module("NumRu");
    cVTTracCore = rb_define_class_under(mNumRu, "VTTracCore", rb_cObject);
    rb_define_alloc_func(cVTTracCore, rb_vttrac_core_alloc);
    rb_define_method(cVTTracCore, "initialize", rb_vttrac_core_initialize, 5);
    rb_define_method(cVTTracCore, "set_ixyhw_from_v",
                     rb_vttrac_core_set_ixyhw_from_v, 2);
    rb_define_method(cVTTracCore, "set_ixyhw_directly",
                     rb_vttrac_core_set_ixyhw_directly, 2);
    rb_define_method(cVTTracCore, "set_basic", rb_vttrac_core_set_basic, 4);
    rb_define_method(cVTTracCore, "set_optional",
                     rb_vttrac_core_set_optional, 10);
    rb_define_method(cVTTracCore, "trac", rb_vttrac_core_trac, 7);

    rb_define_method(cVTTracCore, "t", rb_vttrac_core_t, 0);
    rb_define_method(cVTTracCore, "z", rb_vttrac_core_z, 0);
    rb_define_method(cVTTracCore, "rmiss", rb_vttrac_core_get_rmiss, 0);
    rb_define_method(cVTTracCore, "rmiss=", rb_vttrac_core_set_rmiss, 1);
    rb_define_method(cVTTracCore, "imiss", rb_vttrac_core_get_imiss, 0);
    rb_define_method(cVTTracCore, "imiss=", rb_vttrac_core_set_imiss, 1);
    rb_define_method(cVTTracCore, "nsx", rb_vttrac_core_get_nsx, 0);
    rb_define_method(cVTTracCore, "nsx=", rb_vttrac_core_set_nsx, 1);
    rb_define_method(cVTTracCore, "nsy", rb_vttrac_core_get_nsy, 0);
    rb_define_method(cVTTracCore, "nsy=", rb_vttrac_core_set_nsy, 1);
    rb_define_method(cVTTracCore, "vxch", rb_vttrac_core_get_vxch, 0);
    rb_define_method(cVTTracCore, "vxch=", rb_vttrac_core_set_vxch, 1);
    rb_define_method(cVTTracCore, "vych", rb_vttrac_core_get_vych, 0);
    rb_define_method(cVTTracCore, "vych=", rb_vttrac_core_set_vych, 1);
    rb_define_method(cVTTracCore, "itstep", rb_vttrac_core_get_itstep, 0);
    rb_define_method(cVTTracCore, "itstep=", rb_vttrac_core_set_itstep, 1);
    rb_define_method(cVTTracCore, "ntrac", rb_vttrac_core_get_ntrac, 0);
    rb_define_method(cVTTracCore, "ntrac=", rb_vttrac_core_set_ntrac, 1);
    rb_define_method(cVTTracCore, "subgrid", rb_vttrac_core_get_subgrid, 0);
    rb_define_method(cVTTracCore, "subgrid=", rb_vttrac_core_set_subgrid, 1);
    rb_define_method(cVTTracCore, "subgrid_gaus", rb_vttrac_core_get_subgrid_gaus, 0);
    rb_define_method(cVTTracCore, "subgrid_gaus=", rb_vttrac_core_set_subgrid_gaus, 1);
    rb_define_method(cVTTracCore, "score_method", rb_vttrac_core_get_score_method, 0);
    rb_define_method(cVTTracCore, "score_method=", rb_vttrac_core_set_score_method, 1);
    rb_define_method(cVTTracCore, "score_th0", rb_vttrac_core_get_score_th0, 0);
    rb_define_method(cVTTracCore, "score_th0=", rb_vttrac_core_set_score_th0, 1);
    rb_define_method(cVTTracCore, "score_th1", rb_vttrac_core_get_score_th1, 0);
    rb_define_method(cVTTracCore, "score_th1=", rb_vttrac_core_set_score_th1, 1);
    rb_define_method(cVTTracCore, "peak_inside_th", rb_vttrac_core_get_peak_inside_th, 0);
    rb_define_method(cVTTracCore, "peak_inside_th=", rb_vttrac_core_set_peak_inside_th, 1);
    rb_define_method(cVTTracCore, "min_contrast", rb_vttrac_core_get_min_contrast, 0);
    rb_define_method(cVTTracCore, "min_contrast=", rb_vttrac_core_set_min_contrast, 1);
    rb_define_method(cVTTracCore, "use_init_temp", rb_vttrac_core_get_use_init_temp, 0);
    rb_define_method(cVTTracCore, "use_init_temp=", rb_vttrac_core_set_use_init_temp, 1);
}
