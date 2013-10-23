#ifndef STATS_LES
#define STATS_LES

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

class cstats_les
{
  public:
    cstats_les(cgrid *, cfields *, cmpi *);
    ~cstats_les();

    int readinifile(cinput *);
    int init();
    int create(int);
    int exec(int, double);

  private:
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    bool allocated;
    bool initialized;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var, *iter_var;
    NcVar  *u_var , *v_var, *w_var, *s_var;
    NcVar  *evisc_var;
    NcVar  *u2_var, *v2_var, *w2_var, *s2_var;
    NcVar  *u2_patch_var, *v2_patch_var, *w2_patch_var, *s2_patch_var;
    NcVar  *u2_nopatch_var, *v2_nopatch_var, *w2_nopatch_var, *s2_nopatch_var;
    NcVar  *u3_var, *v3_var, *w3_var, *s3_var;
    NcVar  *ugrad_var, *vgrad_var, *sgrad_var;
    NcVar  *wu_var, *wv_var, *ws_var;
    NcVar  *udiff_var, *vdiff_var, *sdiff_var;
    NcVar  *uflux_var, *vflux_var, *sflux_var;

    double *u , *v , *w , *s;
    double *uabs, *vabs;
    double *evisc;
    double *u2, *v2, *w2, *s2;
    int *nfilter;
    double *u2_patch, *v2_patch, *w2_patch, *s2_patch;
    double *u2_nopatch, *v2_nopatch, *w2_nopatch, *s2_nopatch;
    double *u3, *v3, *w3, *s3;
    double *wu , *wv , *ws;
    double *ugrad, *vgrad, *sgrad;
    double *udiff, *vdiff, *sdiff;
    double *uflux, *vflux, *sflux;

    int calcfilter       (double *, int *, double *, double *);
    int calcmoment_filter(double *, double *, double *, double *, double *, int *, double, int);

    int calcmean     (double *, double *, double);
    int calcmoment   (double *, double *, double *, double, int);
    int calcdiff     (double *, double *, double *, double *, double *, double *, double);
    int calcgrad     (double *, double *, double *);
    int calcflux     (double *, double *, double *, double *, int, int);
    int calctkebudget(double *, double *, double *, double *, double *,
                      double *, double *,
                      double *, double *,
                      double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *,
                      double *, double *,
                      double *, double *, double *,
                      double *, double *,
                      double *, double *, double);

    int nstats;
};
#endif
