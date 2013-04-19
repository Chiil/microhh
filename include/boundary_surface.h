#ifndef BOUNDARY_SURFACE
#define BOUNDARY_SURFACE

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "boundary.h"

class cboundary_surface : public cboundary
{
  public:
    cboundary_surface(cgrid *, cfields *, cmpi *);

    int readinifile(cinput *);
    int init();
    int setvalues();
    int exec();

    double *obuk;
    double *ustar;

  private:
    int setbc_patch(double *, double, double, double);

    // patch type
    int    patch_dim;
    double patch_xh;
    double patch_xr;
    double patch_xi;
    double patch_facr;
    double patch_facl;

    // surface scheme
    int bcvalues();
    int stability(double *, double *, double *,
                  double *, double *, double *,
                  double *, double *, double *,
                  double *, double *);
    int surfm(double *, double *,
              double *, double *, double *, double *,
              double *, double *, double *, double *,
              double, int);
    int surfs(double *, double *, double *,
              double *, double *, double *,
              double, int);
    double calcobuk(double, double, double, double);
    inline double fm(double, double, double);
    inline double fh(double, double, double);
    inline double psim(double);
    inline double psih(double);
    inline double phim(double);
    inline double phih(double);
    double ustarin;
    double z0m;
    double z0h;

    typedef std::map<std::string, int> bcbotmap;
    int surfmbcbot;
    bcbotmap surfsbcbot;
};
#endif
