/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "immersed_boundary.h"

namespace
{
    std::string swib;
    int mblocks;
    int nblocks;
    int iblock;
    int jblock;
    int kblock;

    void set_no_penetration(double* const restrict ut, double* const restrict vt, double* const restrict wt,
                            double* const restrict u, double* const restrict v, double* const restrict w,
                            const int istart, const int iend,
                            const int jstart, const int jend,
                            const int kstart, const int kend,
                            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const int ibc_kstart = kstart;
        const int ibc_kend   = kstart + kblock;

        const int istep = (iend-istart)/mblocks;
        const int jstep = (jend-jstart)/nblocks;

        for (int n=0; n<nblocks; ++n)
        {
            for (int m=0; m<mblocks; ++m)
            {
                const int ibc_istart = istart + m*istep + istep/2 - iblock/2;
                const int ibc_iend   = ibc_istart + iblock;

                const int ibc_jstart = jstart + n*jstep + jstep/2 - jblock/2;
                const int ibc_jend   = ibc_jstart + jblock;

                // Set the u ghost cells, no flow in the block.
                for (int k=ibc_kstart; k<ibc_kend; ++k)
                    for (int j=ibc_jstart; j<ibc_jend; ++j)
                    {
                        const int ijk_istart  = ibc_istart + j*jj + k*kk;
                        const int ijk_iend    = ibc_iend   + j*jj + k*kk;
                        u [ijk_istart] = 0.;
                        u [ijk_iend  ] = 0.;
                        ut[ijk_istart] = 0.;
                        ut[ijk_iend  ] = 0.;
                    }

                // Set the u ghost cells, no flow in the block.
                for (int k=ibc_kstart; k<ibc_kend; ++k)
                    for (int i=ibc_istart; i<ibc_iend; ++i)
                    {
                        const int ijk_jstart  = i + ibc_jstart*jj + k*kk;
                        const int ijk_jend    = i + ibc_jend  *jj + k*kk;
                        v [ijk_jstart] = 0.;
                        v [ijk_jend  ] = 0.;
                        vt[ijk_jstart] = 0.;
                        vt[ijk_jend  ] = 0.;
                    }

                // Set the w ghost cells, no flow in the block.
                for (int j=ibc_jstart; j<ibc_jend; ++j)
                    for (int i=ibc_istart; i<ibc_iend; ++i)
                    {
                        const int ijk_kstart = i + j*jj + ibc_kstart*kk;
                        const int ijk_kend   = i + j*jj + ibc_kend  *kk;
                        w [ijk_kstart] = 0.;
                        w [ijk_kend  ] = 0.;
                        wt[ijk_kstart] = 0.;
                        wt[ijk_kend  ] = 0.;
                    }
            }
        }
    }

    void set_no_slip(double* const restrict ut, double* const restrict vt, double* const restrict wt,
                     const double* const restrict u, const double* const restrict v, const double* const restrict w,
                     const double* const restrict rhoref, const double* const restrict rhorefh,
                     const double* const restrict dzi, const double* const restrict dzhi,
                     const double dxi, const double dyi,
                     const double visc,
                     const int istart, const int iend,
                     const int jstart, const int jend,
                     const int kstart, const int kend,
                     const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const double dxidxi = dxi*dxi;
        const double dyidyi = dyi*dyi;

        const int ibc_kstart = kstart;
        const int ibc_kend   = kstart + kblock;

        const int istep = (iend-istart)/nblocks;
        const int jstep = (jend-jstart)/nblocks;

        for (int n=0; n<nblocks; ++n)
        {
            for (int m=0; m<mblocks; ++m)
            {
                const int ibc_istart = istart + m*istep + istep/2 - iblock/2;
                const int ibc_iend   = ibc_istart + iblock;

                const int ibc_jstart = jstart + n*jstep + jstep/2 - jblock/2;
                const int ibc_jend   = ibc_jstart + jblock;

                // Set the w no slip at the vertical walls, by reverting the advection 
                // and diffusion towards the wall and adding the proper diffusion
                for (int k=ibc_kstart; k<ibc_kend+1; ++k)
                    for (int j=ibc_jstart; j<ibc_jend; ++j)
                    {
                        int ijk = ibc_istart-1 + j*jj + k*kk;
                        wt[ijk] +=
                                + ( interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk], w[ijk+ii]) ) * dxi
                                - visc * ( (w[ijk+ii] - w[ijk]) ) * dxidxi
                                + visc * ( -2.*w[ijk] ) * dxidxi;

                        ijk = ibc_iend + j*jj + k*kk;
                        wt[ijk] +=
                                - ( interp2(u[ijk-kk], u[ijk]) * interp2(w[ijk-ii], w[ijk]) ) * dxi
                                + visc * ( (w[ijk] - w[ijk-ii]) ) * dxidxi
                                - visc * ( 2.*w[ijk] ) * dxidxi;
                    }

                for (int k=ibc_kstart; k<ibc_kend+1; ++k)
                    for (int i=ibc_istart; i<ibc_iend; ++i)
                    {
                        int ijk = i + (ibc_jstart-1)*jj + k*kk;
                        wt[ijk] +=
                                + ( interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk], w[ijk+jj]) ) * dyi
                                - visc * ( (w[ijk+jj] - w[ijk]) ) * dyidyi
                                + visc * ( -2.*w[ijk] ) * dyidyi;

                        ijk = i + ibc_jend*jj + k*kk;
                        wt[ijk] +=
                                - ( interp2(v[ijk-kk], v[ijk]) * interp2(w[ijk-jj], w[ijk]) ) * dyi
                                + visc * ( (w[ijk] - w[ijk-jj]) ) * dyidyi
                                - visc * ( 2.*w[ijk] ) * dyidyi;
                    }

                // Set the u no slip at the horizontal walls
                for (int j=ibc_jstart; j<ibc_jend; ++j)
                    for (int i=ibc_istart; i<ibc_iend+1; ++i)
                    {
                        int k = ibc_kstart-1;
                        int ijk = i + j*jj + k*kk;
                        ut[ijk] +=
                                + ( rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk], u[ijk+kk]) ) / rhoref[k] * dzi[k]
                                - visc * ( (u[ijk+kk] - u[ijk]) * dzhi[k+1]) * dzi[k]
                                + visc * ( -2.*u[ijk] * dzhi[k+1] ) * dzi[k];

                        k = ibc_kend;
                        ijk = i + j*jj + k*kk;
                        ut[ijk] +=
                                - ( rhorefh[k] * interp2(w[ijk-ii], w[ijk]) * interp2(u[ijk-kk], u[ijk]) ) / rhoref[k] * dzi[k]
                                + visc * ( (u[ijk] - u[ijk-kk]) * dzhi[k] ) * dzi[k]
                                - visc * ( 2.*u[ijk] * dzhi[k] ) * dzi[k];
                    }

                // Set the v no slip at the horizontal walls
                for (int j=ibc_jstart; j<ibc_jend+1; ++j)
                    for (int i=ibc_istart; i<ibc_iend; ++i)
                    {
                        int k = ibc_kstart-1;
                        int ijk = i + j*jj + k*kk;
                        vt[ijk] +=
                                + ( rhorefh[k+1] * interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk], v[ijk+kk]) ) / rhoref[k] * dzi[k]
                                - visc * ( (v[ijk+kk] - v[ijk]) * dzhi[k+1]) * dzi[k]
                                + visc * ( -2.*v[ijk] * dzhi[k+1] ) * dzi[k];

                        k = ibc_kend;
                        ijk = i + j*jj + k*kk;
                        vt[ijk] +=
                                - ( rhorefh[k] * interp2(w[ijk-ii], w[ijk]) * interp2(v[ijk-kk], v[ijk]) ) / rhoref[k] * dzi[k]
                                + visc * ( (v[ijk] - v[ijk-kk]) * dzhi[k] ) * dzi[k]
                                - visc * ( 2.*v[ijk] * dzhi[k] ) * dzi[k];
                    }
            }
        }
    }

    void set_scalar(double* const restrict st, double* const restrict s,
                    const double* const restrict u, const double* const restrict v, const double* const restrict w,
                    const double* const restrict rhoref, const double* const restrict rhorefh,
                    const double* const restrict dzi, const double* const restrict dzhi,
                    const double dxi, const double dyi,
                    const double visc,
                    const int istart, const int iend,
                    const int jstart, const int jend,
                    const int kstart, const int kend,
                    const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const double dxidxi = dxi*dxi;
        const double dyidyi = dyi*dyi;

        const int ibc_kstart = kstart;
        const int ibc_kend   = kstart + kblock;

        const int istep = (iend-istart)/mblocks;
        const int jstep = (jend-jstart)/nblocks;

        for (int n=0; n<nblocks; ++n)
        {
            for (int m=0; m<mblocks; ++m)
            {
                const int ibc_istart = istart + m*istep + istep/2 - iblock/2;
                const int ibc_iend   = ibc_istart + iblock;

                const int ibc_jstart = jstart + n*jstep + jstep/2 - jblock/2;
                const int ibc_jend   = ibc_jstart + jblock;

                // Set no flow through the object at the vertical wall and a neumann BC.
                for (int k=ibc_kstart; k<ibc_kend; ++k)
                    for (int j=ibc_jstart; j<ibc_jend; ++j)
                    {
                        int ijk = ibc_istart-1 + j*jj + k*kk;
                        st[ijk] +=
                                 + ( u[ijk+ii] * interp2(s[ijk], s[ijk+ii]) ) * dxi
                                 - visc * ( s[ijk+ii] - s[ijk] ) * dxidxi;

                        ijk = ibc_iend + j*jj + k*kk;
                        st[ijk] +=
                                 - ( u[ijk] * interp2(s[ijk-ii], s[ijk]) ) * dxi
                                 + visc * ( (s[ijk] - s[ijk-ii]) ) * dxidxi;
                    }

                // Set no flow through the object at the vertical wall and a neumann BC.
                for (int k=ibc_kstart; k<ibc_kend; ++k)
                    for (int i=ibc_istart; i<ibc_iend; ++i)
                    {
                        int ijk = i + (ibc_jstart-1)*jj + k*kk;
                        st[ijk] +=
                                 + ( v[ijk+jj] * interp2(s[ijk], s[ijk+jj]) ) * dyi
                                 - visc * ( s[ijk+jj] - s[ijk] ) * dyidyi;

                        ijk = i + ibc_jend*jj + k*kk;
                        st[ijk] +=
                                 - ( v[ijk] * interp2(s[ijk-jj], s[ijk]) ) * dyi
                                 + visc * ( (s[ijk] - s[ijk-jj]) ) * dyidyi;
                    }

                // Set no flow through the object at the horizontal wall
                for (int j=ibc_jstart; j<ibc_jend; ++j)
                    for (int i=ibc_istart; i<ibc_iend; ++i)
                    {
                        int k = ibc_kstart-1;
                        int ijk = i + j*jj + k*kk;
                        st[ijk] +=
                                 + ( rhorefh[k+1] * w[ijk+kk] * interp2(s[ijk], s[ijk+kk]) ) / rhoref[k] * dzi[k]
                                 - visc * (s[ijk+kk] - s[ijk]) * dzhi[k+1] * dzi[k];

                        k = ibc_kend;
                        ijk = i + j*jj + k*kk;
                        st[ijk] +=
                                 - ( rhorefh[k] * w[ijk] * interp2(s[ijk-kk], s[ijk]) ) / rhoref[k] * dzi[k]
                                 + visc * (s[ijk] - s[ijk-kk]) * dzhi[k] * dzi[k];
                    }
            }
        }
    }

}

Immersed_boundary::Immersed_boundary(Master& masterin, Grid& gridin, Input& input) :
    master(masterin),
    grid(gridin)
{
    int nerror = 0;
    nerror += input.get_item(&swib, "ib", "swib", "");

    if (swib != "1")
    {
        input.flag_as_used("ib", "mblocks");
        input.flag_as_used("ib", "nblocks");
        input.flag_as_used("ib", "iblock");
        input.flag_as_used("ib", "jblock");
        input.flag_as_used("ib", "kblock");
    }
    else
    {
        nerror += input.get_item(&mblocks, "ib", "mblocks", "");
        nerror += input.get_item(&nblocks, "ib", "nblocks", "");
        nerror += input.get_item(&iblock , "ib", "iblock" , "");
        nerror += input.get_item(&jblock , "ib", "jblock" , "");
        nerror += input.get_item(&kblock , "ib", "kblock" , "");
    }

    if (nerror > 0)
        throw 1;
}

Immersed_boundary::~Immersed_boundary()
{
}

void Immersed_boundary::exec(Fields& fields)
{
    if (swib != "1")
        return;

    set_no_penetration(fields.ut->data, fields.vt->data, fields.wt->data,
                       fields.u->data, fields.v->data, fields.w->data,
                       grid.istart, grid.iend,
                       grid.jstart, grid.jend,
                       grid.kstart, grid.kend,
                       grid.icells, grid.ijcells);

    set_no_slip(fields.ut->data, fields.vt->data, fields.wt->data,
                fields.u->data, fields.v->data, fields.w->data,
                fields.rhoref, fields.rhorefh,
                grid.dzi, grid.dzhi,
                grid.dxi, grid.dyi,
                fields.visc,
                grid.istart, grid.iend,
                grid.jstart, grid.jend,
                grid.kstart, grid.kend,
                grid.icells, grid.ijcells);

    for (FieldMap::const_iterator it = fields.st.begin(); it!=fields.st.end(); it++)
        set_scalar(it->second->data, fields.sp[it->first]->data,
                   fields.u->data, fields.v->data, fields.w->data,
                   fields.rhoref, fields.rhorefh,
                   grid.dzi, grid.dzhi,
                   grid.dxi, grid.dyi,
                   fields.sp[it->first]->visc,
                   grid.istart, grid.iend,
                   grid.jstart, grid.jend,
                   grid.kstart, grid.kend,
                   grid.icells, grid.ijcells);
}
