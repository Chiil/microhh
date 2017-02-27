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

    void set_east_west_face_no_penetration(double* const restrict ut,
                                           double* const restrict u,
                                           const int iface,
                                           const int jstart, const int jend,
                                           const int kstart, const int kend,
                                           const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Enforce no penetration for u.
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
             {
                 const int ijk = iface + j*jj + k*kk;
                 u [ijk] = 0.;
                 ut[ijk] = 0.;
             }
    }

    void set_north_south_face_no_penetration(double* const restrict vt,
                                             double* const restrict v,
                                             const int istart, const int iend,
                                             const int jface,
                                             const int kstart, const int kend,
                                             const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Enforce no penetration for v.
        for (int k=kstart; k<kend; ++k)
            for (int i=istart; i<iend; ++i)
             {
                 const int ijk = i + jface*jj + k*kk;
                 v [ijk] = 0.;
                 vt[ijk] = 0.;
             }
    }

    void set_top_bottom_face_no_penetration(double* const restrict wt,
                                            double* const restrict w,
                                            const int istart, const int iend,
                                            const int jstart, const int jend,
                                            const int kface,
                                            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Enforce no penetration for v.
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
             {
                 const int ijk = i + j*jj + kface*kk;
                 w [ijk] = 0.;
                 wt[ijk] = 0.;
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

void Immersed_boundary::create()
{
    if (swib != "1")
        return;

    // Set the west faces
    const int istep = grid.itot/mblocks;
    const int jstep = grid.jtot/nblocks;

    for (int n=0; n<nblocks; ++n)
    {
        for (int m=0; m<mblocks; ++m)
        {
            // Calculate the absolute grid indices of the blocks.
            // Note that we use c-style ranging, thus end is not part of the
            // range, and end-start is the number of points in the range.
            const int iblock_start = m*istep + istep/2 - iblock/2;
            const int iblock_end = m*istep + istep/2 + iblock/2;
            const int jblock_start = n*jstep + jstep/2 - jblock/2;
            const int jblock_end = n*jstep + jstep/2 + jblock/2;
            // const int iblock_start = m*istep;
            // const int iblock_end = iblock_start + iblock;
            // const int jblock_start = n*jstep;
            // const int jblock_end = jblock_start + jblock;
            const int kblock_start = 0;
            const int kblock_end = kblock_start + kblock;

            // Check the ranges of i and j for the specific MPI process.
            const int imin_abs = master.mpicoordx*grid.imax;
            const int imax_abs = (master.mpicoordx+1)*grid.imax;
            const int jmin_abs = master.mpicoordy*grid.jmax;
            const int jmax_abs = (master.mpicoordy+1)*grid.jmax;

            // Check whether there is an edge in range.
            const bool iblock_start_in_range = (iblock_start >= imin_abs) && (iblock_start <  imax_abs);
            const bool iblock_end_in_range   = (iblock_end-1 >= imin_abs) && (iblock_end-1 <  imax_abs);
            const bool iblock_fully_in_range = (iblock_start <  imin_abs) && (iblock_end-1 >= imax_abs);

            const bool jblock_start_in_range = (jblock_start >= jmin_abs) && (jblock_start <  jmax_abs);
            const bool jblock_end_in_range   = (jblock_end-1 >= jmin_abs) && (jblock_end-1 <  jmax_abs);
            const bool jblock_fully_in_range = (jblock_start <  jmin_abs) && (jblock_end-1 >= jmax_abs);

            if ( (iblock_start_in_range || iblock_end_in_range || iblock_fully_in_range) &&
                 (jblock_start_in_range || jblock_end_in_range || jblock_fully_in_range) )
            {
                const int istart = (iblock_start_in_range ? iblock_start : imin_abs) - master.mpicoordx*grid.imax + grid.igc;
                const int iend   = (iblock_end_in_range   ? iblock_end   : imax_abs) - master.mpicoordx*grid.imax + grid.igc;
                const int jstart = (jblock_start_in_range ? jblock_start : jmin_abs) - master.mpicoordy*grid.jmax + grid.jgc;
                const int jend   = (jblock_end_in_range   ? jblock_end   : jmax_abs) - master.mpicoordy*grid.jmax + grid.jgc;

                std::printf("CvH CHECK: (%d): %d, %d, %d, %d\n", master.mpiid, istart, iend, jstart, jend);
            }
        }
    }

    throw 1;
}

void Immersed_boundary::exec(Fields& fields)
{
    if (swib != "1")
        return;

    /*
    for (East_west_face& face : west_faces)
        set_east_west_face_no_penetration(fields.ut->data,
                                          fields.u->data,
                                          face.i,
                                          face.jstart, face.jend,
                                          face.kstart, face.kend,
                                          grid.icells, grid.ijcells);

    for (East_west_face& face : east_faces)
        set_east_west_face_no_penetration(fields.ut->data,
                                          fields.u->data,
                                          face.i,
                                          face.jstart, face.jend,
                                          face.kstart, face.kend,
                                          grid.icells, grid.ijcells);

    for (North_south_face& face : south_faces)
        set_north_south_face_no_penetration(fields.vt->data,
                                            fields.v->data,
                                            face.istart, face.iend,
                                            face.j,
                                            face.kstart, face.kend,
                                            grid.icells, grid.ijcells);

    for (North_south_face& face : north_faces)
        set_north_south_face_no_penetration(fields.vt->data,
                                            fields.v->data,
                                            face.istart, face.iend,
                                            face.j,
                                            face.kstart, face.kend,
                                            grid.icells, grid.ijcells);

    for (Top_bottom_face& face : bottom_faces)
        set_top_bottom_face_no_penetration(fields.wt->data,
                                           fields.w->data,
                                           face.istart, face.iend,
                                           face.jstart, face.jend,
                                           face.k,
                                           grid.icells, grid.ijcells);

    for (Top_bottom_face& face : top_faces)
        set_top_bottom_face_no_penetration(fields.wt->data,
                                           fields.w->data,
                                           face.istart, face.iend,
                                           face.jstart, face.jend,
                                           face.k,
                                           grid.icells, grid.ijcells);
                                           */
}
