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
    struct East_west_face
    {
        int i;
        int jstart;
        int jend;
        int kstart;
        int kend;
    };

    struct North_south_face
    {
        int istart;
        int iend;
        int j;
        int kstart;
        int kend;
    };

    struct Top_bottom_face
    {
        int istart;
        int iend;
        int jstart;
        int jend;
        int k;
    };

    std::vector<East_west_face> west_faces;
    std::vector<East_west_face> east_faces;

    std::vector<North_south_face> south_faces;
    std::vector<North_south_face> north_faces;

    std::vector<Top_bottom_face> bottom_faces;
    std::vector<Top_bottom_face> top_faces;

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
        for (int m=0; m<mblocks; ++m)
        {
            // Calculate the absolute grid indices of the faces.
            const int iface_start = m*istep + istep/2 - iblock/2;
            const int iface_end = m*istep + istep/2 + iblock/2;
            const int jface_start = n*jstep + jstep/2 - jblock/2;
            const int jface_end = n*jstep + jstep/2 + jblock/2;
            const int kface_start = 0;
            const int kface_end = kblock;

            // Check the ranges of i and j for the specific MPI process.
            const int imin_abs = master.mpicoordx*grid.imax;
            const int imax_abs = imin_abs + grid.imax;
            const int jmin_abs = master.mpicoordy*grid.jmax;
            const int jmax_abs = jmin_abs + grid.jmax;

            // Check whether there is an edge in range.
            const bool iface_start_in_range = (iface_start >= imin_abs) && (iface_start <  imax_abs);
            const bool iface_end_in_range   = (iface_end   >= imin_abs) && (iface_end   <  imax_abs);
            const bool iface_full_in_range  = (iface_start <  imin_abs) && (iface_end   >= imax_abs);

            const bool jface_start_in_range = (jface_start >= jmin_abs) && (jface_start <  jmax_abs);
            const bool jface_end_in_range   = (jface_end   >= jmin_abs) && (jface_end   <  jmax_abs);
            const bool jface_full_in_range  = (jface_start <  jmin_abs) && (jface_end   >= jmax_abs);

            // EAST-WEST FACES
            if ( (jface_start_in_range || jface_end_in_range || jface_full_in_range) )
            {
                if ( !(iface_start < imin_abs || iface_start >= imax_abs) )
                {
                    // Store the part of the face that is in range and add ghost cells.
                    East_west_face west_face;

                    west_face.i = iface_start%grid.imax + grid.igc;

                    west_face.jstart = (jface_start_in_range ? jface_start%grid.jmax : 0) + grid.jgc;
                    west_face.jend   = (jface_end_in_range ? jface_end%grid.jmax : grid.jmax) + grid.jgc;
                    west_face.kstart = kface_start + grid.kgc;
                    west_face.kend   = kface_end   + grid.kgc;

                    west_faces.push_back(west_face);
                }

                if ( !(iface_end < imin_abs || iface_end >= imax_abs) )
                {
                    // Store the part of the face that is in range and add ghost cells.
                    East_west_face east_face;

                    east_face.i = iface_end%grid.imax + grid.igc;

                    east_face.jstart = (jface_start_in_range ? jface_start%grid.jmax : 0) + grid.jgc;
                    east_face.jend   = (jface_end_in_range ? jface_end%grid.jmax : grid.jmax) + grid.jgc;
                    east_face.kstart = kface_start + grid.kgc;
                    east_face.kend   = kface_end + grid.kgc;

                    east_faces.push_back(east_face);
                }
            }

            // NORTH-SOUTH FACES
            if ( (iface_start_in_range || iface_end_in_range || iface_full_in_range) )
            {
                if ( !(jface_start < jmin_abs || jface_start >= jmax_abs) )
                {
                    // Store the part of the face that is in range and add ghost cells.
                    North_south_face south_face;

                    south_face.j = jface_start%grid.jmax + grid.jgc;

                    south_face.istart = (iface_start_in_range ? iface_start%grid.imax : 0) + grid.igc;
                    south_face.iend   = (iface_end_in_range ? iface_end%grid.imax : grid.imax) + grid.igc;
                    south_face.kstart = kface_start + grid.kgc;
                    south_face.kend   = kface_end + grid.kgc;

                    south_faces.push_back(south_face);
                }

                if ( !(jface_end < jmin_abs || jface_end >= jmax_abs) )
                {
                    // Store the part of the face that is in range and add ghost cells.
                    North_south_face north_face;

                    north_face.j = jface_end%grid.jmax + grid.jgc;

                    north_face.istart = (iface_start_in_range ? iface_start%grid.imax : 0) + grid.igc;
                    north_face.iend   = (iface_end_in_range ? iface_end%grid.imax : grid.imax) + grid.igc;
                    north_face.kstart = kface_start + grid.kgc;
                    north_face.kend   = kface_end + grid.kgc;

                    north_faces.push_back(north_face);
                }
            }

            // TOP-DOWN FACES
            if (false)
            {
                if ( !(jface_start < jmin_abs || jface_start >= jmax_abs) )
                {
                    // Store the part of the face that is in range and add ghost cells.
                    Top_bottom_face bottom_face;

                    bottom_face.k = kface_start + grid.kgc;

                    bottom_face.istart = (iface_start_in_range ? iface_start%grid.imax : 0) + grid.igc;
                    bottom_face.iend   = (iface_end_in_range ? iface_end%grid.imax : grid.imax) + grid.igc;
                    bottom_face.jstart = (jface_start_in_range ? jface_start%grid.jmax : 0) + grid.jgc;
                    bottom_face.jend   = (jface_end_in_range ? jface_end%grid.jmax : grid.jmax) + grid.jgc;

                    bottom_faces.push_back(bottom_face);
                }

                if ( !(jface_start < jmin_abs || jface_start >= jmax_abs) )
                {
                    // Store the part of the face that is in range and add ghost cells.
                    Top_bottom_face top_face;

                    top_face.k = kface_end + grid.kgc;

                    top_face.istart = (iface_start_in_range ? iface_start%grid.imax : 0) + grid.igc;
                    top_face.iend   = (iface_end_in_range ? iface_end%grid.imax : grid.imax) + grid.igc;
                    top_face.jstart = (jface_start_in_range ? jface_start%grid.jmax : 0) + grid.jgc;
                    top_face.jend   = (jface_end_in_range ? jface_end%grid.jmax : grid.jmax) + grid.jgc;

                    top_faces.push_back(top_face);
                }
            }
        }
}

void Immersed_boundary::exec(Fields& fields)
{
    if (swib != "1")
        return;

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
}
