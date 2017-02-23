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
    struct West_face
    {
        int i;
        int jstart;
        int jend;
        int kstart;
        int kend;
    };

    std::vector<West_face> west_faces;

    std::string swib;
    int mblocks;
    int nblocks;
    int iblock;
    int jblock;
    int kblock;

    void set_west_face(double* const restrict ut, double* const restrict vt, double* const restrict wt,
                       double* const restrict u, double* const restrict v, double* const restrict w,
                       const int iface,
                       const int jstart, const int jend,
                       const int kstart, const int kend,
                       const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Enforce no penetration for u
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
             {
                 const int ijk = iface + j*jj + k*kk;
                 u [ijk] = 0.;
                 u [ijk] = 0.;
                 ut[ijk] = 0.;
                 ut[ijk] = 0.;
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
            const int iface = m*istep + istep/2 - iblock/2;
            const int jface_start = n*jstep + jstep/2 - jblock/2;
            const int jface_end = n*jstep + jstep/2 + jblock/2;
            const int kface_start = 0;
            const int kface_end = kblock;

            // Check the ranges of i and j.
            const int imin_abs = master.mpicoordx*grid.imax;
            const int imax_abs = imin_abs + grid.imax;
            const int jmin_abs = master.mpicoordy*grid.jmax;
            const int jmax_abs = jmin_abs + grid.jmax;

            if (iface < imin_abs || iface >= imax_abs)
                continue;

            const bool jface_start_in_range = (jface_start >= jmin_abs) && (jface_start < jmax_abs);
            const bool jface_end_in_range   = (jface_end   >= jmin_abs) && (jface_end   < jmax_abs);

            if ( !(jface_start_in_range || jface_end_in_range) )
                continue;

            // Store the part of the face that is in range and add ghost cells.
            West_face west_face;

            west_face.i = iface % grid.imax + grid.igc;

            west_face.jstart = (jface_start_in_range ? jface_start%grid.jmax : 0) + grid.jgc;
            west_face.jend   = (jface_end_in_range ? jface_end%grid.jmax : grid.jmax) + grid.jgc;
            west_face.kstart = kface_start + grid.kgc;
            west_face.kend   = kface_end   + grid.kgc;

            west_faces.push_back(west_face);
        }

    // for (West_face& wf : west_faces)
    //     std::printf("%d, %d, %d, %d, %d, %d\n", master.mpiid, wf.i, wf.jstart, wf.jend, wf.kstart, wf.kend);
}

void Immersed_boundary::exec(Fields& fields)
{
    if (swib != "1")
        return;

    for (West_face& face : west_faces)
        set_west_face(fields.ut->data, fields.vt->data, fields.wt->data,
                      fields.u->data, fields.v->data, fields.w->data,
                      face.i,
                      face.jstart, face.jend,
                      face.kstart, face.kend,
                      grid.icells, grid.ijcells);

    /*
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
                   */
}
