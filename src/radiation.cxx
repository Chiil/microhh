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

#include "radiation.h"
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"

Radiation::Radiation(Master* masterin, Input* input, Grid* gridin, Fields* fieldsin)
{
    master = masterin;
    grid   = gridin;
    fields = fieldsin;

    rad_tend = 0;

    // Read input options out of ini file
    input->get_item(&swradiation, "radiation", "swradiation", "", "0");
}

Radiation::~Radiation()
{
    delete[] rad_tend;
}

void Radiation::init()
{
    rad_tend = new double[grid->kcells];
}

void Radiation::create(Input *inputin)
{
    int nerror = 0;

    // Read profiles from input
    nerror += inputin->get_prof(&rad_tend[grid->kstart], "rad_tend", grid->kmax);
    
    if (nerror)
        throw 1;
}

#ifndef USECUDA
void Radiation::exec()
{
    if (swradiation == "0")
        return;

    calc_radiation_tendency(fields->st["th"]->data, rad_tend);
}
#endif

void Radiation::calc_radiation_tendency(double* restrict tht, double* restrict rad_tend)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            #pragma GCC ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                tht[ijk] += rad_tend[k];
            }
}
