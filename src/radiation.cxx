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
#include <cmath>
#include <stdlib.h>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "model.h"
#include "radiation.h"

namespace
{
    std::string swradiation;

    std::vector<double> radmeanprof;
    std::vector<double> radflexprof;

    void add_mean_radiation(double* const tvart, const double* const radmeanprof,
                            const int istart, const int jstart, const int kstart,
                            const int iend, const int jend, const int kend,
                            const int jj, const int kk)
    {

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    tvart[ijk] += radmeanprof[k];
                }
    }

    double calc_precipitable_water(double* const pw, const double* const qt,
                                   const double* const rhoref, const double* const dz,
                                   const int istart, const int jstart, const int kstart,
                                   const int iend, const int jend, const int kend,
                                   const int jj, const int kk)
    {

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;
                    pw[ij] += rhoref[k] * dz[k] * qt[ijk];
                }

        return 0.;
    }

    void add_flex_radiation(double* const tvart, const double* const radflexprof,
                            const int istart, const int jstart, const int kstart,
                            const int iend, const int jend, const int kend,
                            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    tvart[ijk] += radflexprof[k];
                }
    }
}

Radiation::Radiation(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    int nerror = 0;
    nerror += inputin->get_item(&swradiation, "radiation", "swradiation", "", "0");

    if (nerror)
        throw 1;
}

Radiation::~Radiation()
{
#ifdef USECUDA
    clear_device();
#endif
}

void Radiation::init()
{
    if (swradiation == "1")
    {
        radmeanprof.resize(grid->kmax);
        radflexprof.resize(grid->kmax);
    }
}

void Radiation::create(Input* inputin)
{
    int nerror = 0;

    if (swradiation == "1")
    {
        // Load the flexible radiation profile.
        nerror += inputin->get_prof(radmeanprof.data(), "radmean", grid->kmax);
        nerror += inputin->get_prof(radflexprof.data(), "radflex", grid->kmax);
    }

    if (nerror)
        throw 1;
}

#ifndef USECUDA
void Radiation::exec()
{
    if (swradiation == "1")
    {
        add_mean_radiation(fields->st["thl"]->data, radmeanprof.data(),
                           grid->istart, grid->jstart, grid->kstart,
                           grid->iend, grid->jend, grid->kend,
                           grid->icells, grid->ijcells);

        // Calculate the precipitable water.
        double precipitable_water_mean = calc_precipitable_water(fields->atmp["tmp1"]->databot,
                                                                 fields->sp["qt"]->data,
                                                                 fields->rhoref, grid->dz,
                                                                 grid->istart, grid->jstart, grid->kstart,
                                                                 grid->iend, grid->jend, grid->kend,
                                                                 grid->icells, grid->ijcells);

        add_flex_radiation(fields->st["thl"]->data, 
                           radflexprof.data(),
                           grid->istart, grid->jstart, grid->kstart,
                           grid->iend, grid->jend, grid->kend,
                           grid->icells, grid->ijcells);
    }
}
#endif
