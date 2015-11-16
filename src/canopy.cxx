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
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "canopy.h"

Canopy::Canopy(Model* modelin, Input* inputin)
{
    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    relative_pad = 0;
    pad = 0;

    int nerror = 0;
    nerror += inputin->get_item(&swcanopy, "canopy", "swcanopy", "", "0");
    nerror += inputin->get_item(&pai, "canopy", "pai", "");
    nerror += inputin->get_item(&Cd,  "canopy", "Cd", "");
}

Canopy::~Canopy()
{
    delete[] relative_pad;
    delete[] pad;
}

void Canopy::init()
{
    if (swcanopy == "1")
    {
        relative_pad = new double[grid->kcells];
        pad          = new double[grid->kcells];
    }
}

void Canopy::create(Input* inputin)
{
    if (swcanopy == "1")
    {
        int nerror = 0;
        nerror += inputin->get_prof(&relative_pad[grid->kstart], "relative_pad", grid->kmax);

        for (int k=grid->kstart; k<grid->kend; ++k)
            pad[k] = pai * relative_pad[k];

        if (nerror)
            throw 1;
    }
}

#ifndef USECUDA
void Canopy::exec()
{
    if (swcanopy == "1")
    {
        calc_canopy_drag(fields->ut->data, fields->vt->data, fields->wt->data,
                         fields->u ->data, fields->v ->data, fields->w ->data,
                         pad);
    }
}
#endif

void Canopy::calc_canopy_drag(double* const restrict ut, double* const restrict vt, double* const restrict wt,
                              const double* const restrict u, const double* const restrict v, const double* const restrict w,
                              const double* const pad)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double utrans = grid->utrans;
    const double vtrans = grid->vtrans;

    const double Cd = this->Cd;

    // Calculate drag for u-component.
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                const double u2 = std::pow(u[ijk] + utrans, 2);
                const double v2 = std::pow(0.25*(v[ijk-ii] + v[ijk-ii+jj] + v[ijk] + v[ijk+jj]) + vtrans, 2);
                const double w2 = std::pow(0.25*(w[ijk-ii] + w[ijk-ii+kk] + w[ijk] + w[ijk+kk]), 2);
                const double utot = std::sqrt(u2 + v2 + w2);
                ut[ijk] -= Cd * pad[k] * utot * (u[ijk]+utrans);
            }

    // Calculate drag for v-component.
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                const double u2 = std::pow(0.25*(u[ijk-jj] + u[ijk+ii-jj] + u[ijk] + u[ijk+ii]) + utrans, 2);
                const double v2 = std::pow(v[ijk] + vtrans, 2);
                const double w2 = std::pow(0.25*(w[ijk-jj] + w[ijk-jj+kk] + w[ijk] + w[ijk+kk]), 2);
                const double utot = std::sqrt(u2 + v2 + w2);
                vt[ijk] -= Cd * pad[k] * utot * (v[ijk]+vtrans);
            }

    // Calculate drag for w-component.
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                const double u2 = std::pow(0.25*(u[ijk-kk] + u[ijk+ii-kk] + u[ijk] + u[ijk+ii]) + utrans, 2);
                const double v2 = std::pow(0.25*(v[ijk-kk] + v[ijk+jj-kk] + v[ijk] + v[ijk+jj]) + vtrans, 2);
                const double w2 = std::pow(w[ijk], 2);
                const double utot = std::sqrt(u2 + v2 + w2);
                wt[ijk] -= Cd * pad[k] * utot * w[ijk];
            }
}
