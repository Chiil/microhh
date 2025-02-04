/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#include "grid.h"
#include "tools.h"
#include "math.h"

template<typename TF>
void Grid<TF>::prepare_device()
{
    // Calculate optimal size thread blocks based on grid
    gd.ithread_block = min(256, 16 * ((gd.itot / 16) + (gd.itot % 16 > 0)));
    gd.jthread_block = 256 / gd.ithread_block;

    const int imemsize = gd.icells*sizeof(TF);
    const int jmemsize = gd.jcells*sizeof(TF);
    const int kmemsize = gd.kcells*sizeof(TF);

    gd.x_g.resize(gd.icells);
    gd.y_g.resize(gd.jcells);
    gd.z_g.resize(gd.kcells);
    gd.zh_g.resize(gd.kcells);
    gd.dz_g.resize(gd.kcells);
    gd.dzh_g.resize(gd.kcells);
    gd.dzi_g.resize(gd.kcells);
    gd.dzhi_g.resize(gd.kcells);
    gd.dzi4_g.resize(gd.kcells);
    gd.dzhi4_g.resize(gd.kcells);

    cuda_safe_call(cudaMemcpy(gd.x_g,     gd.x.data(),     imemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.y_g,     gd.y.data(),     jmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.z_g,     gd.z.data(),     kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.zh_g,    gd.zh.data(),    kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dz_g,    gd.dz.data(),    kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzh_g,   gd.dzh.data(),   kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzi_g,   gd.dzi.data(),   kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzhi_g,  gd.dzhi.data(),  kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzi4_g,  gd.dzi4.data(),  kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzhi4_g, gd.dzhi4.data(), kmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Grid<TF>::clear_device()
{
    gd.x_g.free();
    gd.y_g.free();
    gd.z_g.free();
    gd.zh_g.free();
    gd.dz_g.free();
    gd.dzh_g.free();
    gd.dzi_g.free();
    gd.dzhi_g.free();
    gd.dzi4_g.free();
    gd.dzhi4_g.free();
}


#ifdef FLOAT_SINGLE
template class Grid<float>;
#else
template class Grid<double>;
#endif
