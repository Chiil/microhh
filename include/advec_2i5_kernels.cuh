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

#ifndef ADVEC_2I5_KERNELS_CUH
#define ADVEC_2I5_KERNELS_CUH

#include "finite_difference.h"
#include "cuda_tiling.h"

namespace Advec_2i5_kernels
{
    using namespace Finite_difference::O2;
    using namespace Finite_difference::O4;
    using namespace Finite_difference::O6;

    template<typename TF>
    struct advec_u_g
    {
        DEFINE_GRID_KERNEL("advec_2i5::advec_u", 3)

        template <typename Level>
        CUDA_DEVICE
        void operator()(Grid_layout g, const int i, const int j, const int k, const Level level,
                        TF* __restrict__ ut, const TF* __restrict__ u,
                        const TF* __restrict__ v,  const TF* __restrict__ w,
                        const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefh,
                        const TF* __restrict__ dzi, const TF dxi, const TF dyi)
        {
            const int ii1 = 1*g.istride;
            const int ii2 = 2*g.istride;
            const int ii3 = 3*g.istride;
            const int jj1 = 1*g.jstride;
            const int jj2 = 2*g.jstride;
            const int jj3 = 3*g.jstride;
            const int kk1 = 1*g.kstride;
            const int kk2 = 2*g.kstride;
            const int kk3 = 3*g.kstride;

            const int ijk = g(i, j, k);

            ut[ijk] +=
                    // u*du/dx
                    - ( interp2(u[ijk        ], u[ijk+ii1]) * interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3])
                        - interp2(u[ijk-ii1    ], u[ijk    ]) * interp6_ws(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) ) * dxi

                    + ( fabs(interp2(u[ijk        ], u[ijk+ii1])) * interp5_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3])
                        - fabs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp5_ws(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) ) * dxi

                    // v*du/dy
                    - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp6_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3])
                        - interp2(v[ijk-ii1    ], v[ijk    ]) * interp6_ws(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]) ) * dyi

                    + ( fabs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp5_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3])
                        - fabs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp5_ws(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]) ) * dyi;

            if (level.distance_to_start() == 0)
            {
                // w*du/dz -> second order interpolation for fluxtop, fluxbot = 0. as w=0
                ut[ijk] +=
                        - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 1)
            {
                ut[ijk] +=
                        // w*du/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                        - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                            - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(   u[ijk-kk1], u[ijk    ]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp3_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 2)
            {
                ut[ijk] +=
                        // w*du/dz -> fourth order interpolation for fluxbot
                        - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp6_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                            - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp5_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                            - rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp3_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 2)
            {
                ut[ijk] +=
                        // w*du/dz -> fourth order interpolation for fluxtop
                        - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4_ws(u[ijk-kk1   ], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                            - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp6_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp3_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                            - rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp5_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 1)
            {
                ut[ijk] +=
                        // w*du/dz -> second order interpolation for fluxtop, fourth order for fluxbot
                        - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                            - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) * rhorefi[k] * dzi[k]

                        - ( rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp3_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 0)
            {
                ut[ijk] +=
                        // w*du/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                        - ( -rhorefh[k] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) * rhorefi[k] * dzi[k];
            }
            else
            {
                ut[ijk] +=
                        // w*du/dz
                        - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp6_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                            - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp6_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp5_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                            - rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp5_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
        }
    };

    template<typename TF>
    struct advec_v_g
    {
        DEFINE_GRID_KERNEL("advec_2i5::advec_v", 3)

        template <typename Level>
        CUDA_DEVICE
        void operator()(Grid_layout g, const int i, const int j, const int k, const Level level,
                   TF* __restrict__ vt, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefh,
                   const TF* __restrict__ dzi, const TF dxi, const TF dyi)
        {
            const int ii = g.istride;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ii1 = 1*ii;
            const int ii2 = 2*ii;
            const int ii3 = 3*ii;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            const int ijk = i*ii + j*jj + k*kk;

            vt[ijk] +=
                    // u*dv/dx
                    - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp6_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3])
                        - interp2(u[ijk    -jj1], u[ijk    ]) * interp6_ws(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]) ) * dxi

                    + ( fabs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp5_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3])
                        - fabs(interp2(u[ijk    -jj1], u[ijk    ])) * interp5_ws(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]) ) * dxi

                    // v*dv/dy
                    - ( interp2(v[ijk        ], v[ijk+jj1]) * interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3])
                        - interp2(v[ijk-jj1    ], v[ijk    ]) * interp6_ws(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) ) * dyi

                    + ( fabs(interp2(v[ijk        ], v[ijk+jj1])) * interp5_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3])
                        - fabs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp5_ws(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) ) * dyi;

            if (level.distance_to_start() == 0)
            {
                vt[ijk] +=
                        // w*dv/dz -> second order interpolation for fluxtop, fluxbot=0 as w=0
                        - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 1)
            {
                vt[ijk] +=
                        // w*dv/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                        - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                            - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp3_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 2)
            {
                vt[ijk] +=
                        // w*dv/dz -> fourth order interpolation for fluxbot, sixth for fluxtop
                        - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp6_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                            - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp5_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                            - rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp3_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 2)
            {
                vt[ijk] +=
                        // w*dv/dz -> fourth order interpolation for fluxtop
                        - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                            - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp6_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp3_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                            - rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp5_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 1)
            {
                vt[ijk] +=
                        // w*dv/dz -> second order interpolation for fluxtop, fourth order for fluxbot
                        - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                            - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) * rhorefi[k] * dzi[k]

                        - ( rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp3_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 0)
            {
                vt[ijk] +=
                        // w*dv/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                        - ( -rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) * rhorefi[k] * dzi[k];
            }
            else
            {
                vt[ijk] +=
                        - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp6_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                            - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp6_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp5_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                            - rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp5_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
        }
    };



    template<typename TF>
    struct advec_w_g
    {
        DEFINE_GRID_KERNEL("advec_2i5::advec_w", 3)

        template <typename Level>
        CUDA_DEVICE
        void operator()(Grid_layout g, const int i, const int j, const int k, const Level level,
                       TF* __restrict__ wt, const TF* __restrict__ u,
                       const TF* __restrict__ v,  const TF* __restrict__ w,
                       const TF* __restrict__ rhoref, const TF* __restrict__ rhorefhi,
                       const TF* __restrict__ dzhi, const TF dxi, const TF dyi)
        {
            const int ii = g.istride;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ii1 = 1*ii;
            const int ii2 = 2*ii;
            const int ii3 = 3*ii;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            const int ijk = i*ii + j*jj + k*kk;

            wt[ijk] +=
                    // u*dw/dx
                    - ( interp2(u[ijk+ii1-kk], u[ijk+ii1]) * interp6_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3])
                        - interp2(u[ijk    -kk], u[ijk    ]) * interp6_ws(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]) ) * dxi

                    + ( fabs(interp2(u[ijk+ii1-kk], u[ijk+ii1])) * interp5_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3])
                        - fabs(interp2(u[ijk    -kk], u[ijk    ])) * interp5_ws(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]) ) * dxi

                    // v*dw/dy
                    - ( interp2(v[ijk+jj1-kk], v[ijk+jj1]) * interp6_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3])
                        - interp2(v[ijk    -kk], v[ijk    ]) * interp6_ws(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]) ) * dyi

                    + ( fabs(interp2(v[ijk+jj1-kk], v[ijk+jj1])) * interp5_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3])
                        - fabs(interp2(v[ijk    -kk], v[ijk    ])) * interp5_ws(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]) ) * dyi;

            if (level.distance_to_start() == 1)
            {
                wt[ijk] +=
                        // w*dv/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                        - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                            - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) * rhorefhi[k] * dzhi[k]

                        + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp3_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) * rhorefhi[k] * dzhi[k];
            }
            else if (level.distance_to_start() == 2)
            {
                wt[ijk] +=
                        // w*dv/dz -> fourth order interpolation for fluxbot, sixth order for fluxtop
                        - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp6_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                            - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) * rhorefhi[k] * dzhi[k]

                        + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp5_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                            - rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp3_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) * rhorefhi[k] * dzhi[k];
            }
            else if (level.distance_to_end() == 1)
            {
                wt[ijk] +=
                        // w*dv/dz -> sixth order interpolation for fluxbot, fourth order for fluxtop
                        - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                            - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp6_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) * rhorefhi[k] * dzhi[k]

                        + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp3_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                            - rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp5_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) * rhorefhi[k] * dzhi[k];
            }
            else if (level.distance_to_end() == 0)
            {
                wt[ijk] +=
                        // w*dv/dz -> fourth order interpolation for fluxbot, second order for fluxtop
                        - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                            - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) * rhorefhi[k] * dzhi[k]

                        - ( rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp3_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) * rhorefhi[k] * dzhi[k];
            }
            else if (level.distance_to_start() >= 2 && level.distance_to_end() >= 2)
            {
                wt[ijk] +=
                        // w*dw/dz
                        - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp6_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                            - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp6_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) * rhorefhi[k] * dzhi[k]

                        + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp5_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                            - rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp5_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) * rhorefhi[k] * dzhi[k];
            }
        }
    };


    template<typename TF>
    struct advec_s_g
    {
        DEFINE_GRID_KERNEL("advec_2i5::advec_s", 3)

        template <typename Level>
        CUDA_DEVICE
        void operator()(Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* __restrict__ st, const TF* __restrict__ s,
                const TF* __restrict__ u, const TF* __restrict__ v,  const TF* __restrict__ w,
                const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefh,
                const TF* __restrict__ dzi, const TF dxi, const TF dyi)
        {
            const int ii = g.istride;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            const int ijk = i*ii + j*jj + k*kk;

            st[ijk] +=
                    - ( u[ijk+ii1] * interp6_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3])
                        - u[ijk    ] * interp6_ws(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]) ) * dxi

                    + ( fabs(u[ijk+ii1]) * interp5_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3])
                        - fabs(u[ijk    ]) * interp5_ws(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]) ) * dxi

                    - ( v[ijk+jj1] * interp6_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3])
                        - v[ijk    ] * interp6_ws(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]) ) * dyi

                    + ( fabs(v[ijk+jj1]) * interp5_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3])
                        - fabs(v[ijk    ]) * interp5_ws(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]) ) * dyi;

            if (level.distance_to_start() == 0)
            {
                st[ijk] +=
                        // w*ds/dz -> second order interpolation for fluxtop, fluxbot=0 as w=0
                        - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 1)
            {
                st[ijk] +=
                        // w*ds/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                        - ( rhorefh[k+1] * w[ijk+kk1] * interp4_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                            - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(w[ijk+kk1]) * interp3_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 2)
            {
                st[ijk] +=
                        // w*ds/dz -> fourth order interpolation for fluxbot, sixth for fluxtop
                        - ( rhorefh[k+1] * w[ijk+kk1] * interp6_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3])
                            - rhorefh[k  ] * w[ijk    ] * interp4_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(w[ijk+kk1]) * interp5_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3])
                            - rhorefh[k  ] * fabs(w[ijk    ]) * interp3_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 2)
            {
                st[ijk] +=
                        // w*ds/dz -> fourth order interpolation for fluxtop, sixth order for fluxbot
                        - ( rhorefh[k+1] * w[ijk+kk1] * interp4_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                            - rhorefh[k  ] * w[ijk    ] * interp6_ws(s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(w[ijk+kk1]) * interp3_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                            - rhorefh[k  ] * fabs(w[ijk    ]) * interp5_ws(s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 1)
            {
                st[ijk] +=
                        // w*ds/dz -> second order interpolation for fluxtop, fourth order for fluxbot
                        - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                            - rhorefh[k  ] * w[ijk    ] * interp4_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k]

                        + ( -rhorefh[k  ] * fabs(w[ijk    ]) * interp3_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 0)
            {
                st[ijk] +=
                        // w*ds/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                        - (- rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) * rhorefi[k] * dzi[k];
            }
            else
            {
                st[ijk] +=
                        - ( rhorefh[k+1] * w[ijk+kk1] * interp6_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3])
                            - rhorefh[k  ] * w[ijk    ] * interp6_ws(s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]) ) * rhorefi[k] * dzi[k]

                        + ( rhorefh[k+1] * fabs(w[ijk+kk1]) * interp5_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3])
                            - rhorefh[k  ] * fabs(w[ijk    ]) * interp5_ws(s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]) ) * rhorefi[k] * dzi[k];
            }
        }
    };


    // Implementation flux limiter according to Koren, 1993.
    template<typename TF> __device__
    inline TF flux_lim_g(const TF u, TF sm2, TF sm1, TF sp1, TF sp2)
    {
        const TF eps = TF(1.e-12);

        if (u < TF(0.))
        {
            // swap
            TF tmp1 = sm1;
            sm1 = sp1;
            sp1 = tmp1;

            // swap
            TF tmp = sm2;
            sm2 = sp2;
            sp2 = tmp;
        }

        const TF two_r = TF(2.) * (sp1-sm1+eps) / (sm1-sm2+eps);
        const TF phi = max(
                TF(0.),
                min( two_r, min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
        return u*(sm1 + TF(0.5)*phi*(sm1 - sm2));
    }

    template<typename TF> __device__
    inline TF flux_lim_bot_g(const TF u, const TF sm2, const TF sm1, const TF sp1, const TF sp2)
    {
        const TF eps = TF(1.e-12);

        if (u >= TF(0.))
        {
            return u*sm1;
        }
        else
        {
            const TF two_r = TF(2.) * (sm1-sp1+eps) / (sp1-sp2+eps);
            const TF phi = max(
                    TF(0.),
                    min( two_r, min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
            return u*(sp1 + TF(0.5)*phi*(sp1 - sp2));
        }
    }

    template<typename TF> __device__
    inline TF flux_lim_top_g(const TF u, const TF sm2, const TF sm1, const TF sp1, const TF sp2)
    {
        const TF eps = TF(1.e-12);

        if (u >= TF(0.))
        {
            const TF two_r = TF(2.) * (sp1-sm1+eps) / (sm1-sm2+eps);
            const TF phi = max(
                    TF(0.),
                    min( two_r, min( TF(1./3.)*(TF(1.)+two_r), TF(2.)) ) );
            return u*(sm1 + TF(0.5)*phi*(sm1 - sm2));
        }
        else
        {
            return u*sp1;
        }
    }

    template<typename TF>
    struct advec_s_lim_g
    {
        DEFINE_GRID_KERNEL("advec_2i5::advec_s_lim", 3)

        template <typename Level>
        CUDA_DEVICE
        void operator()(Grid_layout g, const int i, const int j, const int k, const Level level,
                        TF* __restrict__ st, const TF* __restrict__ s,
                        const TF* __restrict__ u, const TF* __restrict__ v,  const TF* __restrict__ w,
                        const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefh,
                        const TF* __restrict__ dzi, const TF dxi, const TF dyi)
        {
            const int ii1 = 1*g.istride;
            const int ii2 = 2*g.istride;
            const int jj1 = 1*g.jstride;
            const int jj2 = 2*g.jstride;
            const int kk = g.kstride;
            const int kk1 = 1*g.kstride;
            const int kk2 = 2*g.kstride;

            const int ijk = i + j*jj1 + k*kk1;
            st[ijk] +=
                    - ( flux_lim_g(u[ijk+ii1], s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                        - flux_lim_g(u[ijk    ], s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                    - ( flux_lim_g(v[ijk+jj1], s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                        - flux_lim_g(v[ijk    ], s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi;

            if (level.distance_to_start() >= 2 && level.distance_to_start() >= 2)
            {
                st[ijk] +=
                        - ( rhorefh[k+1] * flux_lim_g(w[ijk+kk], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                            - rhorefh[k  ] * flux_lim_g(w[ijk   ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 0)
            {
                st[ijk] +=
                        - ( rhorefh[k+1] * flux_lim_bot_g(w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_start() == 1)
            {
                st[ijk] +=
                        - ( rhorefh[k+1] * flux_lim_g    (w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                            - rhorefh[k  ] * flux_lim_bot_g(w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 1)
            {
                st[ijk] +=
                        - ( rhorefh[k+1] * flux_lim_top_g(w[ijk+kk1], s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                            - rhorefh[k  ] * flux_lim_g    (w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
            else if (level.distance_to_end() == 0)
            {
                st[ijk] +=
                        - (
                                - rhorefh[k  ] * flux_lim_top_g(w[ijk    ], s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * rhorefi[k] * dzi[k];
            }
        }
    };

    template<typename TF>
    struct calc_cfl_g
    {
        DEFINE_GRID_KERNEL("advec_2i5::calc_cfl_g", 2)

        template <typename Level>
        CUDA_DEVICE
        void operator()(Grid_layout g, const int i, const int j, const int k, const Level level,
                    TF* const __restrict__ tmp1,
                    const TF* __restrict__ u, const TF* __restrict__ v, const TF* __restrict__ w,
                    const TF* __restrict__ dzi, const TF dxi, const TF dyi)
        {
            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*g.jstride;
            const int jj2 = 2*g.jstride;
            const int jj3 = 3*g.jstride;
            const int kk1 = 1*g.kstride;
            const int kk2 = 2*g.kstride;
            const int kk3 = 3*g.kstride;

            const int ijk = i + j*g.jstride+ k*g.kstride;

            if (level.distance_to_start() == 0 || level.distance_to_end() == 0)
                tmp1[ijk] = fabs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                            + fabs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                            + fabs(interp2(w[ijk], w[ijk+kk1]))*dzi[k];
            else if (level.distance_to_start() == 1 || level.distance_to_end() == 1)
                tmp1[ijk] = fabs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                            + fabs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                            + fabs(interp4_ws(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k];
            else
                tmp1[ijk] = fabs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                            + fabs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                            + fabs(interp6_ws(w[ijk-kk2], w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]))*dzi[k];
        }
    };
}


#endif //MICROHHC_ADVEC_2I5_KERNELS_CUH
