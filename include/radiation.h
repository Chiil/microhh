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

#ifndef RADIATION
#define RADIATION

#include <string>
#include "defines.h"

class Master;
class Input;
class Grid;
class Fields;
class Thermo;
class Stats;
struct Mask;

/**
 * Class for radiation.
 */
class Radiation
{
    public:
        Radiation(Master*, Input*, Grid*,
                  Fields*); ///< Constructor of the force class.
        ~Radiation();          ///< Destructor of the force class.

        void init(Thermo*, Stats*); ///< Initialize the arrays that contain the profiles.

        void create(Input*);    ///< Read the profiles of the forces from the input.
        void exec();            ///< Add the tendencies belonging to the large-scale processes.
        void exec_stats(Mask*); ///< Add the tendencies belonging to the large-scale processes.

    private:
        void calc_radiation_tendency(double* restrict, double* restrict); ///< Add the radiation tendency.
        void calc_radiation_stats(double* restrict, double* restrict); ///< Calculate radiation stats.

        Master* master; ///< Pointer to master class.
        Grid* grid;     ///< Pointer to grid class.
        Fields* fields; ///< Pointer to fields class.
        Stats* stats;   ///< Pointer to stats class.

        double* rad_tend;  ///< Pointer to array radiation tendencies.

        std::string swradiation; ///< Switch for radiation scheme.
};
#endif
