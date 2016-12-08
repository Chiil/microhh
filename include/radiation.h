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

class Master;
class Model;
class Grid;
class Fields;

/**
 * Class for the radiation layer in the top of the domain.
 * This class performs the gravity wave damping in the top of the domain to
 * prevent reflection at the top boundary.
 */
class Radiation
{
    public:
        Radiation(Model*, Input*); ///< Constructor of the radiation class.
        ~Radiation();              ///< Destructor of the radiation class.

        void init();         ///< Initialize the arrays that contain the profiles.
        void create(Input*); ///< Read the profiles of the forces from the input.
        void exec();         ///< Add the tendencies created by the damping.

    private:
        Master* master; ///< Pointer to master class.
        Model*  model;  ///< Pointer to model class.
        Grid*   grid;   ///< Pointer to grid class.
        Fields* fields; ///< Pointer to fields class.
};
#endif
