# 
#  MicroHH
#  Copyright (c) 2011-2014 Chiel van Heerwaarden
#  Copyright (c) 2011-2014 Thijs Heus
#  Copyright (c)      2014 Bart van Stratum
# 
#  This file is part of MicroHH
# 
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#
file(GLOB sourcefiles "*.cxx" "*.cu")
include_directories("../include" SYSTEM ${INCLUDE_DIRS})

# Strip off the path from the file names
foreach(file_ ${sourcefiles})
  get_filename_component(filename ${file_} NAME)
  list(APPEND parsedfiles ${filename})
endforeach()

if(USECUDA)
  cuda_add_library(microhhc STATIC ${parsedfiles})
else()
  add_library(microhhc STATIC ${parsedfiles})
endif()

# Invoke the StencilBuilder parser Python script on each file.
foreach(file_ ${sourcefiles})
  get_filename_component(filename ${file_} NAME)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    COMMAND python ${CMAKE_SOURCE_DIR}/stencilbuilder/parser.py ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
  )
endforeach()
