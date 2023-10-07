/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <elaspect/global.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/revision.h>
#include <deal.II/base/vectorization.h>

#include <cstring>
#include <fstream>

template <class Stream>
void print_elaspect_header(Stream &stream)
{
  const int n_tasks = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  stream << "-----------------------------------------------------------------------------\n"
         << "-- This is elASPECT, the elasticity-based ASPECT.\n";

  stream << "--     . using deal.II " << DEAL_II_PACKAGE_VERSION << "\n"
         << "--     .       with "
#ifdef DEAL_II_WITH_64BIT_INDICES
         << "64"
#else
         << "32"
#endif
         << " bit indices and vectorization level ";
  const unsigned int n_vect_bits =
    dealii::VectorizedArray<double>::size() * 8 * sizeof(double);

  stream << DEAL_II_COMPILER_VECTORIZATION_LEVEL
         << " (" << n_vect_bits << " bits)\n";

  stream << "--     . using Trilinos "
         << DEAL_II_TRILINOS_VERSION_MAJOR << '.'
         << DEAL_II_TRILINOS_VERSION_MINOR << '.'
         << DEAL_II_TRILINOS_VERSION_SUBMINOR << '\n';

  stream << "--     . using p4est "
         << DEAL_II_P4EST_VERSION_MAJOR << '.'
         << DEAL_II_P4EST_VERSION_MINOR << '.'
         << DEAL_II_P4EST_VERSION_SUBMINOR << '\n';

#ifdef DEBUG
  stream << "--     . running in DEBUG mode\n"
#else
  stream << "--     . running in OPTIMIZED mode\n"
#endif
         << "--     . running with " << n_tasks << " MPI process" << (n_tasks == 1 ? "\n" : "es\n");

  const int n_threads = dealii::MultithreadInfo::n_threads();
  if (n_threads > 1)
    stream << "--     . using " << n_threads << " threads" << (n_tasks == 1 ? "\n" : "each\n");

  stream << "-----------------------------------------------------------------------------\n"
         << std::endl;
}

template void print_elaspect_header<std::ostream> (std::ostream &stream);
template void print_elaspect_header<std::ofstream> (std::ofstream &stream);
