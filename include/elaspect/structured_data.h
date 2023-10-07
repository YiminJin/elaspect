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


#ifndef _elaspect_structured_data_h
#define _elaspect_structured_data_h

#include <elaspect/simulator_access.h>

namespace elaspect
{
  namespace Utilities
  {
    using namespace dealii;

    template <int dim>
    class StructuredDataLookup
    {
      public:
        StructuredDataLookup(const unsigned int components,
                             const double scale_factor);

        explicit StructuredDataLookup(const double scale_factor);

        void reinit(const std::vector<std::string> &column_names,
                    std::vector<std::vector<double>> &&coordinate_values,
                    std::vector<Table<dim,double>> &&data_table,
                    const MPI_Comm &mpi_communicator = MPI_COMM_SELF,
                    const unsigned int root_process = numbers::invalid_unsigned_int);

        void
        load_file(const std::string &filename,
                  const MPI_Comm &communicator);

        void 
        load_netcdf(const std::string &filename, 
                    const std::vector<std::string> &data_column_names = {});

        double
        get_data(const Point<dim> &position,
                 const unsigned int component) const;

        Tensor<1,dim>
        get_gradients(const Point<dim> &position,
                      const unsigned int component) const;

        double get_maximum_component_value(const unsigned int component) const;

      private:
        unsigned int components;

        std::vector<std::string> data_component_names;

        std::vector<std::unique_ptr<Function<dim>>> data;

        std::array<std::vector<double>, dim> coordinate_values;

        std::vector<double> maximum_component_value;

        TableIndices<dim> table_points;

        const double scale_factor;

        bool coordinate_values_are_equidistant;

        TableIndices<dim>
        compute_table_indices(const TableIndices<dim> &sizes,
                              const std::size_t idx) const;
    };

    template <int dim>
    class AsciiDataBase
    {
      public:
        AsciiDataBase();

        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename,
                            const std::string &subsection_name = "Ascii data model");

        void
        parse_parameters (ParameterHandler &prm,
                          const std::string &subsection_name = "Ascii data model");

        std::string data_directory;

        std::string data_file_name;

        double scale_factor;
    };

    template <int dim>
    class AsciiDataBoundary : public Utilities::AsciiDataBase<dim>, public SimulatorAccess<dim>
    {
      public:
        AsciiDataBoundary();

        virtual
        void
        initialize(const std::set<types::boundary_id> &boundary_ids,
                   const unsigned int components);

        void
        update();

        double
        get_data_component(const types::boundary_id boundary_indicator,
                           const Point<dim> &       position,
                           const unsigned int       component) const;

        double
        get_maximum_component_value (const types::boundary_id boundary_indicator,
                                     const unsigned int       component) const;

        Tensor<1,dim-1>
        vector_gradient(const types::boundary_id boundary_indicator,
                        const Point<dim>        &p,
                        const unsigned int       component) const;

        static
        void
        declare_parameters (ParameterHandler  &prm,
                            const std::string &default_directory,
                            const std::string &default_filename,
                            const std::string &subsection_name = "Ascii data model");

        void
        parse_parameters (ParameterHandler &prm,
                          const std::string &subsection_name = "Ascii data model");

      protected:
        int current_file_number;

        int first_data_file_number;

        double time_weight;

        bool time_dependent;

        std::map<types::boundary_id,
                 std::unique_ptr<Utilities::StructuredDataLookup<dim-1>>> lookups;

        std::map<types::boundary_id,
                 std::unique_ptr<Utilities::StructuredDataLookup<dim-1>>> old_lookups;

        std::string
        create_filename (const int filenumber,
                         const types::boundary_id boundary_id) const;
    };

    template <int dim>
    class AsciiDataInitial : public Utilities::AsciiDataBase<dim>, public SimulatorAccess<dim>
    {
      public:
    };
  }
}

#endif
