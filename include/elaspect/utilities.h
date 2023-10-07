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


#ifndef _elaspect_utilities_h
#define _elaspect_utilities_h

#include <elaspect/global.h>

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

namespace elaspect
{
  namespace GeometryModel
  {
    template <int dim> class Interface;
  }

  namespace Utilities
  {
    using namespace dealii;
    using namespace dealii::Utilities;

    bool fexists(const std::string &filename);

    bool fexists(const std::string &filename,
                 MPI_Comm comm);

    std::string
    read_and_distribute_file_content(const std::string &filename,
                                     const MPI_Comm &comm);

    void
    collect_and_write_file_content(const std::string &filename,
                                   const std::string &file_content,
                                   const MPI_Comm &comm);

    int
    mkdirp(std::string pathname, const mode_t mode = 0755);

    /**
     * Create directory @p pathname, optionally printing a message.
     *
     * @param pathname String that contains path to create. '/' is used as
     * directory separator.
     * @param comm MPI communicator, used to limit creation of directory to
     * processor 0.
     * @param silent Print a nicely formatted message on processor 0 if set
     * to true.
     */
    void create_directory(const std::string &pathname,
                          const MPI_Comm &comm,
                          bool silent);

    std::string
    expand_ELASPECT_SOURCE_DIR (const std::string &location);

    template <typename T>
    std::vector<T>
    possibly_extend_from_1_to_N (const std::vector<T> &values,
                                 const unsigned int N,
                                 const std::string &id_text);

    bool has_unique_entries (const std::vector<std::string> &strings);

    struct ThousandSep : std::numpunct<char>
    {
      protected:
        char do_thousands_sep() const override
        {
          return ',';
        }

        std::string do_grouping() const override
        {
          return "\003";  // groups of 3 digits (this string is in octal format)
        }

    };

    namespace Coordinates
    {
      enum CoordinateSystem
      {
        depth,
        cartesian,
        spherical,
        ellipsoidal,
        invalid
      };

      /**
       * Returns spherical coordinates of a Cartesian point. The returned array
       * is filled with radius, phi and theta (polar angle). If the dimension is
       * set to 2 theta is omitted. Phi is always normalized to [0,2*pi].
       *
       */
      template <int dim>
      std::array<double,dim>
      cartesian_to_spherical_coordinates(const dealii::Point<dim> &position);

      /**
       * Return the Cartesian point of a spherical position defined by radius,
       * phi and theta (polar angle). If the dimension is set to 2 theta is
       * omitted.
       */
      template <int dim>
      dealii::Point<dim>
      spherical_to_cartesian_coordinates(const std::array<double,dim> &scoord);

      /**
       * Given a vector defined in the radius, phi and theta directions, return
       * a vector defined in Cartesian coordinates. If the dimension is set to 2
       * theta is omitted. Position is given as a Point in Cartesian coordinates.
       */
      template <int dim>
      Tensor<1,dim>
      spherical_to_cartesian_vector(const Tensor<1,dim> &spherical_vector,
                                    const dealii::Point<dim> &position);

      /**
       * A function that takes a string representation of the name of a
       * coordinate system (as represented by the CoordinateSystem enum)
       * and returns the corresponding value.
       */
      CoordinateSystem
      string_to_coordinate_system (const std::string &);
    }

    /**
     * Converts an array of size dim to a Point of size dim.
     */
    template <int dim>
    Point<dim> convert_array_to_point(const std::array<double,dim> &array);

    /**
     * Converts a Point of size dim to an array of size dim.
     */
    template <int dim>
    std::array<double,dim> convert_point_to_array(const Point<dim> &point);

    template <int dim>
    class NaturalCoordinate
    {
      public:
        /**
         * Constructor based on providing the geometry model as a pointer.
         */
        NaturalCoordinate(Point<dim> &position,
                          const GeometryModel::Interface<dim> &geometry_model);

        /**
         * Constructor based on providing the coordinates and associated
         * coordinate system.
         */
        NaturalCoordinate(const std::array<double, dim> &coord,
                          const Utilities::Coordinates::CoordinateSystem &coord_system);

        /**
         * Returns the coordinates in the given coordinate system, which may
         * not be Cartesian.
         */
        const std::array<double,dim> &get_coordinates() const;

        /**
         * The coordinate that represents the 'surface' directions in the
         * chosen coordinate system.
         */
        std::array<double,dim-1> get_surface_coordinates() const;

        /**
         * The coordinate that represents the 'depth' direction in the chosen
         * coordinate system.
         */
        double get_depth_coordinate() const;

      private:
        /**
         * An enum which stores the the coordinate system of this natural
         * point
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * An array which stores the coordinates in the coordinates system
         */
        std::array<double,dim> coordinates;
    };

    /**
     * A class that represents a binary operator between two doubles. The type of
     * operation is specified on construction time, and can be checked later
     * by using the operator ==. The operator () executes the operation on two
     * double parameters and returns the result. This class is helpful for
     * user specified operations that are not known at compile time.
     */
    class Operator
    {
      public:
        /**
         * An enum of supported operations.
         */
        enum operation
        {
          uninitialized,
          add,
          subtract,
          minimum,
          maximum,
          replace_if_valid
        };

        /**
         * The default constructor creates an invalid operation that will fail
         * if ever executed.
         */
        Operator();

        /**
         * Construct the selected operator.
         */
        Operator(const operation op);

        /**
         * Execute the selected operation with the given parameters and
         * return the result.
         */
        double operator() (const double x, const double y) const;

        /**
         * Return the comparison result between the current operation and
         * the one provided as argument.
         */
        bool operator== (const operation op) const;

      private:
        /**
         * The selected operation of this object.
         */
        operation op;
    };


    /**
     * Create a vector of operator objects out of a list of strings. Each
     * entry in the list must match one of the allowed operations.
     */
    std::vector<Operator> create_model_operator_list(const std::vector<std::string> &operator_names);

    /**
     * Create a string of model operators for use in declare_parameters
     */
    const std::string get_model_operator_options(); 


    /**
    * Conversion object where one can provide a function that returns
    * a tensor for the velocity at a given point and it returns something
    * that matches the dealii::Function interface with a number of output
    * components equal to the number of components of the finite element
    * in use.
    */
    template <int dim>
    class VectorFunctionFromVelocityFunctionObject : public Function<dim>
    {
      public:
        /**
         * Given a function object that takes a Point and returns a Tensor<1,dim>,
         * convert this into an object that matches the Function@<dim@>
         * interface.
         *
         * @param n_components total number of components of the finite element system.
         * @param function_object The function that will form one component
         *     of the resulting Function object.
         */
        VectorFunctionFromVelocityFunctionObject (const unsigned int n_components,
                                                  const std::function<Tensor<1,dim> (const Point<dim> &)> &function_object);

        /**
         * Return the value of the
         * function at the given
         * point. Returns the value the
         * function given to the constructor
         * produces for this point.
         */
        double value (const Point<dim>   &p,
                      const unsigned int  component = 0) const override;

        /**
         * Return all components of a
         * vector-valued function at a
         * given point.
         *
         * <tt>values</tt> shall have the right
         * size beforehand,
         * i.e. #n_components.
         */
        void vector_value (const Point<dim>   &p,
                           Vector<double>     &values) const override;

      private:
        /**
         * The function object which we call when this class's value() or
         * value_list() functions are called.
         */
        const std::function<Tensor<1,dim> (const Point<dim> &)> function_object;
    };



    /*------------------------ inline functions ------------------------*/

    template <typename T>
    inline
    std::vector<T>
    possibly_extend_from_1_to_N (const std::vector<T> &values,
                                 const unsigned int N,
                                 const std::string &id_text)
    {
      if (values.size() == 1)
        {
          return std::vector<T> (N, values[0]);
        }
      else if (values.size() == N)
        {
          return values;
        }
      else
        {
          // Non-specified behavior
          AssertThrow(false,
                      ExcMessage("Length of " + id_text + " list must be " +
                                 "either one or " + Utilities::to_string(N) +
                                 ". Currently it is " + Utilities::to_string(values.size()) + "."));
        }

      // This should never happen, but return an empty vector so the compiler
      // will be happy
      return std::vector<T> ();
    }
  }
}

#endif
