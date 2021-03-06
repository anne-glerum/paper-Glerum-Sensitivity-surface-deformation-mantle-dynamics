/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_postprocess_point_values_eigen_vectors_h
#define _aspect_postprocess_point_values_eigen_vectors_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>


namespace aspect
{
  namespace Postprocess
  {
    enum struct SymmetricTensorEigenvectorMethod
          {
            /**
             * A hybrid approach that preferentially uses the characteristic equation to
             * compute eigenvalues and an analytical approach based on the cross-product
             * for the eigenvectors. If the computations are deemed too inaccurate then
             * the method falls back to ql_implicit_shifts.
             *
             * This method potentially offers the quickest computation if the pathological
             * case is not encountered.
             */
            hybrid,
            /**
             * The iterative QL algorithm with implicit shifts applied after
             * tridiagonalization of the tensor using the householder method.
             *
             * This method offers a compromise between speed of computation and its
             * robustness. This method is particularly useful when the elements
             * of $T$ have greatly varying magnitudes, which would typically lead to a
             * loss of accuracy when computing the smaller eigenvalues.
             */
            ql_implicit_shifts,
            /**
             * The iterative Jacobi algorithm.
             *
             * This method offers is the most robust of the available options, with
             * reliable results obtained for even the most pathological cases. It is,
             * however, the slowest algorithm of all of those implemented.
             */
            jacobi
          };

          namespace internal
          {
            namespace SymmetricTensorImplementation
            {
              template <int dim>
              void
              tridiagonalize(const dealii::SymmetricTensor<2, dim, double> &A,
                             dealii::Tensor<2, dim, double> &               Q,
                             std::array<double, dim> &                      d,
                             std::array<double, dim - 1> &                  e);

              template <int dim>
              std::array<std::pair<double, Tensor<1, dim, double>>, dim>
              ql_implicit_shifts(const dealii::SymmetricTensor<2, dim, double> &A);

              std::array<std::pair<double, Tensor<1, 2, double>>, 2>
              hybrid(const dealii::SymmetricTensor<2, 2, double> &A);

              std::array<std::pair<double, Tensor<1, 3, double>>, 3>
              hybrid(const dealii::SymmetricTensor<2, 3, double> &A);

              std::array<double, 2>
              eigenvalues(const SymmetricTensor<2, 2, double> &T);

              std::array<double, 3>
              eigenvalues(const SymmetricTensor<2, 3, double> &T);

              template <int dim>
              std::array<std::pair<double, Tensor<1, dim, double>>, dim>
              perform_eigenvector_decomposition(const SymmetricTensor<2, dim, double> &T,
                  const SymmetricTensorEigenvectorMethod method);

              /**
               * A struct that is used to sort arrays of pairs of eign=envalues and
               * eigenvectors. Sorting is performed in descending order of eigenvalue.
               */
              template <int dim>
              struct SortEigenValuesVectors
              {
                using EigValsVecs = std::pair<double, Tensor<1, dim, double>>;
                bool
                operator()(const EigValsVecs &lhs, const EigValsVecs &rhs)
                {
                  return lhs.first > rhs.first;
                }
              };
            }
          }
    /**
     * A postprocessor that evaluates the solution vector at individual
     * points.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class PointValuesEigenVectors : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        PointValuesEigenVectors ();

        /**
         * Evaluate the solution and determine the values at the
         * selected points.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Save the state of this object.
         */
        virtual
        void save (std::map<std::string, std::string> &status_strings) const;

        /**
         * Restore the state of the object.
         */
        virtual
        void load (const std::map<std::string, std::string> &status_strings);

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:

        std::array<std::pair<double, Tensor<1, dim, double> >, dim>
        eigenvectors(const SymmetricTensor<2, dim, double> &,
                     const SymmetricTensorEigenvectorMethod) const;


        /**
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);

        /**
         * Interval between the generation of output in seconds.
         */
        double output_interval;

        /**
         * A time (in seconds) the last output has been produced.
         */
        double last_output_time;

        /**
         * Vector of Points representing the points where the solution is to be evaluated
         * that can be used by VectorTools.
         */
        std::vector<Point<dim> > evaluation_points_cartesian;
        /**
         * The values of the solution at the evaluation points.
         */
        std::vector<std::pair<double, std::vector<Vector<double> > > > point_values;
        /**
         * Whether or not to interpret the evaluation points in the input file
         * as natural coordinates or not.
         */
        bool use_natural_coordinates;

        /**
         * Whether or not to also compute the strain rate eigenvectors.
         * If false, only the solution variables are outputted.
         */
        bool output_eigen_vectors;
    };
  }
}


#endif
