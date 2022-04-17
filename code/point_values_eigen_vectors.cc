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


#include "point_values_eigen_vectors.h"
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/global.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_tools.h>

#include <aspect/gravity_model/interface.h>

#include <math.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
          {
            namespace SymmetricTensorImplementation
            {
              template <int dim>
              void
              tridiagonalize(const dealii::SymmetricTensor<2, dim, double> &A,
                  dealii::Tensor<2, dim, double> &               Q,
                  std::array<double, dim> &                      d,
                  std::array<double, dim - 1> &                  e)
              {
                // Create some intermediate storage
                double h, g, omega_inv, K, f;

                // Initialize the transformation matrix as the
                // identity tensor
                Q = dealii::unit_symmetric_tensor<dim, double>();

                // Make the first row and column to be of the
                // desired form
                h = 0.0;
                for (int i = 1; i < dim; i++)
                  h += A[0][i] * A[0][i];

                g = 0.0;
                if (A[0][1] > 0.0)
                  g = -std::sqrt(h);
                else
                  g = std::sqrt(h);
                e[0] = g;

                std::array<double, dim> u;
                for (int i = 1; i < dim; i++)
                {
                  u[i] = A[0][i];
                  if (i == 1)
                    u[i] -= g;
                }

                std::array<double, dim> q;
                const double omega = h - g * A[0][1];
                if (omega > 0.0)
                {
                  omega_inv = 1.0 / omega;
                  K         = 0.0;
                  for (int i = 1; i < dim; i++)
                  {
                    f = 0.0;
                    for (int j = 1; j < dim; j++)
                      f += A[i][j] * u[j];
                    q[i] = omega_inv * f;
                    K += u[i] * f;
                  }
                  K *= 0.5 * omega_inv * omega_inv;

                  for (int i = 1; i < dim; i++)
                    q[i] = q[i] - K * u[i];

                  d[0] = A[0][0];
                  for (int i = 1; i < dim; i++)
                    d[i] = A[i][i] - 2.0 * q[i] * u[i];

                  // Store inverse Householder transformation
                  // in Q
                  for (int j = 1; j < dim; j++)
                  {
                    f = omega_inv * u[j];
                    for (int i = 1; i < dim; i++)
                      Q[i][j] = Q[i][j] - f * u[i];
                  }

                  // For dim = 3: Calculate updated A[1][2] and
                  // store it in e[1]
                  for (int i = 1; i < dim - 1; i++)
                    e[i] = A[i][i + 1] - q[i] * u[i + 1] - u[i] * q[i + 1];
                }
                else
                {
                  for (int i = 0; i < dim; i++)
                    d[i] = A[i][i];

                  // For dim = 3:
                  for (int i = 1; i < dim - 1; i++)
                    e[i] = A[i][i + 1];
                }
              }

              template <int dim>
              std::array<std::pair<double, Tensor<1, dim, double>>, dim>
              ql_implicit_shifts(const dealii::SymmetricTensor<2, dim, double> &A)
              {
                // Transform A to real tridiagonal form by the Householder method:
                // The orthogonal matrix effecting the transformation
                // this will ultimately store the eigenvectors
                dealii::Tensor<2, dim, double> Q;
                // The diagonal elements of the tridiagonal matrix;
                // this will ultimately store the eigenvalues
                std::array<double, dim> w;
                // The off-diagonal elements of the tridiagonal
                std::array<double, dim - 1> ee;
                tridiagonalize<dim>(A, Q, w, ee);

                // Number of iterations
                const unsigned int max_n_it = 30;

                // Transfer the off-diagonal entries to an auxiliary array
                // The third element is used only as temporary workspace
                std::array<double, dim> e;
                for (unsigned int i = 0; i < dim - 1; ++i)
                  e[i] = ee[i];

                // Create some intermediate storage
                double g, r, p, f, b, s, c, t;

                // Loop over all off-diagonal elements
                for (int l = 0; l < dim - 1; l++)
                {
                  for (unsigned int it = 0; it <= max_n_it; ++it)
                  {
                    // Check for convergence and exit iteration loop
                    // if the off-diagonal element e[l] is zero
                    int m = l;
                    for (; m <= dim - 2; m++)
                    {
                      g = std::abs(w[m]) + std::abs(w[m + 1]);
                      if (std::abs(e[m]) + g == g)
                        break;
                    }
                    if (m == l)
                      break;

                    // Throw if no convergence is achieved within a
                    // stipulated number of iterations
                    if (it == max_n_it)
                    {
                      AssertThrow(
                          false,
                          ExcMessage(
                              "No convergence in iterative QL eigenvector algorithm.")) return std::
                                  array<std::pair<double, Tensor<1, dim, double>>, dim>();
                    }

                    // Calculate the shift..
                    g = (w[l + 1] - w[l]) / (e[l] + e[l]);
                    r = std::sqrt(g * g + 1.0);
                    // .. and then compute g = d_m - k_s for the
                    // plane rotation (Press2007a eq 11.4.22)
                    if (g > 0.0)
                      g = w[m] - w[l] + e[l] / (g + r);
                    else
                      g = w[m] - w[l] + e[l] / (g - r);

                    // Perform plane rotation, as is done in the
                    // standard QL algorithm, followed by Givens
                    // rotations to recover the tridiagonal form
                    s = c = 1.0;
                    p     = 0.0;
                    for (int i = m - 1; i >= l; i--)
                    {
                      f = s * e[i];
                      b = c * e[i];

                      // Branch to recover from underflow
                      if (std::abs(f) > std::abs(g))
                      {
                        c        = g / f;
                        r        = std::sqrt(c * c + 1.0);
                        e[i + 1] = f * r;
                        c *= (s = 1.0 / r);
                      }
                      else
                      {
                        s        = f / g;
                        r        = std::sqrt(s * s + 1.0);
                        e[i + 1] = g * r;
                        s *= (c = 1.0 / r);
                      }

                      g        = w[i + 1] - p;
                      r        = (w[i] - g) * s + 2.0 * c * b;
                      p        = s * r;
                      w[i + 1] = g + p;
                      g        = c * r - b;

                      // Form the eigenvectors
                      for (int k = 0; k < dim; k++)
                      {
                        t           = Q[k][i + 1];
                        Q[k][i + 1] = s * Q[k][i] + c * t;
                        Q[k][i]     = c * Q[k][i] - s * t;
                      }
                    }
                    w[l] -= p;
                    e[l] = g;
                    e[m] = 0.0;
                  }
                }
                // Structure the data to be outputted
                std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs;
                for (unsigned int e = 0; e < dim; ++e)
                {
                  eig_vals_vecs[e].first = w[e];

                  // The column "e" of Q contains the non-normalized
                  // eigenvector associated with the eigenvalue "e"
                  for (unsigned int a = 0; a < dim; ++a)
                  {
                    eig_vals_vecs[e].second[a] = Q[a][e];
                  }

                  // Normalize
                  Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
                  eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
                }
                return eig_vals_vecs;
              }


              std::array<std::pair<double, Tensor<1, 2, double>>, 2>
              hybrid(const dealii::SymmetricTensor<2, 2, double> &A)
              {
                const unsigned int dim = 2;

                // Calculate eigenvalues
                const std::array<double, dim> w = eigenvalues(A);

                std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs;

                double t, u; // Intermediate storage
                t = std::abs(w[0]);
                for (unsigned int i = 1; i < dim; ++i)
                {
                  u = std::abs(w[i]);
                  if (u > t)
                    t = u;
                }

                if (t < 1.0)
                  u = t;
                else
                  u = t * t;

                // Estimated maximum roundoff error
                const double error =
                    256.0 * std::numeric_limits<double>::epsilon() * u * u;

                // Store eigenvalues
                eig_vals_vecs[0].first = w[0];
                eig_vals_vecs[1].first = w[1];

                // Compute eigenvectors
                // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
                // https://math.stackexchange.com/a/1548616
                if (A[1][0] != 0.0)
                {
                  // First eigenvector
                  eig_vals_vecs[0].second[0] = w[0] - A[1][1];
                  eig_vals_vecs[0].second[1] = A[1][0];

                  // Second eigenvector
                  eig_vals_vecs[1].second[0] = w[1] - A[1][1];
                  eig_vals_vecs[1].second[1] = A[1][0];
                }
                else
                {
                  // First eigenvector
                  eig_vals_vecs[0].second[0] = w[0];
                  eig_vals_vecs[0].second[1] = 0.0;

                  // Second eigenvector
                  eig_vals_vecs[1].second[0] = 0.0;
                  eig_vals_vecs[1].second[1] = w[1];
                }
                // Normalize
                eig_vals_vecs[0].second /= eig_vals_vecs[0].second.norm();
                eig_vals_vecs[1].second /= eig_vals_vecs[1].second.norm();

                // If vectors are nearly linearly dependent, or if there might have
                // been large cancelations in the calculation of A[i][i] - w[0], fall
                // back to QL algorithm
                if (eig_vals_vecs[0].second * eig_vals_vecs[1].second > error)
                {
                  return ql_implicit_shifts(A);
                }

                return eig_vals_vecs;
              }


              std::array<std::pair<double, Tensor<1, 3, double>>, 3>
              hybrid(const dealii::SymmetricTensor<2, 3, double> &A)
              {
                const unsigned int dim = 3;
                double norm; // Squared norm or inverse norm of current eigenvector
                double t, u; // Intermediate storage

                // Calculate eigenvalues
                const std::array<double, dim> w = eigenvalues(A);

                t = std::abs(w[0]);
                for (unsigned int i = 1; i < dim; ++i)
                {
                  u = std::abs(w[i]);
                  if (u > t)
                    t = u;
                }

                if (t < 1.0)
                  u = t;
                else
                  u = t * t;

                // Estimated maximum roundoff error
                const double error =
                    256.0 * std::numeric_limits<double>::epsilon() * u * u;

                // Initialize the transformation matrix as the
                // identity tensor
                dealii::Tensor<2, dim, double> Q;
                Q[0][1] = A[0][1] * A[1][2] - A[0][2] * A[1][1];
                Q[1][1] = A[0][2] * A[0][1] - A[1][2] * A[0][0];
                Q[2][1] = A[0][1] * A[0][1];

                // Calculate first eigenvector by the formula
                //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
                Q[0][0] = Q[0][1] + A[0][2] * w[0];
                Q[1][0] = Q[1][1] + A[1][2] * w[0];
                Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1];
                norm    = Q[0][0] * Q[0][0] + Q[1][0] * Q[1][0] + Q[2][0] * Q[2][0];

                // If vectors are nearly linearly dependent, or if there might have
                // been large cancellations in the calculation of A[i][i] - w[0], fall
                // back to QL algorithm
                // Note that this simultaneously ensures that multiple eigenvalues do
                // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
                // i.e. all columns of A - w[0] * I are linearly dependent.
                if (norm <= error)
                {
                  return ql_implicit_shifts(A);
                }
                else // This is the standard branch
                {
                  norm = std::sqrt(1.0 / norm);
                  for (unsigned j = 0; j < dim; j++)
                    Q[j][0] = Q[j][0] * norm;
                }

                // Calculate second eigenvector by the formula
                //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
                Q[0][1] = Q[0][1] + A[0][2] * w[1];
                Q[1][1] = Q[1][1] + A[1][2] * w[1];
                Q[2][1] = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1];
                norm    = Q[0][1] * Q[0][1] + Q[1][1] * Q[1][1] + Q[2][1] * Q[2][1];
                if (norm <= error)
                {
                  return ql_implicit_shifts(A);
                }
                else
                {
                  norm = std::sqrt(1.0 / norm);
                  for (unsigned int j = 0; j < dim; j++)
                    Q[j][1] = Q[j][1] * norm;
                }

                // Calculate third eigenvector according to
                //   v[2] = v[0] x v[1]
                                     Q[0][2] = Q[1][0] * Q[2][1] - Q[2][0] * Q[1][1];
                                     Q[1][2] = Q[2][0] * Q[0][1] - Q[0][0] * Q[2][1];
                                     Q[2][2] = Q[0][0] * Q[1][1] - Q[1][0] * Q[0][1];

                                     // Structure the data to be outputted
                                     std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs;
                                     for (unsigned int e = 0; e < dim; ++e)
                                     {
                                       eig_vals_vecs[e].first = w[e];

                                       // The column "e" of Q contains the non-normalized
                                       // eigenvector associated with the eigenvalue "e"
                                       for (unsigned int a = 0; a < dim; ++a)
                                       {
                                         eig_vals_vecs[e].second[a] = Q[a][e];
                                       }

                                       // Normalize
                                       Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
                                       eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
                                     }
                                     return eig_vals_vecs;
              }



              std::array<double, 2>
              eigenvalues(const SymmetricTensor<2, 2, double> &T)
              {
                const double upp_tri_sq = T[0][1] * T[0][1];
                if (upp_tri_sq == 0.0)
                {
                  // The tensor is diagonal
                  std::array<double, 2> eig_vals = {{T[0][0], T[1][1]}};

                  // Sort from largest to smallest.
                  std::sort(eig_vals.begin(), eig_vals.end(), std::greater<double>());
                  return eig_vals;
                }
                else
                {
                  const double tr_T    = trace(T);
                  const double det_T   = determinant(T);
                  const double descrim = tr_T * tr_T - 4.0 * det_T;
                  Assert(
                      descrim > 0.0,
                      ExcMessage(
                          "The roots of the characteristic polynomial are complex valued."));
                  const double sqrt_desc = std::sqrt(descrim);

                  const std::array<double, 2> eig_vals = {
                      {(0.5 * (tr_T + sqrt_desc)),
                          (0.5 * (tr_T - sqrt_desc))}};
                  Assert(eig_vals[0] >= eig_vals[1],
                      ExcMessage("The eigenvalue ordering is incorrect."));
                  return eig_vals;
                }
              }

              std::array<double, 3>
              eigenvalues(const SymmetricTensor<2, 3, double> &T)
              {
                const double upp_tri_sq =
                    T[0][1] * T[0][1] + T[0][2] * T[0][2] + T[1][2] * T[1][2];
                if (upp_tri_sq == 0.0)
                {
                  // The tensor is diagonal
                  std::array<double, 3> eig_vals = {{T[0][0], T[1][1], T[2][2]}};

                  // Sort from largest to smallest.
                  std::sort(eig_vals.begin(), eig_vals.end(), std::greater<double>());
                  return eig_vals;
                }
                else
                {
                  // Perform an affine change to T, and solve a different
                  // characteristic equation that has a trigonometric solution.
                  // Decompose T = p*B + q*I , and set q = tr(T)/3
                  // and p = (tr((T - q.I)^{2})/6)^{1/2} . Then solve the equation
                  // 0 = det(\lambda*I - B) = \lambda^{3} - 3*\lambda - det(B)
                  // which has the solution
                  // \lambda = 2*cos(1/3 * acos(det(B)/2) +2/3*pi*k ); k = 0,1,2
                  // when substituting  \lambda = 2.cos(theta) and using trig identities.
                  const double tr_T = trace(T);
                  const double q    = tr_T / 3.0;
                  const double tmp1 = (T[0][0] - q) * (T[0][0] - q) +
                      (T[1][1] - q) * (T[1][1] - q) +
                      (T[2][2] - q) * (T[2][2] - q) + 2.0 * upp_tri_sq;
                  const double                        p = std::sqrt(tmp1 / 6.0);
                  const SymmetricTensor<2, 3, double> B =
                      double(1.0 / p) * (T - q * unit_symmetric_tensor<3, double>());
                  const double tmp_2 = determinant(B) / 2.0;

                  // The value of tmp_2 should be within [-1,1], however
                  // floating point errors might place it slightly outside
                  // this range. It is therefore necessary to correct for it.
                  // Note: The three results in the conditional may lead to different
                  //       number types when using Sacado numbers, so we cast them when
                  //       necessary to a consistent result type.
                  const double phi =
                      (tmp_2 <= -1.0 ?
                          numbers::PI / 3.0 :
                          (tmp_2 >= 1.0 ?
                              0.0 :
                              std::acos(tmp_2) / 3.0));

                  // Due to the trigonometric solution, the computed eigenvalues
                  // should be predictably in the order eig1 >= eig2 >= eig3...
                  std::array<double, 3> eig_vals = {
                      {static_cast<double>(q + 2.0 * p * std::cos(phi)),
                          static_cast<double>(0.0),
                          static_cast<double>(q + 2.0 * p *
                              std::cos(phi + (2.0 / 3.0 * numbers::PI)))}};
                  // Use the identity tr(T) = eig1 + eig2 + eig3
                  eig_vals[1] = tr_T - eig_vals[0] - eig_vals[2];

                  // ... however, when equal roots exist then floating point
                  // errors may make this no longer be the case.
                  // Sort from largest to smallest.
                  std::sort(eig_vals.begin(), eig_vals.end(), std::greater<double>());

                  return eig_vals;
                }
              }

                template <int dim>
                std::array<std::pair<double, Tensor<1, dim, double>>, dim>
                perform_eigenvector_decomposition(
                    const SymmetricTensor<2, dim, double> &T,
                    const SymmetricTensorEigenvectorMethod method)
                {
                  switch (method)
                  {
                    case SymmetricTensorEigenvectorMethod::hybrid:
                      return internal::SymmetricTensorImplementation::hybrid(T);
                      break;
                    case SymmetricTensorEigenvectorMethod::ql_implicit_shifts:
                      //return internal::SymmetricTensorImplementation::ql_implicit_shifts(
                      //  T);
                      break;
                    case SymmetricTensorEigenvectorMethod::jacobi:
                      //  return internal::SymmetricTensorImplementation::jacobi(T);
                      break;
                    default:
                      break;
                  }

                  AssertThrow(false, ExcNotImplemented());
                  return std::array<std::pair<double, Tensor<1, dim, double>>, dim>();
                }


              } // namespace SymmetricTensorImplementation
            } // namespace internal



    template <int dim>
    PointValuesEigenVectors<dim>::PointValuesEigenVectors ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      evaluation_points_cartesian (std::vector<Point<dim> >() ),
      point_values (std::vector<std::pair<double, std::vector<Vector<double> > > >() ),
      use_natural_coordinates (false)
    {}

    template <int dim>
    std::pair<std::string,std::string>
    PointValuesEigenVectors<dim>::execute (TableHandler &)
    {
      // if this is the first time we get here, set the next output time
      // to the current time. this makes sure we always produce data during
      // the first time step
      if (std::isnan(last_output_time))
        last_output_time = this->get_time() - output_interval;

      // see if output is requested at this time
      if (this->get_time() < last_output_time + output_interval)
        return std::pair<std::string,std::string>();

      // evaluate the solution at all of our evaluation points

      // the number of output components is the number of solution components,
      // + the components of the first strain rate eigen vector
      // + the components of the second strain rate eigen vector
      // + the components of possibly the third strain rate eigen vector
      const unsigned int n_solution_variables = this->introspection().n_components;
      const unsigned int n_output_variables =  n_solution_variables + (output_eigen_vectors ? dim*dim : 0);
      std::vector<Vector<double> >
      current_point_values (evaluation_points_cartesian.size(),
                            Vector<double> (n_solution_variables));

      std::vector<Vector<double> >
      all_current_point_values (evaluation_points_cartesian.size(),
                            Vector<double> (n_output_variables));

      std::vector<std::vector<Tensor<1,dim> > >
      current_point_gradients (evaluation_points_cartesian.size(),
                            std::vector<Tensor<1,dim> > (n_solution_variables));

      std::vector<std::vector<Tensor<1,dim> > >
      global_current_point_gradients (evaluation_points_cartesian.size(),
                            std::vector<Tensor<1,dim> > (n_solution_variables));

      const typename Introspection<dim>::ComponentIndices &component_indices = this->introspection().component_indices;

      for (unsigned int p=0; p<evaluation_points_cartesian.size(); ++p)
        {
          // try to evaluate the solution at this point. in parallel, the point
          // will be on only one processor's owned cells, so the others are
          // going to throw an exception. make sure at least one processor
          // finds the given point
          bool point_found = false;

          // first check if this process holds cells in the crust
          bool possibly_in_crust = false;
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();
          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
                if(this->get_geometry_model().depth (cell->center()) < 10000.)
                {
                  possibly_in_crust = true;
                  break;
                }


          if (possibly_in_crust)
          try
            {
              VectorTools::point_value(this->get_mapping(),
                                       this->get_dof_handler(),
                                       this->get_solution(),
                                       evaluation_points_cartesian[p],
                                       current_point_values[p]);
              point_found = true;

              // TODO only get the velocity gradients
              if (output_eigen_vectors)
              VectorTools::point_gradient(this->get_mapping(),
                                       this->get_dof_handler(),
                                       this->get_solution(),
                                       evaluation_points_cartesian[p],
                                       current_point_gradients[p]);
            }
          catch (const VectorTools::ExcPointNotAvailableHere &)
            {
              // ignore
            }
           catch (const GridTools::ExcPointNotFound<dim> &)
             {
             }

          // ensure that at least one processor found things
          const int n_procs = Utilities::MPI::sum (point_found ? 1 : 0, this->get_mpi_communicator());
          // actually, if none of the processors found the point, just continue to the next point
          // especially with the free surface, this can happen, or when a point is exactly on the boundary
          // of two cells
          // All variables will stay zero, as initialized.
          if (n_procs == 0)
            continue;
//          AssertThrow (n_procs > 0,
//                       ExcMessage ("While trying to evaluate the solution at point " +
//                                   Utilities::to_string(evaluation_points_cartesian[p][0]) + ", " +
//                                   Utilities::to_string(evaluation_points_cartesian[p][1]) +
//                                   (dim == 3
//                                    ?
//                                    ", " + Utilities::to_string(evaluation_points_cartesian[p][2])
//                                    :
//                                    "") + "), " +
//                                   "no processors reported that the point lies inside the " +
//                                   "set of cells they own. Are you trying to evaluate the " +
//                                   "solution at a point that lies outside of the domain?"
//                                  ));

          // Reduce all collected values into local Vector
          Utilities::MPI::sum (current_point_values[p], this->get_mpi_communicator(),
              current_point_values[p]);

          // Reduce all collected values into local Tensor
          for (unsigned int n=0; n<n_solution_variables; ++n)
            global_current_point_gradients[p][n] = Utilities::MPI::sum (current_point_gradients[p][n], this->get_mpi_communicator());

          // Normalize in cases where points are claimed by multiple processors
          if (n_procs > 1)
          {
            current_point_values[p] /= n_procs;
            for (unsigned int n=0; n<n_solution_variables; ++n)
              global_current_point_gradients[p][n] /= n_procs;
          }

          // copy the point values
          for (unsigned int n=0; n<n_solution_variables; ++n)
            all_current_point_values[p][n] = current_point_values[p][n];

          // If requested, compute the strain rate eigenvectors
          if (output_eigen_vectors)
          {
            // extract the primal variables
            Tensor<2,dim> grad_u;
            for (unsigned int d = 0; d<dim; ++d)
               grad_u[d] = global_current_point_gradients[p][component_indices.velocities[d]];

            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);

            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            std::array<std::pair<double, Tensor< 1, dim, double> >,dim > eigen_vectors = eigenvectors(compressible_strain_rate, SymmetricTensorEigenvectorMethod::hybrid);
            const Tensor<1,dim> S1 = eigen_vectors[0].first * eigen_vectors[0].second;
            const Tensor<1,dim> S2 = eigen_vectors[1].first * eigen_vectors[1].second;
            const Tensor<1,dim> S3 = eigen_vectors[dim-1].first * eigen_vectors[dim-1].second;

            for (unsigned int d=0; d<dim; ++d)
            {
              all_current_point_values[p](n_output_variables+d) = S1[d];
            }
            for (unsigned int d=0; d<dim; ++d)
            {
              all_current_point_values[p](n_output_variables+dim+d) = S2[d];
            }
            if (dim==3)
            for (unsigned int d=0; d<dim; ++d)
            {
              all_current_point_values[p](n_output_variables+2*dim+d) = S3[d];
            }
          }

        }

      // finally push these point values all onto the list we keep
      point_values.push_back (std::make_pair (this->get_time(),
                                              all_current_point_values));

      // now write all of the data to the file of choice. start with a pre-amble that
      // explains the meaning of the various fields
      const std::string filename = (this->get_output_directory() +
                                    "point_values.txt");
      std::ofstream f (filename.c_str());
      f << ("# <time> "
            "<evaluation_point_x> "
            "<evaluation_point_y> ")
        << (dim == 3 ? "<evaluation_point_z> " : "")
        << ("<velocity_x> "
            "<velocity_y> ")
        << (dim == 3 ? "<velocity_z> " : "")
        << "<pressure> <temperature>";
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        f << " <" << this->introspection().name_for_compositional_index(c) << ">";
      if(output_eigen_vectors)
      {
        f << " <e1_x> <e1_y> ";
        f << (dim == 3 ? "<e1_z> " : "");
        f << " <e2_x> <e2_y> ";
        f << (dim == 3 ? "<e2_z> " : "");
        if (dim == 3)
        {
        f << " <e3_x> <e3_y> ";
        f << (dim == 3 ? "<e3_z> " : "");
        }
      }
      f << '\n';

      for (std::vector<std::pair<double, std::vector<Vector<double> > > >::iterator
           time_point = point_values.begin();
           time_point != point_values.end();
           ++time_point)
        {
          Assert (time_point->second.size() == evaluation_points_cartesian.size(),
                  ExcInternalError());
          for (unsigned int i=0; i<evaluation_points_cartesian.size(); ++i)
            {
              f << /* time = */ time_point->first / (this->convert_output_to_years() ? year_in_seconds : 1.)
                << ' '
                << /* location = */ evaluation_points_cartesian[i] << ' ';

              for (unsigned int c=0; c<time_point->second[i].size(); ++c)
                {
                  // output a data element. internally, we store all point
                  // values in the same format in which they were computed,
                  // but we convert velocities to meters per year if so
                  // requested
                  if ((c>=component_indices.velocities[0] && c<=component_indices.velocities[dim-1])
                      &&
                      (this->convert_output_to_years() == true))
                    f << time_point->second[i][c] * year_in_seconds;
                  else
                    f << time_point->second[i][c];

                  f << (c != time_point->second[i].size()-1 ? ' ' : '\n');
                }
            }

          // have an empty line between time steps
          f << '\n';
        }

      AssertThrow (f, ExcMessage("Writing data to <" + filename +
                                 "> did not succeed in the `point values' "
                                 "postprocessor."));

      // Update time
      set_last_output_time (this->get_time());

      // return what should be printed to the screen. note that we had
      // just incremented the number, so use the previous value
      return std::make_pair (std::string ("Writing point values:"),
                             filename);
    }


     template <int dim>
     std::array<std::pair<double, Tensor<1, dim, double>>,dim>
     PointValuesEigenVectors<dim>::eigenvectors(const SymmetricTensor<2, dim, double> &T,
                  const SymmetricTensorEigenvectorMethod method) const
     {
       std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs = internal::SymmetricTensorImplementation::
           perform_eigenvector_decomposition(T, method);

       // Sort in descending order before output.
       std::sort(
           eig_vals_vecs.begin(),
           eig_vals_vecs.end(),
           internal::SymmetricTensorImplementation::SortEigenValuesVectors<dim>());
       return eig_vals_vecs;
      }


    template <int dim>
    void
    PointValuesEigenVectors<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Point values");
        {
          prm.declare_entry ("Time between point values output", "0",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "point values output. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry("Evaluation points", "",
                            // a list of points, separated by semicolons; each point has
                            // exactly 'dim' components/coordinates, separated by commas
                            Patterns::List (Patterns::List (Patterns::Double(), dim, dim, ","),
                                            0, Patterns::List::max_int_value, ";"),
                            "The list of points at which the solution should be evaluated. "
                            "Points need to be separated by semicolons, and coordinates of "
                            "each point need to be separated by commas.");
          prm.declare_entry("Use natural coordinates", "false",
                            Patterns::Bool (),
                            "Whether or not the Evaluation points are specified in "
                            "the natural coordinates of the geometry model, e.g. "
                            "radius, lon, lat for the chunk model. "
                            "Currently, natural coordinates for the spherical shell "
                            "and sphere geometries are not supported. ");
          prm.declare_entry("Output strain rate eigen vectors", "true",
                            Patterns::Bool (),
                            "Whether or not to also output the maximum horizontal "
                            "compressive stress and the stress regime at each point, "
                            "or only the solution variables. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    PointValuesEigenVectors<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Point values");
        {
          output_interval = prm.get_double ("Time between point values output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          const std::vector<std::string> point_list
            = Utilities::split_string_list(prm.get("Evaluation points"), ';');

          std::vector<std::array<double,dim> > evaluation_points;

          for (unsigned int p=0; p<point_list.size(); ++p)
            {
              const std::vector<std::string> coordinates
                = Utilities::split_string_list(point_list[p], ',');
              AssertThrow (coordinates.size() == dim,
                           ExcMessage ("In setting up the list of evaluation points for the <Point values> "
                                       "postprocessor, one of the evaluation points reads <"
                                       + point_list[p] +
                                       ">, but this does not correspond to a list of numbers with "
                                       "as many coordinates as you run your simulation in."));

              std::array<double,dim> point;
              for (unsigned int d=0; d<dim; ++d)
                point[d] = Utilities::string_to_double (coordinates[d]);
              evaluation_points.push_back (point);
            }

          use_natural_coordinates = prm.get_bool("Use natural coordinates");

          if (use_natural_coordinates)
            AssertThrow (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) == 0 &&
                         dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) == 0,
                         ExcMessage ("This postprocessor can not be used if the geometry "
                                     "is a sphere or spherical shell, because these geometries have not implemented natural coordinates."));

          // Convert the vector of coordinate arrays in Cartesian or natural
          // coordinates to a vector of Point<dim> of Cartesian coordinates.
          // For now, do nothing
          evaluation_points_cartesian.resize(evaluation_points.size());
          for (unsigned int p=0; p<evaluation_points.size(); ++p)
            {
              //if (use_natural_coordinates)
              //  evaluation_points_cartesian[p] = this->get_geometry_model().natural_to_cartesian_coordinates(evaluation_points[p]);
              //else
                for (unsigned int i = 0; i < dim; i++)
                  evaluation_points_cartesian[p][i] = evaluation_points[p][i];
            }
          output_eigen_vectors = prm.get_bool("Output strain rate eigen vectors");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    template <class Archive>
    void PointValuesEigenVectors<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &evaluation_points_cartesian
      & point_values
      & last_output_time;
    }


    template <int dim>
    void
    PointValuesEigenVectors<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);
      oa << (*this);

      status_strings["PointValuesEigenVectors"] = os.str();
    }


    template <int dim>
    void
    PointValuesEigenVectors<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("PointValuesEigenVectors") != status_strings.end())
        {
          std::istringstream is (status_strings.find("PointValuesEigenVectors")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }


    template <int dim>
    void
    PointValuesEigenVectors<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PointValuesEigenVectors,
                                  "point values eigen vectors",
                                  "A postprocessor that evaluates the solution (i.e., velocity, pressure, "
                                  "temperature, and compositional fields along with other fields that "
                                  "are treated as primary variables) at the end of every time step or "
                                  "after a user-specified time interval "
                                  "at a given set of points and then writes this data into the file "
                                  "<point\\_values.txt> in the output directory. The points at which "
                                  "the solution should be evaluated are specified in the section "
                                  "\\texttt{Postprocess/Point values} in the input file."
                                  "\n\n"
                                  "In the output file, data is organized as (i) time, (ii) the 2 or 3 "
                                  "coordinates of the evaluation points, and (iii) followed by the "
                                  "values of the solution vector at this point. The time is provided "
                                  "in seconds or, if the "
                                  "global ``Use years in output instead of seconds'' parameter is "
                                  "set, in years. In the latter case, the velocity is also converted "
                                  "to meters/year, instead of meters/second."
                                  "\n\n"
                                  "\\note{Evaluating the solution of a finite element field at "
                                  "arbitrarily chosen points is an expensive process. Using this "
                                  "postprocessor will only be efficient if the number of evaluation "
                                  "points or output times is relatively small. If you need a very large number of "
                                  "evaluation points, you should consider extracting this "
                                  "information from the visualization program you use to display "
                                  "the output of the `visualization' postprocessor.}")
  }
}
