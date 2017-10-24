#ifndef _JACOBISOLVER_H_
#define _JACOBISOLVER_H_

#include "SparseMatrix.hh"
#include <vector>

class JacobiSolver
{
public:
  /// solves until max norm of residual is less than tolerance, or iterates a_iter times.  returns final residual.
  double solve(
               vector<double>& a_phi,
               const SparseMatrix& a_A, 
               const vector<double>& a_rhs, 
               const double& a_tolerance, 
               int a_iter);
};

#endif
