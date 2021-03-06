#ifndef _LINEAR_H_
#define _LINEAR_H_

#include"cpFEM.h"

#ifdef CHOLMOD
 #include <Eigen/CholmodSupport>
#else
// #include <Eigen/IterativeLinearSolvers>
 #include <Eigen/SparseCholesky>
//#include <Eigen/SparseQR>
#endif

#include <iostream>
#include <vector>

#include <unsupported/Eigen/SparseExtra>

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
//using Eigen::ConjugateGradient;

//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");

class linear {
 public:  // constructor
 linear(Triangulation& TT) : T(TT) {}

  void fill_diff_matrices();

  void w_equation( void );
  void p_equation(const FT dt  );

  void u_star( void );
  void reset_p( void );

  void u_add_press_grad( const FT dt ) ;

  void u_add_spring_force( const FT kdt );

  void divergence(const vfield_list::take from , const sfield_list::take to );
  VectorXd divergence(const vfield_list::take from );

  void gradient(const sfield_list::take from ,
		VectorXd& Dx,VectorXd& Dy);

  VectorXd Delta_laplacian(const sfield_list::take from );

  void copy(const sfield_list::take from, sfield_list::take to  );

  void dd2_stats( void ) ;

  void test_gradient( void );
  void test_Poisson( void );
  
  //  Vector_2 values_at_v(const Point& p0, const vfield_list::take v_field) ;

private:

  
  Triangulation& T; // Its triangulation

  typedef   SparseMatrix<double>  SpMat;
  typedef Eigen::Triplet<double> triplet;

  SpMat Delta;
  SpMat DDx, DDy;
  SpMat LL;
  SpMat stiff;
  
  
  VectorXd field_to_vctr(const sfield_list::take sf );
  void vctr_to_field(const VectorXd& vv, const sfield_list::take sf  );
  void vfield_to_vctrs(const vfield_list::take vf , VectorXd& vx, VectorXd& vy ) ;
  void vctrs_to_vfield(const VectorXd& vx, const VectorXd& vy , const vfield_list::take vf ) ;

#define DIRECT_SOLVER

#ifdef DIRECT_SOLVER
#ifdef CHOLMOD
  Eigen::CholmodSupernodalLLT<SpMat> Delta_solver;
#else
  Eigen::SimplicialLDLT<SpMat> Delta_solver;
  Eigen::SimplicialLDLT<SpMat> LL_solver;
  Eigen::SimplicialLDLT<SpMat> stiff_solver;
#endif
#else
  Eigen::BiCGSTAB<SpMat> solver_stiffp1;
#endif


};


#endif
