//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"



void linear::u_star( void ) {
  VectorXd usx, usy;
  
  vfield_to_vctrs( vfield_list::U0 , usx, usy );

  vctrs_to_vfield( usx, usy, vfield_list::Ustar );

  return;
}


void linear::reset_p( void ) {

  VectorXd p = field_to_vctr( sfield_list::p );

  vctr_to_field( 0*p , sfield_list::p );

  return;
}



void linear::copy(const sfield_list::take from, sfield_list::take to  ) {

  VectorXd p = field_to_vctr( from );

  vctr_to_field( p , to );

  return;
}



void linear::test_gradient( void ) {

  // this makes u = grad p, for debugging purposes
  set_pressure( T );
  volumes( T );

  fill_diff_matrices();

  VectorXd gradPx,gradPy;

  gradient( sfield_list::p  ,  gradPx, gradPy);

  VectorXd vol  = field_to_vctr( sfield_list::Dvol );

  VectorXd U_x, U_y;

  U_x = gradPx.array() / vol.array()  ;
  U_y = gradPy.array() / vol.array() ;

  vctrs_to_vfield( U_x, U_y , vfield_list::U );

  return;
}

void linear::test_Poisson( void ) {

  // this makes p0 = lapl^-1 p, for debugging purposes

  set_pressure( T );
  volumes( T );

  fill_diff_matrices();

  VectorXd source  = field_to_vctr( sfield_list::p );
  VectorXd vol  = field_to_vctr( sfield_list::Dvol );
  VectorXd vol_source = source.array() * vol.array() ;

  VectorXd p0 =  LL_solver.solve( vol_source );

  vctr_to_field( p0 , sfield_list::p0  ) ;

  return;
  
}
