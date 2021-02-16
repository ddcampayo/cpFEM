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



// void linear::test_operators( void ) {

//   // this makes u = grad p, for debugging purposes

//   fill_Delta_DD();

//   VectorXd gradPx,gradPy;

//   DD_times_sfield( sfield_list::p  ,  gradPx, gradPy);

//   VectorXd vol  = field_to_vctr( sfield_list::Vvol );

//   VectorXd U_x, U_y;

//   U_x = gradPx.array() / vol.array()  ;
//   U_y = gradPy.array() / vol.array() ;

//   vctrs_to_vfield( U_x, U_y , vfield_list::U );
  
// }
