#define PPE_DIV_SOURCE

//#include"pParticles.h"
#include"linear.h"
#include"fields_enum.h"
#include"simu.h"


// Solve for pressure
// The famous pressure Poisson equation

void linear::p_equation(const FT dt ) {

  cout << "Solving pressure equation " << endl;
  
  //  fill_Delta_DD(); // This may be important -- or not

  
  // A
  // Approximate Laplacian ~ Delta / V
  // VectorXd divUstar  =  DD_scalar_vfield( vfield_list::Ustar );
  // VectorXd p =  Delta_solver.solve( divUstar );
  // // times (-0.5), because the Laplacian is approximated by -2 Delta / V
  // vctr_to_field( -0.5 * p / ddt ,  sfield_list::p ) ;

  // return;
  
  // B
  //  Laplacian as div of grad :
#ifdef PPE_DIV_SOURCE

  VectorXd divUstar  =  divergence( vfield_list::Ustar );

  // Choice of Laplacian matrix:
  VectorXd p =  LL_solver.solve( divUstar  );
  //VectorXd p =  stiff_solver.solve( divUstar  );

  vctr_to_field( p / dt ,  sfield_list::p ) ;


#else
  // C
  // As B, but Dvol source.

  // diagnostics on volumes.-

  VectorXd vol  = field_to_vctr( sfield_list::Dvol ) ;

  VectorXd vol0  = field_to_vctr( sfield_list::Dvol0 ) ;

  VectorXd Dvol = vol.array() - vol0.array()  ;

  int N = vol.size();

  FT Dvol_sigma =  Dvol.array().square().sum() / N ; // / FT( vol.size() );
  FT Dvol_mean  =  vol.array().sum() / N ; // / FT( vol.size() );

  cout << "Pressure  "
       << " rel Dvol std dev: " << sqrt( Dvol_sigma ) / Dvol_mean 
       << endl;

  // C1: LL Laplacian
  //  VectorXd Dp  =  LL_solver.solve( Dvol );

  // C2: Delta Laplacian
  VectorXd Dp =  Delta_solver.solve( Dvol );

  vctr_to_field( -0.5 * Dp / ( ddt * ddt) , sfield_list::p  ) ;


  
#endif

  return;
}


void linear::u_add_press_grad( const FT dt ) {

  VectorXd gradPx,gradPy;

  gradient( sfield_list::p  ,  gradPx, gradPy);

  VectorXd vol  = field_to_vctr( sfield_list::Dvol );
  // perhaps mean vol would be just fine
  
  VectorXd Ustar_x, Ustar_y;

  vfield_to_vctrs(  vfield_list::Ustar , Ustar_x, Ustar_y );

  VectorXd U_x, U_y;

  U_x = Ustar_x.array() - dt * gradPx.array() / vol.array()  ;
  U_y = Ustar_y.array() - dt * gradPy.array() / vol.array() ;
  
  vctrs_to_vfield( U_x, U_y , vfield_list::U );

}
