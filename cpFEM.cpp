// centered pFEM
// Voronoi weigths are used in order to maintain good meshes
// Pressure is used to enforce incompressibility (constant Delanay volumes)

// half-step froggy predictor-corrector

#include"cpFEM.h"

#include"linear.h"
#include"simu.h"

sim_data simu;


int main() {


  // TODO: read parameter file
  
  const int init_iters = 0;
  const FT  init_tol2 = 1e-3;

  const int inner_iters= 10;
  const FT  inner_tol  = 1e-4;

  const  FT total_time = 2 * M_PI * 0.2 ; // one whole turn

  const std::string particle_file("particles.dat");
  const std::string diagram_file("diagram.dat");

  Triangulation T;

  cout << "Creating point cloud" << endl;

//  simu.do_perturb(0.1);
  create( T , 1.0 );
  number( T );

  //  set_vels_rotating( T );
  //  set_vels_Lamb_Oseen( T );

  volumes( T ); 
  linear algebra( T );
  algebra.copy( sfield_list::Dvol,  sfield_list::Dvol0);
  algebra.copy( sfield_list::Vvol,  sfield_list::Vvol0);

  //  set_vels_Gresho( T );
  
  // // testing .-
   // algebra.test_gradient();
   // algebra.test_Poisson();

   // draw( T , particle_file );
   // draw_diagram( T , diagram_file );  
   // return 0;
  

  // Init loop!
  
  // int iter=1;

  // for( ; iter < init_iters ; ++iter) {
  
  //   volumes( T ); 

  //   copy_weights( T ) ;

  //   //    algebra.solve_for_weights();

  //   FT dd = lloyds( T ) ;

  //   cout << " init loop , iter " << iter << " dd = " << dd << endl;
  //   if( dd < init_tol2) break;

  // }

  // volumes( T ); 
  // simu.set_dt( 0 );  
  // draw( T , particle_file     );
  // draw_diagram( T , diagram_file );  
  // return 0;
  //  cout << "Init loop converged in " << iter << " steps " << endl;
  
//  copy_weights( T ) ;

  set_vels_Gresho( T );

  FT d0;
  FT dt=0.001;

  cout << "Time step, dt = ";
  cin >> dt ;
  cout << endl << dt << endl;

  simu.set_dt( dt );

  FT spring_to_dt;
  cout << "Spring period / dt  = ";
  cin >> spring_to_dt;
  cout << endl << spring_to_dt << endl;

  bool springy = (spring_to_dt > 1e-10); // whether introduce springs or not

  if(!springy) cout << "No spring force will be added" << endl;
  
  // 31 dt is the value for G&M first simulation,
  // "Beltrami flow in the square"
  FT spring_period = spring_to_dt * dt;
  //  FT spring_period = 80 * dt;
  FT omega = 2 * M_PI /  spring_period ;

  FT spring = omega*omega; // factor that appears in the spring force

  // half-step (for e.g. leapfrog)
  FT dt2 = dt / 2.0 ;

  // whole step
  // FT dt2 = dt  ;

  //  algebra.solve_for_weights();

  draw( T , particle_file     );
  draw_diagram( T , diagram_file );
  
  std::ofstream log_file;
  log_file.open("main.log");
  log_file << " #  step   time   iters   kin_energy   L2_velocity " << endl;

  // special first iter.-
  // cout << " First iter, free ";
  // algebra.reset_p();
  // FT displ1 = move( T , dt , d0 );  
  // volumes( T ); 
  // simu.next_step();  simu.advance_time();
  // draw( T , particle_file     );
  // draw_diagram( T , diagram_file );
 
  do {
    simu.next_step();
    simu.advance_time( );

    backup( T );
    
    int iter = 1;

    algebra.u_star( );

    FT displ = 0;
 
    // full-step corrector loop

    for ( ; iter <= inner_iters ; iter++) {

      displ = move( T , dt2 , d0 );

      cout
	<< "********" << endl
	<< "Iter  " << iter
	<< " . Moved from previous (rel.): " << displ <<
	" ; from original (rel.): " << d0
	<< endl ;

      volumes( T ); 

      algebra.fill_diff_matrices();

      if(springy) {
        algebra.w_equation(); 
        copy_weights( T ) ;
	algebra.u_add_spring_force( spring*dt );
      }
      
      algebra.p_equation( dt2 ); 

      algebra.u_add_press_grad( dt2 );


      if( displ < inner_tol ) break;

      
    }

    displ = move( T , dt , d0 );

    update_half_velocity( T );

    volumes( T ); 
      
    cout
      << "Whole step  "
      << " : disp " << displ << endl ;

    draw( T , particle_file     );
    draw_diagram( T , diagram_file );

    log_file
      << simu.current_step() << "  "
      << simu.time() << "  "
      << iter-1 << " "
      << kinetic_E(T) << " "
      << L2_vel_Gresho(T) << " "
      << endl ;
    
  } while ( simu.time() < total_time );

  log_file.close();

  return 0;

}
