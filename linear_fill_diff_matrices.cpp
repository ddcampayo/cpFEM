// Function to fill the various matrices involved. All of
// them. May be split in the future


// Glossary:

// Delta_ij = d V_i / d w_j   (inf.  change of volume of cell i due to change in weight of cell j)

// DD_ij =  d V_i / d r_j, aka "L"
// (inf.  change of _Delaunay_ volume of cell i due to change in position of cell j)
//  , a matrix of vectors, thus stored as DD_ij_x and DD_ij_y

// LL = - DD (1/Dvol) DD^t, Laplacian involved in Ralphson-Newton methods ( notice the sign ),
// as the divergence of the gradient

// stiff_ij = int grad f_i grad f_j, appears in the FEM Laplacian

#include"linear.h"

void linear::fill_diff_matrices( void ) {

  //  std::cout << " Filling Delta _and_ capital D matrices" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet>   aa, ax, ay, ss;

  typedef std::map<int,FT> diag_map;
  diag_map   dd; // diagonal
  diag_map   dd_x, dd_y;
  diag_map   ds;

  int N=1;

  F_e_it eit = T.finite_edges_begin();

  for ( ; eit !=T.finite_edges_end(); ++eit) {

    Face_handle f =  eit -> first ;
    int i0 = eit -> second;

    int i3 = i0 ;
    Vertex_handle v3  = f->vertex( i3 );
    Point p3 = v3->point().point() ;
    
    Vertex_handle v1 = f->vertex( (i3+1) % 3);
    //    Point p1 = v1->point().point() ;

    Vertex_handle v3p = T.mirror_vertex( f , i3 );
    Point p3p = v3p->point().point() ;

    Vector_2 v_3p_3 = p3 - p3p ;

    // negative right angle turn
    Vector_2 v_3p_3_perp = Vector_2( v_3p_3.y() , -v_3p_3.x() );

    //           3
    //
    //    1   -------  2  
    //
    //          3' 

    
    //    Triangle tr( v1->point().point() , v3p->point().point() , v3->point().point() );

    //    CGAL::Orientation ori = tr.orientation();
    //    if( ori == CGAL::NEGATIVE ) v33_perp = -v33_perp;

    // // is this check needed?
    // I don't think so, the orientation remains the same even if
    // vertex 1 would be on other side of edge (3 & 3' would reverse)
    // Vector_2 v_1_3p = p3p - p1 ;

    // CGAL::Orientation ori = CGAL::orientation(  v_1_3p , v_3p_3 );
    // if( ori == CGAL::RIGHT_TURN ) v_3p_3_perp = -v_3p_3_perp;

    Vector_2 DDij = v_3p_3_perp / 6.0 ;
    Vector_2 DDji = -DDij;

    Vertex_handle v0 = v3;
    Vertex_handle v2 = f->vertex( (i0+2) % 3); // aka v2

    Point p1 = v1->point().point();
    Point p2 = v2->point().point();
    Point p0 = p3;
    

    FT ll0 = Vector_2( p1 , p2).squared_length();
    FT ll1 = Vector_2( p2 , p0).squared_length();
    FT ll2 = Vector_2( p0 , p1).squared_length();

    FT ar = CGAL::area( p0 , p1 , p2);

    Point p0p = p3p;

    FT ll1p = Vector_2( p0p , p2).squared_length();
    FT ll2p = Vector_2( p1 , p0p).squared_length();

    FT arp = CGAL::area( p0p , p2 , p1);

    FT st_ij =
      -(ll0 - ll1  - ll2 )/ar /8
      -(ll0 - ll1p - ll2p)/arp/8;

    // FT st_ii = ll1 /ar /8 + ll1p /arp /8 ;
    // FT st_jj = ll2 /ar /8 + ll2p /arp /8 ;
 
    
    Vertex_handle vi = v1 ; // f->vertex( (i0+1) % 3);
    Vertex_handle vj = v2 ;

    
    int i = vi->idx();
    int j = vj->idx();

    CGAL::Object o = T.dual(eit);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    FT Aij = std::sqrt( Vor_segment->squared_length() );

    //    cout << " A0 = " << Aij << endl;
    
    Point pi = vi->point().point();
    Point pj = vj->point().point();

    Vector_2 eij = pj - pi;

    // experimental
    //    Point mij = pi + 0.5 * eij;  // mid-point between i and j
    
    FT lij2 =  eij.squared_length() ;
    FT lij = std::sqrt( lij2 );
  
    FT ddelta = - 0.5 * Aij / lij;

    if( (i >= 0 ) && ( j >= 0) ) {
      aa.push_back( triplet( i, j,  ddelta ));
      aa.push_back( triplet( j, i,  ddelta ));

      ax.push_back( triplet( i, j,  DDij.x() ));
      ay.push_back( triplet( i, j,  DDij.y() ));

      ax.push_back( triplet( j, i,  DDji.x() ));
      ay.push_back( triplet( j, i,  DDji.y() ));

      ss.push_back( triplet( i, j,  st_ij ));
      ss.push_back( triplet( j, i,  st_ij ));
      
    }

    // diagonal terms

    if (i >= 0 ) {
      dd[ i ]  -= ddelta;

      // this diagonal is 0, but just in case
      dd_x[ i ] -= DDji.x();
      dd_y[ i ] -= DDji.y();

      //      ds[ i ]  += st_ii;

      ds[ i ]  -= st_ij;
    }

    if (j >= 0 ) {
      dd[ j ]  -= ddelta;

      dd_x[ j ] -= DDij.x();
      dd_y[ j ] -= DDij.y();

      //      ds[ j ]  += st_jj;
      ds[ j ]  -= st_ij;

    }

    if( i+1 > N ) { N = i+1 ; } // keep maximum index

    //    cout << i << "  " << j << "  " << ddelta << endl;
  }

  // Add diagonal terms .-
  
  for( diag_map::const_iterator it = dd.begin(); it != dd.end(); ++it ) {
    int i = it->first;
    FT diag = it->second;
    aa.push_back( triplet( i, i,  diag ));

    //    cout << i << "  " << i << "  " << diag << endl;

  }
  
  for( diag_map::const_iterator it = dd_x.begin(); it != dd_x.end(); ++it ) {
    int i    = it->first;
    FT diagx = it->second;
    ax.push_back( triplet( i, i,  diagx ));
  }

  
  for( diag_map::const_iterator it = dd_y.begin(); it != dd_y.end(); ++it ) {
    int i    = it->first;
    FT diagy = it->second;
    ay.push_back( triplet( i, i,  diagy ));
  }


  for( diag_map::const_iterator it = ds.begin(); it != ds.end(); ++it ) {
    int i    = it->first;
    FT diag = it->second;
    ss.push_back( triplet( i, i,  diag ));
  }


  
  Delta.resize( N , N );
  Delta.setFromTriplets(aa.begin(), aa.end());
  // std::cout << " Filled Delta  matrix" << std::endl;
  // cout << "matrix size " << Delta.rows() << " times " << Delta.cols() << endl;

  DDx.resize( N , N );
  DDx.setFromTriplets(ax.begin(), ax.end());
  // std::cout << " Filled DDx  matrix" << std::endl;
  // cout << "matrix size " << DDx.rows() << " times " << DDx.cols() << endl;

  DDy.resize( N , N );
  DDy.setFromTriplets(ay.begin(), ay.end());
  // std::cout << " Filled DDy  matrix" << std::endl;
  // cout << "matrix size " << DDy.rows() << " times " << DDy.cols() << endl;

  stiff.resize( N , N );
  stiff.setFromTriplets(ss.begin(), ss.end());

  
  // set up solvers .-
  
  VectorXd vol  = field_to_vctr( sfield_list::Dvol ) ;
  VectorXd inv_vol  = 1.0 / vol.array() ;

  LL =
    - DDx * inv_vol.asDiagonal() * DDx.transpose()
    - DDy * inv_vol.asDiagonal() * DDy.transpose() ;
  
  LL_solver.compute( LL );

  if( LL_solver.info() != Eigen::Success ) {
    std::cout << "Failure decomposing LL matrix " << endl;
  }

  Delta_solver.compute( Delta );

  if(Delta_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing Delta " << //minus 1
      " matrix\n";
  }


  stiff_solver.compute( stiff );

  if(stiff_solver.info()!=Eigen::Success) {
    std::cout << "Failure decomposing stiffness " << //minus 1
      " matrix\n";
  }


  
  return;

}

