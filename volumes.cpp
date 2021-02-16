#include"cpFEM.h"
#include"simu.h"

// #define WARNING

typedef vector<Point> vvP;


// Compute Voronoi volumes (i.e. areas, in 2D)

void volumes(Triangulation& T) {

  const FT threshold=1e-10;
  const FT threshold2= threshold*threshold;
  FT totalV = 0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {
    fv->Dvol.reset();
    fv->Vvol.reset();
    fv->centroid.reset( );
  }

  // Volumes, FEM shape functions (Delaunay areas)
  
  for(F_f_it ff=T.finite_faces_begin();
      ff!=T.finite_faces_end();
      ff++)    {

    Vertex_handle v0 = ff->vertex(0);
    Vertex_handle v1 = ff->vertex(1);
    Vertex_handle v2 = ff->vertex(2);

    Point p0 = v0->point().point();
    Point p1 = v1->point().point();
    Point p2 = v2->point().point();

    Triangle tr( p0 ,  p1 , p2 );

    FT area = std::fabs( tr.area() );

    v0->Dvol += area / 3.0;
    v1->Dvol += area / 3.0;
    v2->Dvol += area / 3.0;
  }


  // Volumes, Voronoi cells

  for(F_e_it fe=T.finite_edges_begin();
      fe!=T.finite_edges_end();
      fe++)    {

    Face_handle f =  fe -> first ;
    int i0 = fe -> second;

    Vertex_handle vi = f->vertex( (i0+1) % 3);
    Point pi = vi->point().point();

    Vertex_handle vj = f->vertex( (i0+2) % 3);
    Point pj = vj->point().point();
    
    CGAL::Object o = T.dual(fe);

    const Segment * Vor_segment = CGAL::object_cast<Segment>( &o );

    if (! Vor_segment ) continue;

    Point p1 = Vor_segment->source() ;
    Point p2 = Vor_segment->target() ;

    // Triangle tri( pi ,  p1 , p2 );
    // Triangle trj( pj ,  p1 , p2 );

    // FT ar_i = std::fabs( tri.area() );
    // FT ar_j = std::fabs( trj.area() );

    // vi->vol +=  ar_i ;
    // vj->vol +=  ar_j ;

    // totalV += ar_i;
    // totalV += ar_j;

    {

      // https://pdfs.semanticscholar.org/b561/d4242952bce7bf986ed670c43532739809d4.pdf


      Vector_2 vi1 = p1 - pi ;
      Vector_2 vi2 = p2 - pi ;

      FT ar_i = std::fabs( vi1.x() * vi2.y() - vi1.y() * vi2.x() ) / 2.0 ;

      vi->Vvol +=  ar_i ;

      totalV += ar_i;

      FT x_cm = ( p1.x() + p2.x() + pi.x() ) / 3.0;
      FT y_cm = ( p1.y() + p2.y() + pi.y() ) / 3.0;
 
      vi->centroid =  vi->centroid.val() + ar_i * Vector_2( x_cm , y_cm);

    }

    {
      Vector_2 vj1 = p1 - pj ;
      Vector_2 vj2 = p2 - pj ;

      FT ar_j = std::fabs( vj1.x() * vj2.y() - vj1.y() * vj2.x() ) / 2.0 ;
      vj->Vvol +=  ar_j ;
      totalV += ar_j;
      
      FT x_cm = ( p1.x() + p2.x() + pj.x() ) / 3.0;
      FT y_cm = ( p1.y() + p2.y() + pj.y() ) / 3.0;
 
      vj->centroid =  vj->centroid.val() + ar_j * Vector_2( x_cm , y_cm);

    }

    // // CGAL::ORIGIN needed because points cannot be added, or multiplied
    // Vector_2 tri_ctr_v = CGAL::centroid( tri ) - CGAL::ORIGIN;
    // Vector_2 trj_ctr_v = CGAL::centroid( trj ) - CGAL::ORIGIN;

    // vi->centroid =  vi->centroid.val() + ar_i * tri_ctr_v;
    // vj->centroid =  vj->centroid.val() + ar_j * trj_ctr_v;

  }

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    FT a = fv->Vvol.val();

    Vector_2 ctr_v = fv->centroid.val() - CGAL::ORIGIN;

    Point cc = CGAL::ORIGIN + ctr_v / a;

    fv->centroid = cc;

    Point p  = fv->point().point();

    //    Vector_2 dA =  a *( cc - p );
    //   fv->dd.set( dA );
    //    fv->dd2.set(  dA.squared_length()  );

    Vector_2 cp =  p - cc ;

    fv->dd.set( cp );

    FT cp_l2 = cp.squared_length() ;

    fv->dd2.set( cp_l2  );

  }


  //  cout << "Volumes: total = " << totalV << " ; ";

  FT  inner_V = 0;
  int inner   = 0;
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)
    if (fv->idx() > -1 ) {
      inner_V +=  fv->Vvol.val();
      ++inner;
    }

  //  cout << "inner = " << inner_V << " ; ";

  simu.set_totalV( totalV );
  simu.set_innerV( inner_V );

  //  cout << "mean inner V = " << simu.meanV() << endl;

  
}
