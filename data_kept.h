struct data_kept {
  int idx;
  Point pos;  // some position
  Point r0;
  weight w, w0;
  FT Dvol0;
  FT Vvol0;
  Vector_2 Dr;
  Vector_2 U, U0;
  Vector_2 Ustar;
  FT p, p0;


  Vector_2 dd;
  FT dd2;

 
  data_kept(const F_v_it fv) {
    idx = fv->idx();

    r0 = fv->r0.val();

    Vvol0 = fv->Vvol0.val();
    Dvol0 = fv->Dvol0.val();

    w = fv->w.val();
    w0 = fv->w0.val();

    Dr = fv->Dr.val();

    U = fv->U.val();
    U0 = fv->U0.val();
    Ustar = fv->Ustar.val();

    p = fv->p.val();
    p0 = fv->p0.val();

    dd  =  fv->dd.val();
    dd2 = fv->dd2.val(); ;


  }

  void restore(Vertex_handle fv) {
    fv->idx.set( idx );

    fv->r0.set( r0 );

    fv->Dvol0.set( Dvol0 );
    fv->Vvol0.set( Vvol0 );

    fv->w.set( w );
    fv->w0.set( w0 );

    fv->Dr.set( Dr );

    fv->U.set( U );
    fv->U0.set( U0 );
    fv->Ustar.set( Ustar );

    fv->p.set( p );
    fv->p0.set( p0 );


    fv->dd.set( dd ) ;
    fv->dd2.set( dd2 ); ;



  }

};
