#include"linear.h"
#include"fields_enum.h"

// NOTE: these two are the divergence (but for a 1/V factor)
void linear::divergence(const vfield_list::take from , const sfield_list::take to )
{

  VectorXd vx, vy;
  vfield_to_vctrs( from , vx, vy );

  VectorXd div = DDx * vx + DDy*vy;

  vctr_to_field( div , to );

  return;
}

VectorXd linear::divergence(const vfield_list::take from )
{

  VectorXd vx, vy;

  vfield_to_vctrs( from , vx, vy );

  // cout << "vx cols " << vx.cols() << endl;
  // cout << "vx rows " << vx.rows() << endl;

  // cout << DDx << endl;
  // cout << "vx " << endl;
  // cout << vx << endl;


  return DDx * vx  + DDy * vy ;
}


// NOTE: this is the gradient (but for a 1/V factor)
// it features a minus sign, and transposition !!

void linear::gradient(const sfield_list::take from ,
			     VectorXd& Dx,VectorXd& Dy)
{

  VectorXd p = field_to_vctr( from );

  Dx = -DDx.transpose() * p;
  Dy = -DDy.transpose() * p;

  return;
}


// a FVM-like Laplacian (but for a 1/V factor)

VectorXd linear::Delta_laplacian(const sfield_list::take from )
{
  VectorXd p = field_to_vctr( from );

  return Delta * p;
}
