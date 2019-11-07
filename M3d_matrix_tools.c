#include <stdio.h>
#include <math.h>
/*
( x')          (x)
( y')  =   M * (y)
( 1 )          (1)
instead of (x',y',1) = (x,y,1) * M
*/
int M3d_print_mat (double a[4][4])
{
 int r,c ;
 for (r = 0 ; r < 4 ; r++ ) {
     for (c = 0 ; c < 4 ; c++ ) {
          printf(" %12.4lf ",a[r][c]) ;
     }
     printf("\n") ;
 }
 return 1 ;
}
int M3d_copy_mat (double a[4][4], double b[4][4])
// a = b
{
 int r,c ;
 for (r = 0 ; r < 4 ; r++ ) {
     for (c = 0 ; c < 4 ; c++ ) {
          a[r][c] = b[r][c] ;
     }
 }
 return 1 ;
}
int M3d_make_identity (double a[4][4])
// a = I
{
 int r,c ;
 for (r = 0 ; r < 4 ; r++ ) {
     for (c = 0 ; c < 4 ; c++ ) {
          if (r == c) a[r][c] = 1.0 ;
              else    a[r][c] = 0.0 ;
     }
 }
 return 1 ;
}
int M3d_make_translation (double a[4][4], double dx, double dy, double dz)
{
 M3d_make_identity(a) ;
 a[0][3] =  dx ;  a[1][3] = dy ; a[2][3] = dz;
 return 1 ;
}
int M3d_make_scaling (double a[4][4], double sx, double sy, double sz)
{
 M3d_make_identity(a) ;
 a[0][0] =  sx ;  a[1][1] = sy ; a[2][2] = sz;
 return 1 ;
}
/*int M3d_make_rotation_radians (double a[4][4],  double radians)
{
 double cs,sn ;
 cs =  cos( radians ) ;
 sn =  sin( radians ) ;
 M3d_make_identity(a) ;
 a[0][0] =   cs ;  a[0][1] = -sn ;
 a[1][0] =   sn ;  a[1][1] =  cs ;
 return 1 ;
}
int M3d_make_rotation_degrees (double a[4][4],  double degrees)
{
 double radians, cs,sn ;
 radians = degrees*M_PI/180 ;
 cs =  cos( radians ) ;
 sn =  sin( radians ) ;
 M2d_make_identity(a) ;
 a[0][0] =   cs ;  a[0][1] = -sn ;
 a[1][0] =   sn ;  a[1][1] =  cs ;
 return 1 ;
}
*/
int M3d_make_z_rotation_cs (double a[4][4], double cs, double sn){
// this one assumes cosine and sine are already known
 M3d_make_identity(a) ;
 a[0][0] =   cs ;  a[0][1] = -sn ;
 a[1][0] =   sn ;  a[1][1] =  cs ;
 return 1 ;
}
int M3d_make_x_rotation_cs(double a[4][4], double cs, double sn){
 M3d_make_identity(a);
 a[1][1] = cs; a[1][2] = -sn;
 a[2][1] = sn; a[2][2] = cs;
 return 1;
}
int M3d_make_y_rotation_cs(double a[4][4], double cs, double sn){
 M3d_make_identity(a);
 a[0][0] = cs; a[0][2] = sn;
 a[2][0] = -sn; a[2][2] = cs;
 return 1;
}
int M3d_mat_mult (double res[4][4], double a[4][4], double b[4][4])
// res = a * b
// this is SAFE, i.e. the user can make a call such as
// M2d_mat_mult(p,  p,q) or M2d_mat_mult(p,  q,p) or  M2d_mat_mult(p, p,p)
{
 double sum ;
 int k ;
 int r,c ;
 double tmp[4][4] ;
 for (r = 0 ; r < 4 ; r++ ) {
     for (c = 0 ; c < 4 ; c++ ) {
          sum = 0.0 ;
          for (k = 0 ; k < 4 ; k++) {
                sum = sum + a[r][k]*b[k][c] ;
          }
          tmp[r][c] = sum ;
     }
 }
 M3d_copy_mat (res,tmp) ;
 return 1 ;
}
int M3d_mat_mult_pt (double P[3], double m[4][4], double Q[3])
// P = m*Q
// SAFE, user may make a call like M2d_mat_mult_pt (W, m,W) ;
{
 double u,v,w ;
 u = m[0][0]*Q[0] + m[0][1]*Q[1] + m[0][2]*Q[2] + m[0][3] ;
 v = m[1][0]*Q[0] + m[1][1]*Q[1] + m[1][2]*Q[2] + m[1][3] ;
 w = m[2][0]*Q[0] + m[2][1]*Q[1] + m[2][2]*Q[2] + m[2][3];
 P[0] = u ;
 P[1] = v ;
 P[2] = w ;
 return 1 ;
}
int M3d_mat_mult_points (double *X, double *Y, double *Z,
                        double m[4][4],
                        double *x, double *y, double *z, int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// | 1  1  1 ...|       | 1  1  1 ...|
// SAFE, user may make a call like M2d_mat_mult_points (x,y, m, x,y, n) ;
{
 double u,v,t ;
 int i ;
 for (i = 0 ; i < numpoints ; i++) {
   u = m[0][0]*x[i] + m[0][1]*y[i] + m[0][2]*z[i] + m[0][3] ;
   v = m[1][0]*x[i] + m[1][1]*y[i] + m[1][2]*z[i] + m[1][3] ;
   t = m[2][0]*x[i] + m[2][1]*y[i] + m[2][2]*z[i] + m[2][3] ;
   X[i] = u ;
   Y[i] = v ;
   Z[i] = t ;
 }
 return 1 ;
}