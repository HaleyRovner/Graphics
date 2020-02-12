//color shading
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <FPT.h>
#include "M3d_matrix_tools.c"

int numpoints[10];
double x[10][15000], y[10][15000], z[10][15000], display_x[10][15000], display_y[10][15000];
int numpolys[10];
int psize[10][15000];
int con[10][15000][8];
int numobjects, poly_size, q, s, sign, action, k;
double xp[15000], yp[15000], x_prime[10][15000], y_prime[10][15000];
double half_angle = 40*M_PI/180;
double V[4][4], translated[4][4], rotated[4][4];
double cs, sn, avx, avy, avz;
double a[3], b[3];
int colnum;
double ambient, diffuse_max, specular_power, intensity, diffuse, specular;
double light_source[3] = {100, 200, 0};

typedef struct {
  int objectnum;
  int polynum;
  double distance;
}

THING;
int n = 0;
THING object[150000];

void init_array() {
  //assign all of the values
  n = 0;
  for (int i = 0; i < numobjects; i++) {
    for (int j = 0; j < numpolys[i]; j++) {
      object[n].objectnum = i;
      object[n].polynum = j;
      object[n].distance = z[i][con[i][j][0]];
      n++;
    }
  }
}


void color(double intensity, double rgb[3]){
  double temp_rgb[3], grey_value; 
  grey_value = ambient + diffuse_max;

  if(colnum == 0){
    temp_rgb[0] = 0.7; temp_rgb[1] = 0.7; temp_rgb[2] = 1;
  }
  else if(colnum == 1){
    temp_rgb[0] = 0.9; temp_rgb[1] = 0.4; temp_rgb[2] = 0.8;
  }
  else if(colnum == 2){
    temp_rgb[0] = 0.8; temp_rgb[1] = 0.6; temp_rgb[2] = 0.6;
  }
  else if(colnum == 3){
    temp_rgb[0] = 0.2; temp_rgb[1] = 0.5; temp_rgb[2] = 0.8;
  }
  else{
    temp_rgb[0] = 0.7; temp_rgb[1] = 0.9; temp_rgb[2] = 0.7;
  }

  if(intensity < grey_value){
    grey_value = (intensity/grey_value);
    rgb[0] = grey_value * temp_rgb[0];
    rgb[1] = grey_value * temp_rgb[1];
    rgb[2] = grey_value * temp_rgb[2];
  }
  else if(intensity > grey_value){
    rgb[0] = (1 - temp_rgb[0]) * ((1 - intensity) / (1 - grey_value)) + temp_rgb[0];
    rgb[1] = (1 - temp_rgb[1]) * ((1 - intensity) / (1 - grey_value)) + temp_rgb[1];
    rgb[2] = (1 - temp_rgb[2]) * ((1 - intensity) / (1 - grey_value)) + temp_rgb[2];
  }
  else{
    rgb[0] = grey_value * temp_rgb[0];
    rgb[1] = grey_value * temp_rgb[1];
    rgb[2] = grey_value * temp_rgb[2];
  }

}

double light_calculate(double xx[], double yy[], double zz[]) {
  double Lu[3], Nu[3], Ru[3], Eu[3], A[3], B[3];
  double NdL, NdE, RdE, square_rootN, square_rootL, square_rootE;
  double L[3], E[3], N[3];
 
  //C is the vector made from the points
  A[0] = xx[1] - xx[0]; A[1] = yy[1] - yy[0]; A[2] = zz[1] - zz[0] ;
  B[0] = xx[2] - xx[0]; B[1] = yy[2] - yy[0]; B[2] = zz[2] - zz[0] ;  

  
  //Lu is the unit vector for the light source and L is the original vector

  
  L[0] = (light_source[0] - xx[0]);
  L[1] = (light_source[1] - yy[0]);
  L[2] = (light_source[2] - zz[0]);
  square_rootL = sqrt((L[0]*L[0]) + (L[1]*L[1]) + (L[2]*L[2]));
  Lu[0] = L[0]/square_rootL; Lu[1] = L[1]/square_rootL; Lu[2] = L[2]/square_rootL;

  //Eu is the unit vector to the origin and E is the origin
  E[0] = (0 - xx[0]);
  E[1] = (0 - yy[0]);
  E[2] = (0 - zz[0]);
  square_rootE = sqrt((E[0]*E[0]) + (E[1]*E[1])+ (E[2]*E[2]));
  Eu[0] = E[0]/square_rootE; Eu[1] = E[1]/square_rootE; Eu[2] = E[2]/square_rootE;

  //Nu is the unit vector
  N[0] =   (A[1] * B[2]) - (B[1] * A[2]);
  N[1] = -((A[0] * B[2]) - (B[0] * A[2]));
  N[2] =   (A[0] * B[1]) - (B[0] * A[1]);
  square_rootN = sqrt((N[0]*N[0]) + (N[1]*N[1]) + (N[2]*N[2]));
  Nu[0] = N[0]/square_rootN; Nu[1] = N[1]/square_rootN; Nu[2] = N[2]/square_rootN;

  
  //NdL is the dot product between Nu and Lu and NdE is the dot product between Nu and Eu
  NdL = ((Nu[0] * Lu[0]) + (Nu[1] * Lu[1]) + (Nu[2] * Lu[2]));
  NdE = ((Nu[0] * Eu[0]) + (Nu[1] * Eu[1]) + (Nu[2] * Eu[2]));

  //if the origin is on the other side of the light source then make the intensity ambient
  if ((NdE * NdL) < 0) {
    return ambient ;
  }
  
  
  //if the dot products are less than zero, then you make the vector positive again
  if (NdE < 0 && NdL < 0) {
    Nu[0] *= -1;
    Nu[1] *= -1;
    Nu[2] *= -1;
  }
  //redo the dot products to account for the change from the if statement
  NdL = ((Nu[0] * Lu[0]) + (Nu[1] * Lu[1]) + (Nu[2] * Lu[2]));
  NdE = ((Nu[0] * Eu[0]) + (Nu[1] * Eu[1]) + (Nu[2] * Eu[2]));


    //Ru is the vector for the reflection off of the polygon
    Ru[0] = (((2 * NdL) * Nu[0]) - Lu[0]);
    Ru[1] = (((2 * NdL) * Nu[1]) - Lu[1]);
    Ru[2] = (((2 * NdL) * Nu[2]) - Lu[2]);

    //RdE is the dot product between R and E
    RdE = (Eu[0] * Ru[0]) + (Eu[1] * Ru[1]) + (Eu[2] * Ru[2]);

    if (RdE < 0) {
      RdE = 0;
    }

  
    
    diffuse = diffuse_max * NdL;
    specular = (1 - ambient - diffuse_max) * pow(RdE, specular_power);
    
    return ambient + diffuse + specular;

}

int compare(const void * p,
  const void * q) {
  //used in qsort
  THING * a, * b;

  a = (THING * ) p;
  b = (THING * ) q;

  if ((( * a).distance) < (( * b).distance)) return 1;
  else if ((( * a).distance) > (( * b).distance)) return -1;
  else return 0;
}

void draw_object() {
  double intensity;
  
  init_array();
  qsort(object, n, sizeof(THING), compare);

  int onum, pnum;
  double H, scale_factor, rgb[3];
  
  H = tan(half_angle);
  
  scale_factor = 300 / H;
  int m;
  int h;
  double xx[20],yy[20],zz[20] ;

  for (int i = 0; i < n; i++) {

    onum = object[i].objectnum;
    pnum = object[i].polynum;

    double tempx[1000], tempy[1000];
    m = psize[onum][pnum];

    for (int j = 0; j < m; j++) {
      h = con[onum][pnum][j];

      xx[j] = x[onum][h] ;
      yy[j] = y[onum][h] ;
      zz[j] = z[onum][h] ;
      
      
      tempx[j] = scale_factor*(xx[j]/zz[j]) + 300 ;
      tempy[j] = scale_factor*(yy[j]/zz[j]) + 300 ;
    }
    
    intensity = light_calculate(xx,yy,zz); //calls the calculation function here
    color(intensity, rgb);

    G_rgb(rgb[0], rgb[1], rgb[2]);
    G_fill_polygon(tempx, tempy, m);
    //G_rgb(0, 0, 0); 
    //G_polygon(tempx, tempy, m);
  }

}

int main(int argc, char ** argv) {
  FILE * fp;
  int key = 0;
  numobjects = argc - 1;

  printf("Input ambient, diffuse maximum, and specular power\n");
  scanf("%lf", & ambient);
  scanf("%lf", & diffuse_max);
  scanf("%lf", & specular_power);
  printf("The light source is at %lf, %lf, %lf\n", light_source[0], light_source[1], light_source[2]);

  G_init_graphics(600, 600);
  G_rgb(0, 0, 0);
  G_clear();

  //reading all of the information and storing into the variables
  for (key = 0; key < numobjects; key++) {
    fp = fopen(argv[key + 1], "r");
    if (fp == NULL) {
      printf("Broken File\n");
      exit(0);
    }
    fscanf(fp, "%d", & numpoints[key]);
    for (int i = 0; i < numpoints[key]; i++) {
      fscanf(fp, "%lf %lf %lf", & x[key][i], & y[key][i], & z[key][i]);
    }
    fscanf(fp, "%d", & numpolys[key]);
    for (int j = 0; j < numpolys[key]; j++) {
      fscanf(fp, "%d", & psize[key][j]);
      for (int k = 0; k < psize[key][j]; k++) {
        fscanf(fp, "%d", & con[key][j][k]);
      }
    }
    draw_object(); //draws it initially after all of the information is read in
  }

  q = '0';
  s = '0';
  sign = 1;
  action = 't';
  key = 0;
  colnum = 0;

  //while loop that does all of the rotation, translation, and flipping
  while (1 == 1) {
    M3d_make_identity(V);
    M3d_make_identity(translated), M3d_make_identity(rotated);
    q = G_wait_key();
    if (q == 'q') {
      exit(0);
    } else if (q == 'c') {
      sign = -sign;
    } else if (q == 't') {
      action = q;
    } else if (q == 'r') {
      action = q;
    } else if (q == 'p') {
      if(colnum != 4){
        colnum ++;
      }
      else{
        colnum = 0;
      }
    }
    else if (('0' <= q) && (q <= '9')) {
      k = q - '0';
      if (k < numobjects) {
        key = k;
      }
    } else if ((q == 'x') && (action == 't')) {
      if (sign > 0) {
        M3d_make_translation(V, 0.2, 0., 0.);
      } else {
        M3d_make_translation(V, -0.2, 0., 0.);
      }
    } else if ((q == 'y') && (action == 't')) {
      if (sign > 0) {
        M3d_make_translation(V, 0, 0.2, 0.);
      } else {
        M3d_make_translation(V, 0, -0.2, 0.);
      }
    } else if ((q == 'z') && (action == 't')) {
      if (sign > 0) {
        M3d_make_translation(V, 0., 0., 0.2);

      } else {
        M3d_make_translation(V, 0., 0., -0.2);
      }
    } else if ((q == 'x') && (action == 'r')) {
      avy = 0;
      avz = 0;
      for (int i = 0; i < numpoints[key]; i++) {
        avy += y[key][i];
        avz += z[key][i];
      }
      avy = avy / numpoints[key];
      avz = avz / numpoints[key];
      M3d_make_translation(translated, 0, -avy, -avz);
      M3d_mat_mult_points(x[key], y[key], z[key], translated, x[key], y[key], z[key], numpoints[key] + 1);
      if (sign > 0) {
        cs = cos(M_PI / 180);
        sn = sin(M_PI / 180);
        M3d_make_x_rotation_cs(rotated, cs, sn);
        M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key] + 1);

      } else {
        cs = cos(-M_PI / 180);
        sn = sin(-M_PI / 180);
        M3d_make_x_rotation_cs(rotated, cs, sn);
        M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key] + 1);

      }
      M3d_make_translation(V, 0, avy, avz);

    } else if ((q == 'y') && (action == 'r')) {
      avx = 0;
      avz = 0;
      for (int i = 0; i < numpoints[key]; i++) {
        avx += x[key][i];
        avz += z[key][i];
      }
      avx = avx / numpoints[key];
      avz = avz / numpoints[key];
      M3d_make_translation(translated, -avx, 0, -avz);
      M3d_mat_mult_points(x[key], y[key], z[key], translated, x[key], y[key], z[key], numpoints[key] + 1);
      if (sign > 0) {
        cs = cos(M_PI / 180);
        sn = sin(M_PI / 180);
        M3d_make_y_rotation_cs(rotated, cs, sn);
        M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key] + 1);
      } else {
        cs = cos(-M_PI / 180);
        sn = sin(-M_PI / 180);
        M3d_make_y_rotation_cs(rotated, cs, sn);
        M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key] + 1);
      }
      M3d_make_translation(V, avx, 0, avz);
    } else if ((q == 'z') && (action == 'r')) {
      avx = 0;
      avy = 0;
      for (int i = 0; i < numpoints[key]; i++) {
        avx += x[key][i];
        avy += y[key][i];
      }
      avx = avx / numpoints[key];
      avy = avy / numpoints[key];
      M3d_make_translation(translated, -avx, -avy, 0);
      M3d_mat_mult_points(x[key], y[key], z[key], translated, x[key], y[key], z[key], numpoints[key] + 1);
      if (sign > 0) {
        cs = cos(M_PI / 180);
        sn = sin(M_PI / 180);
        M3d_make_z_rotation_cs(rotated, cs, sn);
        M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key] + 1);
      } else {
        cs = cos(-M_PI / 180);
        sn = sin(-M_PI / 180);
        M3d_make_z_rotation_cs(rotated, cs, sn);
        M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key] + 1);
      }
      M3d_make_translation(V, avx, avy, 0);
    } else {
      printf("no action\n");
    }
    
    M3d_mat_mult_points(x[key], y[key], z[key], V, x[key], y[key], z[key], numpoints[key] + 1);
    G_rgb(0, 0, 0);
    G_clear();
    draw_object();

  }
}
