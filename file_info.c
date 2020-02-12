#include <stdio.h>
#include <string.h>
#include <math.h>
#include <FPT.h>
#include "M2d_matrix_tools.c"

int numpoints[10];
double x[10][1500], y[10][1500];
int numpolys[10];
int psize[10][1000];
int con[10][1000][20];
double red[10][1000], green[10][1000], blue[10][1000];
int numobjects;
int q, n, h,r;
double xp[100], yp[100];
double xcentered[10][1000], ycentered[10][1000];
double xscaled[10][1000], yscaled[10][1000];
double xorigin[10][1000], yorigin[10][1000];
double xmiddle[10], ymiddle[10];
double m[3][3];
double scalefactor[10], xoriginadjust[10], yoriginadjust[10];
double translated[3][3], scaled[3][3], v[3][3];

void center_and_scale(int key){
  double xhigh, yhigh, xlow, ylow, xadjust, yadjust;
  xhigh = x[key][0];
  xlow = x[key][0];
  yhigh = y[key][0];
  ylow = y[key][0];
  for(int i =1; i<numpoints[key]; i++){
    if(x[key][i]>xhigh){
      xhigh = x[key][i]; 
    }
    if(x[key][i]<xlow){
      xlow = x[key][i];
    }
    if(y[key][i]>yhigh){
      yhigh = y[key][i];
    }
    if(y[key][i]<ylow){
      ylow = y[key][i];
    }
    xmiddle[key] = (xhigh + xlow)/2;
    ymiddle[key] = (yhigh + ylow)/2;
    xoriginadjust[key] =-xmiddle[key];
    yoriginadjust[key] =-ymiddle[key];
  }
   
   if((xhigh-xlow)<=(yhigh-ylow)){
     scalefactor[key]=600/(yhigh-ylow);
   }
   else{
     scalefactor[key]=600/(xhigh-xlow);
   }
}


void draw_polygon(int key){
 for(int o = 0; o<numpolys[key]; o++){
      n = psize[key][o];
      for (int p = 0; p< n; p++){
	h = con[key][o][p] ;
	xp[p] = x[key][con[key][o][p]];
	yp[p] = y[key][con[key][o][p]];	
      }
      G_rgb(red[key][o], green[key][o], blue[key][o]);
      G_fill_polygon(xp, yp, n);
    }
}

void spin (int key){
  r=G_wait_key();
  double temp ;
  double angle = M_PI/180;
  int l ;
  while(r==32){
  
    for (l = 0 ; l < numpoints[key] ; l++) {
     xcentered[key][l] = xcentered[key][l] - 300;
     ycentered[key][l] = ycentered[key][l] - 300;
    }

    for (l = 0 ; l < numpoints[key] ; l++) {
      temp              = xcentered[key][l]*cos(angle) - ycentered[key][l]*sin(angle) ;
      ycentered[key][l] = xcentered[key][l]*sin(angle) + ycentered[key][l]*cos(angle) ;
      xcentered[key][l] = temp ;
    }

    for (l = 0 ; l < numpoints[key] ; l++) {
     xcentered[key][l] = xcentered[key][l] + 300;
     ycentered[key][l] = ycentered[key][l] + 300;
    }

    for(int o = 0; o<numpolys[key]; o++){
      n = psize[key][o];
      for (int p = 0; p< n; p++){
	h = con[key][o][p] ;
	xp[p] = xcentered[key][con[key][o][p]];
	yp[p] = ycentered[key][con[key][o][p]];	
      }
      
      G_rgb(red[key][o], green[key][o], blue[key][o]);
      G_fill_polygon(xp, yp, n);
     
    }
    
    r = G_wait_key();
     G_rgb(1,1,1);
     G_clear();
    
  }

}


int main(int argc, char **argv){
  FILE *fp;
  int onum ;

  numobjects=argc-1;

  G_init_graphics (600, 600);

  for(onum=0; onum<numobjects; onum++){ 
    fp=fopen(argv[onum+1], "r");
    
    if(fp==NULL){
      printf("Can't open file %s\n", argv[2]);
      exit(0);
    }
    fscanf(fp, "%d", &numpoints[onum]);
    for(int i = 0; i<numpoints[onum]; i++){
      fscanf(fp, "%lf %lf", &x[onum][i], &y[onum][i]); 
    }
    fscanf(fp, "%d", &numpolys[onum]);

    for(int j=0; j<numpolys[onum]; j++){ 
      fscanf(fp, "%d", &psize[onum][j]); 
      for (int k=0; k<psize[onum][j]; k++){ 
	fscanf(fp, "%d", &con[onum][j][k]); 
      }
    }
    for(int l = 0; l<numpolys[onum]; l++){
      fscanf(fp, "%lf %lf %lf", &red[onum][l], &green[onum][l], &blue[onum][l]);
    }


    center_and_scale(onum);
    M2d_make_translation(translated, xoriginadjust[onum], yoriginadjust[onum]);
    M2d_make_scaling(scaled, scalefactor[onum], scalefactor[onum]);
    M2d_mat_mult(m, scaled, translated);
    M2d_mat_mult_points(x[onum], y[onum], m, x[onum], y[onum], numpoints[onum]);
    M2d_make_translation(v, 300, 300);
    M2d_mat_mult_points(x[onum], y[onum], v, x[onum], y[onum], numpoints[onum]);
   
  }  // end for onum

  double a[3][3],spin[3][3],c[3][3] ;
  double angle = M_PI/180 ;
  
  q = '0' ;
  while(q != 'q'){
    onum = q - 48 ;


    M2d_make_translation(a,-300, -300);
    M2d_mat_mult_points(x[onum], y[onum], a, x[onum], y[onum], numpoints[onum]);
    M2d_make_rotation_radians(spin, angle);
    M2d_mat_mult_points(x[onum], y[onum], spin, x[onum], y[onum], numpoints[onum]);
    M2d_make_translation(c, 300, 300);
    M2d_mat_mult_points(x[onum], y[onum], c, x[onum], y[onum], numpoints[onum]);

    draw_polygon(onum);    

    q = G_wait_key();
    G_rgb(1,1,1);
    G_clear();
    
  }


}



  
