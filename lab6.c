//painter's algorithm
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
int numobjects, poly_size, q, s, sign, action,k;
double xp[15000], yp[15000], x_prime[10][15000], y_prime[10][15000];
double half_angle = 40;
double V[4][4], translated[4][4], rotated[4][4];
double cs, sn, avx, avy, avz;
double a[3], b[3];
double color[9][3];
typedef struct {int objectnum; int polynum; double distance;} THING;
int n = 0;

THING object[150000];

void init_colors(){
//create different colors
	color[0][0] = 0.7; color[0][1]=0.7; color[0][2]= 1; 	
	color[1][0] = 1; color[1][1] = 1; color[1][2] = 1;
	color[2][0] = 0.9; color[2][1] = 0.1; color[2][2] = 0.5;
	color[3][0] = 0; color[3][1] = 0.7; color[3][2] = 0.4;
	color[4][0] = 0.2; color[4][1] = 0.5; color[4][2] = 0.8;
	color[5][0] = 0.8; color[5][1] = 0.6; color[5][2] = 0.1;
	color[6][0] = 0.9; color[6][1] = 0.9; color[6][2] = 0.9;
	color[7][0] = 0.9; color[7][1] = 0.5; color[7][2] = 0.2;
	color[8][0] = 0; color[8][1] = 0.8; color[8][2] = 0.5;
} 

void move_to_screen(int key){
  double H, scale_factor;
 	H = tan(half_angle);
 	scale_factor = 300/H;
 	for(int i =0; i<numpoints[key]; i++){
   		y_prime[key][i] = (y[key][i])/(z[key][i]);
   		x_prime[key][i] = (x[key][i])/(z[key][i]);

   		display_y[key][i] = (y_prime[key][i]*scale_factor)+300;
   		display_x[key][i] = (x_prime[key][i]*scale_factor)+300;
 	}
}

void draw_image(int key){
	move_to_screen(key);
	int h;
 	for(int i = 0; i<numpolys[key]; i++){
   	poly_size = psize[key][i];
   	for(int j = 0; j<poly_size; j++){
     	h = con[key][i][j] ;
     	xp[j] = display_x[key][con[key][i][j]];
     	yp[j] = display_y[key][con[key][i][j]];
   	}
   	G_rgb(color[key][0], color[key][1], color[key][2]);
   	G_polygon(xp,yp,poly_size);
 	}
}


void init_array(){
	//assign 
	for(int i = 0; i<numobjects; i++){
		for(int j = 0; j<numpolys[i]; j++){
			object[n].objectnum = i;
			object[n].polynum = j;
			object[n].distance = z[i][con[i][j][0]];
			n++;
		}
	}
}
void print_array(){
	int ii;
	for(ii = 0; ii<n; ii++){
		printf("%d %d %lf\n", object[ii].objectnum, object[ii].polynum, object[ii].distance);
	}
	printf("\n");
}

int compare(const void *p, const void *q){
	THING *a, *b;

	a = (THING*)p;
	b = (THING*)q;

	if(((*a).distance) < ((*b).distance)) return 1;
	else if(((*a).distance) > ((*b).distance)) return -1;
	else return 0;
}

void draw_object(){
	qsort(object, n, sizeof(THING), compare);
	int h;
	for(int i = 0; i < n; i++){
		double tempx[10000], tempy[10000];
		for(int j = 0; j < psize[object[i].objectnum][object[i].polynum]; j++){
			h = con[object[i].objectnum][object[i].polynum][j];
			y_prime[onum][j] = (y[onum][j])/(z[onum][j]);
			x_prime[onum][j] = (x[onum][j])/(z[onum][j]);

			display_y[onum][j] = (y_prime[onum][j]*scale_factor)+300;
			display_x[onum][j] = (x_prime[onum][i]*scale_factor)+300;
			
			tempx[j] = display_x[object[i].objectnum][h];
			tempy[j] = display_y[object[i].objectnum][h];
		}
		double red = color[object[i].objectnum][0]; double green = color[object[i].objectnum][1]; double blue = color[object[i].objectnum][2];
		G_rgb(red, green, blue);
		G_fill_polygon(tempx, tempy, psize[object[i].objectnum][object[i].polynum]);
		//G_rgb(0,0,0);
		//G_polygon(tempx, tempy, psize[object[i].objectnum][object[i].polynum]);
	}

}

int main(int argc, char **argv){
 	FILE *fp;
 	int key = 0;
    numobjects = argc-1;
    init_colors();

 	G_init_graphics(600,600);
 	G_rgb(0,0,0);
 	G_clear();

 	for(key = 0; key<numobjects; key++){
   		fp = fopen(argv[key+1], "r");
   		if(fp == NULL){
        	printf("Broken File\n");
     	  	exit(0);
   	  	}
   		fscanf(fp, "%d", &numpoints[key]);
   		for(int i = 0; i<numpoints[key]; i++){
     		fscanf(fp, "%lf %lf %lf", &x[key][i], &y[key][i], &z[key][i]);
   		}
   		fscanf(fp, "%d", &numpolys[key]);
   		for(int j = 0; j<numpolys[key]; j++){
     		fscanf(fp, "%d", &psize[key][j]);
     		for(int k = 0; k<psize[key][j]; k++){
			   fscanf(fp, "%d", &con[key][j][k]);
     		}
   		}
		init_array();
		print_array();
   		draw_object();
   		print_array();
 	}

 	q = '0';
 	s = '0';
 	sign = 1;
 	action = 't';
 	key = 0;

 	while(1==1){
   		M3d_make_identity(V); M3d_make_identity(translated), M3d_make_identity(rotated);
   		q = G_wait_key();
   		if(q =='q'){
     		exit(0);}
   		else if(q == 'c'){
    		sign = -sign;
      }
   		else if(q == 't'){
     		action = q;}
   		else if(q == 'r'){
     		action = q;}
   		else if(('0' <= q) && (q <= '9')){
     		k = q - '0';
     		if(k < numobjects){
				  key = k;}}
   		else if((q=='x') && (action == 't')){
        if(sign > 0){
          M3d_make_translation(V, 0.2, 0., 0.);
        }
        else{
          M3d_make_translation(V, -0.2, 0., 0.);
        }
   		}
   		else if((q == 'y') && (action == 't')){
        if(sign > 0){
          M3d_make_translation(V, 0, 0.2, 0.);
        }
        else{
          M3d_make_translation(V, 0, -0.2, 0.);
        }
   		}
   		else if((q == 'z') && (action == 't')){
     		if(sign > 0){
          M3d_make_translation(V, 0., 0., 0.2);

        }
        else{
          M3d_make_translation(V, 0., 0., -0.2);
        }
   		}
   		else if((q == 'x') && (action == 'r')){
        avy = 0; avz = 0;
        for (int i = 0; i<numpoints[key]; i++){
          avy += y[key][i];
          avz += z[key][i];
        }
        avy = avy/numpoints[key];
        avz = avz/numpoints[key];
        M3d_make_translation(translated, 0, -avy, -avz);
        M3d_mat_mult_points(x[key], y[key], z[key], translated, x[key], y[key], z[key], numpoints[key]+1);
        if(sign > 0){
          cs = cos(M_PI/180);
          sn = sin(M_PI/180);
          M3d_make_x_rotation_cs(rotated, cs, sn);
          M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key]+1);

        }
        else{
          cs = cos(-M_PI/180);
          sn = sin(-M_PI/180);
          M3d_make_x_rotation_cs(rotated, cs, sn);
          M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key]+1);

       	}
        M3d_make_translation(V, 0, avy, avz);

        }
   		else if((q == 'y') && (action == 'r')){
        avx = 0; avz = 0;
        for(int i = 0; i<numpoints[key]; i++){
          avx += x[key][i];
          avz += z[key][i];
        }         
        avx = avx/numpoints[key];
        avz = avz/numpoints[key];
        M3d_make_translation(translated, -avx, 0, -avz);
        M3d_mat_mult_points(x[key], y[key], z[key], translated, x[key], y[key], z[key], numpoints[key]+1);
        if(sign>0){
          cs = cos(M_PI/180);
          sn = sin(M_PI/180);
          M3d_make_y_rotation_cs(rotated, cs, sn);
          M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key]+1);
        }
        else{
          cs = cos(-M_PI/180);
          sn = sin(-M_PI/180);
          M3d_make_y_rotation_cs(rotated, cs, sn);
          M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key]+1);
        }
        M3d_make_translation(V, avx, 0, avz);
   		}
   		else if((q == 'z') && (action == 'r')){
        avx = 0; avy = 0;
        for(int i = 0; i<numpoints[key]; i++){
          avx += x[key][i];
          avy += y[key][i];
        }
        avx = avx/numpoints[key];
        avy = avy/numpoints[key];
        M3d_make_translation(translated, -avx, -avy, 0);
        M3d_mat_mult_points(x[key], y[key], z[key], translated, x[key], y[key], z[key], numpoints[key]+1);
        if(sign>0){
          cs = cos(M_PI/180);
          sn = sin(M_PI/180);
          M3d_make_z_rotation_cs(rotated, cs, sn);
          M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key]+1);
        }
        else{
          cs = cos(-M_PI/180);
          sn = sin(-M_PI/180);
          M3d_make_z_rotation_cs(rotated, cs, sn);
          M3d_mat_mult_points(x[key], y[key], z[key], rotated, x[key], y[key], z[key], numpoints[key]+1);
        }
        M3d_make_translation(V, avx, avy, 0);
   		}
   		else{
     		printf("no action\n");
   		}
   		M3d_mat_mult_points(x[key], y[key], z[key], V, x[key], y[key], z[key], numpoints[key]+1);
   		G_rgb(0,0,0);
   		G_clear();
   		draw_object();
   		
 	}
}