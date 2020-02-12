#include <FPT.h>

int main()
{
  double swidth, sheight ;
  double p[2];
  double xcenter, ycenter, xedge, yedge, q;
  int i;
  double x,y;
  double radius;
  double angle;

  swidth = 600; sheight = 600;
  G_init_graphics(swidth, sheight);

  G_wait_click(p);
  xcenter=p[0]; ycenter=p[1];
  G_fill_rectangle(xcenter-2, ycenter-2, 4,4);

  G_wait_click(p);
  xedge=p[0]; yedge=p[1];
  G_fill_rectangle(xedge-2, yedge-2, 4,4);

  radius =sqrt( (xedge- xcenter)*(xedge- xcenter) + (yedge-ycenter)*(yedge-ycenter));

  G_circle(xcenter, ycenter, radius);
  for (i=0; i<16; i++){
    angle = i*2*M_PI/16;
    x=xcenter+radius*cos(angle);
    y=ycenter+radius*sin(angle);
    G_line(xcenter, ycenter, x,y);
  }

  q = G_wait_key();
}
  
  
  
