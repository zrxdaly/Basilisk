#ifndef BASILISK_HEADER_2
#define BASILISK_HEADER_2
#line 1 "./scatter2.h"
#if dimension == 3 
void glPointParameterfv(GLenum pname, const GLfloat * params);
#endif

struct _scatter {
  Particles p;    // particles
  float s, pc[3], coefs[3]; // point size, colour and distance attenuation coefs.
};

trace
void scatter (struct _scatter p){
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#else // Dimension == 3
  if (!p.coefs[0]){ // A guess:
    p.coefs[0] = 0.01;
    p.coefs[1] = 0.2;
    p.coefs[2] = 0.5;
  }
  glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, p.coefs);
#endif
  glEnable (GL_BLEND);
  glEnable (GL_POINT_SMOOTH);
  if (p.pc)
    glColor3f(p.pc[0], p.pc[1], p.pc[2]);
  if (!p.s)
    p.s = 20;
  glPointSize(p.s);

  glBegin (GL_POINTS);
  foreach_particle_in(p.p) {
#if dimension == 2    
    glvertex2d (view, x, y);
#else // dimension == 3
    glvertex3d (view, x, y, z);
#endif
  }
  glEnd();

  view->ni++; 
}
#endif
