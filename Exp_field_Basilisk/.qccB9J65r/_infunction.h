
  coord v = {(xp - x)/Delta, (yp - y)/Delta}, np;
  {
#line 712

    np.x = - interp (point, v, n.x);
#line 712

    np.y = - interp (point, v, n.y);}
  glNormal3d (np.x, np.y, 1.);
  glvertex3d (view, xp, yp, zp);
