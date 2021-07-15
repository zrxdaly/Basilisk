void eddyviscosity(double Cs, vector u, double molv, scalar Evis){	
  double d1v1, d2v1, d3v1, d1v2, d2v2, d3v2, d1v3, d2v3, d3v3; 
  double b11, b12, b13, b22, b23, b33;
  double abeta, bbeta;	
  foreach(){
    d1v1 = (u.x[1,0,0] - u.x[-1,0,0])/2/Delta;
    d2v1 = (u.x[0,1,0] - u.x[0,-1,0])/2/Delta;
    d3v1 = (u.x[0,0,1] - u.x[0,0,-1])/2/Delta;
    d1v2 = (u.y[1,0,0] - u.y[-1,0,0])/2/Delta;
    d2v2= (u.y[0,1,0] - u.y[0,-1,0])/2/Delta;
    d3v2= (u.y[0,0,1] - u.y[0,0,-1])/2/Delta;
    d1v3= (u.z[1,0,0] - u.z[-1,0,0])/2/Delta;
    d2v3= (u.z[0,1,0] - u.z[0,-1,0])/2/Delta;
    d3v3= (u.z[0,0,1] - u.z[0,0,-1])/2/Delta;
    b11 = sq(Delta)*(d1v1*d1v1 + d2v1*d2v1 + d3v1*d3v1);
    b12 = sq(Delta)*(d1v1*d1v2 + d2v1*d2v2 + d3v1*d3v2);
    b13 = sq(Delta)*(d1v1*d1v3 + d2v1*d2v3 + d3v1*d3v3);
    b22 = sq(Delta)*(d1v2*d1v2 + d2v2*d2v2 + d3v2*d3v2);
    b23 = sq(Delta)*(d1v2*d1v3 + d2v2*d2v3 + d3v2*d3v3);
    b33 = sq(Delta)*(d1v3*d1v3 + d2v3*d2v3 + d3v3*d3v3);
    abeta = sq(d1v1) + sq(d2v1) + sq(d3v1)+
            sq(d1v2) + sq(d2v2) + sq(d3v2)+
            sq(d1v3) + sq(d2v3) + sq(d3v3);
    bbeta = b11*b22 - sq(b12) + b11*b33 - sq(b13) + b22*b33 - sq(b23);
    // Accordig to [1], Evis = 0 for abeta = 0.   
    Evis[] =  (abeta > 10E-5 && bbeta > (abeta/10E6))? 2.5*sq(Cs)*sqrt(bbeta/(abeta)) + molv: molv; 
  }	
}
