#include "GaussQuad.hh"

#include <cmath>
#include <iostream>

using namespace std; 

vector<vector<double>> createRef1dQuad(double x1, double x2, int order)
{
   double z, x, w, p1, p2, p3, pp;
   const double pi = 3.141592653589793;

   vector<vector<double> > quad(order);
   for(int i = 0; i < order; i++) quad[i].resize(2);

   // Newton-Raphson
   int m = (order + 1) / 2;
   double tol = 1.0e-12;
   double xm  = 0.5*(x2 + x1);
   double xl  = 0.5*(x2 - x1);
   for(int i = 1; i <= m; i++) {
       z = cos(pi * (i - 0.25) / (order + 0.5));
       do {
           p1 = 1.0;
           p2 = 0.0;
           for(int j = 0; j < order; j++) {
               p3 = p2;
               p2 = p1;
               p1 = ((2.0*j + 1.)*z*p2 - j*p3) / (j + 1.);
           }
           pp = order*(z*p1 - p2) / ((z*z) - 1.);
           z = z - p1 / pp;
       } while (abs(p1 / pp) > tol);
       quad[i - 1][0] = xm - xl*z;
       quad[order - i][0] = xm + xl*z;
       w = 2.0*xl / ((1.0 - (z*z)) * pp * pp);
       quad[i - 1][1] = w;
       quad[order - i][1] = w;
   }
   return quad;
}