#include "LinearSolver.hh"

#include <iostream> 
#include <cmath> 

/* Gauss elimintation */ 
int gauss_elim(int dim, vector<vector<double> > &a, vector<double> &x, vector <double> &b) 
{
	vector<int> mGaussLastRow; 
	vector<int> mGaussLastCol; 
    int ierr = 0;
    if(mGaussLastCol.size() != dim) {
       mGaussLastCol.clear();
       mGaussLastRow.clear();
       mGaussLastCol.resize(dim, dim-1);
       mGaussLastRow.resize(dim, dim-1);
       // if(mGaussPatternOn) {
       //    setGaussPattern(dim, a);
       // }
    }

    // We first do dim reduction steps.
    for(int k = 0; k < dim; k++) {
        double aSubK = a[k][k];
        double bSubK = b[k];
        if(aSubK == 0.0) {              /* check */
            //fprintf(stderr, "gauss: a[%d][%d] = 0\n", k, k);
            //cout << "gauss: a[" << k << "][" << k << "] = 0\n";
            //exit(1);
            return -1;
        }

        for(int i = k+1; i <= mGaussLastCol[k]; i++) {
            double xtemp = a[i][k];
            if(xtemp == 0.0) continue;
            a[i][k] = 0.0;                      /* eq (4.1a) */
            xtemp /= aSubK;
            for(int j = k+1; j <= mGaussLastRow[k]; j++) {
                a[i][j] -= a[k][j]*xtemp;
            }
            b[i] -= bSubK*xtemp;              /* eq (4.1c) */
        }
    }
 
    // Now we perform dim back substitutions.
    for(int i = dim-1; i >= 0; i--) {
        b[i] /= a[i][i];
        for(int k = i - 1; k >= 0; k -= 1) {
            b[k] -= a[k][i] * b[i];
        }
    }

    // Return solution in x and b vectors
    for(int i = 0; i < dim; i++) x[i] = b[i];
    return ierr;
}