/*
  An implementation of SVD from Numerical Recipes in C and Mike Erhdmann's lectures
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Maxfiles.h" 			// Includes .max files
#include <MaxSLiCInterface.h>	// Simple Live CPU interface

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1 = (a),iminarg2 = (b),(iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

int svdcmp(double **a, int nRows, int nCols, double *w, double **v);

// prints an arbitrary size matrix to the standard output
void printMatrix(double **a, int rows, int cols);
void printMatrix(double **a, int rows, int cols) {
    int i,j;

    for(i=0;i<rows;i++) {
        for(j=0;j<cols;j++) {
            printf("%.4lf ",a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// prints an arbitrary size vector to the standard output
void printVector(double *v, int size);
void printVector(double *v, int size) {
    int i;

    for(i=0;i<size;i++) {
        printf("%.4lf ",v[i]);
    }
    printf("\n\n");
}

/*
  Modified from Numerical Recipes in C
  Given a matrix a[nRows][nCols], svdcmp() computes its singular value
  decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
  returns.  The diagonal matrix W is output as a vector w[nCols].
  V (not V transpose) is output as the matrix V[nCols][nCols].
*/
int svdcmp(double **a, int nRows, int nCols, double *w, double **v) {
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

    rv1 = malloc(sizeof(double)*nCols);
    if(rv1 == NULL) {
        printf("svdcmp(): Unable to allocate vector\n");
        return(-1);
    }

    g = scale = anorm = 0.0;
    for(i=0;i<nCols;i++) {
        l = i+1;
        rv1[i] = scale*g;
        g = s = scale = 0.0;
        if(i < nRows) {
            for(k=i;k<nRows;k++) scale += fabs(a[k][i]);
            if(scale) {
                for(k=i;k<nRows;k++) {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s),f);
                h = f * g - s;
                a[i][i] = f - g;
                for(j=l;j<nCols;j++) {
                    for(s=0.0,k=i;k<nRows;k++) s += a[k][i] * a[k][j];
                    f = s / h;
                    for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
                }
                for(k=i;k<nRows;k++) a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if(i < nRows && i != nCols-1) {
            for(k=l;k<nCols;k++) scale += fabs(a[i][k]);
            if(scale)  {
                for(k=l;k<nCols;k++) {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = - SIGN(sqrt(s),f);
                h = f * g - s;
                a[i][l] = f - g;
                for(k=l;k<nCols;k++) rv1[k] = a[i][k] / h;
                for(j=l;j<nRows;j++) {
                    for(s=0.0,k=l;k<nCols;k++) s += a[j][k] * a[i][k];
                    for(k=l;k<nCols;k++) a[j][k] += s * rv1[k];
                }
                for(k=l;k<nCols;k++) a[i][k] *= scale;
            }
        }
        anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));

        printf(".");
        fflush(stdout);
    }

    for(i=nCols-1;i>=0;i--) {
        if(i < nCols-1) {
            if(g) {
                for(j=l;j<nCols;j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for(j=l;j<nCols;j++) {
                    for(s=0.0,k=l;k<nCols;k++) s += a[i][k] * v[k][j];
                    for(k=l;k<nCols;k++) v[k][j] += s * v[k][i];
                }
            }
            for(j=l;j<nCols;j++) v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
        printf(".");
        fflush(stdout);
    }

    for(i=IMIN(nRows,nCols) - 1;i >= 0;i--) {
        l = i + 1;
        g = w[i];
        for(j=l;j<nCols;j++) a[i][j] = 0.0;
        if(g) {
            g = 1.0 / g;
            for(j=l;j<nCols;j++) {
                for(s=0.0,k=l;k<nRows;k++) s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
            }
            for(j=i;j<nRows;j++) a[j][i] *= g;
        }
        else
            for(j=i;j<nRows;j++) a[j][i] = 0.0;
        ++a[i][i];
        printf(".");
        fflush(stdout);
    }

    for(k=nCols-1;k>=0;k--) {
        for(its=0;its<30;its++) {
            flag = 1;
            for(l=k;l>=0;l--) {
                nm = l-1;
                if((fabs(rv1[l]) + anorm) == anorm) {
                    flag =  0;
                    break;
                }
                if((fabs(w[nm]) + anorm) == anorm) break;
            }
            if(flag) {
                c = 0.0;
                s = 1.0;
                for(i=l;i<=k;i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if((fabs(f) + anorm) == anorm) break;
                    g = w[i];
                    double* h1;
                    double* f1 = &f;
                    double* g1 = &g;
                    svd(3, f1, g1, h1);
                    h = *h1;
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for(j=0;j<nRows;j++) {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];
            if(l == k) {
                if(z < 0.0) {
                    w[k] = -z;
                    for(j=0;j<nCols;j++) v[j][k] = -v[j][k];
                }
                break;
            }
            if(its == 29) printf("no convergence in 30 svdcmp iterations\n");
            x = w[l];
            nm = k-1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            double* f1 = &f;
            double nums = 1.0;
            double* nums1 = &nums;
            double* g1;
            svd(3, f1, nums1, g1);
            g = *g1;
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f))) - h)) / x;
            c = s = 1.0;
            for(j=l;j<=nm;j++) {
                i = j+1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                double* z1;
                f1 = &f;
                double* h1 = &h;
                svd(3, f1, h1, z1);
                z = *z1;
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for(jj=0;jj<nCols;jj++) {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                f1 = &f;
                h1 = &h;
                svd(3, f1, h1, z1);
                z = *z1;
                w[j] = z;
                if(z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for(jj=0;jj < nRows;jj++) {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
        printf(".");
        fflush(stdout);
    }
    printf("\n");

    free(rv1);

    return(0);
}

int main(){

    int n1=45, n2=27;

    // Create matrix
    double *M;
    for(size_t i=0;i<n1*n2;i++) M[i] = (i % 10);

    int m = n2;
    int n = n1;
    int k = (m<n?m:n);
    double tU[m*k];
    double tS[k];
    double tVT[k*n];

    { // Compute SVD
        int INFO=0;
        char JOBU  = 'S';
        char JOBVT = 'S';
        int wssize = 3*(m<n?m:n)+(m>n?m:n);
        int wssize1 = 5*(m<n?m:n);
        wssize = (wssize>wssize1?wssize:wssize1);
        double* wsbuf[wssize];
        svdcmp(&M, n1, n2, &tU, &tVT);
    }

    { // Check Error
        double max_err=0;
        for(size_t i0=0;i0<m;i0++)
            for(size_t i1=0;i1<n;i1++){
                double E = M[i1*m+i0];
                for(size_t i2=0;i2<k;i2++) E -= tU[i2 * m+i0] * tS[i2] * tVT[i1*k+i2];
                if(max_err<fabs(E)) max_err=fabs(E);
            }
        printf("%lf\n", max_err);
    }

    return 0;
}
