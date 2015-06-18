/* Programa en C que resuelve el problema de una particula 
 * en un pozo de potencial usando B-splines.
 * El potencial es de la forma V(r) = -lambda si r1<r<r2 y cero fuera.
 * El programa guarda los autovalores en funcion de lambda.
 * Usa el metodo variacional de Rayleight-Ritz usando como base del 
 * espacio los B-splines, como esta no es una base ortonormal el 
 * problema de autovalores queda de la siguiente forma:
 *
 *              H |psi> = e S |psi>
 *
 * donde H es la matriz del hamiltoniano y S es la matriz de solapamiento
 * de los B-splines.
 *
 * Usa la funcion gauleg() del Numerical Recipies C para calcular los 
 * puntos de evaluacion y los pesos para la cuadratura.
 * La funcion KNOTS_PESOS() esta hecha por mi guiandome de la version de 
 * fortran, calcula los knots de los B-splines con una distribucion 
 * uniforme solamente, en el caso de querer otra distribucion es solo 
 * cuestion de modificar el codigo para usar la distribucion que uno quiera.
 * Para el calculo de los autovalores usa la funcion dsygvx_() de lapack.
 * Las funciones para evaluar los B-splines y las derivadas son bsplvb() y 
 * bder() respectivamente, son versiones hechas por mi a partir de las 
 * versiones de fortran.
 *
 * Una observacion importante: este codigo anda solo para l>=25, 
 * para valores mas chicos hace cosas raras no se bien por que.
 *
 * Al final del archivo se encuentran las diferencias con la version 
 * anterior, leerlas para tener en cuenta.
*/

#include <stdio.h> // prinft
#include <stdlib.h> // malloc y free
#include <malloc.h>
#include <math.h> // posibles operaciones matematicas
#include <assert.h>
#include <omp.h> //omp_get_wtime()

// defino algunas constantes para el programa //
#define EPS 3.0e-14  // EPS precision relativa para gauleg //

// defino algunos parametros del programa //

#ifndef Rmin
#define Rmin 0.0 // R minimo donde empieza el intervalo para la integracion //
#endif

#ifndef Rmax
#define Rmax 50.0 // R maximo donde termina el intervalo para la integracion //
#endif

#ifndef l
#define l 50 // numero de intervalos en el que divido al intervalo [Rmin, Rmax] //
#endif

#ifndef kord
#define kord 5 // orden de los B-splines, el grado es kord-1 //
#endif

#ifndef r1
#define r1 5.0 // radio interio del pozo //
#endif

#ifndef r2
#define r2 10.0 // radio exterior del pozo //
#endif

#ifndef me
#define me 1.0 // masa de la particula //
#endif

#ifndef intg
#define intg 50 // grado de integracion por cuadratura //
#endif

#ifndef nev
#define nev 15 // numero de autovalores que vamos a calcular //
#endif

#ifndef lmax
#define lmax 0 // momento angular que vamos a usar //
#endif

#ifndef lambda_in
#define lambda_in 0.0 // lambda inicial para el pozo //
#endif

#ifndef lambda_fin
#define lambda_fin 20.0 // lambda final para el pozo //
#endif

#ifndef numero_puntos_lambda
#define numero_puntos_lambda 200 // numero de puntos para calcular //
#endif

// escribo las funciones del programa //
int dsygvx_(int *itype, char *jobz, char *range, char * uplo, 
    int *n, double *a, int *lda, double *b, int *ldb, 
    double *vl, double *vu, int *il, int *iu, double *abstol, 
    int *m, double *w, double *z__, int *ldz, double *work, 
    int *lwork, int *iwork, int *ifail, int *info);


int idx(unsigned int y, unsigned int x, unsigned int numcolumns){
    return y*numcolumns + x;
}

int cleari(unsigned int N, int * __restrict__ vec) {
    
    for(unsigned int i = 0; i<N; ++i) 
        vec[i] = 0;

    return 0;
}

int cleard(unsigned int N, double * __restrict__ vec) {
    
    for(unsigned int i = 0; i<N; ++i)
        vec[i] = 0.f;

    return 0;
}


void gauleg(double x1, double x2, double x[], double w[], int n) {
/* Given the lower and upper limits of integration x1 and x2, 
 * and given n, this routine returns arrays x[1..n] and w[1..n]
 * of length n, containing the abscissas and weights of the Gauss-
 * Legendre n-point quadrature formula.
*/
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);

    for (i = 1; i<=m; i++) {
        z = cos(3.141592654*(i-0.25)/(n+0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 1;j<=n;j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp = n*(z*p1-p2)/(z*z-1.0);
            z1 = z;
            z = z1-p1/pp;
        } while (fabs(z-z1) > EPS);
        x[i] = xm-xl*z;
        x[n+1-i] = xm+xl*z;
        w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i] = w[i];
    }
}

int KNOTS_PESOS(unsigned int nk, int * __restrict__ k, 
        double * __restrict__ t, double * __restrict__ x, 
        double * __restrict__ w){

    double dr, ri, rf;
    double vx[intg+1], vw[intg+1];
//  double * vx = calloc(intg+1, sizeof(double));
//  double * vw = calloc(intg+1, sizeof(double));

    dr = (Rmax-Rmin)/l;


    for(unsigned int i = 0; i<l; ++i) {
        ri = Rmin + i*dr;
        rf = ri + dr;
        gauleg(ri, rf, vx, vw, intg);

        for(unsigned int j = 0; j<intg; j += 1 ) {
            x[idx(i, j, intg)]  = vx[j+1] ;//  x[idx(i, j+1, intg)] = vx[j+2];

            w[idx(i, j, intg)]  = vw[j+1] ;//  w[idx(i, j+1, intg)] = vw[j+2];
        }
    }

    /// en esta parte controlar bien el tema de los indices de los vectores //

    for(unsigned int i = 0; i<kord; ++i) {
        t[i] = Rmin;
        k[i] = 0;
    }

    for(unsigned int i = kord; i<kord+l; ++i) {
        double pos = (i - kord + 1);
        t[i] = Rmin + dr * pos;
        k[i] = pos;
    }

    for(unsigned int i = kord+l; i<nk; ++i) {
        t[i] = Rmin + dr * l;;
        k[i] = l;
    }

    return 0;
}

int last_ind = -1;
int bsplvb(double * __restrict__ t, unsigned int jhigh, double rr, int left, double * __restrict__ biatx, int ind ) {

    if(ind == last_ind){
        return 0;  
    }
    last_ind = ind;
    unsigned int j, jp1, JMAX;
    double saved, term;
    double * deltar, * deltal;

    //if(1 != index) printf("index no es igual a 1");

    JMAX = 100;

    deltar = (double *) calloc(JMAX, sizeof(double));
    deltal = (double *) calloc(JMAX, sizeof(double));

    biatx[0] = 1.0;

    for(j=0; j<jhigh-1; ++j) {
        jp1 = j+1;

        deltar[j] = t[left+j+1]-rr;
        deltal[j] = rr-t[left-j];

        saved = 0.0;
        for(unsigned int i = 0; i<j+1; ++i) {
            term = biatx[i]/(deltar[i]+deltal[jp1-i-1]);
            biatx[i] = saved + deltar[i]*term;
            saved = deltal[jp1-i-1]*term;
        }

        biatx[jp1] = saved;
    }

    free(deltar);
    free(deltal);

    return 0; 
}

double bder(double rr, double * __restrict__ t, unsigned int korder, 
     unsigned int np, unsigned int indexm, unsigned int left,
     double * __restrict__ Sp, double dm, int ind ) {

    unsigned int i;

    if(t[0]<rr && rr<t[np-1]) {

        if(abs(rr-t[np-1])<1.e-10) {
            if(indexm==np-korder) {
                dm = (korder-1)/(t[np-1]-t[np-1-korder]);
            }
            else {
                dm = -(korder-1)/(t[np-1]-t[np-1-korder]);
            }
        }

        bsplvb(t, korder-1, rr, left, Sp, ind);

        if(indexm-left+korder>=1 || indexm-left+korder<=korder) {
            i = indexm-left+korder;
            if(1==i) {
                dm = (korder-1)*(-Sp[i-1]/(t[indexm+korder]-t[indexm+1]));
            }
            else if(korder==i) {
                dm = (korder-1)*(Sp[i-1-1]/(t[indexm+korder-1]-t[indexm]));
            }
            else {
                dm = (korder-1)*(Sp[i-1-1]/(t[indexm+korder-1]-t[indexm])
                    - Sp[i-1]/(t[indexm+korder]-t[indexm+1]));
            }
        }

    }

    return dm;
}

void calculo_matrices(unsigned int nk, unsigned int nb,
              int * __restrict__ k, double * __restrict__ t,
              double * __restrict__ x, double * __restrict__ w,
              double * __restrict__ s, double * __restrict__ v0,
              double * __restrict__ ke) {

    double ma, rr, _rr2;
    double * Sp;
    FILE * file;

    ma = 0.5*lmax*(lmax+1);

    Sp = (double *) calloc(kord, sizeof(double));
    // ojo con los limites en los for's //
    for(unsigned int i = kord-1; i<kord+l-1; ++i) {
        for(unsigned int j = 0; j<intg; ++j) {
            rr = x[idx(k[i], j, intg)];
            _rr2= 1.0/(rr*rr);

            bsplvb(t, kord, rr, i, Sp, idx(k[i], j, intg));

            for(unsigned int m = kord-1 == i, im = i-kord + m ; m<kord && im<nb; ++m, ++im) {

                for(unsigned int n = kord-1 == i, in = i-kord+n; n<kord && in<nb; ++n, ++in) {

                    s[idx(im, in, nb)] = s[idx(im, in, nb)] + Sp[m]*Sp[n]*w[idx(k[i], j, intg)];

                    ke[idx(im, in, nb)] = ke[idx(im, in, nb)] + ma*Sp[m]*Sp[n]*w[idx(k[i], j, intg)]*_rr2;

                    if(r1<rr && rr<r2) v0[idx(im, in, nb)] = v0[idx(im, in, nb)] + Sp[m]*Sp[n]*w[idx(k[i] , j, intg)];
                }
            }   
        }
    }

    for(unsigned int m = 1; m<=kord-1; ++m) {

        for(unsigned int n = m; n<=kord-1; ++n) {

            for(unsigned int j = 0; j<intg; ++j) {
                
                double bm = 0, bn = 0;

                rr = x[idx(k[kord-1], j, intg)];

                bm = bder(rr, t, kord, nk, m, kord-1, Sp, bm, idx(k[kord-1], j, intg));
                bn = bder(rr, t, kord, nk, n, kord-1, Sp, bn, idx(k[kord-1], j, intg));

                ke[idx(m-1, n-1, nb)] = ke[idx(m-1, n-1, nb)] + 0.5*w[idx(k[kord-1], j, intg)]*bm*bn/me;

            }
        }
    }

    for(unsigned int i = kord; i<kord+l-1; ++i) {
        // ojo con los indices en esta parte //
        for(unsigned int m = i-kord+1; m<=i; ++m) {
//          if(m>0 && m<nb+1) {
            if(m>0 && m<nb) {
                for(unsigned int n = m; n<=i; ++n) {
//                  if(n<nb+1) {
                    if(n<nb) {
                        for(unsigned int j = 0; j<intg; ++j) {

                            double bm = 0, bn = 0;

                            rr = x[idx(k[i], j, intg)];

                            bm = bder(rr, t, kord, nk, m, i, Sp, bm, idx(k[i], j, intg));
                            bn = bder(rr, t, kord, nk, n, i, Sp, bn, idx(k[i], j, intg));

                            ke[idx(m-1, n-1, nb)] = ke[idx(m-1, n-1, nb)] + 0.5*w[idx(k[i], j, intg)]*bm*bn/me;

                        }
                    }
                }
            }
        }
    }

    for(unsigned int i = 0; i<nb; ++i) {
        for(unsigned int j = i+1; j<nb; ++j) {
            s[idx(j, i, nb)] = s[idx(i, j, nb)];
            v0[idx(j, i, nb)] = v0[idx(i, j, nb)];
            ke[idx(j, i, nb)] = ke[idx(i, j, nb)];
        }
    }

    file = fopen("matrices.dat", "w");
    for(unsigned int i = 0; i<nb; ++i)
        for(unsigned int j=0; j<nb; ++j)
            fprintf(file, "%i   %i  %.12f   %.12f   %.12f\n", i, j, s[idx(i,j,nb)], v0[idx(i,j,nb)], ke[idx(i,j,nb)]);

    fclose(file);

    free(Sp);
}

/*
void eigenvalues(int n, int m, double * __restrict__ a,
         double * __restrict__ b, double * __restrict__ w, 
         double * __restrict__ z) {

    int itype, lda, ldb, ldz;
    int il, iu, lwork, info;
    char jobz, range, uplo;
    double abstol, vl, vu;
    double * work;
    int * iwork, * ifail;

    // le doy los valores correspondientes a las distintas variables //
    itype = 1; lda = n; ldb = n; ldz = n;
    vl = 0.0; vu = 0.0;
    jobz = 'V'; range = 'I'; uplo = 'U';
    il = 1; iu = m; lwork = 9*n;

    // le doy memoria a las matrices que neceista dsygvx //
    work = (double *) malloc(lwork*sizeof(double));
    iwork = (int *) malloc(5*n*sizeof(int));
    ifail = (int *) malloc(n*sizeof(int));

    dsygvx_( &itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, 
         &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, 
         &lwork, iwork, ifail, &info);

    free(work); 
    free(iwork); 
    free(ifail);
} 

void hamiltoniano_autovalores(unsigned int nb, double * __restrict__ s,
                  double * __restrict__ v0, double * __restrict__ ke,
                  FILE * archivo) {

    double * h, * auval, * auvec, * s_copy;
    double lambda, delta;
    
    // doy memoria a las matrices y vectores //
    h = (double *) calloc( nb*nb, sizeof(double));
    s_copy = (double *) calloc( nb*nb, sizeof(double));
    auval = (double *) calloc( nev, sizeof(double));
    auvec = (double *) calloc( nev*nb, sizeof(double));
    
    lambda = 1.0;
    delta = (lambda_fin-lambda_in)/numero_puntos_lambda;

    // abro el archivo para guardar los datos //
    fprintf(archivo, "# Autovalores calculados\n");
    fprintf(archivo, "# Lambda  auval[0]   auval[1] ....\n");

    for(unsigned int j = 0; j<=numero_puntos_lambda; ++j) {

        lambda = lambda_in + delta*j;
    
        for(unsigned int m = 0; m<nb; ++m) {
            for(unsigned int n = 0; n<nb; ++n){
                h[idx(n, m, nb)] = ke[idx(n, m, nb)] - lambda*v0[idx(n, m, nb)];
                s_copy[idx(n, m, nb)] = s[idx(n, m, nb)];
            }
        }

        eigenvalues( nb, nev, h, s_copy, auval, auvec );

        fprintf(archivo, "%.5f   ", lambda);
        for(unsigned int i = 0; i<nev; ++i) {
            fprintf(archivo, "%.15f   ", auval[i]);
        }
        fprintf(archivo, "\n");

    }

    free(h); free(s_copy); free(auval); free(auvec);
}*/

int main(void) {

    // defino algunas variables que voy a usar //
    unsigned int nk, nb;
    int *k;
    double *t, *x, *w;
    double *s, *v0, *ke;
    double t_in, t_fin, t_n;
    FILE * archivo;
    
    nk = l+2*kord-1; // numero de knots //
    nb = nk-kord-2; // tamaño de la base //

    // controlo algunos parametros //
    assert(intg>kord);
    assert(nev>0);

    archivo = fopen("autovalores_pozo-1.dat", "w");
    // imprimo los parametros //
    fprintf(archivo, "# Rmin=%.12f y Rmax=%.12f\n", Rmin, Rmax);
    fprintf(archivo, "# Numero de intervalos l=%i\n", l);
    fprintf(archivo, "# Orden los B-splines kord=%i\n", kord);
    fprintf(archivo, "# Radios del pozo r1=%.12f y r2=%.12f\n", r1, r2);
    fprintf(archivo, "# Masa de la particula me=%.12f\n", me);
    fprintf(archivo, "# Grado de integracion de la cuadratura intg=%i\n", intg);
    fprintf(archivo, "# Numero de autovalores nev=%i\n", nev);
    fprintf(archivo, "# Momento angular que usamos lmax=%i\n", lmax);
    fprintf(archivo, "# Numero de knots nk=%i\n", nk);
    fprintf(archivo, "# Tamaño de la base nb=%i\n", nb);
    fprintf(archivo, "# Valores inicial y final del pozo, lambda_in=%.12f, lambda_fin=%.12f\n", lambda_in, lambda_fin);
    fprintf(archivo, "# Numero de puntos lambda = %i\n", numero_puntos_lambda);

    // doy memoria a las matrices que voy a necesitar //
    k = (int *) malloc( nk*sizeof(int));
    t = (double *) malloc( nk*sizeof(double));
    x = (double *) malloc( l*intg*sizeof(double));
    w = (double *) malloc( l*intg*sizeof(double));
    s = (double *) malloc( nb*nb*sizeof(double));
    v0 = (double *) malloc( nb*nb*sizeof(double));
    ke = (double *) malloc( nb*nb*sizeof(double));

    
    cleari(nk, k);
    cleard(nk, t);
    cleard(l*intg, x);
    cleard(l*intg, w);
    cleard(nb*nb, s);
    cleard(nb*nb, v0);
    cleard(nb*nb, ke);

    t_in = omp_get_wtime();
    // primero calculos los knost y los pesos para hacer la cuadratura //
    KNOTS_PESOS(nk, k, t, x, w);

    // calculo las matrices que necesito para resolver el problema //
    calculo_matrices(nk, nb, k, t, x, w, s, v0, ke);
    t_fin = omp_get_wtime();

    // armo el hamiltoniano y calculo los autovalores //
//  hamiltoniano_autovalores(nb, s, v0, ke, archivo);
    
    t_n = (t_fin-t_in)/(nb*nb*l*intg),
    printf("%i     %i     %i     %.12f  %.12f\n", l, nk, nb, t_fin-t_in, t_n); 
    FILE *file;
    file = fopen("matrices.dat", "w");
    for(unsigned int i = 0; i<nb; ++i)
        for(unsigned int j=0; j<nb; ++j)
            fprintf(file, "%i\t%i\t%.12f\t%.12f\t%.12f\n", i, j, s[idx(i,j,nb)], v0[idx(i,j,nb)], ke[idx(i,j,nb)]);
    fclose(file);
    
    // libero la memoria //
    fclose(archivo);
    free(ke);
    free(v0);
    free(s);
    free(w);
    free(x);
    free(t);
    free(k);
    return 0;
}

