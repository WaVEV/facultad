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
#include <string.h>
#include <omp.h> //omp_get_wtime()

// defino algunas constantes para el programa //
#define EPS 3.0e-14  // EPS precision relativa para gauleg //

// defino algunos parametros del programa //

#ifndef R_MIN
#define R_MIN 0.0 // R minimo donde empieza el intervalo para la integracion //
#endif

#ifndef R_MAX
#define R_MAX 50.0 // R maximo donde termina el intervalo para la integracion //
#endif

#ifndef L_INTERVALS
#define L_INTERVALS 510 // numero de intervalos en el que divido al intervalo [R_MIN, R_MAX] //
#endif

#ifndef KORD
#define KORD 5 // orden de los B-splines, el grado es kord-1 //
#endif

#ifndef RADIO_1
#define RADIO_1 5.0 // radio interio del pozo //
#endif

#ifndef RADIO_2
#define RADIO_2 10.0 // radio exterior del pozo //
#endif

#ifndef ME
#define ME 1.0 // masa de la particula //
#endif

#ifndef INT_G
#define INT_G 500 // grado de integracion por cuadratura //
#endif

#ifndef NEV
#define NEV 15 // numero de autovalores que vamos a calcular //
#endif

#ifndef L_MAX
#define L_MAX 0 // momento angular que vamos a usar //
#endif

#ifndef LAMBDA_IN
#define LAMBDA_IN 0.0 // lambda inicial para el pozo //
#endif

#ifndef LAMBDA_FIN
#define LAMBDA_FIN 20.0 // lambda final para el pozo //
#endif

#ifndef NUMEROS_PUNTO_LAMBDA
#define NUMEROS_PUNTO_LAMBDA 200 // numero de puntos para calcular //
#endif

#ifndef BASE_KORD
#define BASE_KORD 0
#endif

 #ifndef JMAX
 #define JMAX 100
 #endif

const unsigned int nk = L_INTERVALS+2*KORD-1;
int k[L_INTERVALS+2*KORD-1];
double  t[L_INTERVALS+2*KORD-1];
double x[L_INTERVALS*INT_G], 
        w[L_INTERVALS*INT_G];

unsigned int nb = (L_INTERVALS+2*KORD-1)-KORD-2; // tamaño de la base //
    
double  s[((L_INTERVALS+2*KORD-1)-KORD-2)*((L_INTERVALS+2*KORD-1)-KORD-2)],
        v0[((L_INTERVALS+2*KORD-1)-KORD-2)*((L_INTERVALS+2*KORD-1)-KORD-2)],
        ke[((L_INTERVALS+2*KORD-1)-KORD-2)*((L_INTERVALS+2*KORD-1)-KORD-2)];

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
    memset(vec, 0, sizeof(double) * N);
    //for(unsigned int i = 0; i<N; ++i)
    //    vec[i] = 0.f;

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

int knots_pesos(double * __restrict__ x, double * __restrict__ w){

    double dr, ri, rf;
    double vx[INT_G+1], vw[INT_G+1];
//  double * vx = calloc(INT_G+1, sizeof(double));
//  double * vw = calloc(INT_G+1, sizeof(double));

    dr = (R_MAX-R_MIN)/L_INTERVALS;


    for(unsigned int i = 0; i<L_INTERVALS; ++i) {
        ri = R_MIN + i*dr;
        rf = ri + dr;
        gauleg(ri, rf, vx, vw, INT_G);

        for(unsigned int j = 0; j<INT_G; j += 1 ) {
            x[idx(i, j, INT_G)]  = vx[j+1] ;//  x[idx(i, j+1, INT_G)] = vx[j+2];

            w[idx(i, j, INT_G)]  = vw[j+1] ;//  w[idx(i, j+1, INT_G)] = vw[j+2];
        }
    }

    return 0;
}

double ti(int i){
    double dr = (R_MAX-R_MIN)/L_INTERVALS;
    int pos = (i - KORD + 1);
    return R_MIN + (KORD - 1 < i) * dr * (pos - (pos - L_INTERVALS) * (KORD + L_INTERVALS - 1 < i));
}

int last_ind = -1;
int bsplvb(unsigned int jhigh, double rr, int left, double * __restrict__ biatx, int ind ) {

    if(ind == last_ind){
        return 0;  
    }
    //printf("%d\n", ind);
    last_ind = ind;
    unsigned int j;
    double saved, term;
    

    //if(1 != index) printf("index no es igual a 1");

    double deltar[JMAX];
    double deltal[JMAX];

    biatx[0] = 1.0;


    for(j=0; j<jhigh-1; ++j) {
        deltar[j] = ti(left+j+1) - rr;
        deltal[j] = rr-ti(left-j);

        saved = 0.0;
        for(unsigned int i = 0; i<j+1; ++i) {
            term = biatx[i]/(deltar[i]+deltal[j-i]);
            biatx[i] = saved + deltar[i]*term;
            saved = deltal[j-i]*term;
        }

        biatx[j+1] = saved;
    }

    return 0; 
}

double bder(double rr, unsigned int indexm, unsigned int left,
     double * __restrict__ Sp, double dm, int ind ) {

    unsigned int i;
    
    if(ti(0)<rr && rr<ti(nk-1)) {

        if(abs(rr-ti(nk-1))<1.e-10) {
            if(indexm==nk-KORD) {
                dm = (KORD-1)/(ti(nk-1)-ti(nk-1-KORD));
            }
            else {
                dm = -(KORD-1)/(ti(nk-1)-ti(nk-1-KORD));
            }
        }

        bsplvb(KORD-1, rr, left, Sp, ind);

        if(indexm-left+KORD>=1 || indexm-left+KORD<=KORD) {
            i = indexm-left+KORD;
            if(1==i) {
                dm = (KORD-1)*(-Sp[i-1]/(ti(indexm+KORD)-ti(indexm+1)));
            }
            else if(KORD==i) {
                dm = (KORD-1)*(Sp[i-1-1]/(ti(indexm+KORD-1)-ti(indexm)));
            }
            else {
                dm = (KORD-1)*(Sp[i-1-1]/(ti(indexm+KORD-1)-ti(indexm))
                    - Sp[i-1]/(ti(indexm+KORD)-ti(indexm+1)));
            }
        }

    }

    return dm;
}

void calculo_matrices(double * __restrict__ x, double * __restrict__ w,
              double * __restrict__ s, double * __restrict__ v0,
              double * __restrict__ ke) {

    double ma, rr, _rr2;
    double Sp[KORD];

    ma = 0.5*L_MAX*(L_MAX+1);

    for(unsigned int i = KORD-1, basek=0; i<KORD+L_INTERVALS-1; ++i, ++basek) {

        //int  = i - (KORD-1);
        for(unsigned int j = 0; j<INT_G; ++j) {
            rr = x[idx(basek, j, INT_G)];
            _rr2= 1.0/(rr*rr);
            
            bsplvb(KORD, rr, i, Sp, idx(basek, j, INT_G));
            double wikj = w[idx(basek, j, INT_G)];

            for(unsigned int m = (KORD-1 == i), im = i-KORD + m ; m<KORD && im<nb; ++m, ++im) {
                double sp_m = Sp[m];
                for(unsigned int n = (KORD-1 == i), in = i-KORD+n; n<KORD && in<nb; ++n, ++in) {

                    s[idx(im, in, nb)] += sp_m * Sp[n] * wikj;

                    ke[idx(im, in, nb)] += ma*sp_m * Sp[n] * wikj * _rr2;

                    if(RADIO_1<rr && rr<RADIO_2) v0[idx(im, in, nb)] += sp_m * Sp[n] * wikj;
                }
            }   
        }
    }

    for(unsigned int m = 1; m<=KORD-1; ++m) {

        for(unsigned int n = m; n<=KORD-1; ++n) {

            for(unsigned int j = 0; j<INT_G; ++j) {

                double bm = 0, bn = 0;

                rr = x[idx(BASE_KORD, j, INT_G)];

                bm = bder(rr, m, KORD-1, Sp, bm, idx(BASE_KORD, j, INT_G));
                bn = bder(rr, n, KORD-1, Sp, bn, idx(BASE_KORD, j, INT_G));

                ke[idx(m-1, n-1, nb)] = ke[idx(m-1, n-1, nb)] + 0.5*w[idx(BASE_KORD, j, INT_G)]*bm*bn/ME;

            }
        }
    }

    for(unsigned int i = KORD, basek=1; i<KORD+L_INTERVALS-1; ++i, ++basek) {
        // ojo con los indices en esta parte //
        //int basek = i - (KORD-1);
        for(unsigned int m = i-KORD+1; m<=i && m<nb ; ++m) {
            for(unsigned int n = m; n<=i && n<nb; ++n) {

                for(unsigned int j = 0; j<INT_G; ++j) {

                    double bm = 0, bn = 0;

                    rr = x[idx(basek, j, INT_G)];

                    bm = bder(rr, m, i, Sp, bm, idx(basek, 0, INT_G) + j);
                    bn = bder(rr, n, i, Sp, bn, idx(basek, 0, INT_G) + j);

                    ke[idx(m-1, n-1, nb)] += 0.5*w[idx(basek, 0, INT_G) + j]*bm*bn/ME;

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
    auval = (double *) calloc( NEV, sizeof(double));
    auvec = (double *) calloc( NEV*nb, sizeof(double));
    
    lambda = 1.0;
    delta = (LAMBDA_FIN-LAMBDA_IN)/NUMEROS_PUNTO_LAMBDA;

    // abro el archivo para guardar los datos //
    fprintf(archivo, "# Autovalores calculados\n");
    fprintf(archivo, "# Lambda  auval[0]   auval[1] ....\n");

    for(unsigned int j = 0; j<=NUMEROS_PUNTO_LAMBDA; ++j) {

        lambda = LAMBDA_IN + delta*j;
    
        for(unsigned int m = 0; m<nb; ++m) {
            for(unsigned int n = 0; n<nb; ++n){
                h[idx(n, m, nb)] = ke[idx(n, m, nb)] - lambda*v0[idx(n, m, nb)];
                s_copy[idx(n, m, nb)] = s[idx(n, m, nb)];
            }
        }

        eigenvalues( nb, NEV, h, s_copy, auval, auvec );

        fprintf(archivo, "%.5f   ", lambda);
        for(unsigned int i = 0; i<NEV; ++i) {
            fprintf(archivo, "%.15f   ", auval[i]);
        }
        fprintf(archivo, "\n");

    }

    free(h); free(s_copy); free(auval); free(auvec);
}*/

int main(void) {

    // defino algunas variables que voy a usar //
    double t_in, t_fin, t_n;
    FILE * archivo;

    // controlo algunos parametros //
    assert(INT_G>KORD);
    assert(NEV>0);

    archivo = fopen("autovalores_pozo-1.dat", "w");
    // imprimo los parametros //
    fprintf(archivo, "# Rmin=%.12f y R_MAX=%.12f\n", R_MIN, R_MAX);
    fprintf(archivo, "# Numero de intervalos l=%i\n", L_INTERVALS);
    fprintf(archivo, "# Orden los B-splines kord=%i\n", KORD);
    fprintf(archivo, "# Radios del pozo RADIO_1=%.12f y RADIO_2=%.12f\n", RADIO_1, RADIO_2);
    fprintf(archivo, "# Masa de la particula me=%.12f\n", ME);
    fprintf(archivo, "# Grado de integracion de la cuadratura INT_G=%i\n", INT_G);
    fprintf(archivo, "# Numero de autovalores NEV=%i\n", NEV);
    fprintf(archivo, "# Momento angular que usamos L_MAX=%i\n", L_MAX);
    fprintf(archivo, "# Numero de knots nk=%i\n", nk);
    fprintf(archivo, "# Tamaño de la base nb=%i\n", nb);
    fprintf(archivo, "# Valores inicial y final del pozo, LAMBDA_IN=%.12f, LAMBDA_FIN=%.12f\n", LAMBDA_IN, LAMBDA_FIN);
    fprintf(archivo, "# Numero de puntos lambda = %i\n", NUMEROS_PUNTO_LAMBDA);

    // doy memoria a las matrices que voy a necesitar //
    //k = (int *) malloc( nk*sizeof(int));
    //t = (double *) malloc( nk*sizeof(double));
   // x = (double *) malloc( L_INTERVALS*INT_G*sizeof(double));
   // w = (double *) malloc( L_INTERVALS*INT_G*sizeof(double));
    //s = (double *) malloc( nb*nb*sizeof(double));
  //  v0 = (double *) malloc( nb*nb*sizeof(double));
//    ke = (double *) malloc( nb*nb*sizeof(double));


    
    cleari(nk, k);
    cleard(nk, t);
    cleard(L_INTERVALS*INT_G, x);
    cleard(L_INTERVALS*INT_G, w);
    cleard(nb*nb, s);
    cleard(nb*nb, v0);
    cleard(nb*nb, ke);

    t_in = omp_get_wtime();
    // primero calculos los knost y los pesos para hacer la cuadratura //
    knots_pesos(x, w);

    // calculo las matrices que necesito para resolver el problema //
    calculo_matrices(x, w, s, v0, ke);
    t_fin = omp_get_wtime();

    // armo el hamiltoniano y calculo los autovalores //
//  hamiltoniano_autovalores(nb, s, v0, ke, archivo);
    
    t_n = (t_fin-t_in)/(nb*nb*L_INTERVALS*INT_G),
    printf("%i     %i     %i     %.12f  %.12f\n", L_INTERVALS, nk, nb, t_fin-t_in, t_n); 
    FILE *file;
    file = fopen("matrices.dat", "w");
    for(unsigned int i = 0; i<nb; ++i)
        for(unsigned int j=0; j<nb; ++j)
            fprintf(file, "%i\t%i\t%.12f\t%.12f\t%.12f\n", i, j, s[idx(i,j,nb)], v0[idx(i,j,nb)], ke[idx(i,j,nb)]);
    fclose(file);
    
    // libero la memoria //
    fclose(archivo);
//    free(ke);
//    free(v0);
//    free(s);
//    free(w);
//    free(x);
//    free(t);
 //   free(k);
    return 0;
}

