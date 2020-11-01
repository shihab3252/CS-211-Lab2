#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 66;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    int i, maxI;
    double max;
    double *tmpr = (double*) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++){
        // pivoting
        maxI = i;
        max = fabs(A[i*n + i]);
        
        int j;
        for (j = i+1; j < n; j++){
            if (fabs(A[j*n + i]) > max)
            {
                maxI = j;
                max = fabs(A[j*n + i]);
            }
        }
        if (max == 0)
        {
            printf("LU factorization failed: coefficient matrix is singular.\n");
            return -1;
        }
        else
        {
            if (maxI != i)
            {
                // pivoting is done here
                int temp = ipiv[i];
                ipiv[i] = ipiv[maxI];
                ipiv[maxI] = temp;
                // swap rows
                memcpy(tmpr, A + i*n, n * sizeof(double));
                memcpy(A + i*n, A + maxI*n, n * sizeof(double));
                memcpy(A + maxI*n, tmpr, n * sizeof(double));
            }
        }

        // factorization is done here
        for (j = i+1; j < n; j++){
            A[j*n + i] = A[j*n + i] / A[i*n + i];
            int k;
            for (k = i+1; k < n; k++){
                A[j*n + k] -= A[j*n +i] * A[i*n + k];
            }
        }
    }
    free(tmpr);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    double *y = (double*) malloc(n * sizeof(double));
    int i, j;
    double sum;
    if (UPLO == 'L')
    {
        y[0] = B[ipiv[0]];
        for (i = 1; i < n; i++){
            sum = 0.0;
            for (j = 0; j < i; j++){
                sum += y[j] * A[i*n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    }
    else if (UPLO == 'U')
    {
        y[n - 1] = B[n - 1] / A[(n-1)*n + n-1];
        for (i = n-2; i >= 0; i--){
            sum = 0;
            for (j = i+1; j < n; j++){
                sum += y[j] * A[i*n + j];
            }
            y[i] = (B[i] - sum) / A[i*n + i];
        }
    }

    memcpy(B, y, sizeof(double) * n);
    free(y);
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    int i1 = i, j1 = j, k1 = k;
    int ni = i + b > n ? n : i + b;
    int nj = j + b > n ? n : j + b;
    int nk = k + b > n ? n : k + b;

    for (i1 = i; i1 < ni; i1 += 3){
        for (j1 = j; j1 < nj; j1 += 3){
            int t = i1 * n + j1;
            int tt = t + n;
            int ttt = tt + n;
            register double c00 = C[t];
            register double c01 = C[t + 1];
            register double c02 = C[t + 2];
            register double c10 = C[tt];
            register double c11 = C[tt + 1];
            register double c12 = C[tt + 2];
            register double c20 = C[ttt];
            register double c21 = C[ttt + 1];
            register double c22 = C[ttt + 2];

            for (k1 = k; k1 < nk; k1 += 3){
		int l;
                for (l = 0; l < 3; l++){
                    int ta = i1 * n + k1 + l;
                    int tta = ta + n;
                    int ttta = tta + n;
                    int tb = k1 * n + j1 + l * n;
                    register double a0 = A[ta];
                    register double a1 = A[tta];
                    register double a2 = A[ttta];
                    register double b0 = B[tb];
                    register double b1 = B[tb + 1];
                    register double b2 = B[tb + 2];

                    c00 -= a0 * b0;
                    c01 -= a0 * b1;
                    c02 -= a0 * b2;
                    c10 -= a1 * b0;
                    c11 -= a1 * b1;
                    c12 -= a1 * b2;
                    c20 -= a2 * b0;
                    c21 -= a2 * b1;
                    c22 -= a2 * b2;
                }
            }
            C[t] = c00;
            C[t + 1] = c01;
            C[t + 2] = c02;
            C[tt] = c10;
            C[tt + 1] = c11;
            C[tt + 2] = c12;
            C[ttt] = c20;
            C[ttt + 1] = c21;
            C[ttt + 2] = c22;
        }
    }
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, i, j, k, maxI;
    double max, sum;
    double *tmpr = (double*) malloc(sizeof(double) * n);

    for (ib = 0; ib < n; ib += b){
        for (i = ib; i < ib+b && i < n; i++){
            // pivoting
            maxI = i;
            max = fabs(A[i*n + i]);
            
            int j;
            for (j = i+1; j < n; j++){
                if (fabs(A[j*n + i]) > max)
                {
                    maxI = j;
                    max = fabs(A[j*n + i]);
                }
            }
            if (max == 0)
            {
                printf("LU factorization failed: coefficient matrix is singular.\n");
                return -1;
            }
            else
            {
                if (maxI != i)
                {
                    // save pivoting information
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[maxI];
                    ipiv[maxI] = temp;
                    // swap rows
                    memcpy(tmpr, A + i*n, n * sizeof(double));
                    memcpy(A + i*n, A + maxI*n, n * sizeof(double));
                    memcpy(A + maxI*n, tmpr, n * sizeof(double));
                }
            }

            // factorization
            for (j = i+1; j < n; j++){
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                int k;
                for (k = i+1; k < ib+b && k < n; k++){
                    A[j*n + k] -= A[j*n +i] * A[i*n + k];
                }
            }
        }

        // update A(ib:end, end+1:n)
        for (i = ib; i < ib+b && i < n; i++){
            for (j = ib+b; j < n; j++){
                sum = 0;
                for (k = ib; k < i; k++){
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }

        // update A(end+1:n, end+1:n)
        for (i = ib+b; i < n; i += b){
            for (j = ib+b; j < n; j += b){
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    return 0;
}
void swap(double* A, double* tmpr, int n, int r1, int r2)
{
    memcpy(tmpr, A + r1 * n, n * sizeof(double));
    memcpy(A + r1 * n, A + r2 * n, n * sizeof(double));
    memcpy(A + r2 * n, tmpr, n * sizeof(double));
}
void transpose(double* A, int m, int n)
{
    int i, j;
    double *tmp = (double*)malloc(sizeof(double) * n * m);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            tmp[i * n + j] = A[j * n + i];
        }
    }
    memcpy(A, tmp, sizeof(double) * n * m);
    free(tmp);
}
int mydgetrf_non_squrare_naive(double* A, int pos, int* ipiv, int n, int bm, int bn, int b)
{
    /* add your code here */
    int i, j, k, i1, j1;
    int bn2 = bm - bn;
    double* tmpr = (double*)malloc(sizeof(double) * n);
    double* LLT = (double*)malloc(sizeof(double) * bn * bn);
    double* AUR = (double*)malloc(sizeof(double) * bn * bn2);
    double* AURD = (double*)malloc(sizeof(double) * bn * bn2);
    double* ALLD = (double*)malloc(sizeof(double) * bn2 * bn);
    double* LL = (double*)malloc(sizeof(double) * bn * bn);
    int* ipivl = (int*)malloc(sizeof(int) * bn);

    for (i = 0; i < bn; i++)
    {
        int maxidx = i;
        double max = fabs(A[i * n + i]);
        for (j = i + 1; j < bm; j++)
        {
            double tmp = fabs(A[j * n + i]);
            if (tmp - max > 1e-6)
            {
                maxidx = j;
                max = tmp;
            }
        }

        //too small pivot is also unacceptable
        if (fabs(max - 0.0) < 1e-3)
            return -1;

        if (maxidx != i)
        {
            int newMaxidx = pos + maxidx;
            int newI      = pos + i;
            ipiv[newMaxidx] = ipiv[newMaxidx] ^ ipiv[newI];
            ipiv[newI] = ipiv[newMaxidx] ^ ipiv[newI];
            ipiv[newMaxidx] = ipiv[newMaxidx] ^ ipiv[newI];

            swap(A-pos, tmpr, n, i, maxidx);
        }

        for (j = i + 1; j < bm; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            double A_j = A[j * n + i];
            for (k = i + 1; k < bn; k++)
            {
                A[j * n + k] -= A_j * A[i * n + k];
            }
        }
    }

    if (bn2 > 0)
    {
        memset(LLT, 0, bn * bn * sizeof(double));
        memset(LL, 0, bn * bn * sizeof(double));
        memset(ipivl, 0, bn * sizeof(int));
        memset(AUR, 0, bn * bn2 * sizeof(double));
        memset(AURD, 0, bn * bn2 * sizeof(double));
        memset(ALLD, 0, bn * bn2 * sizeof(double));

        for (i = 0; i < bn; i++)
        {
            LLT[i * bn + i] = 1;
            ipivl[i] = i;
            
            LL[i * bn + i] = 1;
            for (j = 0; j < i; j++)
            {
                LL[i * bn + j] = A[i * n + j];
            }
        }
            
        for (i = 0; i < bn; i++)
        {
            memcpy(AUR + i * bn2, A + i * n + bn, bn2 * sizeof(double));
        }

        //get LL inverse and store in LLT
        for (i = 0; i < bn; i++)
        {
            mydtrsv('L', LL, LLT + i * bn, bn, ipivl);
        }

        transpose(LLT, bn, bn);

        //A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n)
        mydgemm(LLT, AUR, AURD, bn, bn, bn2, b);

        for (i = 0; i < bn; i++)
        {
            memcpy(A + i * n + bn, AURD + i * bn2, bn2 * sizeof(double));
        }
        //A(end+1:n , end+1:n )-= A(end+1:n , ib:end) * A(ib:end , end+1:n)    
        mydgemm(A + bn * n, A + bn, A + bn * n + bn, bn2, bn, bn2, n, b);
    }

    free(LLT);
    free(LL);
    free(ipivl);
    free(AUR);
    free(AURD);
    free(tmpr);
    return 0;
}
int mydgetrf_block_naive(double *A, int *ipiv, int n, int b) 
{
    int i, j, k;

    double* Aptr = A;
    for (i = 0; i < n - b; i += b)
    {
        mydgetrf_non_squrare_naive(Aptr, i, ipiv, n, n - i, b, b);
        Aptr += b * n + b;
    }
    int blocksize = n % b > 0 ? n % b : b;
    int bias = n - blocksize;
    mydgetrf_non_squrare_naive(Aptr, bias, ipiv, n, blocksize, blocksize, blocksize);
    return 0;
}


