#include "../include/lu_functions.h"
#include "../include/for_you_to_do.h"

void block_naive(double *A, double *B, int n, int b)
{
    int ipiv[n], i;
    for (i = 0; i < n; i++) {
        ipiv[i] = i;
    }

    int success = mydgetrf_block_naive(A, ipiv, n, b);

    if (success) 
    {
        printf("LU factoration failed: coefficient matrix is singular in naive.\n");
        return;
    }

    mydtrsv('L', A, B, n, ipiv);
    mydtrsv('U', A, B, n, ipiv);

}