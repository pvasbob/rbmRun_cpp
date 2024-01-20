#include <iostream>

#include "DiagMat.h"
#include "MsgToScreen.h"
#include "MultiDimArraySetToValue.h"
#include "MultiDimArrayPrint.h"
#include "MultiDimArrayAllocate.h"

// external C library.
extern "C"
{
    // Declaration of zgeev from LAPACK
    extern void zgeev_(
        char *jobvl,
        char *jobvr,
        int *n,
        std::complex<double> *a,
        int *lda,
        std::complex<double> *w,
        std::complex<double> *vl,
        int *ldvl,
        std::complex<double> *vr,
        int *ldvr,
        std::complex<double> *work,
        int *lwork,
        double *rwork,
        int *info);
}
//
extern "C"
{
    // LAPACK zgesv function for solving linear equations
    void zgesv_(
        const int *n,
        const int *nrhs,
        std::complex<double> *a,
        const int *lda,
        int *ipiv,
        std::complex<double> *b,
        const int *ldb,
        int *info);
}
//
void diagGenComplexMat(int &dim,
                       std::complex<double> **mat,
                       std::complex<double> **pmat,
                       std::complex<double> *eigenvalue,
                       std::complex<double> **pmatinv,
                       int &info)
{
    msgToScreen("Inside diagGenComplexMat");

    std::complex<double> *a = 0;
    a = new std::complex<double>[dim * dim];
    //
    copy2dTo1dColMajor(mat, a, dim, dim, dim * dim);
    // for (int i = 0; i < dim; i++)
    //{
    //     for (int j = 0; j < dim; j++)
    //     {
    //         a[i * dim + j] = mat[i][j];
    //     }
    // }

    char jobvl = 'N';
    char jobvr = 'V';
    int n = dim;
    int m = dim;
    int lda = dim;
    int ldvl = dim;
    int ldvr = dim;
    int lwork_query = -1; // Query optimal workspace size
    // int info;

    // Variables for workspace query
    std::complex<double> work_query;
    double rwork_query;

    // Call zgeev with query to get optimal workspace size
    // zgeev_()
    zgeev_(&jobvl, &jobvr, &n, a, &lda, nullptr, nullptr, &ldvl, nullptr, &ldvr, &work_query, &lwork_query, &rwork_query, &info);

    // Check for successful query
    if (info != 0)
    {
        std::cerr << "Error: zgeev workspace query failed with info = " << info << std::endl;
        exit(1);
    }

    // Get optimal workspace size
    int lwork = static_cast<int>(work_query.real()) + 1; // Add 1 as a safety margin

    // Allocate arrays for zgeev
    std::complex<double> *w = new std::complex<double>[n];
    std::complex<double> *vr = new std::complex<double>[ldvr * n];
    std::complex<double> *vl = new std::complex<double>[ldvl * n];
    std::complex<double> *work = new std::complex<double>[lwork];
    double *rwork = new double[2 * n];
    //
    set1dArrayToValue<std::complex<double>>(w, n, 0.0);
    set1dArrayToValue<std::complex<double>>(vr, ldvr * n, 0.0);
    set1dArrayToValue<std::complex<double>>(vl, ldvl * n, 0.0);
    set1dArrayToValue<std::complex<double>>(work, lwork, 0.0);
    set1dArrayToValue<double>(rwork, 2 * n, 0.0);

    // Call zgeev with allocated workspace
    zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

    // Check for successful computation
    if (info != 0)
    {
        std::cerr << "Error: zgeev failed with info = " << info << std::endl;
        exit(1);
    }

    // Print eigenvalues
    // std::cout << "Eigenvalues:" << std::endl;
    msgToScreen("eigen value:");
    print1d<std::complex<double>>(w, n);
    // for (int i = 0; i < dim; ++i)
    // {
    // std::cout << w[i] << std::endl;
    // }
    //
    // In fortran this is 2d to 2d, but in cpp vr is 1d
    copy1dTo2dColMajor(vr, pmat, ldvr * n, dim, dim);
    copy1dTo1d(w, eigenvalue, n);
    //
    msgToScreen("pmat:");
    print2d(pmat, dim, dim);
    msgToScreen("eigenvalue:");
    print1d(eigenvalue, n);
    //
    //
    int ldb = dim;
    int nrhs = dim;
    // std::complex<double> ** b =
    std::complex<double> **bmat;
    std::complex<double> *b;
    int *ipiv;

    allocate2dArray<std::complex<double>>(bmat, ldb, nrhs);
    allocate1dArray<std::complex<double>>(b, ldb * nrhs);
    allocate1dArray<int>(ipiv, n);
    //
    set2dArrayToValue<std::complex<double>>(bmat, ldb, nrhs, 0.0);
    set1dArrayToValue<int>(ipiv, n, 0.0);
    set1dArrayToValue<std::complex<double>>(b, ldb * nrhs, 0.0);
    //
    setDiagToValue(bmat, ldb, std::complex<double>(1.0, 0.0));
    copy2dTo1dColMajor(bmat, b, ldb, nrhs, ldb * nrhs);
    //
    copy1dTo1d(vr, a, ldvr * n);
    //
    // msgToScreen("before zgesv:");
    // print1d(b, ldb * nrhs);
    zgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    // msgToScreen("after zgesv:");
    // print1d(b, ldb * nrhs);
    // zgesv_()
    // msgToScreen("After zgesv_:");
    //
    copy1dTo2dColMajor(b, pmatinv, ldb * nrhs, dim, dim);
    //
    deallocate1dArray<std::complex<double>>(a);
    deallocate1dArray<std::complex<double>>(w);
    deallocate1dArray<std::complex<double>>(vl);
    deallocate1dArray<std::complex<double>>(vr);
    deallocate1dArray<std::complex<double>>(work);
    deallocate1dArray<double>(rwork);
    deallocate1dArray<std::complex<double>>(b);
    deallocate1dArray<int>(ipiv);
    deallocate2dArray<std::complex<double>>(bmat, ldb);

    msgToScreen("pmatinv:");
    print2d(pmatinv, dim, dim);
    //
    std::complex<double> **tempmat;
    allocate2dArray(tempmat, dim, dim);
    set2dArrayToValue<std::complex<double>>(tempmat, dim, dim, 0.0);
    // int dim = nrbm;
    for (int i = 1; i <= dim; i++)
    {
        for (int j = 1; j <= dim; j++)
        {
            for (int k = 1; k <= dim; k++)
            {
                tempmat[-1 + i][-1 + j] += pmat[-1 + i][-1 + k] * eigenvalue[-1 + k] * pmatinv[-1 + k][-1 + j];
            }
        }
    }
    //
    msgToScreen("tempmat:");
    print2d(tempmat, dim, dim);
}
