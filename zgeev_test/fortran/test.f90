
program DiagonalizeComplexMatrix
    implicit none



        !   Implicit None
      Integer, Parameter:: pr =Kind(1.0d0)
      Integer:: dim = 3
      Complex(pr), dimension(:,:) , Allocatable :: mat, pmat, pmatinv
      Complex(pr), dimension(:)   , Allocatable :: eigenvalue
      Integer :: info
   
      !zgesvd/zgesv variables
      Character :: jobvl, jobvr
      Integer :: m, n, lda, ldb, nrhs, ldvl, ldvr, lwork
      Integer, Allocatable :: ipiv(:)
      Complex(pr), Allocatable :: a(:,:), w(:), vl(:,:), vr(:,:), work(:), b(:,:)
      Real(pr), Allocatable :: rwork(:)
   
      Complex(pr), Allocatable :: tempmat(:,:)
      Integer :: i, j, k
       

    Allocate(mat(dim, dim))
    Allocate(pmat(dim, dim))
    Allocate(pmatinv(dim, dim))
    !
    Allocate(eigenvalue(dim))
    !
    Allocate(tempmat(dim, dim))




      jobvl = 'n'
      jobvr = 'v'
      m = dim
      n = dim
      lda = dim
      ldvl = dim
      ldvr = dim
      lwork = 4*n+10 ! >= 2*N


      ! Initialize complex matrix A using two loops
      !do i = 1, dim
      !    do j = 1, dim
      !        mat(i, j) = cmplx(1.0 * (i + 1) + j, 2.0 * (i + 1) - j)
      !      !   mat(i, j) = mat(i, j) + 1.0
      !    end do
      !end do
      mat(1,1) = cmplx(0.0, 0.0)
      mat(1,2) = cmplx(1.0, -1.0)
      mat(1,3) = cmplx(0.0, 0.0)
      !
      mat(2,1) = cmplx(1.0, 1.0)
      mat(2,2) = cmplx(0.0, 0.0)
      mat(2,3) = cmplx(1.0, -1.0)
      !
      mat(3,1) = cmplx(0.0, 0.0)
      mat(3,2) = cmplx(1.0, 1.0)
      mat(3,3) = cmplx(0.0, 0.0)

      Allocate(a(lda,n), w(n),vl(ldvl,n),vr(ldvr,n), work(lwork),rwork(2*n))
      a(:,:) = mat(:,:)
      w = 0.0_pr; vl = 0.0_pr; vr = 0.0_pr; work = 0.0_pr; rwork = 0.0_pr  
!      print *, "zgeev"
      Call zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
!      print *, "zgeev done"
    !   If (info .ne. 0) Then

    write(*,*) w
end program DiagonalizeComplexMatrix
