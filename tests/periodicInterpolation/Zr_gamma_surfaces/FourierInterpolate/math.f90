MODULE Math

CONTAINS

  SUBROUTINE Mat2Inv(A, B)
    ! return matrix B which is the inverse of matrix A

    implicit none

    REAL(kind(0.d0)), dimension(2,2), intent(in) :: A
    REAL(kind(0.d0)), dimension(2,2), intent(out) :: B

    REAL(kind(0.d0)) :: invdet

    invdet=1.d0/( A(1,1)*A(2,2)-A(1,2)*A(2,1) )
    B(1,1)=A(2,2)*invdet
    B(2,2)=A(1,1)*invdet
    B(1,2)=-A(1,2)*invdet
    B(2,1)=-A(2,1)*invdet

  END SUBROUTINE Mat2Inv

END MODULE Math
