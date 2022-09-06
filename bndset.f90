subroutine bndset(nr,ntheta,alpha,alphainv,c)
    implicit none
    integer :: nr,ntheta
    double precision :: c(-1:nr,0:ntheta-1)
    double precision :: alphainv,alpha
    integer :: j

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j)
do j = 0,ntheta-1
  c(-1,j) = -c(   0,j)+2.0d0
  c(nr,j) = -c(nr-1,j)
enddo
!$OMP  END PARALLEL DO

endsubroutine bndset
