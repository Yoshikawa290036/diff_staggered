subroutine update(nr,ntheta,c,advis,dt)
    implicit none
    integer :: nr,ntheta
    double precision :: c(-1:nr,0:ntheta-1)
    double precision :: advis(-1:nr,0:ntheta-1)
    double precision :: dt
    integer :: i,j

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do j = 0,ntheta-1
        do i = 0,nr-1
            c(i,j) = c(i,j)+dt*advis(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

endsubroutine update
