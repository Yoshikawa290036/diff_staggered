subroutine cal_sh(nr,ntheta,drmininv,dtheta,thetas,c,sh,drinv)
    implicit none

    integer :: nr,ntheta
    double precision :: drmininv,dtheta
    double precision :: sh
    double precision :: thetas(0:ntheta)
    double precision :: c(-1:nr,0:ntheta-1)
    double precision :: drinv(0:nr)

    integer :: j
    sh = 0.0d0

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j) &
!$OMP& REDUCTION(+:sh)
    do j = 0,ntheta-1
        sh = sh+((-cos(thetas(j))+cos(thetas(j+1)))*(C(0,j)-C(-1,j))*drinv(0))
    enddo

    ! do j = 0,ntheta-1
    !     sh = sh+((-cos(thetas(j))+cos(thetas(j+1))))
    ! enddo

!$OMP  END PARALLEL DO
    ! write(*,*) sh
    ! stop
endsubroutine cal_sh
