subroutine cal_thetas(ntheta,dtheta,thetas,dsinv)
    implicit none

    integer :: ntheta
    double precision :: dtheta
    double precision :: thetas(0:ntheta)
    double precision :: dsinv(0:ntheta-1)

    integer :: j

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j)
    do j = 0,ntheta
        thetas(j) = dtheta*dble(j)
    enddo
!$OMP  END PARALLEL DO


!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j)
    do j = 0,ntheta-1
        dsinv(j) = 1.0d0/(cos(thetas(j))-cos(thetas(j+1)))
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_thetas
