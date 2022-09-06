! C の定義点における x,y

subroutine cal_xs_ys(nr,ntheta,rs,thetas,xs,ys)
    implicit none

    integer :: nr,ntheta
    double precision :: rs(0:nr)
    double precision :: thetas(0:ntheta)
    double precision :: xs(0:nr-1,0:ntheta-1)
    double precision :: ys(0:nr-1,0:ntheta-1)

    integer :: i,j
    double precision :: r,theta

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,r,theta)
    do j = 0,ntheta-1
        ! theta = (thetas(j)+thetas(j+1))*0.5d0
        do i = 0,nr-1
            theta = (thetas(j)+thetas(j+1))*0.5d0

            r = (rs(i)+rs(i+1))*0.5d0
            xs(i,j) = r*cos(theta)
            ys(i,j) = r*sin(theta)
        enddo
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_xs_ys
