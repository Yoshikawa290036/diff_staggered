
subroutine cal_vel(nr,ntheta,lambda,xs,ys,rs,thetas,ur,uw)
    implicit none
    integer :: nr,ntheta
    double precision :: lambda
    double precision :: xs(0:nr-1,0:ntheta-1)
    double precision :: ys(0:nr-1,0:ntheta-1)
    double precision :: rs(0:nr)
    double precision :: thetas(0:ntheta)
    double precision :: ur(0:nr,0:ntheta-1)
    double precision :: uw(0:nr-1,0:ntheta)
    double precision :: nur(0:nr,0:ntheta-1)  ! natural ur
    double precision :: nuw(0:nr-1,0:ntheta)  ! natural uw
    double precision :: ux(0:nr-1,0:ntheta-1)
    double precision :: uy(0:nr-1,0:ntheta-1)
    integer :: i,j
    double precision :: theta,dcos
    double precision :: rinv,rrrinv,tmp
    double precision :: nurt,nuwt

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,theta,dcos)
    do j = 0,ntheta-1
        ! theta = 0.5d0*(thetas(j)+thetas(j+1))
        ! dcos = cos(thetas(j))-cos(thetas(j+1))
        do i = 0,nr
            theta = 0.5d0*(thetas(j)+thetas(j+1))
            dcos = cos(thetas(j))-cos(thetas(j+1))
            nur(i,j) = cos(theta)*(1.0d0-(2.0d0*lambda/rs(i))+(2.0d0*lambda-1.0d0)/rs(i)**3)
            ur(i,j) = (rs(i)**2)*dcos*nur(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,rinv,rrrinv,tmp)
    do i = 0,nr-1
        ! rinv = 2.0d0/(rs(i)+rs(i+1))
        ! rrrinv = 8.0d0/((rs(i)+rs(i+1))**3)
        ! tmp = -(rs(i)**2)+(rs(i+1)**2)
        do j = 0,ntheta
            rinv = 2.0d0/(rs(i)+rs(i+1))
            rrrinv = 8.0d0/((rs(i)+rs(i+1))**3)
            tmp = -(rs(i)**2)+(rs(i+1)**2)
            nuw(i,j) = sin(thetas(j))*(-1.0d0+(lambda*rinv)+(2.0d0*lambda-1.0d0)*0.5d0*rrrinv)
            uw(i,j) = 0.5d0*tmp*sin(thetas(j))*nuw(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,theta,nuwt,nurt)
    do j = 0,ntheta-1
        ! theta = 0.5d0*(thetas(j)+thetas(j+1))
        do i = 0,nr-1
            theta = 0.5d0*(thetas(j)+thetas(j+1))

            nuwt = 0.5d0*(nuw(i,j)+nuw(i,j+1))
            nurt = 0.5d0*(nur(i,j)+nur(i+1,j))
            ux(i,j) = nurt*cos(theta)-(0.5d0*(rs(i)+rs(i+1)))*sin(theta)*nuwt
            uy(i,j) = nurt*sin(theta)+(0.5d0*(rs(i)+rs(i+1)))*cos(theta)*nuwt
        enddo
    enddo
!$OMP  END PARALLEL DO

    open(10,file='xyuv')

    do j = 0,ntheta-1,2
        do i = 0,nr-1,2
            write(10,'(20e20.10)') xs(i,j),-ys(i,j),ux(i,j),-uy(i,j)
        enddo
        write(10,'()')
    enddo

    do j = 0,ntheta-1,2
        do i = 0,nr-1,2
            write(10,'(20e20.10)') xs(i,j),ys(i,j),ux(i,j),uy(i,j)
        enddo
        write(10,'()')
    enddo

    close(10)
    write(*,*) 'output :  xyuv'
endsubroutine cal_vel
