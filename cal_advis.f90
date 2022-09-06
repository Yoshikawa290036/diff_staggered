subroutine cal_advis(nr,ntheta,peinv,rs,thetas,drinv,dr3inv,dsinv,dthetainv,ur,uw,c,advis)
    implicit none
    integer :: nr,ntheta
    double precision :: peinv,dthetainv
    double precision :: rs(0:nr)
    double precision :: thetas(0:ntheta)
    double precision :: dr3inv(0:nr-1)
    double precision :: drinv(0:nr)
    double precision :: dsinv(0:ntheta-1)
    double precision :: ur(0:nr,0:ntheta-1)
    double precision :: uw(0:nr-1,0:ntheta)
    double precision :: c(-1:nr,0:ntheta-1)
    double precision :: advis(-1:nr,0:ntheta-1)

    integer :: i,j
    double precision :: adv,vis
    double precision :: cp1,cm1,aw,ae,as,an,dr


!        write(*,*)peinv,dthetainv
!        stop

!    do j=0,ntheta-1
!    do i=0,nr
!    write(201,'(20e20.10)')dble(i),dble(j),rs(i),0.5d0*(thetas(j)+thetas(j+1)),ur(i,j)*dsinv(j)/rs(i)**2
!    enddo
!    write(201,'()')
!    enddo
!
!    do j=0,ntheta
!    do i=0,nr-1
!    write(202,'(20e20.10)')dble(i),dble(j),(rs(i)+rs(i+1))*0.5d0,thetas(j),uw(i,j)/sin(thetas(j))*2.0d0/(-(rs(i)**2)+(rs(i+1)**2))
!    enddo
!    write(202,'()')
!    enddo
!
!    stop

!    do i=0,nr-1
!    write(203,'(20e20.10)')dble(i),(rs(i)+rs(i+1))*0.5d0,drinv(i)/dr3inv(i)
!    enddo
!    stop

!    do j=0,ntheta-1
!    write(204,'(20e20.10)')dble(j),sin(0.5d0*(thetas(j)+thetas(j+1))),dthetainv/dsinv(j)
!    enddo
!    stop

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,adv,vis,cm1,cp1,aw,ae,as,an,dr)
do j = 0,ntheta-1
do i = 0,nr-1

dr = rs(i+1)-rs(i)
if(j.ge.1)then
  cm1=c(i,j-1)
else
  cm1=c(i,j)
endif
if(j.le.ntheta-2)then
  cp1=c(i,j+1)
else
  cp1=c(i,j)
endif

adv = dr3inv(i)*dsinv(j)*0.5d0* ( &
      ( (-c(i-1,j)+c(i,j))*ur(i,j) + (-c(i,j)+c(i+1,j))*ur(i+1,j) ) &
     +( (-cm1     +c(i,j))*uw(i,j) + (-c(i,j)+cp1     )*uw(i,j+1) ) &
     )

aw=dr3inv(i)         *   rs(i  )**2*     drinv(i  )
ae=dr3inv(i)         *   rs(i+1)**2*     drinv(i+1)
as=dr3inv(i)*dr      *dsinv(j)   *sin(thetas(j  ))*dthetainv
an=dr3inv(i)*dr      *dsinv(j)   *sin(thetas(j+1))*dthetainv

vis = 2.0d0*peinv * (       &
     + aw*(c(i-1,j  )-c(i,j)) &
     + ae*(c(i+1,j  )-c(i,j)) &
     + as*(cm1       -c(i,j)) &
     + an*(cp1       -c(i,j)) &
     )

advis(i,j) = -adv+vis

enddo
enddo
!$OMP  END PARALLEL DO

endsubroutine cal_advis
