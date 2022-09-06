subroutine cal_dt(pe,drmin,dtheta,dt)
    implicit none
    double precision :: pe,drmin,dtheta,dt
    ! dt = (drmin**2)*pe/8.0d0
    ! dt=1.d-5
    dt =drmin*0.01d0
    dt = min(dt,(drmin**2)*pe/8.0d0)
    ! dt = drmin*0.1d0
    ! dt = min(dt,(dtheta**2)*pe/8.0d0)
endsubroutine cal_dt

! subroutine cal_dt(nr,ntheta,advis,c,time,dt)
!     implicit none
!     integer :: nr,ntheta
!     double precision :: c(-1:nr,0:ntheta-1)
!     double precision :: advis(-1:nr,0:ntheta-1)
!     double precision :: time,dt
!     integer :: i,j
!     double precision :: mmm,tmp,cfl

!     mmm = 1.d-2

! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i,j,tmp) &
! !$OMP& REDUCTION(min:mmm)

!     do j = 0,ntheta-1
!         do i = 0,nr-1
!             tmp = (1.0d0-c(i,j))/advis(i,j)
!             if(tmp>0) then
!                 mmm = min(mmm,tmp)
!             endif
!         enddo
!     enddo
! !$OMP  END PARALLEL DO

!     if(mmm>0.01) then
!         cfl = 0.9
!     else if(mmm<1.d-5) then
!         cfl = 0.1
!     else
!         cfl = 0.25
!     endif

!     dt = cfl*mmm
!     ! dt=1.d-10
!     time = time+dt
! endsubroutine cal_dt
