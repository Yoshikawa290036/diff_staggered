! calculate rs

!       rs(0) ------- rs(1)  ....       rs(nr)
!   ^   wall    ^     fluid              wall   ^
!  C(-1)       C(0)                            C(nr)

subroutine cal_rs(nr,drmin,alpha,rs,dr3inv,drinv)
    implicit none

    integer :: nr
    double precision :: drmin,alpha
    double precision :: rs(0:nr)
    double precision :: dr3inv(0:nr-1)
    double precision :: drinv(0:nr)

    integer :: i
    double precision :: alinv

    rs(0) = 1.0d0
    alinv = 1.0d0/(1.0d0-alpha)

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i)
    do i = 1,nr
        rs(i) = 1.0d0+drmin*(1.0d0-alpha**i)*alinv
    enddo
!$OMP  END PARALLEL DO

    ! do i = 0, nr
    !     write (*, *) rs(i)
    ! end do
    ! stop

    dr3inv(0) = 3.0d0/(-rs(0)**3+rs(1)**3)
    drinv(0) = 2.0d0/(drmin+drmin/alpha)
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i)
    do i = 1,nr-1
        dr3inv(i) = 3.0d0/(-rs(i)**3+rs(i+1)**3)
        drinv(i) = 2.0d0/(-rs(i-1)+rs(i+1))
    enddo
!$OMP  END PARALLEL DO

    ! drinv(nr) = 2.0d0/(-rs(nr-1)+(1.0d0+drmin*(1.0d0-alpha**(nr+1))*alinv))
    drinv(nr) = 2.0d0/(drmin*alpha**(nr-1)+drmin*alpha**(nr))
    ! drinv(nr) = 2.0d0/(drmin*alpha**(nr-1)+ )
endsubroutine cal_rs
