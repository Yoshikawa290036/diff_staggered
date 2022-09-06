program main
    implicit none

    integer :: i,j
    integer :: nr,ntheta
    integer :: max_step,nstep
    double precision :: PI,Pe,peinv,sh,alpha,alphainv
    double precision :: lambda
    double precision :: drmin,drmininv
    double precision :: dtheta,dthetainv
    double precision :: time,dt

    double precision,dimension(:),allocatable :: rs,thetas,dr3inv,drinv,dsinv
    double precision,dimension(:,:),allocatable :: ur,uw,c,advis,xs,ys

    integer :: nmkxyc,nmkstd
    character(32) fname

    PI = atan(1.0d0)*4.0d0
    read(*,*) Pe

! ================ system set ================ !
    max_step = 10000000

!    nr = 512
!    ntheta = 256
    nr = 512/2
    ntheta = 256/2

    nmkxyc = 100000
    nmkstd = 10000

    drmin = 0.001d0

!    alpha = 1.02d0
    alpha = 1.04d0

    lambda = 3.0d0/4.0d0
! ================ system set ================ !

    peinv = 1.0d0/Pe
    alphainv = 1.0d0/alpha
    dtheta = PI/dble(ntheta)
    dthetainv = dble(ntheta)/PI
    drmininv = 1.0d0/drmin
    nstep = 0
    time = 0.0d0
    dt = 0.0d0

    include'allocate.h'

    call cal_rs(nr,drmin,alpha,rs,dr3inv,drinv)
    call cal_thetas(ntheta,dtheta,thetas,dsinv)
    call cal_xs_ys(nr,ntheta,rs,thetas,xs,ys)

    write(*,'("max step                 ",1i9)') max_step
    write(*,'("nmkxyc                   ",1i9)') nmkxyc
    write(*,'("nmkstd                   ",1i9)') nmkstd
    write(*,*)
    write(*,'("Nr                       ",1i9)') nr
    write(*,'("Ntheta                   ",1i9)') ntheta
    write(*,'("Pe                       ",20e20.10)') Pe
    write(*,'("Rmax                     ",20e20.10)') rs(nr)
    write(*,'("drmin                    ",20e20.10)') drmin
    write(*,'("alpha                    ",20e20.10)') alpha
    write(*,'("lambda                   ",20e20.10)') lambda
    write(*,*)
    write(*,*)

    ! do i = 0,nr-1
    !     tmp = 1.0d0/(rs(nr)-1.0d0)*(rs(nr)/rs(i)-1.0d0)
    !     do j = 0, ntheta-1
    !         c(i,j)=tmp
    !     end do
    ! enddo

    call cal_vel(nr,ntheta,lambda,xs,ys,rs,thetas,ur,uw)

!    do j=0,ntheta-1
!    do i=0,nr
!    write(101,'(20e20.10)')dble(i),dble(j),rs(i),0.5d0*(thetas(j)+thetas(j+1)),ur(i,j)*dsinv(j)/rs(i)**2
!    enddo
!    write(101,'()')
!    enddo
!
!    do j=0,ntheta
!    do i=0,nr-1
!    write(102,'(20e20.10)')dble(i),dble(j),(rs(i)+rs(i+1))*0.5d0,thetas(j),uw(i,j)/sin(thetas(j))*2.0d0/(-(rs(i)**2)+(rs(i+1)**2))
!    enddo
!    write(102,'()')
!    enddo
!
!    stop


    call cal_sh(nr,ntheta,drmininv,dtheta,thetas,c,sh,drinv)
    call cal_dt(pe,drmin,dtheta,dt)

    open(12,file='tdtsh')
    write(12,'(20e20.10)') time,dt,sh

    do nstep = 1,max_step
        call bndset(nr,ntheta,alpha,alphainv,c)
!        write(*,*)peinv,dthetainv
        call cal_advis(nr,ntheta,peinv,rs,thetas,drinv,dr3inv,dsinv,dthetainv,ur,uw,c,advis)

!do j = 0,ntheta-1
!do i = 0,nr-1
!write(301,'(20e20.10)')dble(i),dble(j),(rs(i)+rs(i+1))*0.5d0,0.5d0*(thetas(j)+thetas(j+1)),advis(i,j)*(-rs(i)+rs(i+1))
!enddo
!enddo
!
!do j = 0,ntheta-1
!do i = 0,nr-1
!write(302,'(20e20.10)')dble(i),dble(j),(rs(i)+rs(i+1))*0.5d0,0.5d0*(thetas(j)+thetas(j+1)),advis(i,j)*(rs(i)+rs(i+1))*0.5*dtheta
!enddo
!enddo
!
!stop



        call update(nr,ntheta,c,advis,dt)

!        write(*,*)pe,peinv

        if(mod(nstep,nmkstd)==0) then
            call cal_sh(nr,ntheta,drmininv,dtheta,thetas,c,sh,drinv)
            time = dt*dble(nstep)
            write(12,'(20e20.10)') time,dt,sh
            write(*,*) '---------------------------------------'
            write(*,'("nstep        ",1i9.9)') nstep
            write(*,'("time         ",20e20.10)') time
            write(*,'("dt           ",20e20.10)') dt
            write(*,'("Sh           ",20e20.10)') Sh
            ! call flush (6)
        endif

        if(mod(nstep,nmkxyc)==0) then
            include'mk_xyc.h'
        endif
        call flush (6)
    enddo
    close(12)

endprogram main
