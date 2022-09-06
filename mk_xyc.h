write(fname,'("xyc",i7.7)') nstep
open(11,file=fname)

do j = 0,ntheta-1
    do i = 0,nr-1
        write(11,'(20e20.10)') xs(i,j),ys(i,j),c(i,j)
    enddo
    write(11,'()')
enddo

do j = ntheta-1,0,-1
    do i = 0,nr-1
        write(11,'(20e20.10)') xs(i,j),-ys(i,j),c(i,j)
    enddo
    write(11,'()')
enddo


close(11)
write(*,*) 'output  : ',fname
