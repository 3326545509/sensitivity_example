program main
    implicit none
    real(kind=8), external :: radian
!100 format(f20.15)
	integer,parameter::maxnfreq = 20,maxix = 4000,maxiy = 4000
	real(kind=8),parameter::radius=6371.0d0
    real(kind=8)::phsens(maxix,maxiy),      ampsens(maxix,maxiy)
	real(kind=8)::avgphsens(maxix,maxiy),   avgampsens(maxix,maxiy)
	real(kind=8)::wgttemp(maxix,maxiy)
	real(kind=8)::freq(maxnfreq),amplitude(maxnfreq)
    ! data uesed to compute
	real(kind=8):: slo,sla,rlo,rla,dist,phvel,scalelen
    real(kind=8):: pi,sumamp,kk,lamda,period,delta,delta1,delta2
    ! variable quantity grid
    integer::nfreq,ifreq,ix,iy
    real(kind=8)::xbeg,xend,ybeg,yend,dx,dy,x,y,ds
    integer::xout,yout,nx,ny
	character(len=70)::outfn1,outfn2,spcfn

	write(*,*) 'please input velocity'
	read(*,*) phvel
	write(*,*) 'please input the spectral file'
	read(*,*) spcfn
	write(*,*) 'please input outfn of amplitude sens'
	read(*,*) outfn1
	write(*,*) 'please input outfn of phase sens'
	read(*,*) outfn2
    ! source and station
	write(*,*) 'please input dist form source to receiver(km)'
	read(*,*) dist

	slo=0.0d0
	sla=0.0d0
	rlo=dist/radius*180.0d0/(4.0d0* atan(1.0d0))
	rla=0.0d0


	open(10,file = spcfn)   !spetral file 
	open(20,file = outfn1)  !unsmooth kernel outputed
    !open(80,file = 'check')
	open(30,file = outfn2)

	xbeg = -1.50d0
	xend = nint(rlo*100.0d0)/100.0d0+1.50d0
	ybeg = -2.0d0
	yend =  2.0d0
	dx   =  0.01d0
	dy   =  0.01d0
	
    !====check============
    xbeg=-15
    xend=15
    ybeg=-15
    yend=15
    dx = 0.1d0
    dy = 0.1d0
    !=====================

	nx = nint((xend-xbeg)/dx + 1.0)
	ny = nint((yend-ybeg)/dy + 1.0)

    if(nx>maxix .or. ny>maxiy) then
        write(*,*)"!!!!Error, nx>maxix"
        stop
    end if

    pi = 4.0d0*atan(1.0d0)
    ! read spectral file. normalize the amplitude.
    sumamp=0.0d0
    read(10,*)nfreq
    do ifreq=1,nfreq,1
        read(10,*)freq(ifreq),amplitude(ifreq)
        sumamp = sumamp + amplitude(ifreq)
    end do
    close(10)
    do ifreq=1,nfreq,1
        amplitude(ifreq) = amplitude(ifreq)/sumamp
    end do

    ! output into check file
    ! write(80,*)'freq and amplitude before normalization'
    ! do ifreq =1,nfreq,1
    !    write(80,*)freq(ifreq),amplitude(ifreq)
    ! end do

    ! initialization
    phsens=0.0d0
    avgphsens=0.0d0
    ampsens=0.0d0
    avgampsens=0.0d0
    !write(80,*)'ampsens init'
    !write(80,*)ampsens

    ! compute different freq, add them together.
    call lola_delta(slo,sla,rlo,rla,delta) ! epicentral distance
    do ifreq=1,nfreq,1
        period = 1.0d0/freq(ifreq)
        lamda  = phvel*period/radius
        kk     = 2.0d0*pi/lamda
        do ix=1,nx,1
            x = (ix-1)*dx + xbeg
            do iy=1,ny,1
                y = (iy-1)*dy + ybeg
                !ds=radius**2*sin(radian(90.0d0-y))*radian(dx)*radian(dy), kilometers
                ds=sin(radian(90.0d0-y))*radian(dx)*radian(dy)
                call lola_delta(slo,sla,x,y,delta1)
                call lola_delta(rlo,rla,x,y,delta2)
                ! when denominator<1e-8, ampsens+0
                if (abs(sqrt(8*pi*kk*abs( sin(delta1)*sin(delta2)/sin(delta) )))<1e-8) then
                    ampsens(ix,iy) = ampsens(ix,iy) + 0
                    phsens(ix,iy)  = phsens(ix,iy) + 0
                    !write(80,*)abs(sqrt(8*pi*kk*abs( sin(delta1)*sin(delta2)/sin(delta) ))),'is zero'
                ! set kernel near source and receiver to 0
                else if (radian(sqrt((x-slo)**2+(y-sla)**2)).le.lamda .or. radian(sqrt((x-rlo)**2+(y-rla)**2)).le.lamda) then
                    ampsens(ix,iy) = ampsens(ix,iy) + 0
                    phsens(ix,iy)  = phsens(ix,iy) + 0
                    !write(80,*)x,y,'is near source or receiver'
                else if (kk*(delta1+delta2-delta)+pi/4>4*pi) then
                    ampsens(ix,iy) = ampsens(ix,iy) + 0
                    phsens(ix,iy)  = phsens(ix,iy) + 0
                else
                    ! write(80,*)'fenmu:',abs(sqrt(8*pi*kk*abs( sin(delta1)*sin(delta2)/sin(delta) )))
                    ampsens(ix,iy) = ampsens(ix,iy) + &
                        amplitude(ifreq)*(-2)*kk**2*cos(kk*(delta1+delta2-delta)+pi/4) &
                        /sqrt(8*pi*kk*abs( sin(delta1)*sin(delta2)/sin(delta) ))/6371**2

                    phsens(ix,iy) = phsens(ix,iy) + &
                        amplitude(ifreq)*ds*(-2)*kk**2*sin(kk*(delta1+delta2-delta)+pi/4) &
                        /sqrt(8*pi*kk*abs( sin(delta1)*sin(delta2)/sin(delta) ))
                end if
                ! write(80,100)x,y,delta1**2+delta2**2
            end do
        end do
    end do

    do ix=1,nx,1
        x=(ix-1)*dx + xbeg
        do iy=1,ny,1
            y = (iy-1)*dy + ybeg
            write(20,*)x,y,ampsens(ix,iy)
            write(30,*)x,y,phsens(ix,iy)
        end do
    end do 
    !=====================
    close(20)
    close(30)
    close(80)
end program main

function radian(degree)
    implicit none
    real(kind=8) :: radian
    real(kind=8) :: degree,pi
    pi = 4.0d0*atan(1.0d0)
    radian = degree *pi/180.0d0
    return
end function radian

subroutine lola_delta(lo1,la1,lo2,la2,delta)
    implicit none
    real(kind=8),external::radian
    real(kind=8),intent(in)::   lo1,la1,lo2,la2
    real(kind=8),intent(out)::  delta

    real(kind=8):: lo1_temp, la1_temp, lo2_temp, la2_temp, aa_temp, d_lo, d_la
    lo1_temp = radian(lo1)
    la1_temp = radian(la1)
    lo2_temp = radian(lo2)
    la2_temp = radian(la2)

    d_lo  = lo2_temp - lo1_temp
    d_la  = la2_temp - la1_temp

    aa_temp = sin(d_la/2)**2 + cos(la1_temp)*cos(la2_temp)*sin(d_lo/2)**2
    delta = 2*asin(sqrt(aa_temp))
end subroutine lola_delta
