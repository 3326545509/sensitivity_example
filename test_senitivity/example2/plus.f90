program main
    implicit none
    
    integer:: imax_v, imax_k
    common imax_v,imax_k

    character *70 :: file_velocity, file_kernel,out_kernel,out_dlnA
    integer::   stat
    real::  slo,sla,rlo,rla,dist,dlnA
    real,allocatable:: velocity(:,:),kernel(:,:)
    real,dimension(3,3)::romax
    integer:: j

    real::t1,t2
    real,external::radian
    call cpu_time(t1)

    !if(1>2)then
    write(*,*) 'please input velocity file name'
    read(*,*) file_velocity
    write(*,*) 'please input kernel file name'
    read(*,*) file_kernel
    write(*,*) 'please input slo,sla,rlo,rla,dist'
    read(*,*) slo,sla,rlo,rla,dist
    write(*,*) 'please input rotation kernel outfile name'
    read(*,*) out_kernel
    write(*,*) 'please input outfile dln(obseve)'
    read(*,*) out_dlnA
    call count_row(file_velocity,imax_v)
    call count_row(file_kernel,imax_k)    
    allocate(velocity(imax_v,3))
    allocate(kernel(imax_k,3))
    !read velocity and kernel files
    call readfile(file_velocity,imax_v,velocity)
    call readfile(file_kernel,imax_k,kernel)
    !kernel rotation
    call zero_array(kernel,imax_k,dist)
    call rotate_matrix(dist,slo,sla,rlo,rla,romax)
    call rotate(kernel,romax)
    open(unit=10,file=out_kernel)
    do j=1,imax_k,1
        write(10,*) kernel(j,1),kernel(j,2),kernel(j,3)
    end do
    close(10)
    !将kernel和phvel相乘相加
    call plus(kernel,velocity,dlnA)
    
    open(unit=10,file=out_dlnA)
    write(10,*)file_kernel,dlnA
    close(10)

    call cpu_time(t2)
    ! write(*,*)'cpu time:',t2-t1
end program main

subroutine plus(kernel,velocity,dlnA)
    implicit none
    integer:: imax_v,imax_k
    common imax_v,imax_k

    real,intent(in):: kernel(imax_k,3), velocity(imax_v,3)
    real,intent(out):: dlnA
    ! ke：单个点的上的kernel值
    real :: dv,ke,vlo,vla,klo,kla
    integer :: i,j
    real:: lomin_kernel_area,lomax_kernel_area,lamin_kernel_area,lamax_kernel_area
    
    dlnA=0

    ! open(20,file='k_v')
    lomin_kernel_area=minval(kernel(:,1))
    lomax_kernel_area=maxval(kernel(:,1))
    lamin_kernel_area=minval(kernel(:,2))
    lamax_kernel_area=maxval(kernel(:,2))
    do i=1,imax_v,1
        vlo =   velocity(i,1)
        vla =   velocity(i,2)
        if (vlo<lomin_kernel_area .or. vlo>lomax_kernel_area .or. vla<lamin_kernel_area .or. vla>lamax_kernel_area)cycle
        dv  =   (velocity(i,3)-3)/3
        do j=1,imax_k,1
            klo =   kernel(j,1)
            kla =   kernel(j,2)
            ke  =   kernel(j,3)
            if (vlo-klo<=0.25 .and. vlo-klo>-0.25 .and. vla-kla<=0.25 .and. vla-kla>-0.25)then
                dlnA = dlnA + ke*dv
                ! write(20,*)ke,dv,ke*dv
            end if
        end do
        !write(*,*)ke,dv
    end do
    ! close(20)
    write(*,*)'dln(obseve)=',dlnA
end subroutine plus

subroutine rotate_matrix(dist,slo,sla,rlo,rla,romax)
    implicit none
    real,external::radian
    real ,intent(in) :: slo,sla,rlo,rla,dist
    real ,intent(out):: romax(3,3)

    real,dimension(3,3):: matrix_z, matrix_y, matrix_x, temp
    real,dimension(3,1):: source, receiver, vector0, vector1, check1, vector_check!vector shi duo yu de
    real :: pi, R, theta, fi, x0, y0, z0, x1, y1, z1, belta,sin_belta,cos_belta

    pi  =  4.0d0* atan(1.0d0)  
    R   =  6371.0d0

    call sphere_to_rectangular(r, pi/2, dist/r, vector0(1,1), vector0(2,1), vector0(3,1))

    theta = pi/2 - radian(rla)
    fi  =   radian(rlo)
    call sphere_to_rectangular(R,theta,fi,receiver(1,1),receiver(2,1),receiver(3,1))
    ! source 是用来check的, 以确保确实把它转到了x轴上
    theta = pi/2 - radian(sla)
    fi  =   radian(slo)
    call sphere_to_rectangular(R,theta,fi,source(1,1),source(2,1),source(3,1))

    call rotate_matrix_z(-radian(slo), matrix_z)
    call rotate_matrix_y(radian(sla), matrix_y )

    vector1 = matmul(matrix_y, matmul(matrix_z,receiver))
    y0 = vector0(2,1)
    z0 = vector0(3,1)   !z0 ~= 0
    y1 = vector1(2,1)
    z1 = vector1(3,1)

    cos_belta=(y0/z1) / ( y1/z1 + z1/y1 )
    sin_belta=y0 / (-y1**2/z1 -z1)

    call rotate_matrix_x(cos_belta, sin_belta, matrix_x)

    temp = matmul(matrix_x, matmul(matrix_y,matrix_z))
    call inv(temp, romax)
end subroutine rotate_matrix

subroutine rotate_matrix_z(belta,matrix_z)
    implicit none
    real,intent(in)::belta
    real,dimension(3,3)::matrix_z

    matrix_z=reshape((/ cos(belta), sin(belta), 0.0, -sin(belta), cos(belta), 0.0, 0.0, 0.0, 1.0 /),(/3,3/))
end subroutine rotate_matrix_z

subroutine rotate_matrix_y(belta,matrix_y)
    implicit none
    real,intent(in)::belta
    real,dimension(3,3)::matrix_y

    matrix_y=reshape((/ cos(belta), 0.0, -sin(belta), 0.0, 1.0, 0.0, sin(belta), 0.0, cos(belta) /),(/3,3/))
end subroutine rotate_matrix_y

subroutine rotate_matrix_x(cos_belta,sin_belta,matrix_x)
    implicit none
    real,intent(in)::cos_belta,sin_belta
    real,dimension(3,3)::matrix_x

    matrix_x=reshape((/ 1.0, 0.0, 0.0, 0.0, cos_belta, sin_belta, 0.0, -sin_belta, cos_belta /),(/3,3/))
end subroutine rotate_matrix_x

subroutine inv(a,e)
    implicit none
    real,dimension(3,3),intent(in)::a
    real,dimension(3,3),intent(out)::e  ! output inverse(a)
    real,dimension(3,3)::m  !temp copy a 
    real::temp
    integer::i,j,k,n    ! n is dimension  of matrix which can be changed to suit different matrix

    m=a
    !unit matrix
    e=reshape((/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/),(/3,3/))

    n=3
    do i=1,n,1
        temp=m(i,i)
        do j=1,n,1
            m(i,j)=m(i,j)/temp
            e(i,j)=e(i,j)/temp
        end do
        
        do k=i+1,n,1
            temp=m(k,i)
            do j=1,n,1
                m(k,j)=m(k,j)-m(i,j)*temp
                e(k,j)=e(k,j)-e(i,j)*temp
            end do
        end do
    end do
    
    do i=1,n-1,1
        do k=i+1,n,1
            temp=m(i,k)
            do j=1,n,1
                m(i,j)=m(i,j)-m(k,j)*temp
                e(i,j)=e(i,j)-e(k,j)*temp
            end do
        end do
    end do

end subroutine inv

subroutine rotate(kernel,romax)
    implicit none
    integer :: imax_v,imax_k
    common imax_v,imax_k
    real,dimension(imax_k,3)::  kernel
    real,dimension(3,3),intent(in):: romax


    real :: x0,y0,z0,theta0,fi0,theta1,fi1,R,pi
    real,dimension(3,1)::vector0,vector1
    integer :: i 

    pi  =   4.* atan(1.)  
    R   =   6371.0

    do i=1,imax_k,1
        fi0     =   kernel(i,1)*pi/180.0
        theta0  =   pi/2 - kernel(i,2)*pi/180.0
        call sphere_to_rectangular(R,theta0,fi0,x0,y0,z0)

        vector0 = reshape((/x0,y0,z0/),(/3,1/))
        vector1 = matmul(romax,vector0)

        call rectangular_to_sphere(vector1(1,1),vector1(2,1),vector1(3,1),r,theta1,fi1)

        kernel(i,1) =   fi1*180.0/pi
        kernel(i,2) =   (pi/2-theta1)*180.0/pi
    end do
end subroutine rotate

subroutine zero_array(array,imax,dist)
    implicit none
    integer,intent(in)::    imax
    real,intent(in)::dist
    real,dimension(imax,3)::array
    integer::i

    do i=1,imax,1
        if( abs(array(i,3))>100 ) then
            array(i,3)=0
        end if 
        ! or. array(i,1)>arcdist+0.5 .or. array(i,1)<-0.5 ) then
    end do
end subroutine zero_array

subroutine readfile(filein,imax,array)
    implicit none
    character(len=*),intent(in)::  filein
    integer,intent(in)::    imax
    real,dimension(imax,3),intent(out)::  array
    integer::   stat, i

    open(unit=40,file=filein)
    do i =1,imax,1
        read(40,*,iostat=stat) array(i,1),array(i,2),array(i,3)
    end do
    close(unit=40)
end subroutine readfile

subroutine count_row(filein,imax)
    implicit none
    character(len=*),intent(in)::    filein
    integer,intent(out)::   imax
    integer::   stat
    !temp
    real::      a,b,c

    open(unit=40,file = filein)
    imax=0
    do while(.true.)
        read(40,*,iostat=stat)a,b,c
        if (stat/=0) exit
        imax=imax+1
    end do
    close(unit=40)
end subroutine count_row

subroutine rectangular_to_sphere(x,y,z,r,theta,fi)
    implicit none
    real,intent(in)::   x,y,z
    real,intent(out)::  r,theta,fi

    real :: pi
    pi = 4.0*atan(1.0)

    r   =   sqrt(x**2+y**2+z**2)
    theta = atan(sqrt(x**2+y**2)/z)
    ! theta in coordinate is 0 to pi, arctan is -pi/2 to pi/2
    if (theta<0) theta=pi+theta
    ! longitude is 0 to 360
    fi  =   atan(y/x)
    if ( x>0 .and. y>0 ) fi = fi
    if ( x>0 .and. y<0 ) fi = fi + 2*pi
    if ( x<0 .and. y>0 ) fi = fi + pi
    if ( x<0 .and. y<0 ) fi = fi + pi

end subroutine rectangular_to_sphere

subroutine sphere_to_rectangular(r,theta,fi,x,y,z)
    implicit none
    real,intent(in):: r,theta,fi
    real,intent(out)::x,y,z

    x   =   r*sin(theta)*cos(fi)
    y   =   r*sin(theta)*sin(fi)
    z   =   r*cos(theta)
end subroutine sphere_to_rectangular

function radian(degree)
    implicit none
    real :: degree
    real :: radian
    radian = degree *3.141593d0/180.0d0
    return
end function radian
