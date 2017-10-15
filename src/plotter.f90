module plotter
contains
  function save_data(q,x,y,n_x,n_y,t,id)
    integer :: n_x,id
    real :: t,gamma=1.4
    real,dimension(n_x,n_y,4) :: q
    real,dimension(n_x,n_y) :: u,p,rho,E,a,M,v
    real, dimension(n_x) :: x
    real, dimension(n_y) :: y
    character(len=5) :: charI
    character(len=1024) :: fname

    rho = q(:,:,1)
    u = q(:,:,2)/rho
    v = q(:,:,3)/rho
    E = q(:,:,4)/rho
    p = (gamma-1)*rho*(E-0.5*(u*u+v*v))
    a = SQRT(gamma*p/rho)
    M = u/a

    write(charI,'(I5)') id
    write(fname,'(a,i4.4,a,f6.4,a)') "./data/frame",id,"t_",t,".dat"


    open(unit = 100,file = fname)
    do j=1,n_y
      do i=1,n_x
        write(100,*)x(i),y(j),u(i,j),v(i,j),p(i,j),rho(i,j)
      end do
    end do
    close(100)

  end function save_data


  function plot_data(q,x,y,n_x,n_y,t,id,val)
    integer :: n_x,id,check,val
    real :: t
    real,dimension(n_x,n_y,4) :: q
    real,dimension(n_x,n_y) :: u,p,rho,E,a,M
    real, dimension(n_x) :: x
    real, dimension(n_y) :: y
    character(len=5) :: charI
    character(len=1024) :: property

    check = save_data(q,x,y,n_x,n_y,t,id)

    select case (val)
      case(1)
        property='u velocity'
      case(2)
        property='v velocity'
      case(3)
        property='pressure'
      case(4)
        property='density'
      case default
        property='y'
    end select

    ! Create the gnuplot file
    open(unit = 100,file='./plots/1plot.py')

    write(100,'(a)') 'import numpy as np'
    write(100,'(a)') 'import matplotlib.pyplot as plt'
    write(100,'(a,i4.4,a,f6.4,a)') "fname = ""./data/frame",id,"t_",t,".dat"""
    write(100,'(a)') 'data = np.loadtxt(fname)'
    write(100,'(a,i4)')'nx = ',n_x
    write(100,'(a,i4)')'ny = ',n_y
    write(100,'(a,i1,a)')'Z = np.reshape(data[:,',val+1,'],(ny,nx))'
    write(100,'(a)') 'X = np.reshape(data[:,0],(ny,nx))'
    write(100,'(a)') 'Y = np.reshape(data[:,1],(ny,nx))'
    write(100,'(3a)') "plt.title('",trim(property),"',loc='left')"
    write(100,'(a,i4.4,a,f6.4,a)') "plt.title('Frame ",id,", time = ",t,"')"
    write(100,'(a)') "plt.xlabel('x')"
    write(100,'(a)') "plt.xlabel('y')"
    write(100,'(a)') 'plt.contourf(Z,100)'
    write(100,'(a,a,i4.4,a)') 'plt.savefig(''./plots/',trim(property),id,'.png'')'

    close(100)

    call system('python ./plots/1plot.py')
  end function plot_data



  function get_video(val)
    integer :: val
    character(len=1024) :: property,command
    select case (val)
      case(1)
        property='velocity'
      case(2)
        property='pressure'
      case(3)
        property='density'
      case(4)
        property='Mach_number'
      case default
        property='y'
    end select
    write(command,'(a,a,a,a,a)') "avconv -i ""./plots/",trim(property),"%04d.png"" -r 30 ./plots/",trim(property),".mp4"
    call system(command)
  end function get_video
end module plotter
