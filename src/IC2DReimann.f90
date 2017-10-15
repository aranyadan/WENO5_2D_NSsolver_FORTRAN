subroutine IC2DReimann(Prim,q,n_x,n_y,x,y,case_id,tend)
  implicit none
  integer :: n_x,n_y,mid_x,mid_y,i,j, case_id
  real :: delx,del_y,gamma=1.4,tend,cfl
  real, dimension(4) :: p,u,v,rho
  real,dimension(n_x,n_y,4) :: Prim,q
  real, dimension(n_x) :: x
  real, dimension(n_y) :: y
  real,dimension(n_x,n_y) :: E


  select case (case_id)
  case (1)  ! config 1
      p =   (/ 1.0, 0.4,     0.0439,  0.15    /)
      rho = (/ 1.0, 0.5197,  0.1072,  0.2579  /)
      u =   (/ 0.0, -0.7259, -0.7259, 0.0     /)
      v =   (/ 0.0, 0.0,     -1.4045, -1.4045 /)
      tend = 0.2
    case (2)
      p =   (/ 1.0, 0.4,     1.0,     0.4     /)
      rho = (/ 1.0, 0.5197,  1.0,     0.5197  /)
      u =   (/ 0.0, -0.7259, -0.7259, 0.0     /)
      v =   (/ 0.0, 0.0,     -0.7259, -0.7259 /)
      tend = 0.2
    case (3)
      p =   (/ 1.5, 0.3,    0.029, 0.3    /)
      rho = (/ 1.5, 0.5323, 0.138, 0.5323 /)
      u =   (/ 0.0, 1.206,  1.206, 0.0    /)
      v =   (/ 0.0, 0.0,    1.206, 1.206  /)
      tend = 0.3
    case (4)
      p =   (/ 1.1, 0.35,   1.1,    0.35   /)
      rho = (/ 1.1, 0.5065, 1.1,    0.5065 /)
      u =   (/ 0.0, 0.8939, 0.8939, 0.0    /)
      v =   (/ 0.0, 0.0,    0.8939, 0.8939 /)
      tend = 0.25
    case (5)
      p =   (/ 1.0,    1.0,   1.0,  1.0  /)
      rho = (/ 1.0,    2.0,   1.0,  3.0  /)
      u =   (/ -0.75, -0.75,  0.75, 0.75 /)
      v =   (/ -0.5,   0.5,   0.5, -0.5  /)
      tend = 0.23
    case (6)
      p =   (/  1.0,  1.0,   1.0,   1.0  /)
      rho = (/  1.0,  2.0,   1.0,   3.0  /)
      u =   (/  0.75, 0.75, -0.75, -0.75 /)
      v =   (/ -0.5,  0.5,   0.5,  -0.5  /)
      tend = 0.3
    case (7)
      p =   (/ 1.0,  0.4,    0.4,  0.4    /)
      rho = (/ 1.0,  0.5197, 0.8,  0.5197 /)
      u =   (/ 0.1, -0.6259, 0.1,  0.1    /)
      v =   (/ 0.1,  0.1,    0.1, -0.6259 /)
      tend = 0.25
    case (8)
      p =   (/ 0.4,     1.0,    1.0,  1.0    /)
      rho = (/ 0.5197,  1.0,    0.8,  1.0    /)
      u =   (/ 0.1,    -0.6259, 0.1,  0.1    /)
      v =   (/ 0.1,     0.1,    0.1, -0.6259 /)
      tend = 0.25
    end select


  mid_x =(n_x+1)/2
  mid_y =(n_y+1)/2




  ! Set the primitive variables

  ! Quadrant 1
  Prim(mid_x+1:n_x, mid_y+1:n_y, 1) = u(1)
  Prim(mid_x+1:n_x, mid_y+1:n_y, 2) = v(1)
  Prim(mid_x+1:n_x, mid_y+1:n_y, 3) = p(1)
  Prim(mid_x+1:n_x, mid_y+1:n_y, 4) = rho(1)

  ! Quadrant 2
  Prim(1:mid_x, mid_y+1:n_y, 1) = u(2)
  Prim(1:mid_x, mid_y+1:n_y, 2) = v(2)
  Prim(1:mid_x, mid_y+1:n_y, 3) = p(2)
  Prim(1:mid_x, mid_y+1:n_y, 4) = rho(2)

  ! Quadrant 3
  Prim(1:mid_x, 1:mid_y, 1) = u(3)
  Prim(1:mid_x, 1:mid_y, 2) = v(3)
  Prim(1:mid_x, 1:mid_y, 3) = p(3)
  Prim(1:mid_x, 1:mid_y, 4) = rho(3)

  ! Quadrant 4
  Prim(mid_x+1:n_x, 1:mid_y, 1) = u(4)
  Prim(mid_x+1:n_x, 1:mid_y, 2) = v(4)
  Prim(mid_x+1:n_x, 1:mid_y, 3) = p(4)
  Prim(mid_x+1:n_x, 1:mid_y, 4) = rho(4)


  do j = 1,n_y
    do i = 1,n_x
      q(i,j,1) = Prim(i,j,4)
      q(i,j,2) = Prim(i,j,1) * Prim(i,j,4)
      q(i,j,3) = Prim(i,j,2) * Prim(i,j,4)
      E(i,j) = Prim(i,j,3) / ( (gamma-1)*Prim(i,j,4) ) &
               + 0.5 * (Prim(i,j,1)**2 + Prim(i,j,2)**2)
      q(i,j,4) = E(i,j) * Prim(i,j,4)
    end do
  end do
  ! Check IC setup
  ! open(unit=48,file='data')
  ! do i=1,n_x
  !   write(48,*)x(i),F(i,3)
  ! end do
  ! close(48)
  ! call system('gnuplot -p plotter.plt')
  return
end subroutine IC2DReimann
