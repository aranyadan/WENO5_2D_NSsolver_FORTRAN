! Main file
program main
  use flux
  use plotter
  use transform
  implicit none
  ! Plot values: 1=velocity, 2=pressure, 3=density 4=Mach number
  integer,parameter :: n_x = 101, n_y =101, SAVE=1,PLOT=1,PLOTVAL=4,VIDEO=0
  real, parameter :: startx = 0, endx = 1, starty = 0, endy = 1, gamma = 1.4
  real :: delx,dely,dt,cfl,tend,lambda_0,lambda,t,dt_0
  integer :: I,id=0,check,case_id
  real,dimension(n_x,n_y) :: a_0,p,rho,u,v,E,a                             ! Stores x coordinate of the points, primitive values
  real,dimension(n_x) :: x
  real,dimension(n_y) :: y
  real, dimension(n_x,n_y,4) :: Prim,Prim_0,q_0,q,qo,dF,hp,hn,g,f,dL,dG       ! Stores primitive values and flux values

  delx = abs(endx-startx)/(n_x-1)
  dely = abs(endy-starty)/(n_y-1)
  cfl = 0.475
  tend = 0.3

  x = (/ (startx + (I-1)*delx,I = 1,n_x) /)
  y = (/ (starty + (I-1)*dely,I = 1,n_y) /)

  print*,'Setting initial conditions'
  case_id = 3
  call IC2DReimann(Prim_0,q_0,n_x,n_y,x,y,case_id,cfl,tend)
  print*,'Initial conditions Set'

  q = q_0
  a_0 = SQRT(gamma*Prim_0(:,:,3)/Prim_0(:,:,4))
  lambda_0 = MAXVAL(MAXVAL( ABS( Prim_0(:,:,1) )+a_0,1))
  dt_0 = cfl * delx/lambda_0

  !! Solver Loop
  Prim = Prim_0
  q = q_0
  t=0
  dt = dt_0
  lambda = lambda_0
  if(PLOT==1) then
    check=plot_data(q,x,y,n_x,n_y,t,id,PLOTVAL)
    id=id+1
  else if(SAVE==1)then
    check=save_data(q,x,y,n_x,n_y,t,id)
    id=id+1
  end if

  print*,'Starting time stepping'

  do while (t < tend)
    ! Starting RK
    qo = q

    ! RK 1st step
    call build_flux(q,n_x,n_y,f,g)
    call WENO52d(lambda,f,q,n_x,n_y,hp,hn,1)
    dF = ((hp - turn(hp,n_x,n_y,1,1)) + (hn - turn(hn,n_x,n_y,1,1)))/delx
    call WENO52d(lambda,g,q,n_x,n_y,hp,hn,2)
    dG = ((hp - turn(hp,n_x,n_y,1,2)) + (hn - turn(hn,n_x,n_y,1,2)))/dely
    dL = dF + dG

    q = qo - dt*dL

    q(1,:,:)=q(2,:,:)
    q(n_x,:,:)=q(n_x-1,:,:)
    q(:,1,:)=q(:,2,:)
    q(:,n_y,:)=q(:,n_y-1,:)



    ! RK 2nd step
    call build_flux(q,n_x,n_y,f,g)
    call WENO52d(lambda,f,q,n_x,n_y,hp,hn,1)
    dF = ((hp - turn(hp,n_x,n_y,1,1)) + (hn - turn(hn,n_x,n_y,1,1)))/delx
    call WENO52d(lambda,g,q,n_x,n_y,hp,hn,2)
    dG = ((hp - turn(hp,n_x,n_y,1,2)) + (hn - turn(hn,n_x,n_y,1,2)))/dely
    dL = dF + dG

    q = 0.75*qo + 0.25*( q - dt*dL)

    q(1,:,:)=q(2,:,:)
    q(n_x,:,:)=q(n_x-1,:,:)
    q(:,1,:)=q(:,2,:)
    q(:,n_y,:)=q(:,n_y-1,:)

    ! RK 3rd step
    call build_flux(q,n_x,n_y,f,g)
    call WENO52d(lambda,f,q,n_x,n_y,hp,hn,1)
    dF = ((hp - turn(hp,n_x,n_y,1,1)) + (hn - turn(hn,n_x,n_y,1,1)))/delx
    call WENO52d(lambda,g,q,n_x,n_y,hp,hn,2)
    dG = ((hp - turn(hp,n_x,n_y,1,2)) + (hn - turn(hn,n_x,n_y,1,2)))/dely
    dL = dF + dG

    q = (qo + 2.0*( q - dt*dL))/3.0

    q(1,:,:)=q(2,:,:)
    q(n_x,:,:)=q(n_x-1,:,:)
    q(:,1,:)=q(:,2,:)
    q(:,n_y,:)=q(:,n_y-1,:)
    ! Extract primitive values
    call primitives(q,n_x,n_y,rho,u,v,E,p,a)

    lambda = MAXVAL(MAXVAL(ABS(u)+a,1))
    dt = cfl*delx/lambda
    if(t+dt>tend) then
      dt = tend-t
    end if
    t=t+dt
    if(PLOT==1) then
      check=plot_data(q,x,y,n_x,n_y,t,id,PLOTVAL)
      id=id+1
    else if(SAVE==1)then
      check=save_data(q,x,y,n_x,n_y,t,id)
      id=id+1
    end if
    if(MOD(id,20)==0) then
      write( *, '(a,i4,a,f6.4)')'Writing file id: ', id,'    Time = ',t
    end if

  end do
if(VIDEO==1) then
  check=get_video(PLOTVAL)
end if
end program main
