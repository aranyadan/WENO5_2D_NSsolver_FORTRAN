subroutine viscous_fluxes(vel,p,rho,n_x,n_y,delx,dely,Suth,Re,Pr,tauxx,tauyy,tauxy,q_x,q_y)
  implicit none
  real, dimension(n_x,n_y,3) :: vel     ! 1-xvel, 2-yvel, 3-Temp
  real, dimension(n_x,n_y) :: p,rho,tauxx,tauyy,tauxy,q_x,q_y
  real, dimension(-1:n_x+3,-1:n_y+3,3) :: velnew
  real :: delx,dely,Suth,R_gas_const = 287.0, T,meu,lambda,Re,Pr,Therm_const
  real :: u_x,u_y,v_x,v_y,T_x,T_y
  real, dimension(1,3) :: temp_der,temp_der2
  real, dimension(6,3) :: temp,temp2
  real, parameter :: gamma=1.4
  integer :: n_x,n_y,i,j

  Therm_const = gamma / ((gamma-1.0)*Pr*Re)

  call pad_fluxes(vel,velnew,n_x,n_y,3)
  ! Calculate the derivatives in x and y directions
  do j=1,n_y
    do i =1,n_x

      temp(1:6,:) = velnew(i-2:i+3,j,:)
      call get_first_deriv(temp,delx,temp_der)
      u_x = temp_der(1,1); v_x = temp_der(1,2); T_x = temp_der(1,3)

      temp(1:6,:) = velnew(i,j-2:j+3,:)
      call get_first_deriv(temp,dely,temp_der)
      u_y = temp_der(1,1); v_y = temp_der(1,2); T_y = temp_der(1,3)

      T = vel(i,j,3)
      meu = (T**(3.0/2.0)) * ( (1+Suth)/(T+Suth) )
      lambda = (-2.0/3.0)*meu

      tauxx(i,j) = (lambda*(u_x+v_y) + 2*meu*u_x)/Re
      tauyy(i,j) = (lambda*(u_x+v_y) + 2*meu*v_y)/Re
      tauxy(i,j) = (meu*(u_y+v_x))/Re
      q_x(i,j) = T_x * Therm_const
      q_y(i,j) = T_y * Therm_const

    end do
  end do


end subroutine viscous_fluxes


subroutine get_first_deriv(vel,delta,deriv)
  implicit none
  real, dimension(6,3) :: vel
  real :: delta
  real, dimension(1,3) :: deriv

  deriv(1,:) = ( (1.0/20.0)*vel(1,:) + (-0.5)*vel(2,:) + (-1.0/3.0)*vel(3,:) &
               +	vel(4,:) + (-0.25)*vel(5,:) + (1.0/30.0)*vel(6,:) )/delta

end subroutine get_first_deriv
