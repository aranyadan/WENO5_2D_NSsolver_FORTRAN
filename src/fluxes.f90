module flux
contains
  function turn(F,n_x,n_y,times,direction)
    integer :: n_x,times,direction
    real,dimension(n_x,n_y,4) :: F,F_new,turn

    if(times == 0) then
      turn = F
    else
      if(direction==1) then
        do i = 1,n_x
          if (mod(i+times,n_x)>0) then
            F_new(mod(i+times,n_x),:,:) = F(i,:,:)
          else
            F_new(n_x,:,:) = F(i,:,:)
          end if
        end do
        turn = F_new
      else
        do i = 1,n_y
          if (mod(i+times,n_y)>0) then
            F_new(:,mod(i+times,n_y),:) = F(:,i,:)
          else
            F_new(:,n_y,:) = F(:,i,:)
          end if
        end do
        turn = F_new
      end if
    end if

    return
  end function turn

  function turnx(F,n_x,times)
    integer :: n_x,times
    real,dimension(n_x,4) :: F,F_new,turn
    if(times == 0) then
      turn = F
    else
      do i = 1,n_x
        if (mod(i+times,n_x)>0) then
          F_new(mod(i+times,n_x),:) = F(i,:)
        else
          F_new(n_x,:) = F(i,:)
        end if
      end do
      turn = F_new
    end if
    return
  end function turnx

  subroutine build_flux(q,n_x,n_y,F,G)
    integer :: n_x,n_y
    real :: gamma=1.4
    real, dimension(n_x,n_y,4) :: q,F, G
    real, dimension(n_x,n_y) :: u,v,p,rho,E,a

    call primitives(q,n_x,n_y,rho,u,v,E,p,a)

    F(:,:,1) = rho*u
    F(:,:,2) = rho*u*u + p
    F(:,:,3) = rho*u*v
    F(:,:,4) = rho*u*E + p*u

    G(:,:,1) = rho*v
    G(:,:,2) = rho*u*v
    G(:,:,3) = rho*v*v +p
    G(:,:,4) = rho*v*E + p*v

  end subroutine build_flux

  ! function get_deriv(lambda,hp,hn,q,n_x,dx)
  !   integer :: n_x
  !   real :: lambda,dx
  !   real, dimension(n_x,3) :: hp,hn,res1,res2,res,get_deriv,q
  !   res1 = 0.5 * ( build_flux(hp,n_x) + build_flux(hn,n_x) - lambda*(hn - hp) )
  !   res2 = 0.5 * ( build_flux((turn(hp,n_x,1)),n_x) + build_flux((turn(hn,n_x,1)),n_x) - lambda*( turn(hn,n_x,1) - turn(hp,n_x,1)) )
  !   res = (res1-res2)/dx
  !
  !   get_deriv= res
  !
  ! end function get_deriv
end module flux
