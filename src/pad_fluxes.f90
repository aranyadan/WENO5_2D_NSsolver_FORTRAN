subroutine pad_fluxes(F,Fnew,n_x,n_y,vars)
  implicit none
  real, dimension(n_x,n_y,vars) :: F
  real, dimension(-1:n_x+3,-1:n_y+3,vars) :: Fnew
  integer :: n_x,n_y,vars

  Fnew(1:n_x,1:n_y,:) = F(1:n_x,1:n_y,:)

  Fnew(-1:0,:,:) = spread(Fnew(1,:,:),1,2)
  Fnew(:,-1:0,:) = spread(Fnew(:,1,:),2,2)
  Fnew(n_x+1:n_x+3,:,:) = spread(Fnew(n_x,:,:),1,3)
  Fnew(:,n_y+1:n_y+3,:) = spread(Fnew(:,n_y,:),2,3)
end subroutine pad_fluxes
