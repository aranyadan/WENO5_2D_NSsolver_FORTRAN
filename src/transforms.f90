module transform
  implicit none
  integer :: unitR=0

contains
  function RFinv(u,v,a)
    real,dimension(4,4) :: RFinv
    real :: u,v,a,q2
    real :: gamma = 1.4

    q2=0.5*(u*u+v*v)

    RFinv(1,1) = ((gamma-1.0)*q2 + a*u)/(2.0*a*a)
    RFinv(2,1) = (a*a - (gamma-1.0)*q2)/(a*a)
    RFinv(3,1) = ((gamma-1.0)*q2 - a*u)/(2.0*a*a)
    RFinv(4,1) = v

    RFinv(1,2) = ((1.0-gamma)*u - a)/(2.0*a*a)
    RFinv(2,2) = ((gamma-1.0)*u)/(a*a)
    RFinv(3,2) = ((1.0-gamma)*u + a)/(2.0*a*a)
    RFinv(4,2) = 0.0

    RFinv(1,3) = ((1.0-gamma)*v)/(2.0*a*a)
    RFinv(2,3) = ((gamma-1.0)*v)/(a*a)
    RFinv(3,3) = ((1.0-gamma)*v)/(2.0*a*a)
    RFinv(4,3) = -1.0

    RFinv(1,4) = (gamma-1.0)/(2.0*a*a)
    RFinv(2,4) = (1.0-gamma)/(a*a)
    RFinv(3,4) = (gamma-1.0)/(2.0*a*a)
    RFinv(4,4) = 0.0

    if(unitR==1) then
      RFinv(1,:) = (/1,0,0,0/) !
      RFinv(2,:) = (/0,1,0,0/) !
      RFinv(3,:) = (/0,0,1,0/) !
      RFinv(4,:) = (/0,0,0,1/) !
    endif

  end function RFinv

  function RFcalc(u,v,a)
    real,dimension(4,4) :: RFcalc
    real :: gamma = 1.4
    real :: u,v,a,q2

    q2=0.5*(u*u+v*v)

    RFcalc(1,1) = 1.0
    RFcalc(2,1) = u-a
    RFcalc(3,1) = v
    RFcalc(4,1) = a*a/(gamma - 1) + q2 - u*a

    RFcalc(1,2) = 1.0
    RFcalc(2,2) = u
    RFcalc(3,2) = v
    RFcalc(4,2) = q2

    RFcalc(1,3) = 1.0
    RFcalc(2,3) = u+a
    RFcalc(3,3) = v
    RFcalc(4,3) = a*a/(gamma - 1) + q2 + u*a

    RFcalc(1,4) = 0.0
    RFcalc(2,4) = 0.0
    RFcalc(3,4) = -1.0
    RFcalc(4,4) = -v

    if(unitR==1) then
      RFcalc(1,:) = (/1,0,0,0/) !
      RFcalc(2,:) = (/0,1,0,0/) !
      RFcalc(3,:) = (/0,0,1,0/) !
      RFcalc(4,:) = (/0,0,0,1/) !
    endif

  end function RFcalc


  function RGinv(u,v,a)
    real,dimension(4,4) :: RGinv
    real :: u,v,a,q2
    real :: gamma = 1.4

    q2 = 0.5*(u*u+v*v)

    RGinv(1,1) = ((gamma-1.0)*q2 + a*v)/(2.0*a*a)
    RGinv(2,1) = (a*a - (gamma-1.0)*q2)/(a*a)
    RGinv(3,1) = ((gamma-1.0)*q2 - a*v)/(2.0*a*a)
    RGinv(4,1) = -u

    RGinv(1,2) = ((1.0-gamma)*u)/(2.0*a*a)
    RGinv(2,2) = ((gamma-1.0)*u)/(a*a)
    RGinv(3,2) = ((1.0-gamma)*u)/(2.0*a*a)
    RGinv(4,2) = 1.0

    RGinv(1,3) = ((1.0-gamma)*v - a)/(2.0*a*a)
    RGinv(2,3) = ((gamma-1.0)*v)/(a*a)
    RGinv(3,3) = ((1.0-gamma)*v + a)/(2.0*a*a)
    RGinv(4,3) = 0.0


    RGinv(1,4) = (gamma-1.0)/(2.0*a*a)
    RGinv(2,4) = (1.0-gamma)/(a*a)
    RGinv(3,4) = (gamma-1.0)/(2.0*a*a)
    RGinv(4,4) = 0.0

    if(unitR==1) then
      RGinv(1,:) = (/1,0,0,0/) !
      RGinv(2,:) = (/0,1,0,0/) !
      RGinv(3,:) = (/0,0,1,0/) !
      RGinv(4,:) = (/0,0,0,1/) !
    endif

  end function RGinv


  function RGcalc(u,v,a)
    real,dimension(4,4) :: RGcalc
    real :: gamma = 1.4
    real :: u,v,a,q2

    q2 = 0.5*(u*u+v*v)

    RGcalc(1,1) = 1.0
    RGcalc(2,1) = u
    RGcalc(3,1) = v-a
    RGcalc(4,1) = a*a/(gamma - 1) + q2 - v*a

    RGcalc(1,2) = 1.0
    RGcalc(2,2) = u
    RGcalc(3,2) = v
    RGcalc(4,2) = q2

    RGcalc(1,3) = 1.0
    RGcalc(2,3) = u
    RGcalc(3,3) = v+a
    RGcalc(4,3) = a*a/(gamma - 1) + q2 + v*a

    RGcalc(1,4) = 0.0
    RGcalc(2,4) = 1.0
    RGcalc(3,4) = 0.0
    RGcalc(4,4) = u

    if(unitR==1) then
      RGcalc(1,:) = (/1,0,0,0/) !
      RGcalc(2,:) = (/0,1,0,0/) !
      RGcalc(3,:) = (/0,0,1,0/) !
      RGcalc(4,:) = (/0,0,0,1/) !
    endif

  end function RGcalc
end module transform
