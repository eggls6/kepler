module kepler_m

  use orbdyn_m, only:rv2kema,kema2rv,semia2n,epoch2ma,com2ke,pi,pix2,deg2rad,rad2deg,dist3d,lighttraveltime
  use trafo3_m, only:icrf2he,lat2cart,cart2lat,he2icrf,icrf2radec,rvhe2icrf,rvicrf2dradec,norm,angv1v2,angv1v2n
  
  contains
!*******************************************************
  subroutine two_body_rv(rv_in,epoch,t_obs,mass1,mass2,rv_out)
  implicit none

  !input
  real*8,dimension(1:6),intent(in)::rv_in           !heliocentric state x,y,z,vx,vy,vz [au,au/day]
  real*8,intent(in)::epoch                                 !epoch of initial conditions [JD]
  real*8,intent(in)::t_obs                                 !time to propagate to [JD]
  real*8,intent(in)::mass1                                    !mass of the main body, e.g. Sun [Msun]
  real*8,intent(in)::mass2                                    !mass of the second body, e.g. asteroid [Msun]

  !ouput
  real*8,dimension(1:6),intent(out)::rv_out       !heliocentric state x,y,z,vx,vy,vz [au,au/day]
    
  !local
  real*8::meanmotion                            !mean motion [rad/day]
  real*8,dimension(1:6)::ele                      !heliocentric Keplerian elements (a,e,i,w,Om,M) [au,-,4 x rad]
  real*8::ma_obs

  !transform heliocentric, ecliptic, Cartesian states to heliocentric Keplerian orbital elements
  !open(21,file='tst.txt',status='replace')
  !write(21,*)'entered',rv_in
  !ele=1.d0
  call rv2kema(rv_in,ele,mass1,mass2)
  !write(21,*)'rv2kema',meanmotion
  ele(3:6)=ele(3:6)*deg2rad  
  !write(21,*)'rv2kema',ele
  call semia2n(ele(1),mass1,mass2,meanmotion)
  !calculate new mean anomaly
  call epoch2ma(meanmotion,epoch,ele(6),t_obs,ma_obs)
  ele(6)=ma_obs
  !transform heliocentric Keplerian orbital elements to heliocentric, ecliptic, Cartesian states
  call kema2rv(ele,rv_out)

  close(21)
  return
  end subroutine
!*******************************************************
  subroutine two_body_ke(ke_in,epoch,t_obs,mass1,mass2,rv_out)
  implicit none

  !input
  real*8,dimension(6),intent(in)::ke_in         !heliocentric Keplerian elements (a,e,i,w,Om,M) [au,-,4 x rad]
  real*8,intent(in)::epoch                                 !epoch of initial conditions [JD]
  real*8,intent(in)::t_obs                                 !time to propagate to [JD]
  real*8,intent(in)::mass1                                    !mass of the main body, e.g. Sun [Msun]
  real*8,intent(in)::mass2                                    !mass of the second body, e.g. asteroid [Msun]

  !ouput
  real*8,dimension(6),intent(out)::rv_out       !heliocentric state x,y,z,vx,vy,vz [au,au/day]
  
  !local
  real*8::meanmotion                            !mean motion [rad/day]
  real*8,dimension(6)::ele                      !heliocentric Keplerian elements (a,e,i,w,Om,M) [au,-,4 x rad]

  call semia2n(ke_in(1),mass1,mass2,meanmotion)
  call epoch2ma(meanmotion,epoch,ke_in(6),t_obs,ele(6))
  ele(1:5)=ke_in(1:5)  
  call kema2rv(ele,rv_out)

  return
  end subroutine

end module




 
