module orbdyn_m
!//////////////////////////////////////////////////////////////////////////////////
! ORBITAL DYNAMICS SUBROUTINES
!
! written by S.Eggl 20170101
!//////////////////////////////////////////////////////////////////////////////////

 use kind_m, only:wp,ip
 use const_m, only:gk,deg2rad,rad2deg,pi,pix2,jd_mjd
 use trafo3_m

 public::dist3d     !calculate Euclidean distances between objects, solar elongation and phase angle
 public::lighttraveltime !1-way light travel time

 public::com2ke     !Cometetary to Keplerian orbital elements
 public::ke2com     !Keplerian to cometary orbital elements
 
 public::keta2rv    !Keplerian oribtal elements to heliocentric ecliptic state vectors (position angle : True Anomaly)
 public::kema2rv    !Keplerian oribtal elements to heliocentric ecliptic state vectors (position angle : Mean Anomaly)
 
 public::rv2kema    !calculate Keplerian orbital elements from rv vectors
 
 public::ma2ea      !Mean Anomaly (ma) to Eccentric Anomaly (ea)
 public::ea2ma      !Eccentric Anomaly (Ea) to Mean Anomaly (ma) 
 public::ma2ta      !Mean Anomaly (ma) to True Anomaly (ta) 
 public::ta2ma      !True Anomaly (ta) to Mean Anomaly (ma)
 
 public::t2ma       !calculate mean anomaly after time t
 public::epoch2ma   !calculate mean anomaly at epoch
 
 public::semia2p    !semimajor axis to period
 public::semia2n    !semimajor axis to mean motion
 
 public::visvivav2  !calculate orbital speed squared (v^2) from vis-viva equation
 public::visvivar   !calculate orbital focal distance (r) from vis-viva equation
 
 public::tofma1ta2  !time of flight given mixed pair of mean anomaly and true anomaly angles
 
 public::oe2xy2     !planar orbital element to x,y coordinate conversion
 
 public::orb_dist_ta !calculate distance between two points on orbits given true anomalies
 public::orb_dist_ma !calculate distance between two points on orbits given mean anomalies
 public::orb_dist_tm !calculate distance between two points on orbits given true and mean anomalies
 
 public::wrap360
 
! interface
!  function lighttraveltime(d)
!   use kind_m, only:wp
!   use const_m,only:caupd
!   real(kind=wp),intent(in)::d
!   real(kind=wp)::lighttraveltime
!  end function
!end interface

contains
!***********************************************
subroutine dist3d(ra,re,rea,dea,dsa,dse,phase,solarelon)
        implicit none
        real*8,intent(in)::ra(3),re(3) !heliocentric positions of asteroid and Earth
        real*8,intent(out)::rea(3) !Earth-asteroid vector
        real*8,intent(out)::dea,dsa,dse !distance Earth-asteroid, Sun-asteroid and Sun-Earth
        real*8,intent(out)::phase  !phase angle
        real*8,intent(out)::solarelon !solar elongation
        !vector Earth-asteroid
        rea=ra-re
        !distance Earth-asteroid
        dea=sqrt(dot_product(rea,rea))
        !distance asteroid - Sun
        dsa=sqrt(dot_product(ra,ra))
        !distance Sun-Earth
        dse=sqrt(dot_product(re,re))
        !phase angle is angle between vectors asteroid -> Earth and asteroid -> Sun (careful with the signs)
        phase=acos(dot_product(rea,ra)/dea/dsa)
        !write(*,*)xya,xye,d,phase
        solarelon=acos(dot_product(-re,rea)/dse/dea)
        return
end subroutine
!******************************************************
function lighttraveltime(d)
    use kind_m, only:wp
    use const_m,only:caupd
!1-way light travel time

 implicit none
 real(kind=wp),intent(in)::d !distance between objects (non-relativistic) [au]
 real(kind=wp)::lighttraveltime !light travel time between objects [days]

 lighttraveltime=d/caupd
 return
end function
!******************************************************
subroutine wrap_angle(a,l,u)
! wraps the angle a to lie between 0 and w deg
 implicit none
 real(kind=wp),intent(inout)::a
 real(kind=wp),intent(in)::l,u
 
 if(a.lt.l) then
  do while(a.lt.l)
   a=a+u
  end do
 end if
 
 if (a.gt.u) a=mod(a,u)
return
end subroutine
!*********************************************** 
SUBROUTINE com2ke(com,epoch,m1,m2,ke,n)
!Transform cometetary to Keplerian orbital elements
!com(1:6)={q,e,i,node,w,Tp} all angles in deg
!ke(1:6)={a,e,i,w,node,M} all angles in deg

implicit none
real(kind=wp),intent(in),dimension(1:6)::com
real(kind=wp),intent(in)::epoch,m1,m2 !Epoch of orbital elements [JD], mass of host(star)[Msun], mass of orbiting body[Msun]
real(kind=wp),intent(out),dimension(1:6)::ke !Keplerian orbital elements {a,e,i,w,node,M}

real(kind=wp),intent(out),optional::n !mean motion [rad/day]

real(kind=wp)::mm,depoch

!pericenter distance to semimajor axis
if(com(2).ge.0._wp) then
 ke(1)=com(1)/(1._wp-com(2))
else
 write(*,*)'Warning: cometary eccentricity negative'
 ke(1)=com(1)
end if

call semia2n(ke(1),m1,m2,mm)

!mean anomaly
depoch=epoch-com(6)

if(abs(depoch).gt.jd_mjd) then
 write(*,*)'Warning, epoch of osculating elements and time of pericenter passage may not be in the same units (jd / mjd)'
end if

ke(6)=(epoch-com(6))*mm*rad2deg

call wrap_angle(ke(6),0._wp,360._wp)

!e, i
ke(2:3)=com(2:3)
!w
ke(4)=com(5)
!node
ke(5)=com(4)

if(present(n)) then
  n=mm
end if

return
end subroutine
!*********************************************** 
SUBROUTINE ke2com(ke,epoch,m1,m2,com,n)
!Transform Keplerian to cometary orbital elements
!ke(1:6)={a,e,i,w,node,M} all angles in deg
!com(1:6)={q,e,i,node,w,Tp} all angles in deg

implicit none
real(kind=wp),intent(in)::epoch,m1,m2 !Epoch of orbital elements [JD], mass of host(star)[Msun], mass of orbiting body[Msun]
real(kind=wp),intent(in),dimension(1:6)::ke
real(kind=wp),intent(out),dimension(1:6)::com

real(kind=wp),intent(out),optional::n !mean motion [rad/day]

real(kind=wp)::period,mm 


com(1)=ke(1)*(1._wp-ke(2))
com(2:3)=ke(2:3)
com(4)=ke(5)
com(5)=ke(4)

call semia2n(ke(1),m1,m2,mm)
period=pix2/mm

com(6)=epoch-ke(6)*deg2rad/mm+period

if(present(n)) then
 n=mm
end if

return
end subroutine

!*********************************************** 

SUBROUTINE ta2r(a,e,ta,r,deg)
!calculate distance to focus from true anomaly
implicit none
real(kind=wp),intent(in)::a,e,ta   !semimajor axis [au], eccentricity, true anomaly
logical,intent(in),optional::deg
real(kind=wp),intent(out)::r       !distance to focus [au] (e.g. sun)

real(kind=wp)::f

if(present(deg)) then
 if(deg) then
  f=ta*deg2rad
 else
 f=ta
 end if
else
f=ta
end if

 r=(a*(1._wp-e*e))/(1._wp+e*cos(f))
 
 return
 end subroutine

!*********************************************** 
SUBROUTINE t2ma(n,t,ma0,ma,deg)
!calculate mean anomaly after time t
implicit none
  real(kind=wp),intent(in)::n !mean motion [1/Gaussian Day]
  real(kind=wp),intent(in)::t !time [Gaussian Days]
  real(kind=wp),intent(in)::ma0 !mean anomaly at time t=0 
  logical,intent(in),optional::deg !degrees or radians
  
  real(kind=wp),intent(out)::ma  !mean anomaly at time t 
  
  real(kind=wp)::m0
  
  if(present(deg)) then
   if(deg) then
    ma=n*t+ma0*deg2rad
   else
   ma=n*t+ma0
  end if
  else
  ma=n*t+ma0
  end if

  call wrap_angle(ma,0._wp,pix2)

  if(present(deg)) then
   if(deg) then
    ma=ma*rad2deg
   end if
  end if
   
return
end subroutine
!*********************************************** 
SUBROUTINE epoch2ma(n,epoch0,ma0,epoch,ma,deg)
!calculate mean anomaly at epoch 
implicit none
  real(kind=wp),intent(in)::n !mean motion [rad/Gaussian Day]
  real(kind=wp),intent(in)::epoch0 !inital epoch [JD or MJD]
  real(kind=wp),intent(in)::ma0 !mean anomaly at epoch0
  real(kind=wp),intent(in)::epoch !final epoch to progagate ma to [JD or MJD, same as epoch0!]
  logical,intent(in),optional::deg !degrees or radians
  
  real(kind=wp),intent(out)::ma  !mean anomaly at time t 
  
  real(kind=wp)::m0,t
  
  t=epoch-epoch0
  
  if(present(deg)) then
   if(deg) then
    ma=n*t*rad2deg+ma0
   else
    ma=(n*t+ma0)*rad2deg
   end if
  else
   ma=(n*t+ma0)*rad2deg
  end if
  
  call wrap_angle(ma,0._wp,360._wp)

  if(present(deg)) then
   if(deg) then
   else
    ma=ma*deg2rad
   end if
  else
    ma=ma*deg2rad
  end if
   
return
end subroutine

!*********************************************** 
SUBROUTINE tofma1ta2(n,e,ma1,ta2,tof)
!calculates time of flight from a mixed pair of angles mean anomaly at time 1 and true anomaly at time 2
!angles are in [rad] !
  implicit none
      real(kind=wp),intent(in)::n,e !mean motion [1/Gaussian Days], eccentricity
      real(kind=wp),intent(in)::ma1,ta2 !mean anomaly at time 1, true anomaly at time 2  
      real(kind=wp),intent(out)::tof !time of flight [Gaussian Days] 
      
      real(kind=wp)::ma2,dma
      
      call ta2ma(e,ta2,ma2)
      dma=(ma2-ma1)
      if(dma<0._wp) dma=dma+pix2
      tof=dma/n      
    return
end subroutine

!***********************************************     
subroutine semia2p(a,m1,m2,p)
      !calculates period from semimajor axis and masses
      use const_m, only: wp, gk2, pix2
      implicit none
      real(kind=wp),intent(in)::m1,m2 !masses of the two bodies [Msun]
      real(kind=wp),intent(in)::a !semimajor axis in [au]
      real(kind=wp),intent(out)::p !period [Gaussian Days]
      
      real(kind=wp)::mu,n !GM and mean motion
    
      mu=gk2*(m1+m2)

      n=sqrt(mu/a**3)
      
      p=pix2/n

      return
      end subroutine
!***********************************************     
subroutine semia2n(a,m1,m2,n)
!calculate mean motion n from semimajor axis and masses 
use const_m, only:wp,gk2
implicit none
      real(kind=wp),intent(in)::m1,m2 !masses of the two bodies
      real(kind=wp),intent(in)::a !semimajor axis in [au]
      real(kind=wp),intent(out)::n !mean motion [rad/Gaussian Day]
      real(kind=wp)::mu
    
      mu=gk2*(m1+m2)

      n=sqrt(mu/a**3)
    return
end subroutine
      
      
 !***********************************************     
     
subroutine visvivav2(m1,m2,r,a,v2)
!calculates square of orbital velocity from vis viva equation
use const_m, only:wp,gk2
implicit none
   real(kind=wp),intent(in)::m1,m2 !masses of the two bodies [Msun]
   real(kind=wp),intent(in)::r,a  ! orbital distance from focus, semimajor axis [au]
   real(kind=wp),intent(out)::v2  ! orbital velocity^2 [(au/D)^2]
   real(kind=wp)::mu
    
      mu=gk2*(m1+m2)
      v2=mu*(2._wp/r+1._wp/a)

   return
end subroutine

!***********************************************

subroutine visvivar(m1,m2,v2,a,r)
!calculates orbital distance to focus from vis viva equation
use const_m, only:wp,gk2
implicit none
   real(kind=wp),intent(in)::m1,m2 !masses of the two bodies
   real(kind=wp),intent(in)::v2,a  ! orbital velocity^2 [(au/D)^2], semimajor axis [au]
   real(kind=wp),intent(out)::r  ! orbital distance from focus [au] 
   real(kind=wp)::mu
    
      mu=gk2*(m1+m2)
      r=2._wp/(v2/mu+1._wp/a)

   return
end subroutine


!***********************************************
SUBROUTINE ea2ma(e,ea,ma,deg)
!------------------------------------------------------------------------------------------------------------
! converts eccentric anomaly (ea) to mean anomaly (ma)
!
! written by Siegfried EGGL  20170131
!
! standard input of angles is [rad], with the option of having [deg]
! 
!------------------------------------------------------------------------------------------------------------
implicit none
real(kind=wp),intent(in)::e,ea
logical,intent(in),optional::deg

real(kind=wp),intent(out)::ma

real(kind=wp)::ea1

if(e.eq.0._wp) then
 ma=ea
else
if(present(deg)) then
 if(deg) then
  ea1=ea*deg2rad
 else
  ea1=ea
 end if
else
 ea1=ea
end if

 ma=ea1-e*sin(ea1)

if(present(deg)) then
   if(deg) then
    ma=ma*rad2deg
   end if
end if

end if !e=0

return
end subroutine
!***********************************************

subroutine ma2ea(ecc,ma,ea,deg)
!-------------------------------------------------------------------------------------------------
!solves Kepler’s equation by applying Newton Raphson Method
!with an accuracy limit of "deat"
!
!----------------------------------------------------------------------------------------------------------
!Input:
!
!m[real].................mean longitude (radians!)
!ecc[real]................Eccentricity (numerical Eccentricity <1!)
!
!----------------------------------------------------------------------------------------------------------
!Output:
!
!ea[real]................Eccentric Anomaly (radians!)
!
!
!written by Siegfried Eggl  20061222
!modified 20111026
!dependencies: none
!-------------------------------------------------------------------------------------------------
implicit none
        integer(kind=ip)::i,iter
        real(kind=wp),intent(in)::ma,ecc
        logical,intent(in),optional::deg
        real(kind=wp),intent(out)::ea
        real(kind=wp)::ea0,dea,deat,m


if(ecc.eq.0._wp) then
 ea=ma
else

ea=0._wp
ea0=1.3421_wp
deat=1.E-8
dea=ABS(ea-ea0)
iter=12

if(present(deg)) then
 if(deg) then
  m=ma*deg2rad
 else 
  m=ma
 end if
else
m=ma
end if

if(ma.eq.0._wp) then
 ea=0._wp
 dea=deat
!elseif (abs(ma-pi).lt.deat) then
! ea=pi
! dea=deat
else

do i=1,iter
   if (dea>deat)then
      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._wp-ecc*COS(ea0))
      dea=ABS(ea-ea0)
      ea0=ea
   else 
      continue
   end if !deat
end do !iter

if (ea.le.1.E-14) then
   ea=0._wp
end if
!if precision is not achieved try with different initial condition
if(dea>deat) then

ea=0._wp
ea0=3.121_wp
dea=ABS(ea-ea0)
 do i=1,iter
   if (dea>deat)then
      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._wp-ecc*COS(ea0))
      dea=ABS(ea-ea0)
      ea0=ea
   else 
      continue
   end if
 end do !iter
end if


if(dea>deat) then

ea=0._wp
ea0=6.1_wp
dea=ABS(ea-ea0)
 do i=1,iter
   if (dea>deat)then
      ea=ea0-(ea0-ecc*SIN(ea0)-m)/(1._wp-ecc*COS(ea0))
      dea=ABS(ea-ea0)
      ea0=ea
   else 
      continue
   end if
 end do !iter
end if !deat
end if !0, pi

if(dea>deat) then
   write(unit=*,fmt=*)'convergence-error in subroutine ma2ea'
   write(unit=*,fmt=*)'target precision:',deat
   write(unit=*,fmt=*)'achieved precision:',dea
end if


if (present(deg)) then
 if(deg) then
  ea=ea*rad2deg
 end if
end if

end if !e=0

if (ea.lt.0._wp) ea=ea+pix2

return

end subroutine

!***********************************************
SUBROUTINE ma2ta(e,ma,ta,deg)
!------------------------------------------------------------------------------------------------------------
! converts mean anomaly (ma) to true anomaly (ta)
!
! written by Siegfried EGGL  20170131
!
! standard input of angles is [rad], with the option of having [deg]
! requires: subroutine ma2ea 
!------------------------------------------------------------------------------------------------------------

implicit none
real(kind=wp),intent(in)::e,ma
logical,intent(in),optional::deg

real(kind=wp),intent(out)::ta

real(kind=wp)::ea,rsinta,rcosta


if(e.eq.0._wp) then
 ta=ma
else

call ma2ea(e,ma,ea,deg)

if(present(deg)) then
 if(deg) then
  ea=ea*deg2rad
 end if
end if

rcosta=(cos(ea)-e)    !*a
rsinta=sqrt(1._wp-e*e)*sin(ea) !*a

ta=atan2(rsinta,rcosta)

if(ta.lt.0._wp) ta=ta+pix2

if(present(deg)) then
 if(deg) then
  ta=ta*rad2deg
 end if
end if

end if !e=0

return
end subroutine

!***************************************************

SUBROUTINE ta2ma(e,ta,ma,deg)
!------------------------------------------------------------------------------------------------------------
! converts true anomaly (ta) to mean anomaly (ma)
!
! written by Siegfried EGGL  20170131
!
! standard input of angles is [rad], with the option of having [deg]
!------------------------------------------------------------------------------------------------------------

implicit none
real(kind=wp),intent(in)::e !orbital eccentricity
real(kind=wp),intent(in)::ta !true anomaly
logical,optional::deg

real(kind=wp),intent(out)::ma

real(kind=wp)::ea,tarad
real(kind=wp)::costa,sinta,ecostap1,cosE,sinE


if(e.eq.0._wp) then
 ma=ta
else

if(present(deg)) then
 if(deg) then
   tarad=ta*deg2rad
 else
   tarad=ta
 end if
 else
   tarad=ta
end if

costa=cos(tarad)
sinta=sin(tarad)

ecostap1=e*costa+1._wp

   !calculation of eccentric anomaly eanom (=E)
               cosE=(e+Costa)/(ecostap1)
               sinE=Sqrt(1._wp-e*e)*sinta/(ecostap1)
               ea=ATAN2(sinE,cosE)  
               
               if(ea.lt.0._wp) ea=ea+pix2
            
               !mean anomaly 'ma' via Kepler's equation
               ma=(ea-e*sinE) 
              
if(present(deg)) then
 if (deg) then
  ma=ma*rad2deg
 end if
end if

end if ! e=0

return

end subroutine

!*******************************************
SUBROUTINE keta2rv(ele,rv,deg,m0,m1)
!------------------------------------------------------------------------------------------------------------
!transforms Keplerian orbital elements ele (a,e,i,w,Omega,true anomaly) to Cartesian positions and velocities
!
! written by Siegfried EGGL  20161117
!
! standard input of angles is [rad], with the option of having [deg]
!------------------------------------------------------------------------------------------------------------
        
implicit none   
        
real(kind=wp),dimension(1:6),intent(in)::ele
logical,optional,intent(in)::deg
real(kind=wp),optional,intent(in)::m0,m1

real(kind=wp),dimension(1:6),intent(out)::rv
        

real(kind=wp)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,w,om,ta       
real(kind=wp)::cosom,sinom,cosi,sini,px,py,pz,qx,qy,qz,ah,cappaq
real(kind=wp)::cosex,sinex,rh,bh,cosw,sinw
real(kind=wp)::costa,sq1me2,denom
        
a=ele(1)
e=ele(2)

if(present(deg)) then
 if(deg) then
   incl=ele(3)*deg2rad
   w=ele(4)*deg2rad
   om=ele(5)*deg2rad
   ta=ele(6)*deg2rad
  else
   incl=ele(3)
   w=ele(4)
   om=ele(5)
   ta=ele(6)
  end if
        else   
incl=ele(3)
w=ele(4)
om=ele(5)
ta=ele(6)   
end if

IF(present(m0)) then
 cappaq=m0   
else
 cappaq=1._wp
end if

IF(present(m1)) then
 cappaq=m0+m1
end if

 IF (E .EQ.0._wp)  then
    w  = 0._wp
 end if
 IF (INCL .EQ.0._wp) then
    om  = 0._wp
 end if
 COSW = COS(W)
 SINW = SIN(W)
 COSom = COS(om)
 SINom = SIN(om)
 COSI = COS(INCL)
 SINI = SIN(INCL) 
 PX= COSW*COSom-SINW*SINom*COSI
 PY= COSW*SINom+SINW*COSom*COSI
 PZ= SINW*SINI
 QX=-SINW*COSom-COSW*SINom*COSI
 QY=-SINW*SINom+COSW*COSom*COSI
 QZ= COSW*SINI

 costa=cos(ta)
 sq1me2=sqrt(1._wp-e*e)
 denom=(1._wp+e*costa)
 
 cosex=(e+costa)/denom
 sinex=sq1me2*sin(ta)/denom

 AH = A *sq1me2
 XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
 XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
 XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
 RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )

 if(A.eq.0._wp.or.RH.eq.0._wp) then
    AH=1._wp
 else
    AH = SQRT(cappaq )/(SQRT(A )*RH)
 end if
  
 BH = A *sq1me2*COSEX
 VH1  = AH*(-A *PX*SINEX+BH*QX)
 VH2  = AH*(-A *PY*SINEX+BH*QY)
 VH3  = AH*(-A *PZ*SINEX+BH*QZ)

   
 rv(1)=  XH1
 rv(2)=  XH2
 rv(3)=  XH3
 rv(4)=VH1*gk
 rv(5)=VH2*gk
 rv(6)=VH3*gk 

        RETURN
END  subroutine
!****************************************

SUBROUTINE kema2rv(ele,rv,mass_primary,mass_secondary)
!------------------------------------------------------------------------------------------------------------
!transforms Keplerian orbital elements ele (a,e,i,w,Omega, mean anomaly) to Cartesian positions and velocities
!
! written by Siegfried EGGL  20161208
!
! standard input of angles is [deg]
!------------------------------------------------------------------------------------------------------------
        
        implicit none   
        
real(kind=wp),dimension(1:6),intent(in)::ele
real(kind=wp),intent(in)::mass_primary,mass_secondary

real(kind=wp),dimension(1:6),intent(out)::rv
        

real(kind=wp)::xh1,xh2,xh3,vh1,vh2,vh3,a,e,incl,w,om,ma,ea,ta       
real(kind=wp)::cosom,sinom,cosi,sini,px,py,pz,qx,qy,qz,ah,cappaq
real(kind=wp)::cosex,sinex,rh,bh,cosw,sinw
real(kind=wp)::sq1me2,denom
        

a=ele(1)
e=ele(2)
incl=ele(3)*deg2rad
w=ele(4)*deg2rad
om=ele(5)*deg2rad
ma=ele(6)*deg2rad

cappaq=mass_primary+mass_secondary

 IF (E .EQ.0._wp)  then
    w  = 0._wp
 end if
 IF (INCL .EQ.0._wp) then
    om  = 0._wp
 end if
 COSW = COS(W)
 SINW = SIN(W)
 COSom = COS(om)
 SINom = SIN(om)
 COSI = COS(INCL)
 SINI = SIN(INCL) 
 PX= COSW*COSom-SINW*SINom*COSI
 PY= COSW*SINom+SINW*COSom*COSI
 PZ= SINW*SINI
 QX=-SINW*COSom-COSW*SINom*COSI
 QY=-SINW*SINom+COSW*COSom*COSI
 QZ= COSW*SINI

 call ma2ea(e,ma,ea)
 
 sq1me2=sqrt(1._wp-e*e)
 
 cosex=cos(ea)
 sinex=sin(ea)

 AH = A *sq1me2
 XH1  = A *PX*(COSEX-E )+AH*QX*SINEX
 XH2  = A *PY*(COSEX-E )+AH*QY*SINEX
 XH3  = A *PZ*(COSEX-E )+AH*QZ*SINEX
 RH = SQRT(XH1 *XH1 +XH2 *XH2 +XH3 *XH3 )

 if(A.eq.0._wp.or.RH.eq.0._wp) then
    AH=1._wp
 else
    AH = SQRT(cappaq )/(SQRT(A )*RH)
 end if
  
 BH = A *sq1me2*COSEX
 VH1  = AH*(-A *PX*SINEX+BH*QX)
 VH2  = AH*(-A *PY*SINEX+BH*QY)
 VH3  = AH*(-A *PZ*SINEX+BH*QZ)

 rv(1)=  XH1
 rv(2)=  XH2
 rv(3)=  XH3
 rv(4)=VH1*gk
 rv(5)=VH2*gk
 rv(6)=VH3*gk 

        RETURN
END  subroutine


!******************************************************      
      
subroutine orb_dist_ta(ke1,ke2,dr,dv,deg,m0,m1,m2)
!calculates the distance between two particles on heliocentric orbits with Keplerian elements ke1 and ke2
!ke1=(a,e,i,w,OM,true anomaly)
!
implicit none

real(kind=wp),dimension(1:6),intent(in)::ke1,ke2
logical,optional,intent(in)::deg
real(kind=wp),optional,intent(in)::m0,m1,m2

real(kind=wp),intent(out)::dr,dv

real(kind=wp),dimension(1:6)::rv1,rv2
real(kind=wp),dimension(1:3)::dp,du

call keta2rv(ke1,rv1,deg,m0,m1)
call keta2rv(ke2,rv2,deg,m0,m2) 

dp(1:3)=rv2(1:3)-rv1(1:3)
du(1:3)=rv2(4:6)-rv1(4:6)

dr=sqrt(dot_product(dp,dp))
dv=sqrt(dot_product(du,du))

return
end subroutine
!******************************************************   

subroutine orb_dist_ma(ke1,ke2,dr,dv,deg,m0,m1,m2)
!calculates the distance between two particles on heliocentric orbits with Keplerian elements ke1 and ke2
!ke1=(a,e,i,w,OM,Mean anomaly)
!
implicit none

real(kind=wp),dimension(1:6),intent(in)::ke1,ke2
logical,optional,intent(in)::deg
real(kind=wp),optional,intent(in)::m0,m1,m2

real(kind=wp),intent(out)::dr,dv

real(kind=wp),dimension(1:6)::ele1,ele2,rv1,rv2
real(kind=wp),dimension(1:3)::dp,du

call kema2rv(ke1,rv1,deg,m0,m1)
call kema2rv(ke2,rv2,deg,m0,m2) 

dp(1:3)=rv2(1:3)-rv1(1:3)
du(1:3)=rv2(4:6)-rv1(4:6)


dr=sqrt(dot_product(dp,dp))
dv=sqrt(dot_product(du,du))

return
end subroutine
 
!******************************************************************************
subroutine orb_dist_tm(ke1,ke2,dr,dv2,deg,m0,m1,m2)
!calculates the distance between two particles on heliocentric orbits with Keplerian elements ke1 and ke2
!ke1=(a,e,i,w,OM,true anomaly)
!ke2=(a,e,i,w,OM,Mean anomaly)
!MIXED ANGLES! OUTPUT SQUARED VELOCITIES!!!
implicit none

real(kind=wp),dimension(1:6),intent(in)::ke1,ke2
logical,optional,intent(in)::deg
real(kind=wp),optional,intent(in)::m0,m1,m2

real(kind=wp),intent(out)::dr,dv2

real(kind=wp),dimension(1:6)::rv1,rv2
real(kind=wp),dimension(1:3)::dp,du

call keta2rv(ke1,rv1,deg,m0,m1)
call kema2rv(ke2,rv2,deg,m0,m2) 

dp(1:3)=rv2(1:3)-rv1(1:3)
du(1:3)=rv2(4:6)-rv1(4:6)


dr=sqrt(dot_product(dp,dp))
dv2=dot_product(du,du)

return
end subroutine

!*******************************************
subroutine oe2xy2(a,e,n,m0,w,t,x)
!---------------------------------------------------------------------------
!converts orbital elements to xy coordinates in the planar two body problem 
!
! written by Siegfried Eggl  20140304
!
! dependencies: nr
!----------------------------------------------------------------------------------------------------------
!Input:
!
!real::
!a[au] semimajor axis
!e[] eccentricity
!n[rad/D] mean motion
!m0        [rad] initial mean longitude
!w[rad] argument of pericenter
!t[D] time-steps
!----------------------------------------------------------------------------------------------------------
!Output
!
!real::
!x(1:2)    [au] x,y coordinates wrt focus
!

!--------------------------------------------------------------------------
implicit none
real(kind=wp),intent(in)::a,e,n,m0,w,t
real(kind=wp),intent(out)::x(1:2)
!local
real(kind=wp)::rotmat(1:2,1:2),ea,mm,ml,rcosf,rsinf
real(kind=wp),dimension(1:3)::x2

!calculate mean longitude as a function of time (Ml=n*t, n=2pi/T)
!since phi is the angle between r1 and r2 rather than r0 and r2 the supplementary angle to phi has to be added to the mean longitude ml
ml=mod(n*t+m0,pix2)

!find mean longitude as a function of time (M=n*t+l0-w, n=2pi/T)
mm=mod(ml-w+pix2,pix2)

!eccentric anomaly
call ma2ea(e,mm,ea)
!true anomaly -> x,y

rcosf=a*(Cos(ea)-e)     !x
rsinf=a*sqrt(1._wp-e*e)*Sin(ea)  !y

!rotate orbit to correct pericenter
rotmat(1,1:2)=(/cos(w),-sin(w)/)
rotmat(2,1:2)=(/sin(w),cos(w)/)

!final xy coordinates
x=matmul(rotmat,(/rcosf,rsinf/))

return
end subroutine

!************************************************************************************************************

!...Calculates Keplerian Orbital Elements from heliocentric position and velocity vectors

SUBROUTINE rv2kema (rv,ele,mass_primary_msun,mass_secondary_msun)
use const_m, only: gk,rad2deg,deg2rad,pi,pix2

implicit none

real(kind=wp),dimension(1:6),intent(in)::rv     
real(kind=wp),intent(in)::mass_primary_msun,mass_secondary_msun

real(kind=wp),dimension(1:6),intent(out)::ele

real(kind=wp)::mm

integer(kind=ip)::j
integer(kind=ip),parameter::gd=3

real(kind=wp)::e2,sq1me2,eanom,Lprojxy,U(1:3,1:3)
real(kind=wp):: tanom,Hnorm,sinE,cosE,costa,ectp1
real(kind=wp),dimension(1:3)::LRL,L,r,v,H,node,rU
real(kind=wp)::cappaq,a,e,incl,om,gom,mmm,L2,rnorm,nnorm

!determine mass coefficients for orbit (units: solar masses)
cappaq=mass_primary_msun+mass_secondary_msun

!write(*,*)'cappaq',cappaq
r=rv(1:3)
v=rv(4:6)/gk
         
!calculation of unit mass angular momentum vector L and its projection onto xy-plane
call crossp3d(r,v,L)
L2=Dot_Product(L,L)
Lprojxy=SQRT( L(1)*L(1)+L(2)*L(2) )

 ! inclination "incl"
incl=ATAN2(Lprojxy,L(3))*rad2deg

!argument of the ascending node "gom"
       gom=Atan2(L(1),-L(2))*rad2deg
       if (gom.lt.0._wp) then
           gom = gom+360._wp
        end if

! Laplace-Runge-Lenz vector
        call crossp3d(v,L,LRL)
        rnorm=Sqrt(Dot_Product(r(:),r(:)))

        LRL(:)=LRL(:)/cappaq-r(:)/rnorm
         
! eccentricity "e" = |LRL|
        e=Sqrt(Dot_Product(LRL,LRL))

        if(e.lt.1.d-14)then 
           e=0._wp
           e2=0._wp
           om=0._wp
           sq1me2=1._wp
         else
           e2=e*e
           sq1me2=sqrt(1._wp-e2)
         end if
         
    ! semi-major axis "a"
        a=L2/(cappaq*(1._wp-e2))    
   
   !calculation of vector "node" pointing to ascending node    
    ! attention: will change node for retrograde orbits!!!!
         node(1)=-L(2)*L(3)
         node(2)=L(1)*L(3)
         node(3)=0._wp

         nnorm=Sqrt(Dot_Product(node,node))
 
       !coplanar, fix nodes
       if(incl.lt.1.d-14) then
          gom=0._wp
         ! argument of pericenter "om" = angle of LRL to x-axis
          om=Acos(LRL(1)/Sqrt(Dot_Product(LRL,LRL)))
              if(LRL(2).lt.0._wp) then
                 om=360._wp-om*rad2deg
              else
                 om=om*rad2deg
              end if
  
        else

          call crossp3d(L,node,H)
          Hnorm=Sqrt(Dot_Product(H,H))
          om=Atan2(Dot_Product(LRL,H)*nnorm,Dot_Product(LRL,node)*Hnorm)*rad2deg
              if(om.lt.0._wp) om=om+360._wp 
  
        end if  

        !coplanar, circular
        if(e.lt.1.2d-14.and.incl.le.1.d-14) then
              gom=0._wp
              om=0._wp

              !mean anomaly 'mmm'
              mmm=Atan2(r(2),r(1))*rad2deg
           
        !non-coplanar, circular
         elseif (e.lt.1.d-14.and.incl.gt.1.d-14) then             
     
              !transformation of positionvector r into coordinate-system U spanned by Angular Momentum vector, vector Sun-ascending Node, 
              !and their crossproduct H, corresponding to motion of rU in the orbital plane. 
              

              call crossp3d(L,node,H)
              Hnorm=Sqrt(Dot_Product(H,H))
              ! calculation of the base - transformation Matrix
              U(:,1)=node(:)/nnorm
              U(:,2)=H/Hnorm
              U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
              call gauss(gd,U,r,rU)

              !true anomaly in orbital plane
              tanom=Atan2(rU(2),rU(1))

              !mean anomaly 'mmm'
              mmm=tanom*rad2deg     
            
         !coplanar, non-circular
         elseif (incl.lt.1.d-14.and.e.gt.1.d-14) then
              !calculation of true anomaly via vector "H" at right angle to LRL and L (alternative)
                      call crossp3d(L,LRL,H)
                      Hnorm=Sqrt(Dot_Product(H,H))
                      tanom= Atan2(Dot_Product(H,r)*e,Dot_Product(LRL,r)*Hnorm) 
     

              !calculation of eccentric anomaly eanom (=E)
               costa=cos(tanom)
               ectp1=(1._wp+e*costa)
               cosE=(e+costa)/ectp1
               sinE=sq1me2*Sin(tanom)/(1._wp+e*costa)
               eanom=ATAN2(sinE,cosE)  

              !mean anomaly 'mmm' via Kepler's equation
              mmm=(eanom-e*sinE)*rad2deg              

         ! general case
         else
  
               call crossp3d(L,LRL,H)
               Hnorm=Sqrt(Dot_Product(H,H))
               ! calculation of the base - transformation Matrix
               U(:,1)=LRL(:)/e
               U(:,2)=H/Hnorm
               U(:,3)=L(:)/Sqrt(L2)
              
              !coordinate transformation
               call gauss(gd,U,r,rU)

              !true anomaly in orbital plane
               tanom=Atan2(rU(2),rU(1))


              !calculation of eccentric anomaly eanom (=E)
               costa=cos(tanom)
               ectp1=(1._wp+e*costa)
               cosE=(e+costa)/ectp1
               sinE=sq1me2*Sin(tanom)/ectp1
               eanom=Atan2(sinE,cosE)  
                       
               !mean anomaly 'mmm' via Kepler's equation
               mmm=(eanom-e*sinE)*rad2deg
         end if
         
          if(mmm.lt.0._wp) mmm=mmm+360._wp
          if(e.lt.0._wp.or.e.gt.1._wp) e=1._wp
          if(a.lt.0._wp.or.a.gt.huge(a)) a=0._wp 
          if(incl.gt.180._wp.or.incl.lt.0._wp) om=-100._wp

          
          ele(1)=a
          ele(2)=e
          ele(3)=incl
          ele(4)=om
          ele(5)=gom
          ele(6)=mmm

          do j=4,6
            if(ele(j).gt.360._wp.or.ele(j).lt.0._wp.or.ele(j).ne.ele(j)) then
              ele(j)=-100._wp
            end if
          end do
      
          
!mean motion
          mm=gk*sqrt(cappaq)*a**(-1.5_wp)
        
          
       RETURN
     END subroutine
!*************************************************************************************************
subroutine cosmic_v(m0,m1,u2,r,vinf,vimp)
!calculates impact velocity of two body encounter based on the true v_infinity (velocity at infinity, unperturbed case) 
!no correction for evaluation of v_infinity at a finite distance from gravitational attractor 
implicit none
real(kind=wp),intent(in)::m0,m1 !mass of center of hyperbolic orbit and secondary body  [Msun]
real(kind=wp),intent(in)::u2,r  !relative velocitites, relative distance [au/d], [au]

real(kind=wp)::vimp,vinf

real(kind=wp)::a,mu,vesc2,vinf2

mu=gk2*(m0+m1)

a=1._wp/(u2/mu-2._wp/r)

vesc2=ve2earth
vinf=sqrt(mu/a)*vau2vkm
vimp=sqrt(vesc2+vinf2)

return
end subroutine 

!############################################################
end module
