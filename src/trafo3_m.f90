module trafo3_m
use kind_m
use const_m
   implicit none
      
  public::rotx
  public::roty
  public::rotz
  
  public::rotx2
  public::roty2
  public::rotz2
  
  public::rotax
  public::rotax2
  
  public::jd2greg
  
  public::he2icrf
  public::icrf2he
  
  public::icrf2rtn
  
  public::icrf2radec
  public::radec2icrf
  
  public::rvicrf2dradec
  public::rvhe2icrf
  public::rvicrf2he
    
  public::lambet2icrf
  
  public::cart2pol
  public::pol2cart
  
  public::cart2lat
  public::lat2cart
  
  public::v2basis
  
  public::crossp3d
  
  public::invert3x3
  
  public::gauss
!////////////////////////////////////////////////////////////////////////////
! functions
  public::angv1v2
  public::angv1v2n
  
  public::norm


  contains

!/////////////////////////////////////////////////////////////////////////////////
!                              Transformations and Rotations
!
!                                         by Siegfried Eggl  20170131 
!--------------------------------------------------------------------------------
!*********************************************************************************************************
subroutine rvicrf2dradec(rvicrf,rradec,drradec)
use const_m
!calculate equatorial coordinates and derivatives from ICRF position and velocities 
!INPUT:     rvicrf=(x,y,z,vx,vy,vz)_ICRF, 
!OUTPUT:    rradec=(r,RA [rad],DEC [rad]), drradec=(dr/dt,dRA/dt,dDEC/dt) 
!attention: some applicatons require dRA/dt*cos(dec)!

implicit none
real(kind=wp),dimension(1:6),intent(in)::rvicrf
real(kind=wp),dimension(1:3),intent(out)::rradec,drradec

real(kind=wp)::r2,x2py2,rdotv

!r.r
r2=dot_product(rvicrf(1:3),rvicrf(1:3))
!r.v
rdotv=dot_product(rvicrf(1:3),rvicrf(4:6))
!x^2+y^2
x2py2=dot_product(rvicrf(1:2),rvicrf(1:2))

!r
rradec(1)=sqrt(r2)
!RA [rad]
rradec(2)=atan2(rvicrf(2),rvicrf(1))
if(rradec(2).lt.0._wp) then
 rradec(2)=rradec(2)+pix2
end if
!DEC [rad]
rradec(3)=asin(rvicrf(3)/rradec(1))

!dr/dt
drradec(1)=rdotv/rradec(1)
!dRA/dt=(-y dx + x dy)/(x^2+y^2)
drradec(2)=(rvicrf(1)*rvicrf(5)-rvicrf(2)*rvicrf(4))/x2py2
!dDEC/dt=[(x^2+y^2)dz-z(x dx + y dy)]/{sqrt[(x^2+y^2)/r^2]*r^3}
drradec(3)=(rvicrf(6)*x2py2-rvicrf(3)*(dot_product(rvicrf(1:2),rvicrf(4:5))))/(sqrt(x2py2/r2)*r2**1.5d0)
   
return
end subroutine
!*********************************************************************************************************


subroutine icrf2rtn(rvicrf,rvtn,bt)
!rotates vectors from ICRF positions and velocities to RTN (radial transversal normal) frame
implicit none
real(kind=wp),dimension(1:6),intent(in)::rvicrf
real(kind=wp),dimension(1:6),intent(out)::rvtn
real(kind=wp),dimension(1:3,1:3),intent(out),optional::bt

real(kind=wp),dimension(1:3)::rn,vn
real(kind=wp)::t(1:3),n(1:3)
real(kind=wp),dimension(1:3,1:3)::rtnb

rn=rvicrf(1:3)/norm(rvicrf(1:3))
vn=rvicrf(4:6)/norm(rvicrf(4:6))

call crossp3d(rn,vn,n)
call crossp3d(n,rn,t)

rtnb(1,1:3)=rn
rtnb(2,1:3)=t
rtnb(3,1:3)=n

!write(*,*)'rvicrf',rvicrf
!write(*,*)'rn,t,n',rn,t,n

call invert3x3(rtnb,bt)

rvtn(1:3)=matmul(bt(1:3,1:3),rvicrf(1:3))
rvtn(4:6)=matmul(bt(1:3,1:3),rvicrf(4:6))

   
return
end subroutine
 !*********************************************************************************************************
subroutine icrf2radec(x,ra,dec)
    !transform icrf vector to Right Ascension and Declination [rad]
    implicit none
    real(kind=wp),intent(in),dimension(3)::x
    real(kind=wp),intent(out)::ra,dec

    real(kind=wp)::p(3)

    call cart2lat(x,p)

    dec=p(2)
    ra=p(3)

return
end subroutine

 !*********************************************************************************************************
subroutine radec2icrf(ra,dec,x)
    !transform Right Ascension and Declination [rad] to icrf unit vector
    implicit none
     real(kind=wp),intent(in)::ra,dec
    real(kind=wp),intent(out),dimension(3)::x
   
    real(kind=wp)::p(3)

    p(:)=(/1._wp,dec,ra/)
    
    call lat2cart(p,x)
    
return
end subroutine

 !*********************************************************************************************************
subroutine jd2greg(jd,a,mon,day)
!calculates Gregorian year, month, day from JD following Richards
implicit none
integer,intent(out)::a,mon,day
real(kind=wp),intent(in)::jd

integer,parameter::y=4716,n=12,m=2
real(kind=wp),parameter::j=1401,r=4,p=1461
real(kind=wp),parameter::v=3,u=5,s=153,w=2,B=274277,C=-38
real(kind=wp)::f,e,g,h

f=jd+j
f=f+real(int((int((4._wp*jd+B)/146097._wp)*3._wp)/4._wp))+C
e=r*f+v
g=real(int(mod(e,p)/r))
h=u*g+w
day=int(mod(h,s)/u)+1
mon=mod(int(real(h)/real(s))+m,n)+1
a=int(e/p)-y+int(real(n+m-mon)/real(n))

return
end subroutine 
!****************************************************************************************************
subroutine lambet2icrf(lambda,beta,xicrf,deg)
!changes geocentric elliptic longitude (lambda) and latitude (beta) into icrf orientaion vector (xicrf(1:3))
implicit none
real(kind=wp),intent(in)::lambda,beta
logical,intent(in)::deg !input in degree? .true./.false. (if .false. input is in radians)
real(kind=wp),dimension(1:3),intent(out)::xicrf

real(kind=wp),dimension(1:3)::elat,ecar


!ecliptic latitudinal coordinates
if(deg) then
 elat(:)=(/1._wp,lambda*deg2rad,beta*deg2rad/)
else
 elat(:)=(/1._wp,lambda,beta/)
end if 

call lat2cart(elat,ecar)

call he2icrf(ecar)

xicrf=ecar

return

end subroutine
!*****************************************************************************************************

       subroutine cart2pol(x,r)
!----------------------------------------------------

          implicit none
          real(kind=wp),intent(in)::x(1:3)
          real(kind=wp),intent(out)::r(1:3)
        
          r(1)=sqrt(Dot_Product(x,x))
          r(2)=atan2(x(2),x(1))
          r(3)=acos(x(3)/r(1))
          
     
          return
        end subroutine

!*************************************************************        
        subroutine pol2cart(r,x)
!----------------------------------------------------

        implicit none
        real(kind=wp),intent(in):: r(1:3)
        real(kind=wp),intent(out)::x(1:3)
        real(kind=wp)::sr3
        
        sr3=sin(r(3))

        x(1) = r(1)*sr3*cos(r(2))
        x(2) = r(1)*sr3*sin(r(2))
        x(3) = r(1)*cos(r(3))

        return
        end subroutine
!*******************************************************************
      subroutine cart2lat(x,p)
!------------------------------------------------------------------------------
!     Changes 3D Cartesian to latitudinal coordinates.
!     Attention: Latitudes {p(2)} are counted from the equator upwards!!!
!     
!     written by Siegfried Eggl 20140611
!
!     dependencies: none
!------------------------------------------------------------------------------
!     Input:
!     real::
!     x(1:3)			...[distance unit = DU] Cartesian coordinates
!------------------------------------------------------------------------------
!     Output:
!     real::
!     p(1:3)  ...[mixed] polar coordinate vector,
!     p(1)	  ...[DU] polar distance r
!	  p(2)    ...[rad] latitude (declination)
!     p(3)    ...[rad] longitude
!------------------------------------------------------------------------------
 implicit none
 real(kind=wp),intent(in)::x(1:3) 
 real(kind=wp),intent(out)::p(1:3)
        
 p(1)=sqrt(Dot_Product(x,x)) !r
 p(2)=asin(x(3)/p(1))        !latitude (declination)
 p(3)=atan2(x(2),x(1))       !longitude
 
 if(p(3).lt.0._wp) then
    p(3)=p(3)+pix2
 end if
     
 return
end subroutine 
       

!*************************************************************        
        subroutine lat2cart(p,x)
!------------------------------------------------------------------------------
!     Changes 3D latitudinal to Cartesian coordinates.
!     Attention: Latitudes {p(2)} are counted from the equator upwards!!!
!     
!     written by Siegfried Eggl 20140611
!
!     dependencies: none
!------------------------------------------------------------------------------
!     Input:
!     real::
!     p(1:3)  ...[mixed] polar coordinate vector,
!     p(1)	  ...[DU] polar distance r
!	  p(2)    ...[rad] latitude (declination)
!     p(3)    ...[rad] longitude
!------------------------------------------------------------------------------
!     Output:
!     real::
!     x(1:3)			...[distance unit = WDU] Cartesian coordinates
!------------------------------------------------------------------------------
        implicit none
        !input
        real(kind=wp),intent(in)::p(1:3)
        !output
        real(kind=wp),intent(out)::x(1:3)
        
        real(kind=wp)::cp2
        
        cp2=cos(p(2))
       
        x(1) = p(1)*cp2*cos(p(3))
        x(2) = p(1)*cp2*sin(p(3))
        x(3) = p(1)*sin(p(2))

        return
        end subroutine
!************************************************************
        subroutine rotx(theta,u)
!----------------------------------------------------
!       rotate vector u around x-axis for an angle of theta radians 
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::theta
        real(kind=wp)::st,ct,x(1:2)

        st=sin(theta)
        ct=cos(theta)
        x(1:2)=u(2:3)
        
        u(2) = x(1)*ct-x(2)*st
        u(3) = x(2)*ct+x(1)*st

        return
        end subroutine
        
        
!*****************************************************************
        subroutine roty(theta,u)
!----------------------------------------------------
!       rotate vector u around y-axis for an angle of theta radians 
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::theta
        real(kind=wp)::st,ct,x(1:2)

        st=sin(theta)
        ct=cos(theta)
        x(1)=u(1)
        x(2)=u(3)
        
        u(1) = x(1)*ct+x(2)*st
        u(3) = x(2)*ct-x(1)*st

        return
        end subroutine
        
        
!*****************************************************************

        subroutine rotz(theta,u)
!----------------------------------------------------
!       rotate vector u around y-axis for an angle of theta radians 
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::theta
        real(kind=wp)::st,ct,x(1:2)

        st=sin(theta)
        ct=cos(theta)
        x(1:2)=u(1:2)
        
        u(1) = x(1)*ct-x(2)*st
        u(2) = x(2)*ct+x(1)*st

        return
        end subroutine
        
        
!*****************************************************************
        subroutine rotx2(ct,st,u)
!----------------------------------------------------
!       rotate vector u around x-axis for an angle of theta radians 
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::ct,st !cos(theta), sin(theta)
        real(kind=wp)::x(1:2)

        x(1:2)=u(2:3)
        
        u(2) = x(1)*ct-x(2)*st
        u(3) = x(2)*ct+x(1)*st

        return
        end subroutine
        
        
!*****************************************************************
        subroutine roty2(ct,st,u)
!----------------------------------------------------
!       rotate vector u around y-axis for an angle of theta radians 
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::ct, st !cos(theta), sin(theta)
        real(kind=wp)::x(1:2)


        x(1)=u(1)
        x(2)=u(3)
        
        u(1) = x(1)*ct+x(2)*st
        u(3) = x(2)*ct-x(1)*st

        return
        end subroutine
        
        
!*****************************************************************

        subroutine rotz2(ct,st,u)
!----------------------------------------------------
!       rotate vector u around y-axis for an angle of theta radians 
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::ct, st !cos(theta), sin(theta)
        real(kind=wp)::x(1:2)
        
        x(1:2)=u(1:2)
        
        u(1) = x(1)*ct-x(2)*st
        u(2) = x(2)*ct+x(1)*st

        return
        end subroutine
!*****************************************************************
        subroutine rotax(theta,a,u)
!----------------------------------------------------
!       rotate vector u around arbitray axis a (||a||=1) for an angle of theta radians
!       input variables are cos(theta) and sin(theta), axis a=(ax,ay,az)
!       ouptut variables are u=(ux,uy,uz)
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::a(1:3),theta 
        
        integer::i
        real(kind=wp)::ct,st !cos(theta), sin(theta)
        real(kind=wp)::rot(1:3,1:3),ct1,axy,axz,ayz,axs,ays,azs
        
        ct=cos(theta)
        st=sin(theta)
        ct1=1._wp-ct
        
        axy=a(1)*a(2)*ct1
        axz=a(1)*a(3)*ct1
        ayz=a(2)*a(3)*ct1
        
        axs=a(1)*st
        ays=a(2)*st
        azs=a(3)*st
        
        do i=1,3
         rot(i,i)=ct+u(i)*u(i)*ct1
        end do
         rot(1,2)=axy-azs
         rot(1,3)=axz+ays
         rot(2,1)=axy+azs
         rot(2,3)=ayz-axs
         rot(3,1)=axz-ays
         rot(3,2)=ayz+axs
        
         u=matmul(rot,a)

        return
        end subroutine
!*****************************************************************
        subroutine rotax2(ct,st,a,u)
!----------------------------------------------------
!       rotate vector u around arbitray axis a (||a||=1) for an angle of theta radians
!       input variables are cos(theta) and sin(theta), axis a=(ax,ay,az)
!       ouptut variables are u=(ux,uy,uz)
!--------------------------------------------------
        implicit none
        real(kind=wp),intent(inout):: u(1:3)
        real(kind=wp),intent(in)::a(1:3),ct, st !cos(theta), sin(theta)
        
        integer::i
        real(kind=wp)::rot(1:3,1:3),ct1,axy,axz,ayz,axs,ays,azs
        
        ct1=1._wp-ct
        axy=a(1)*a(2)*ct1
        axz=a(1)*a(3)*ct1
        ayz=a(2)*a(3)*ct1
        
        axs=a(1)*st
        ays=a(2)*st
        azs=a(3)*st
        
        do i=1,3
         rot(i,i)=ct+u(i)*u(i)*ct1
        end do
         rot(1,2)=axy-azs
         rot(1,3)=axz+ays
         rot(2,1)=axy+azs
         rot(2,3)=ayz-axs
         rot(3,1)=axz-ays
         rot(3,2)=ayz+axs
        
         u=matmul(rot,a)

        return
        end subroutine
!*****************************************************************
subroutine icrf2he(r)
!rotates vectors in ICRF system of JPL Ephemeris to ecliptic system
implicit none
real(kind=wp),dimension(1:3),intent(inout)::r

!Earth's obliquity on J2000 in degrees
!!old  real(kind=wp),parameter::obl=23.439280444444444444444444444444444445_wp
!real(kind=wp),parameter::obl=23.439291111111_wp

!cos(obliquity)
real(kind=wp),parameter::co=0.91748206206918259120186576183186843991279602050781_wp
!sin(obliquity)
real(kind=wp),parameter::so=0.39777715593191187437582811980973929166793823242188_wp

real(kind=wp)::x(1:2)

        x(1:2)=r(2:3)
        
        r(2) = x(1)*co+x(2)*so
        r(3) = x(2)*co-x(1)*so
        
 
!Rotation around x-axis
!erotj(1,:)=(/1._wp,0._wp,0._wp/)
!erotj(2,:)=(/0._wp,co,so/)
!erotj(3,:)=(/0._wp,-so,co/)
   
!rotate positions
!r(1:3)=matmul(erotj(1:3,1:3),r(1:3))
   
return
end subroutine
 
!*********************************************************************
subroutine he2icrf(r)
!rotates heliocentric ecliptic r vectors into ICRS system of JPL Ephemeris at reference J2000
implicit none
real(kind=wp),dimension(1:3),intent(inout)::r

!Earth's obliquity on J2000 in degrees
!!old  real(kind=wp),parameter::obl=23.439280444444444444444444444444444445_wp
!real(kind=wp),parameter::obl=23.439291111111_wp

!cos(obliquity)
real(kind=wp),parameter::co=0.91748206206918259120186576183186843991279602050781_wp
!sin(obliquity)
real(kind=wp),parameter::so=0.39777715593191187437582811980973929166793823242188_wp

real(kind=wp)::x(1:2)

        x(1:2)=r(2:3)
        
        r(2) = x(1)*co-x(2)*so
        r(3) = x(1)*so+x(2)*co

 

!Earth's obliquity on J2000 in degrees
!!  real(kind=wp),parameter::obl=23.439280444444444444444444444444444445_wp
!real(kind=wp),parameter::obl=23.439291111111_wp

!  !rotation matrix (no correction of coordinate center) ecliptic -> J2000 (~ICRF)
! jrote(1,:)=(/1._wp,0._wp,0._wp/)
! jrote(2,:)=(/0._wp,cos(orad),-sin(orad)/)
! jrote(3,:)=(/0._wp,sin(orad),cos(orad)/)
!
!!rotate positions
!   r(1:3)=matmul(jrote(:,:),r(1:3))


return
end subroutine
!*********************************************************************
subroutine rvhe2icrf(rv)
!rotates both, heliocentric ecliptic position vectors and velocities into ICRF system 
!rotation matrix M is 
!       1     0       0
! M=    0 cos(obl) -sin(obl)
!       0 sin(obl)  cos(obl)
!fortunately M=M^-1T, so that we can use the same matrix to transform velocities


implicit none
real(kind=wp),dimension(1:6),intent(inout)::rv

!Earth's obliquity on J2000 in degrees
!!old  real(kind=wp),parameter::obl=23.439280444444444444444444444444444445_wp
!real(kind=wp),parameter::obl=23.439291111111_wp

!cos(obliquity)
real(kind=wp),parameter::co=0.91748206206918259120186576183186843991279602050781_wp
!sin(obliquity)
real(kind=wp),parameter::so=0.39777715593191187437582811980973929166793823242188_wp

real(kind=wp)::x(1:2)

        x(1:2)=rv(2:3)
        
        rv(2) = x(1)*co-x(2)*so
        rv(3) = x(1)*so+x(2)*co

        x(1:2)=rv(5:6)
        
        rv(5) = x(1)*co-x(2)*so
        rv(6) = x(1)*so+x(2)*co

return
end subroutine
!****************************************************************
subroutine rvicrf2he(rv)
!rotates both, ICRF position vectors and velocities into heliocentric ecliptic system 
!rotation matrix M is 
!       1     0       0
! M=    0  cos(obl)  sin(obl)
!       0 -sin(obl)  cos(obl)
!fortunately M=M^-1T, so that we can use the same matrix to transform velocities


implicit none
real(kind=wp),dimension(1:6),intent(inout)::rv

!Earth's obliquity on J2000 in degrees
!!old  real(kind=wp),parameter::obl=23.439280444444444444444444444444444445_wp
!real(kind=wp),parameter::obl=23.439291111111_wp

!cos(obliquity)
real(kind=wp),parameter::co=0.91748206206918259120186576183186843991279602050781_wp
!sin(obliquity)
real(kind=wp),parameter::so=0.39777715593191187437582811980973929166793823242188_wp

real(kind=wp)::x(1:2)

        x(1:2)=rv(2:3)
        
        rv(2) = x(1)*co+x(2)*so
        rv(3) = x(2)*co-x(1)*so

        x(1:2)=rv(5:6)
        
        rv(5) = x(1)*co+x(2)*so
        rv(6) = x(2)*co-x(1)*so

return
end subroutine
!****************************************************************
subroutine invert3x3(m,minv)
 !inverts 3x3 matrix
 implicit none
 real(kind=wp),dimension(3,3),intent(in)::m
 real(kind=wp),dimension(3,3),intent(out)::minv
 real(kind=wp),dimension(3)::x0,x1,x2,dum
 real(kind=wp)::deta
 
 x0(1:3)=m(1:3,1)
 x1(1:3)=m(1:3,2)
 x2(1:3)=m(1:3,3)
 
 call crossp3d(x1,x2,dum)
 
 deta=Dot_Product(x0,dum)
 
 if (deta.eq.0._wp) then
   write(*,*)'Error in subroutine invert, module global_m: could not invert martix' 
 end if
 
 call crossp3d(x1,x2,minv(1,1:3))
 call crossp3d(x2,x0,minv(2,1:3))
 call crossp3d(x0,x1,minv(3,1:3))
 
 minv=minv/deta
  
 return
 end subroutine
 
!************************************************************ 
subroutine crossp3d(a,b,c)
!--------------------------------
!  3 dimensional crossproduct
!  a,b .... input vectors
!  c   .... output vector
!  c = a x b
!------------------------------- 
implicit none
real(kind=wp),dimension(1:3),intent(in)::a,b
real(kind=wp),dimension(1:3),intent(out)::c

 c(1)=a(2)*b(3)-b(2)*a(3)
 c(2)=a(3)*b(1)-b(3)*a(1)
 c(3)=a(1)*b(2)-b(1)*a(2)

return
end subroutine crossp3d
!***************************************************************

subroutine v2basis(a,v,va)
use kind_m, only: wp
!change basis of vector v 
implicit none
real(kind=wp),dimension(:),intent(in)::v
real(kind=wp),dimension(:,:),intent(in)::a
real(kind=wp),dimension(:),intent(out)::va

integer::d

d=SIZE(v)

!change basis by solving linear equation A.va = v corresponding to va = A^-1.v
call gauss(d,a,v,va)

return
end subroutine

!***************************************************************
subroutine gauss(d,a,b,c)
!----------------------------------------------------------
! Gauß'sches Eliminationsverfahren zur Lösung eines 
! Systems Linearer Gleichungen der Form 
!                                    A.c=b
! mit Teilpivotisierung
!
! d....[Int] Dimension des Gleichungssystems
! a...[Real] d x d Koeffizientenmatrix des Gleichungssystems
! b...[Real] d Vektor der rechten Seite der Gleichung
! c...[Real] d  Vektor der Unbekannten
! dependencies: none
!----------------------------------------------------------
implicit none
  integer(kind=ip),intent(in)::d
  integer(kind=ip)::i,j,k,pivot(1:d),amax(1:1),idum
  real(kind=wp)::a(1:d,1:d+1),b(1:d),dum(1:d+1)
  real(kind=wp)::c(1:d)

  do i=1,d
     pivot(i)=i
  end do

  dum(:)=0._wp
  c(:)=0._wp
  a(:,d+1)=b(:)
 
 do i=1,d
!-----------------------------------------------------------------------
!  Partial Pivoting (nur Zeilen werden vertauscht)
!-----------------------------------------------------------------------
    amax(:)=maxloc(abs(a(:,i))) !Welche Zeile hat den größten Wert?
    idum=amax(1)
   if (idum.gt.i) then
      pivot(i)=idum
       dum(:)=a(i,:)                                !Sichere die Zeile die ersetzt wird
      a(i,:)=a(pivot(i),:)                       !Ersetze Zeile i mit jener mit größtem Wert 
       a(pivot(i),:)=dum(:)                 !Beende die Vertauschung, der gesicherte Wert ersetzt die ersetzende Zeile
     end if
!------------------------------------------------------
!   Gauss Elimination
!-----------------------------------------------------     
      a(i,:)=a(i,:)/a(i,i)
      if(i+1.le.d) then
         do j=i+1,d
            a(j,:)=a(j,:)-a(j,i)*a(i,:)         
        end do  
      end if
  end do
  
!-----------------------------------------------
!   Rückwärtssubstitution
!-----------------------------------------------

do i=0,d-1
   j=d-i  
   dum(1)=0._wp
   if (j+1.le.d) then
      do k=j+1,d
         dum(1)=dum(1)+a(j,k)*c(k)
      end do
   end if
   c(j)=(a(j,d+1)-dum(1))
end do

!-----------------------------------------------------------------------
!  Back Pivoting enfällt, da nur Zeilen vertauscht wurden, und sich
!  somit zwar die Reihenfolge der Gleichungen, nicht aber die Koeffizienten
!  der Variablen geändert haben.
!-----------------------------------------------------------------------

return
end subroutine gauss

!******************************************************

function norm(a)
!use const_m, only: wp
implicit none
real(kind=wp),dimension(:),intent(in)::a
real(kind=wp)::norm

norm=sqrt(dot_product(a,a))
return
end function
!******************************************************
function angv1v2(d,v1,v2)
!use const_m, only:wp
!calculate angle between two vectors [rad]
!verified
implicit none
    integer,intent(in)::d !dimension
    real(kind=wp),intent(in),dimension(d)::v1,v2 !two non unit vectors
    real(kind=wp)::angv1v2 !angle between v1 and v2 in rad
    
    angv1v2=acos(dot_product(v1,v2)/sqrt(dot_product(v1,v1)*dot_product(v2,v2)))
    return
end function
!***********************************************
function angv1v2n(d,v1,v2)
use const_m, only: wp
!calculate angle between two normalized vectors [rad]
!verified
implicit none
    integer,intent(in)::d !dimension
    real(kind=wp),intent(in),dimension(d)::v1,v2 !two non unit vectors
    real(kind=wp)::angv1v2n !angle between v1 and v2 in rad
    
    angv1v2n=acos(dot_product(v1,v2))
    return
end function
!******************************************************
end module
