program kepler_tst

      use kepler_m

      implicit none
      real*8,dimension(6)::oe0,oe
      real*8,dimension(6)::rv0,rv
      real*8::epoch,m1,m2,t_obs
      real*8::twopi

      twopi=8.d0*atan(1.d0)

      !Hungaria
      oe0=(/1.944371546951209,.07369922825254498,22.51048102001419,123.6597366228673,175.3768049807603,238.8604173008058/)

      oe0(3:6)=oe0(3:6)*deg2rad
        
      epoch=2455668.5d0
      t_obs=2455669.5d0
      
      m1=1.d0
      m2=0.d0

      write(*,*)'--------------------------------------------------------'
      write(*,*)'testing orbital element to heliocentric state conversion'
      write(*,*)'Keplerian Orbital elements',oe0
      call kema2rv(oe0,rv,m1,m2)
      call rv2kema(rv,oe,m1,m2)

      write(*,*)'difference in orbital elements after conversion'
      write(*,*)oe-oe0
      
      write(*,*)'--------------------------------------------------------'
      write(*,*)'testing two body propagation'
      call two_body_ke(oe0,epoch,t_obs,m1,m2,rv)
      write(*,*)'initial state at epoch',epoch
      write(*,*)rv

      call two_body_rv(rv,t_obs,epoch,m1,m2,rv0)
      write(*,*)'final state at epoch'
      write(*,*)rv0

      write(*,*)'--------------------------------------------------------'



end program
