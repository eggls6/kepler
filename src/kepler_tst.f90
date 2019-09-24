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

      call two_body_ke(oe0,epoch,t_obs,m1,m2,rv)

      write(*,*)rv

      call two_body_rv(rv,t_obs,epoch,m1,m2,rv0)
      write(*,*)rv0
end program
