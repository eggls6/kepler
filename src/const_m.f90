module const_m
!-------------------------------------------------------------------------------------------------------------------------------------------
! This module contains constants commonly used in dynamical astronomy
!
! written by Siegfried EGGL     20161117
!
! dependencies: none
!-------------------------------------------------------------------------------------------------------------------------------------------

use kind_m
implicit none

      PUBLIC


      !PI
      real(kind=wp),parameter::pi=3.1415926535897932384626433832795028841971693993751_wp                !PI
      real(kind=wp),parameter::pix2=6.2831853071795864769252867665590057683943387987502_wp              !2*PI
      real(kind=wp),parameter::pib2=1.5707963267948966192313216916397514420985846996876_wp              !PI/2
     
      !GRAD / RAD
      real(kind=wp),parameter::rad2deg=57.29577951308232286464772187173366546630859375_wp               ![deg] radians to degree conversion
      real(kind=wp),parameter::deg2rad=0.017453292519943295769236907684886127134428718885417_wp         ![] degree to radians conversion

      !GRAVITATIONAL CONSTANTS
      real(kind=wp),parameter::gk=0.01720209895_wp         ![au^(3/2) Msun^(-1/2) D^(-1)] Gaussian gravitational constant
      real(kind=wp),parameter::gk2=0.0002959122082855911487452027497369044795050285756588_wp         !square of Gaussian gravitational constant
      real(kind=wp),parameter::gsi=6.6740831E-11         ![m^3 kg^-1 s^-2] Newtonian gravitational constant

      !EPOCH
      real(kind=wp),parameter::J2000=2451545._wp         !Julian Date of epoch J2000
      real(kind=wp),parameter::jd_mjd=2400000.5_wp       !Julian Date / MJD offset

      !DAY
      real(kind=wp),parameter::spd=86400._wp               ![s] seconds per mean solar day [D]

      !YEAR
      real(kind=wp),parameter::dpgy=365.2568983263280983919685240835_wp  ![D] Gaussian year in days (2 PI/k)
      real(kind=wp),parameter::spgy=31558196.015394747257232666015625_wp ![s] Gaussian year in seconds
    
      real(kind=wp),parameter::dpsidy=365.25636_wp                       ![D] Sidereal year in SI days, Epoch J2000 
      real(kind=wp),parameter::spsidy=31558149.503999996930360794067383_wp                       ![s] Sidereal year in SI days, Epoch J2000 
      
      !SPEED OF LIGHT
      real(kind=wp),parameter::cmps=299792458._wp          ![m/s] speed of light in vacuum
      real(kind=wp),parameter::caupd=173.144632674240341430049738846719264984130859375_wp    ![au/D] speed of light in vacuum
      real(kind=wp),parameter::caupgy=63242.178284130088286474347114563_wp                      ![au / Gaussian year] speed of light in vacuum
      
      !ASTRONOMICAL UNIT
      real(kind=wp),parameter::au=149597870700._wp  ![m] astronomical unit
      real(kind=wp),parameter::aukm=149597870.7_wp   ![km] astronomical unit
      real(kind=wp),parameter::vau2vkm=1731.456836805555440150783397257328033447265625 ![au/D] to [km/s]
      real(kind=wp),parameter::vau2vkm2=2997942.777720700018107891082763671875 ![au/D]^2 to [km/s]^2

      !SUN
      real(kind=wp),parameter::rsun=696342.E3                                                            ![m] mean radius of the Sun (SOHO)
      real(kind=wp),parameter::rsunau=0.00465475876589465387357281416172_wp                              ![au] mean radius of the Sun (SOHO) in au

      real(kind=wp),parameter::teffsun=5777._wp  ![K] effective temperature of the Sun
      real(kind=wp),parameter::zsun=0.02_wp      ![]  metalicity of the sun

      real(kind=wp),parameter::msun=1.98855E30           ![kg] mass of the sun
      real(kind=wp),parameter::gmsun=1.327124400E11      ![km^3/s^2]
      
      !EARTH MOON
      real(kind=wp),parameter::ldm=385000600._wp         ![m] time averaged lunar distance, i.e. distance between the Earth and the moon
      real(kind=wp),parameter::ldkm=385000.6_wp          ![km] time averaged lunar distance, i.e. distance between the Earth and the moon
      real(kind=wp),parameter::ldau=0.0025735700528256248015290807984456478152424097061157_wp ![au] time averaged lunar distance, i.e. distance between the Earth and the moon

      !EARTH
      real(kind=wp),parameter::mekg=5972362487305559718494208._wp           ![kg] Earth mass 
      real(kind=wp),parameter::memsun=0.0000030034896161241036080503073735226138296638964675367_wp ![msun] Earth mass 
      real(kind=wp),parameter::gmearth=398600.435436_wp       ![km^3 s^-2] Earth gravitational parameter  mu=G*M_Earth JPL DE431
      real(kind=wp),parameter::sq2gme=892.86105910830269749567378312349_wp     ![km^(3/2) s^-1] Square root of 2 times Earth gravitational parameter sqrt(2*mu)
      
      real(kind=wp),parameter::rearthkm=6371._wp           ![km] Earth mean radius WGS-84
      real(kind=wp),parameter::rearth=6371000._wp          ![m] Earth mean radius WGS-84
      real(kind=wp),parameter::rep=6356752.3_wp            ![m] Earth polar radius (semi minor axis b) WGS-84
      real(kind=wp),parameter::ree=6378137._wp             ![m] Earth equatorial radius (semi major axis a) WGS-84

      real(kind=wp),parameter::rhillau=0.010003875848834508505147411483449_wp ![au] Earth's Hill radius in au assuming a circular orbit around the sun at 1 au
        
      real(kind=wp),parameter::ve2earth=125.1296297083660391535886446945369243621826171875_wp  ![km^2/s^2] vesc^2 Earth
      real(kind=wp),parameter::veearth=11.186135602090923057971849630121141672134399414062_wp   ![km/s] vesc Earth

      !JUPITER
      real(kind=wp),parameter::gmjup=126712764.8_wp             ![km^3/s^2] GM of Jupiter system JPL DE431
      real(kind=wp),parameter::mjupkg=1898579368902374059958861824._wp          ![kg] mass of Jupiter system  (1.898 x 10^27 kg)                  
      real(kind=wp),parameter::mjupmsun=0.00095479191563363149681392672007973_wp   ![Msun] mass of Jupiter system 
      
      real(kind=wp),parameter::rjup=69911000._wp                   ![m] mean radius of Juptiter
      real(kind=wp),parameter::rjupkm=69911._wp                   ![km] mean radius of Juptiter
      real(kind=wp),parameter::rjupau=0.00046732617030490928525937599502527_wp ![au] mean radius of Jupiter
      
           

end module          
