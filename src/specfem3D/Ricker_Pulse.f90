program Ricker_Pulse_Marcel
implicit none

integer, parameter :: Nt = 20000
integer            :: it
double precision   :: dt,t

dt=0.003d0

do it = 1, Nt
  t=dble(it-1)*dt
  write(*,*) t,Pt_Ricker_Pulse_M(t)

enddo


contains

double precision  function Pt_Ricker_Pulse_M(t)
implicit none

integer          :: Nt ! Number points
double precision :: t


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% INPUT DATA SPECIFIED BY USER %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double precision, parameter :: AMP = 1.d0           !% peak amplitude of the pulse [a.u]
double precision, parameter :: freq = 0.2d0        !% center frequency of the pulse [Hz]
double precision, parameter :: t0 = 1/freq       !% time when peak of the pulse is reached [s]
double precision, parameter :: PI = 3.141592653589793d0

!variable declaration


!%%%%%%%%%%%%%%%%%%%%%%%
!% WAVEFORM GENERATION %
!%%%%%%%%%%%%%%%%%%%%%%%

!% time history of the pulse
Pt_Ricker_Pulse_M = AMP*(1.d0-2.d0*((t-t0)*freq*pi)**2)*exp(-((t-t0)*freq*pi)**2);

end function

end program
