PROGRAM PCAttempt2
IMPLICIT NONE
!-------------------------Declaring Variables----------------------------!

!Integers for loops and recording time
INTEGER :: i,j,n,x,y,z,tcount

DOUBLE PRECISION :: dt,t

!One Dimension Variables/Constants
DOUBLE PRECISION :: G,AU,AbsD,Yr,KE_i,PE_i,PE,KE,E_i,E 

!Arrays for kinematic calculations
DOUBLE PRECISION :: COM(0:2),COV(0:2),ai(0:2,0:9),da(0:2,0:9),m(0:9),D(0:2)

!Position, velocity, and acceleration vectors, with time history
DOUBLE PRECISION :: r(0:2,0:9,-8:1)
DOUBLE PRECISION :: v(0:2,0:9,-8:1)
DOUBLE PRECISION :: a(0:2,0:9,-8:1)

!Variables for Randomising Initial Positions and Velocities
DOUBLE PRECISION :: Abs(0:9),AbsVel(0:9),Val,RN
INTEGER :: Parity(0:1),ParityVel(0:1)

!-------------------------------Set Up------------------------------------!

!Open data files
OPEN(3,file="../Data/DataStable.txt",STATUS = 'REPLACE')	!This is for the position vectors against time

!Initialise Arrays, Constants etc..
!Constants
AU 	= 1.496e11	!AU in metres
G 	= 6.67e-11	!G in SI
Yr	= 3.154e7	!YR in seconds

n=6			!number of planets in simulation

!1000 seconds timestep
t = 0.			!global time value
dt = 10		    !10 second timesteps
tcount = 0 		!Timestep counter

!Initialising Body parameters Sun to Neptune
!Masses in SI
m(0) = 1.99e30	!Sun
m(1) = 5.97e24	!Earth
m(2) = 0.642e24	!Mars
m(3) = 1.898e27	!etc..
m(4) = 568e24
m(5) = 86.8e24
m(6) = 102e24

!initialise vectors as 0
v = 0.		
a = 0.		
r = 0.

!Absolute Values of Orbital Radii in AU 
abs(1) = 1. 	!Earth
abs(2) = 1.524 	!Mars
abs(3) = 5.203	!etc...
abs(4) = 9.582
abs(5) = 19.20 
abs(6) = 30.05

!Absolute Values of Orbital Velocities
absVel(1) = 29.78e3	!Earth
absVel(2) = 24.1e3	!etc...
absVel(3) = 13.07e3
absVel(4) = 9.536e3
absVel(5) = 6.687e3
absVel(6) = 5.372e3

!Producing random positions and velocities in orbit

do j=1,n 		!Required For Earth to Final planet in system
	!Parity
	do i=0,1	!z axis ignored; all bodies on x-y plane initially
		!Get random number to determine whether in positive or negative plane, for x and y axes
		call random_number(val)
		if (val > 0.5) then
			parity(i) = 1
		else 
			parity(i) = -1
		end if
	end do
	!Get random number to determine multiplicative factor to orbital radius
	call random_number(RN)
	
	!Apply Random number along one axis - produces corresponding distance on other axis - using Trigonometry
	r(0,j,0) = parity(0)*RN*abs(j)*AU
	r(1,j,0) = parity(1)*abs(j)*SQRT(1 - RN)*AU
	
	!Algorithm to work out what parities velocities should take in corresponding positions - uses parity vector used for position
	if (parity(0)==parity(1)) then 
		if (parity(0) == 1) then
			ParityVel = [-1,1]
		else
			parityVel = [1,-1]
		end if
	else
		if (parity(0) == 1) then
			ParityVel = [1,1]
		else
			ParityVel = [-1,-1]
		end if
	end if
	
	!Same multiplicative factor (RN) required for velocity vector, but on opposite axis
	v(1,j,0) = parityVel(1)*RN*absvel(j)
	v(0,j,0) = parityvel(0)*absvel(j)*SQRT(1 - RN)	!Trigonometry to calculate the other axis
	
end do

!Get CoM - Then centre system on CoM
do i=0,2
	do j=0,n
		COM(i) = COM(i) + m(j)*r(i,j,0)
	end do
end do 
COM = COM/sum(m(0:n))
do i=0,n
	r(0:2,i,0) = r(0:2,i,0) - COM(0:2)
end do

!Get CoV - Then centre velocities on CoV
do i=0,2
	do j=0,n
		COV(i) = COV(i) + m(j)*v(i,j,0)
	end do
end do 
COV = COV/sum(m(0:n))
do i=0,n
	v(0:2,i,0) = v(0:2,i,0) - COV(0:2)
end do

!-------------------Bootstrap-----------------!

!Use second order method for 8 timesteps, so time history can be filled for acceleration and velocity

!initial acceleration
do i=0,n	!getting acceleration
	do j=0,n
		if (i==j) then
			cycle
		end if 
		D = r(0:2,i,0) - r(0:2,j,0)	!Get distance between two objects		
		AbsD = SQRT(SUM(D**2))	!Get absolute distance
		a(0:2,i,0) = a(0:2,i,0) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
	end do
end do  

do tcount=0,100000

	do i=-8,-1,1
		r(0:2,0:n,i) = r(0:2,0:n,i+1)
		v(0:2,0:n,i) = v(0:2,0:n,i+1)
		a(0:2,0:n,i) = a(0:2,0:n,i+1)
	end do

	ai = a(0:2,0:n,0)

	!Change in r
	
	r(0:2,0:n,0) = r(0:2,0:n,0) + v(0:2,0:n,0)*(dt) + 0.5*a(0:2,0:n,0)*((dt)**2) 
	
	do i=0,n	!getting acceleration
		do j=0,n
			if (i==j) then
				cycle
			end if 
			D = r(0:2,i,0) - r(0:2,j,0)	!Get distance between two objects		
			AbsD = SQRT(SUM(D**2))	!Get absolute distance
			a(0:2,i,0) = a(0:2,i,0) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
		end do
	end do
	
	da(0:2,0:n)  = a(0:2,0:n,0) - ai(0:2,0:n) !Change in acceleration
	v(0:2,0:n,0) = v(0:2,0:n,0) + a(0:2,0:n,0)*dt + 0.5*da(0:2,0:n)*dt
	
	write(3,*) (t/Yr), r(0:1,0:n,0)/AU
	
	t = t + dt
end do

!----------------------Predictor-----------------------!

do tcount=0,0
	exit
	
	r(0:2,0:n,1) = r(0:2,0:n,0) + (dt/24)*(-(9*v(0:2,0:n,-3))+(37*v(0:2,0:n,-2))-(59*v(0:2,0:n,-1))+(55*v(0:2,0:n,0)))
	v(0:2,0:n,1) = v(0:2,0:n,0) + (dt/24)*(-(9*a(0:2,0:n,-3))+(37*a(0:2,0:n,-2))-(59*a(0:2,0:n,-1))+(55*a(0:2,0:n,0)))
	
	do i=0,n	!getting acceleration
		do j=0,n
			if (i==j) then
				cycle
			end if 
			D = r(0:2,i,1) - r(0:2,j,1)	!Get distance between two objects		
			AbsD = SQRT(SUM(D**2))	!Get absolute distance
			a(0:2,i,1) = a(0:2,i,1) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
		end do
	end do
	
	
	if (a(0,1,1) == a(0,1,0)) then
		write (6,*) tcount
		exit
	end if
	
	write(3,*) (t/Yr), r(0:1,0:n,0)
	
	do i=-8,0,1
		r(0:2,0:n,i) = r(0:2,0:n,i+1)
		v(0:2,0:n,i) = v(0:2,0:n,i+1)
		a(0:2,0:n,i) = a(0:2,0:n,i+1)
	end do

	t = t + dt

end do
	
	
write(6,*) " " 
write(6,*) tcount
write(6,*) " " 	
write(6,*) r(0:2,1,-8:1)/AU
write(6,*) " " 	 
write(6,*) v(0:2,1,-8:1) 
write(6,*) " " 
write(6,*) a(0:2,1,-8:1) 
write(6,*) " " 	

END PROGRAM PCAttempt2