PROGRAM PC Attempt
IMPLICIT NONE
!-------------------------Declaring Variables----------------------------!

!Integers for counting
INTEGER :: i,j,n,x,y,z,t,dt,tcount,fifty

!One Dimension Variables/Constants
DOUBLE PRECISION :: G,AU,AbsD,Yr,KE_i,PE_i,PE,KE,E_i,E 

!Arrays for kinematic calculations
DOUBLE PRECISION :: COM(0:2),COV(0:2),a(0:2,0:9),ai(0:2,0:9),da(0:2,0:9),m(0:9),v(0:2,0:9),D(0:2),r(0:2,0:9)

!For getting eccentricity
DOUBLE PRECISION :: Apihelion(1:9),Perihelion(1:9),ecc(1:9)

!Time History of velocity and acceleration vectors (Will delete original vectors later)
DOUBLE PRECISION :: vHist(0:2,0:9,0:9)
DOUBLE PRECISION :: aHist(0:2,0:9,0:9)

!Variables for Randomising Initial Positions and Velocities
DOUBLE PRECISION :: Abs(0:9),AbsVel(0:9),Val,RN
INTEGER :: Parity(0:1),ParityVel(0:1)


!-------------------------------Set Up------------------------------------!


!Open data files
OPEN(3,file="../Data/DataStable.txt",STATUS = 'REPLACE')	!This is for the position vectors against time
OPEN(4,file="../Data/EnergyStable.txt",STATUS = 'REPLACE')	!This is for System energy against time
OPEN(1,file="../Data/Eccentricities.txt",STATUS = 'REPLACE')!This is for orbital eccentricity against planetary mass.. just kidding, it's against time

!Initialise Arrays, Constants etc..
!Constants
AU = 1.496e11	!AU in metres
G = 6.67e-11	!G in SI
Yr = 3.154e7	!YR in seconds
n=7 			!number of planets in simulation

!1000 seconds timestep,
dt = 1000		!1000 second timesteps
tcount = 1000 	!Timestep counter, initialised at 1000 so first values are written
Fifty = 0		!50Yr Counter, for command line output (Every 50 years) 

!Initialising Body parameters Sun to Neptune
!Masses in SI
m(0) = 1.99e30	!Sun
m(1) = 5.97e24	!Earth
m(2) = 0.642e24	!Mars
m(3) = 1.898e27	!etc..
m(4) = 568e24
m(5) = 86.8e24
m(6) = 102e24

!Positions initially set as 0
r(0:2,0:n) = 0.

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
	r(0,j) = parity(0)*RN*abs(j)
	r(1,j) = parity(1)*(SQRT(abs(j)**2 - r(0,j)**2))
	
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
	v(1,j) = parityVel(1)*RN*absvel(j)
	v(0,j) = parityvel(0)*(SQRT(absvel(j)**2 - v(1,j)**2))	!Trigonometry to calculate the other axis
	
end do

!Initialising Max/Min Api/Perihelion for each body
do i=1,n !Only required for planets, not the sun
	Apihelion(i) = abs(i)	!Give both values the absolute distance initially
	Perihelion(i) = abs(i)
end do

!Get CoM - Then centre system on CoM
do i=0,2
	do j=0,n
		COM(i) = COM(i) + m(j)*r(i,j)
	end do
end do 
COM = COM/sum(m(:))
do i=0,n
	r(0:2,i) = r(0:2,i) - COM(0:2)
end do

!Get CoV - Then centre velocities on CoV
do i=0,2
	do j=0,n
		COV(i) = COV(i) + m(j)*v(i,j)
	end do
end do 
COV = COV/sum(m(0:2))
do i=0,n
	v(0:2,i) = v(0:2,i) - COV(0:2)
end do

!For Energy Conservation
!Initial KE
KE_i=0.
do i=1,n
	KE_i = KE_i + 0.5*m(i)*sum(v(0:2,i)**2)
end do

!Initial PE
PE_i = 0
do i=0,n
	do j=0,n
		if (i==j) then
			cycle
		end if 
		D = r(0:2,i) - r(0:2,j)	
		AbsD = SQRT(SUM(D**2))
		PE_i = PE_i - (G*m(i)*m(j))/(AbsD*AU)
	end do
end do
PE_i = 0.5*PE_i	

!Total Initial Energy
E_i = PE_i + KE_i

!-------------------------------Bootstrap------------------------------------!

END PROGRAM PC Attempt