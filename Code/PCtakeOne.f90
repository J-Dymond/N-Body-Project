PROGRAM PCAttempt
IMPLICIT NONE
!-------------------------Declaring Variables----------------------------!

!Integers for loops and recording time
INTEGER :: i,j,n,x,y,z,t,tcount,fifty

DOUBLE PRECISION :: dt

!One Dimension Variables/Constants
DOUBLE PRECISION :: G,AU,AbsD,Yr,KE_i,PE_i,PE,KE,E_i,E 

!Arrays for kinematic calculations
DOUBLE PRECISION :: COM(0:2),COV(0:2),a(0:2,0:9),ai(0:2,0:9),da(0:2,0:9),m(0:9),v(0:2,0:9),D(0:2),r(0:2,0:9)

!For getting eccentricity
DOUBLE PRECISION :: Apihelion(1:9),Perihelion(1:9),ecc(1:9)

!Time History of position, velocity, and acceleration vectors (Will delete original vectors later)
DOUBLE PRECISION :: rHist(0:2,0:9,0:8)
DOUBLE PRECISION :: vHist(0:2,0:9,0:8)
DOUBLE PRECISION :: aHist(0:2,0:9,0:8)

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
AU 	= 1.496e11	!AU in metres
G 	= 6.67e-11	!G in SI
Yr	= 3.154e7	!YR in seconds

n=6			!number of planets in simulation

!1000 seconds timestep
t = 0			!global time value
dt = 10		    !10 second timesteps
tcount = 0 		!Timestep counter
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

!Positions,Velocity, and accelerations initially set as 0 - hardcode to length of array rather than using n
r = 0.
v = 0.
a = 0.

vhist = 0.		!initialise time histories as 0
ahist = 0.		
rhist = 0.

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
COM = COM/sum(m(0:n))
do i=0,n
	r(0:2,i) = r(0:2,i) - COM(0:2)
end do

!Get CoV - Then centre velocities on CoV
do i=0,2
	do j=0,n
		COV(i) = COV(i) + m(j)*v(i,j)
	end do
end do 
COV = COV/sum(m(0:n))
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


write(6,*) " " 
!-------------------------------Bootstrap------------------------------------!

!Use second order method for 8 timesteps, so time history can be filled for acceleration and velocity

!initial acceleration
do i=0,n	!getting acceleration
	do j=0,n
		if (i==j) then
			cycle
		end if 
		D = r(0:2,i) - r(0:2,j)	!Get distance between two objects		
		AbsD = SQRT(SUM(D**2))	!Get absolute distance
		a(0:2,i) = a(0:2,i) - ((G*m(j))/((AbsD*AU)**3))*(D*AU) !Make calculation
	end do
end do

do
	t = 0
	do
	
		do i=8,2,-1
			rhist(0:2,0:n,i) = rhist(0:2,0:n,i-1)
			vhist(0:2,0:n,i) = vhist(0:2,0:n,i-1)
			ahist(0:2,0:n,i) = ahist(0:2,0:n,i-1)
		end do
		
		ai = a	!Previous acceleration
		a = 0.	!Current acceleration
	
		do i=0,n	!getting acceleration
			do j=0,n
				if (i==j) then
					cycle
				end if 
				D = r(0:2,i) - r(0:2,j)	!Get distance between two objects		
				AbsD = SQRT(SUM(D**2))	!Get absolute distance
				a(0:2,i) = a(0:2,i) - ((G*m(j))/((AbsD*AU)**3))*(D*AU) !Make calculation
			end do
		end do	
	
		da(0:2,0:n) = a(0:2,0:n) - ai(0:2,0:n) !Change in acceleration

		!Change in r, and change in v, using the da in v calculation
		r(0:2,0:n) = r(0:2,0:n) + v(0:2,0:n)*(dt/AU) + 0.5*a(0:2,0:n)*((dt/AU)**2)
		v(0:2,0:n) = v(0:2,0:n) + a(0:2,0:n)*dt + 0.5*da(0:2,0:n)*dt
		
		rhist(0:2,0:n,1) = r(0:2,0:n)
		vhist(0:2,0:n,1) = v(0:2,0:n)
		ahist(0:2,0:n,1) = a(0:2,0:n)
	
		tcount = tcount + 1
	
		t = t + dt
	
		if (tcount >= 100000) then
			write(3,*) (t/Yr + 50*Fifty), r(0:1,0:n)
			tcount=0
		
			!Get KE
			KE=0.
			do i=1,n
				KE = KE + 0.5*m(i)*sum(v(0:2,i)**2)
			end do
	
			!Get PE
			PE = 0
			do i=0,n
				do j=0,n
					if (i==j) then
						cycle
					end if 
					D = r(0:2,i) - r(0:2,j)	
					AbsD = SQRT(SUM(D**2))
					PE = PE - (G*m(i)*m(j))/(AbsD*AU)
				end do
			end do
			PE = 0.5*PE	
	
			!Total Energy
			E = KE + PE
			write(4,*) (t/Yr + 50*Fifty), (E-E_i)/E_i,(KE-KE_i)/KE_i,(PE-PE_i)/PE_i
			!write(6,*) Perihelion(7)
			write(1,*) (t/Yr + 50*Fifty), ecc(1:n)
		end if 
	
		if (t>=50*Yr) then
			exit
		end if

	end do
	
	fifty = fifty + 1
	write(6,*) fifty, (fifty)/(yr)
	
	
	if (fifty>=1)then
		exit
	end if
	
end do

write(6,*) " " 	
write(6,*) r(0:2,1)
write(6,*) " " 	
write(6,*) rhist(0:2,1,0:8)
write(6,*) " " 	 
write(6,*) v(0:2,1)
write(6,*) " " 	
write(6,*) vhist(0:2,1,0:8) 
write(6,*) " " 
write(6,*) a(0:2,1)
write(6,*) " " 	
write(6,*) ahist(0:2,1,0:8) 
write(6,*) " " 	

!-----------------------------Predictor-Corrector Algorithm------------------------------------
tcount = 0
dt=10

do
	t = 0
	do
		!if (tcount<10) then
		!	write(6,*) a(0:2,1)
		!end if	
	
		do i=8,2,-1
			rhist(0:2,0:n,i) = rhist(0:2,0:n,i-1)
			vhist(0:2,0:n,i) = vhist(0:2,0:n,i-1)
			ahist(0:2,0:n,i) = ahist(0:2,0:n,i-1)
		end do
		
		do i=0,n	!getting acceleration
			do j=0,n
				if (i==j) then
					cycle
				end if 
				D = r(0:2,i) - r(0:2,j)	!Get distance between two objects		
				AbsD = SQRT(SUM(D**2))	!Get absolute distance
				a(0:2,i) = a(0:2,i) - ((G*m(j))/((AbsD*AU)**3))*(D*AU) !Make calculation
			end do
		end do

		r(0:2,0:n) = r(0:2,0:n) + (dt/(AU*24))*(-(9*vhist(0:2,0:n,3))+(37*vhist(0:2,0:n,2))-(59*vhist(0:2,0:n,1))+(55*v(0:2,0:n)))
		v(0:2,0:n) = v(0:2,0:n) + (dt/24)*(-(9*ahist(0:2,0:n,3))+(37*ahist(0:2,0:n,2))-(59*ahist(0:2,0:n,1))+(55*a(0:2,0:n)))
		
		ahist(0:2,0:n,1) = a(0:2,0:n)
		rhist(0:2,0:n,1) = r(0:2,0:n)
		vhist(0:2,0:n,1) = v(0:2,0:n)
		

		tcount = tcount + 1
	
		t = t + dt

		if (tcount >= 100000) then
			write(3,*) (t/Yr + 50*Fifty), rhist(0:1,0:n,1)
			tcount=0

			!Get KE
			KE=0.
			do i=1,n
				KE = KE + 0.5*m(i)*sum(v(0:2,i)**2)
			end do

			!Get PE
			PE = 0
			do i=0,n
				do j=0,n
					if (i==j) then
						cycle
					end if 
					D = r(0:2,i) - r(0:2,j)	
					AbsD = SQRT(SUM(D**2))
					PE = PE - (G*m(i)*m(j))/(AbsD*AU)
				end do
			end do
			PE = 0.5*PE	

			!Total Energy
			E = KE + PE
			write(4,*) (t/Yr + 50*Fifty), (E-E_i)/E_i,(KE-KE_i)/KE_i,(PE-PE_i)/PE_i
			!write(6,*) Perihelion(7)
			write(1,*) (t/Yr + 50*Fifty), ecc(1:n)
		end if

		if (t>=50*Yr) then
			exit
		end if

	end do

	fifty = fifty + 1
	write(6,*) fifty

	if (fifty>=1) then
		exit
	end if

end do

write(6,*) " " 	
write(6,*) r(0:2,1)
write(6,*) " " 	
write(6,*) rhist(0:2,1,0:8)
write(6,*) " " 	 
write(6,*) v(0:2,1)
write(6,*) " " 	
write(6,*) vhist(0:2,1,0:8) 
write(6,*) " " 
write(6,*) a(0:2,1)
write(6,*) " " 	
write(6,*) ahist(0:2,1,0:8) 

END PROGRAM PCAttempt