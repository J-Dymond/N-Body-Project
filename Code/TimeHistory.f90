PROGRAM StableSystem
IMPLICIT NONE

!Declaring Variables
INTEGER :: i,j,n,x,y,z,t,dt,tcount,fifty

DOUBLE PRECISION :: G,AU,AbsD,Yr,KE_i,PE_i,PE,KE,E_i,E 

!Initial Conditions for Perturber
DOUBLE PRECISION :: CA,MP,VP

DOUBLE PRECISION :: COM(0:2),COV(0:2),a(0:2,0:9),ai(0:2,0:9),da(0:2,0:9),m(0:9),v(0:2,0:9),D(0:2)

!RandomisingPositions
DOUBLE PRECISION :: Abs(0:9),AbsVel(0:9),Val,RN
INTEGER :: Parity(0:1),ParityVel(0:1)

!For getting eccentricity
DOUBLE PRECISION :: Apihelion(1:9),Perihelion(1:9),ecc(1:9)

!Position vector
DOUBLE PRECISION :: r(0:2,0:9)

!Time History of velocity and acceleration vectors (Will delete original vectors later)
DOUBLE PRECISION :: vHist(0:2,0:9,0:9)
DOUBLE PRECISION :: aHist(0:2,0:9,0:9)

 
!Open data files
OPEN(3,file="../Data/DataStable.txt",STATUS = 'REPLACE')		!orbit data	
!OPEN(2,file="../Data/EnergyStableF.txt",STATUS = 'REPLACE')	!First Derivative energies
OPEN(4,file="../Data/EnergyStable.txt",STATUS = 'REPLACE')		!Second Derivative Energies	
OPEN(1,file="../Data/Eccentricities.txt",STATUS = 'REPLACE')	!Eccentricities

!Initialising Constants
AU = 1.496e11
G = 6.67e-11
Yr = 3.154e7
n=7


!1000 seconds timestep, initialise counter for data writing
dt = 1000
tcount = 1000
Fifty = 0

!Initialising planets Sun to Neptune
!Masses
m(0) = 1.99e30
m(1) = 5.97e24
m(2) = 0.642e24
m(3) = 1.898e27
m(4) = 568e24
m(5) = 86.8e24
m(6) = 102e24

!Initialising Positions 
!Sun
r(0:2,0) = 0.

!Planets
r(1:2,:) = 0.

vHist = 0.

!Orbital radii in AU 
abs(1) = 1.
abs(2) = 1.524 
abs(3) = 5.203
abs(4) = 9.582
abs(5) = 19.20 
abs(6) = 30.05

!Orbital velocities - moving parallel to y axis initially
absVel(1) = 29.78e3
absVel(2) = 24.1e3
absVel(3) = 13.07e3
absVel(4) = 9.536e3
absVel(5) = 6.687e3
absVel(6) = 5.372e3

!Producing random positions and velocities in orbit
do j=1,n
	!Parity
	do i=0,1
		call random_number(val)
		if (val > 0.5) then
			parity(i) = 1
		else 
			parity(i) = -1
		end if
	end do
	
	call random_number(RN)
	
	!Set Values for distance
	r(0,j) = parity(0)*RN*abs(j)
	r(1,j) = parity(1)*(SQRT(abs(j)**2 - r(0,j)**2))
	
	!Algorithm to work out what parities velocities should take in corresponding positions
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
	
	v(1,j) = parityVel(1)*RN*absvel(j)
	v(0,j) = parityvel(0)*(SQRT(absvel(j)**2 - v(1,j)**2))
	
end do

!Initialising Api/Perihelion
do i=1,n
	Apihelion(i) = abs(i)
	Perihelion(i) = abs(i)
end do

!Get COM - Put in COM 
do i=0,2
	do j=0,n
		COM(i) = COM(i) + m(j)*r(i,j)
	end do
end do 
COM = COM/sum(m(:))

do i=0,n
	r(0:2,i) = r(0:2,i) - COM(0:2)
end do

!Get COV - Put in COV
do i=0,2
	do j=0,n
		COV(i) = COV(i) + m(j)*v(i,j)
	end do
end do 
COV = COV/sum(m(0:2))

do i=0,n
	v(0:2,i) = v(0:2,i) - COV(0:2)
end do

!Initialise Perturber parameters - In solar masses and AU - Do this after getting COM of the solar system
MP = 5
CA = 160.
VP = 30

!Setting mass
m(7) = m(0)*MP 
!Intended closest approach value - assign to x axis
r(0,7) = CA
r(1,7) = -15*r(0,7)
!Velocity parallel to y-axis like that of planets
v(1,7) = VP*1e3
Perihelion(7) = SQRT(SUM((r(0:2,7))**2))

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


!Makes output more readable
print *, ''

do
	t = 0
	do
		!vHist(0:2,0:9,0:9)
		!Fill time history vectors
		do i=9,1,-1
			vHist(0:2,0:n,i) = vHist(0:2,0:n,i-1)
			!aHist(0:2,0:n,i) = aHist(0:2,0:n,i-1)
		end do
		vHist(0:2,0:n,0) = v(0:2,0:n)
		!aHist(0:2,0:n,0) = a(0:2,0:n)

		!Calculating individual accelerations
		!To get direction multiply by the distance vector which stores a plus/minus
		!Return values to zero as to not increment acceleration
		
		!Store previous values for acceleration to use in the 2nd derivative position calculation
		ai = a
		a = 0.
	
		do i=0,n
			do j=0,n
				if (i==j) then
					cycle
				end if 
				D = r(0:2,i) - r(0:2,j)	
				AbsD = SQRT(SUM(D**2))
				
				!Getting Apihelion and Perihelion, only needs to be done when calculating a from Sun
				if(j==0) then
					if (AbsD > Apihelion(i)) then
						Apihelion(i) = AbsD
					end if
					if (AbsD < Perihelion(i)) then
						Perihelion(i) = AbsD
					end if
					ecc(i) = (Apihelion(i) - Perihelion(i))/(Apihelion(i) + Perihelion(i))
				end if
				
				a(0:2,i) = a(0:2,i) - ((G*m(j))/((AbsD*AU)**3))*(D*AU)
			end do
		end do	
	
		!Calculate the velocities - with second derivatives
		da(0:2,0:n) = a(0:2,0:n) - ai(0:2,0:n)
		
		do i=0,n
			r(0:2,i) = r(0:2,i) + v(0:2,i)*(dt/AU) + 0.5*a(0:2,i)*((dt/AU)**2)
			v(0:2,i) = v(0:2,i) + a(0:2,i)*dt + 0.5*da(0:2,i)*dt
		end do
		
		!And write radii every 1000th time-steps
		if (tcount >= 1000) then
			write(3,*) (t/Yr + 50*Fifty), r(0:1,0:n), Perihelion(7)
			tcount = 0

		
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
			write(4,*) (t/Yr + 50*Fifty), (E-E_i)/E_i,KE/KE_i,PE/PE_i
			!write(6,*) Perihelion(7)
			write(1,*) (t/Yr + 50*Fifty), ecc(1:n)
	
		end if
		tcount = tcount + 1
	
		!Increment time
		t = t + dt

	
		!Check time has not exceded time limit
		if (t>=50*Yr) exit
	

	end do 

	Fifty = Fifty + 1
	print *, Fifty, Perihelion(7)
	
	if (Fifty>=20) exit
end do
	

!Close data files
close(3)	!orbit data
close(4)	!Second Derivative energies
!close(2)  	!First Derivative Energies
close(1)	!Eccentricities

print *, ''
!print *, a(:,1)
!print *, v(:,1)
!print *, r(:,1)

END PROGRAM StableSystem