PROGRAM PCAttempt2
IMPLICIT NONE
!-------------------------Declaring Variables----------------------------!

!Integers for loops and recording time
INTEGER :: i,j,n,x,y,z,tcount,tchange

DOUBLE PRECISION :: dt,t,tWrite,tPrint

!One Dimension Variables/Constants
DOUBLE PRECISION :: G,AU,AbsD,Yr,KE_i,PE_i,PE,KE,E_i,E 

!For getting error and difference between predicted values
!first axis - position or velocity, second axis - position axis, third axis - body
DOUBLE PRECISION :: Err, Small(1:2), RealErr(1:2,1:3,1:9), AbsDifference(1:2,1:3,1:9)

!Arrays for kinematic calculations
DOUBLE PRECISION :: COM(1:3),COV(1:3),ai(1:3,1:9),da(1:3,1:9),m(1:9),D(1:3)

!Position, velocity, and acceleration vectors, with time history
DOUBLE PRECISION :: r(1:3,1:9,-8:0) 
DOUBLE PRECISION :: v(1:3,1:9,-8:0)
DOUBLE PRECISION :: a(1:3,1:9,-8:0)

!When halving the timestep it necessary to have placeholder values
DOUBLE PRECISION :: a_place(1:3,1:9,-4:-1)
DOUBLE PRECISION :: v_place(1:3,1:9,-4:-1)
DOUBLE PRECISION :: r_place(1:3,1:9,-4:-1)

DOUBLE PRECISION :: rp(1:3,1:9)
DOUBLE PRECISION :: vp(1:3,1:9)

!Variables for Randomising Initial Positions and Velocities
DOUBLE PRECISION :: Abs(1:9),AbsVel(1:9),Val,RN
INTEGER :: Parity(1:2),ParityVel(1:2)

!-------------------------------Set Up------------------------------------!

!Open data files
OPEN(3,file="../Data/DataStable.txt",STATUS = 'REPLACE') !This is for the position vectors against time
OPEN(4,file="../Data/EnergyStable.txt",STATUS = 'REPLACE')

!Initialise Arrays, Constants etc..
!Constants
AU = 1.496e11 !AU in metres
G = 6.67e-11  !G in SI
Yr = 3.154e7  !YR in seconds

!error calculation
RealErr = 0.
Small = 1e-4
Err = 1e-11

n=8  !number of bodies in simulation

!1000 seconds timestep
t = 0.       !global time value
dt = 0.01   !10 second timesteps
tcount = 0   !Timestep counter
tchange = 8  !Counter for when timestep is permitted to change, initialise at 8
tPrint = 1e4 !Time counter for print outs: Every 10,000 years 

!Initialising Body parameters Sun to Neptune
!Masses in SI

m(1) = 1.99e30  !Sun
m(2) = 5.97e24  !Earth
m(3) = 0.642e24 !Mars
m(4) = 1.898e27 !etc..
m(5) = 568e24
m(6) = 86.8e24
m(7) = 102.0e24
m(8) = 4.87e24  !Venus - put in last

!initialise vectors as 0
v = 0.
a = 0.
r = 0.

!Absolute Values of Orbital Radii in AU 
abs(2) = 1.      !Earth
abs(3) = 1.524   !Mars
abs(4) = 5.203   !etc...
abs(5) = 9.582
abs(6) = 19.20 
abs(7) = 30.05
abs(8) = 0.723

!Absolute Values of Orbital Velocities
absVel(2) = 29.78e3 !Earth
absVel(3) = 24.1e3  !etc...
absVel(4) = 13.07e3
absVel(5) = 9.536e3
absVel(6) = 6.687e3
absVel(7) = 5.372e3
absVel(8) = 35.0e3

!Producing random positions and velocities in orbit 


do j=2,n 		!Required For Earth to Final planet in system
	
	!Parity
	do i=1,2	!z axis ignored; all bodies on x-y plane initially
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
	r(1,j,-8) = parity(1)*RN*abs(j)*AU
	r(2,j,-8) = parity(2)*abs(j)*SQRT(1 - RN*RN)*AU
	
	!Algorithm to work out what parities velocities should take in corresponding positions - uses parity vector used for position
	if (parity(1)==parity(2)) then 
		if (parity(1) == 1) then
			ParityVel = [-1,1]
		else
			parityVel = [1,-1]
		end if
	else
		if (parity(1) == 1) then
			ParityVel = [1,1]
		else
			ParityVel = [-1,-1]
		end if
	end if
	
	!Same multiplicative factor (RN) required for velocity vector, but on opposite axis
	v(2,j,-8) = parityVel(2)*RN*absvel(j)
	v(1,j,-8) = parityvel(1)*absvel(j)*SQRT(1 - RN*RN)	!Trigonometry to calculate the other axis
	
end do


!Get CoM - Then centre system on CoM
do i=1,3
    do j=1,n
        COM(i) = COM(i) + m(j)*r(i,j,-8)
    end do
end do 
COM = COM/sum(m(1:n))
do i=1,n
    r(1:3,i,-8) = r(1:3,i,-8) - COM(1:3)
end do

!Get CoV - Then centre velocities on CoV
do i=1,3
    do j=1,n
        COV(i) = COV(i) + m(j)*v(i,j,-8)
    end do
end do 
COV = COV/sum(m(1:n))
do i=1,n
    v(1:3,i,-8) = v(1:3,i,-8) - COV(1:3)
end do


!For Energy Conservation
!Initial KE
KE_i=0.
do i=1,n
	KE_i = KE_i + 0.5*m(i)*sum(v(1:3,i,-8)**2)
end do


!Initial PE
PE_i = 0
do i=1,n
    do j=1,n
		if (i==j) then
			cycle
		end if 
		D = r(1:3,i,-8) - r(1:3,j,-8)	
		AbsD = SQRT(SUM(D**2))
		PE_i = PE_i - (G*m(i)*m(j))/(AbsD)
	end do
end do
PE_i = 0.5*PE_i	

!Total Initial Energy
E_i = PE_i + KE_i


!Makes output more readable
print *, ''

write(3,*) (t/Yr), r(1:2,1:n,-8)/AU

!-------------------Bootstrap-----------------!

!Use second order method for 8 timesteps, so time history can be filled for acceleration and velocity

!initial acceleration
do i=1,n  !getting acceleration
    do j=1,n
        if (i==j) then
            cycle
        end if 
        D = r(1:3,i,-8) - r(1:3,j,-8) !Get distance between two objects		
        AbsD = SQRT(SUM(D**2)) !Get absolute distance
        a(1:3,i,-8) = a(1:3,i,-8) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
    end do
end do  

do tcount=-7,0,1

    ai = a(1:3,1:n,tcount-1)

    !Change in r
	
    do i=1,n
        r(1:3,i,tcount) = r(1:3,i,tcount-1) + v(1:3,i,tcount-1)*(dt) + 0.5*a(1:3,i,tcount-1)*((dt)**2) 
    end do
	
    do i=1,n    !getting acceleration
        do j=1,n
            if (i==j) then
                cycle
            end if 
            D = r(1:3,i,tcount) - r(1:3,j,tcount) !Get distance between two objects		
            AbsD = SQRT(SUM(D**2)) !Get absolute distance
            a(1:3,i,tcount) = a(1:3,i,tcount) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
        end do
    end do
	
    do i=1,n
        da(1:3,i)  = a(1:3,i,tcount) - ai(1:3,i) !Change in acceleration
        v(1:3,i,tcount) = v(1:3,i,tcount-1) + a(1:3,i,tcount-1)*dt + 0.5*da(1:3,i)*dt
    end do
	
    t = t + dt
end do

write(6,*) " "
write(6,*) tcount
write(6,*) " "
write(6,*) r(1:3,1,-8:0)/AU
write(6,*) " " 
write(6,*) v(1:3,1,-8:0) 
write(6,*) " "
write(6,*) a(1:3,1,-8:0) 
write(6,*) " "

!---------------------- Predictor-Corrector -----------------------!
tWrite=0
!do x=0,1
!do tcount=0,1000000
do

	if (tPrint >= 1e4*Yr) then
		write(6,*) " "
 		write(6,*) " Time: "
 		write(6,*) t/yr
		write(6,*) " "
 		write(6,*) " Timestep: "
 		write(6,*) dt 
 		write(6,*) " "
 		write(6,*) " Current Error: "
 		write(6,*) (maxval(RealErr))
 		tPrint = 0
	end if

	!-------- Predicted Values --------!

    do i=1,n
        rp(1:3,i) = r(1:3,i,0) + (dt/24)*( -(9*v(1:3,i,-4)) + (37*v(1:3,i,-3)) - (59*v(1:3,i,-2)) + (55*v(1:3,i,-1)))
        vp(1:3,i) = v(1:3,i,-1) + (dt/24)*( -(9*a(1:3,i,-4)) + (37*a(1:3,i,-3)) - (59*a(1:3,i,-2)) + (55*a(1:3,i,-1)))
    end do
	
	a(1:3,1:n,0)=0.
    do i=1,n !getting acceleration
        do j=1,n
            if (i==j) then
                cycle
            end if 
            D = rp(1:3,i) - rp(1:3,j) !Get distance between two objects		
            AbsD = SQRT(SUM(D**2)) !Get absolute distance
	
            a(1:3,i,0) = a(1:3,i,0) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
				
        end do
    end do
    
    
    !-------- corrector --------!
    
    do i=1,n
    	v(1:3,i,0) = v(1:3,i,-1) + (dt/24)*( a(1:3,i,-3) - (5*a(1:3,i,-2)) + (19*a(1:3,i,-1)) + (9*a(1:3,i,0)) )
    	r(1:3,i,0) = r(1:3,i,0) +  (dt/24)*( v(1:3,i,-3) - (5*v(1:3,i,-2)) + (19*v(1:3,i,-1)) + (9*v(1:3,i,0)) )
    end do
    
    a(1:3,1:n,0)=0.
    do i=1,n !getting acceleration
        do j=1,n
            if (i==j) then
                cycle
            end if 
            D = r(1:3,i,0) - r(1:3,j,0) !Get distance between two objects		
            AbsD = SQRT(SUM(D**2)) !Get absolute distance
	
            a(1:3,i,0) = a(1:3,i,0) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
				
        end do
    end do
    
	do i=1,n
		do j=1,3
			AbsDifference(1,j,i) =  r(j,i,0) - rp(j,i)
			AbsDifference(2,j,i) =  v(j,i,0) - vp(j,i)
		end do
	end do
	
	RealErr = 0.
	do i=1,n
		RealErr(1,1:3,i) = sqrt(AbsDifference(1,1:3,i)**2)/(sqrt(r(1:3,i,0)**2) + small(1))
		RealErr(2,1:3,i) = sqrt(AbsDifference(2,1:3,i)**2)/(sqrt(v(1:3,i,0)**2) + small(1))
	end do
    
    if ( maxval(RealErr(1:2,1:3,1:n)) > err ) then
    
    	dt = dt*0.5
    
    	!write(6,*) "Halving:" , tcount, dt

		do i=1,n
			
			v_place(1:3,i,-1) = (1./128.)*( -5.*v(1:3,i,-4) + 28.*v(1:3,i,-3) - 70.*v(1:3,i,-2) + 140.*v(1:3,i,-1) + 35.*v(1:3,i,0) )
			a_place(1:3,i,-1) = (1./128.)*( -5.*a(1:3,i,-4) + 28.*a(1:3,i,-3) - 70.*a(1:3,i,-2) + 140.*a(1:3,i,-1) + 35.*a(1:3,i,0) )

			v_place(1:3,i,-2) = v(1:3,i,-1)
			a_place(1:3,i,-2) = a(1:3,i,-1)
			
			v_place(1:3,i,-3) = (1./128.)*( 3.*v(1:3,i,-4) - 20.*v(1:3,i,-3) + 90.*v(1:3,i,-2) + 60.*v(1:3,i,-1) - 5.*v(1:3,i,0) )
			a_place(1:3,i,-3) = (1./128.)*( 3.*a(1:3,i,-4) - 20.*a(1:3,i,-3) + 90.*a(1:3,i,-2) + 60.*a(1:3,i,-1) - 5.*a(1:3,i,0) )
			
			v_place(1:3,i,-4) = v(1:3,i,-2)
			a_place(1:3,i,-4) = a(1:3,i,-2)
			
			! if (i == 2) then
! 				write(6,*) v(1,2,-4:0)
! 				write(6,*) a(1,2,-4:0)
! 			end if
			
			v(1:3,i,-1) = v_place(1:3,i,-1)
			a(1:3,i,-1) = a_place(1:3,i,-1)
			
			v(1:3,i,-2) = v_place(1:3,i,-2)
			a(1:3,i,-2) = a_place(1:3,i,-2)
			
			v(1:3,i,-3) = v_place(1:3,i,-3)
			a(1:3,i,-3) = a_place(1:3,i,-3)
			
			v(1:3,i,-4) = v_place(1:3,i,-4)
			a(1:3,i,-4) = a_place(1:3,i,-4)
			
		 end do
		 
		 ! write(6,*) " "
! 	     write(6,*) maxval(RealErr(1:2,1:3,2:n)), dt
! 		 write(6,*) " "

	! write(6,*) " "
! 	write(6,*) v(1,2,-4:0)
! 	write(6,*) a(1,2,-4:0)
! 	
! 	pause
		 
	end if
    
    if (maxval(RealErr(1:2,1:3,1:n)) < err*0.001) then
        
		if (tchange > 8) then
			!double timestep

		dt = dt*2
			
		!write(6,*) "Doubling:" , tcount, dt

		do i=1,n
			
			! if (i == 2) then
! 				write(6,*) v(1,2,-4:0)
! 				write(6,*) a(1,2,-4:0)
! 			end if
	
			v(1:3,i,-1) = v(1:3,i,-2)
			a(1:3,i,-1) = a(1:3,i,-2)
	
			v(1:3,i,-2) = v(1:3,i,-4)
			a(1:3,i,-2) = a(1:3,i,-4)
	
			v(1:3,i,-3) = v(1:3,i,-6)
			a(1:3,i,-3) = a(1:3,i,-6)
	
			v(1:3,i,-4) = v(1:3,i,-8)
			a(1:3,i,-4) = a(1:3,i,-8)
		
		end do
		
		tchange = 4
		
		! write(6,*) " "
! 		write(6,*) v(1,2,-4:0)
! 		write(6,*) a(1,2,-4:0)
! 	
! 		pause
		
		end if
    
    end if
    
    do i=-8,-1,1
        do j=1,n
			r(1:3,j,i) = r(1:3,j,i+1)
			v(1:3,j,i) = v(1:3,j,i+1)
			a(1:3,j,i) = a(1:3,j,i+1)
        end do
    end do
	
    if (tWrite >= 26*604800) then !Write data once every 2 months (week = 604800s in simulation) based on simulation time, since timestep is variable
        ! For Energy Conservation
		! KE
		KE=0.
		do i=1,n
			KE = KE + 0.5*m(i)*sum(v(1:3,i,0)**2)
		end do


		! PE
		PE = 0.
		do i=1,n
			do j=1,n
				if (i==j) then
					cycle
				end if 
				D = r(1:3,i,0) - r(1:3,j,0)	
				AbsD = SQRT(SUM(D**2))
				PE = PE - (G*m(i)*m(j))/(AbsD)
			end do
		end do
		PE = 0.5*PE	
		
		!Total Energy
		E = PE + KE	
			
        write(3,*) (t/Yr), r(1:2,1:n,0)/AU
        write(4,*) (t/Yr), (E-E_i)/E_i,(KE-KE_i)/KE_i,(PE-PE_i)/PE_i
        
        tWrite = 0
    end if

	tcount = tcount + 1
	tChange = tchange + 1 
    t = t + dt
    tWrite = tWrite + dt
    tPrint = tPrint + dt

if (t/yr >= 1000000) then
	exit
end if

end do
!end do
	
	
write(6,*) " "
write(6,*) tcount, dt
write(6,*) " "
write(6,*) r(1:3,2,-8:0)/AU
write(6,*) " "
write(6,*) v(1:3,2,-8:0) 
write(6,*) " "
write(6,*) a(1:3,2,-8:0) 
write(6,*) " "

END PROGRAM PCAttempt2