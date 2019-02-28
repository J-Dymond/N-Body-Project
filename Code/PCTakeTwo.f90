PROGRAM PCAttempt2
IMPLICIT NONE
!-------------------------Declaring Variables----------------------------!

!Integers for loops and recording time
INTEGER :: i,j,n,x,y,z,tcount,tcountWrite

DOUBLE PRECISION :: dt,t

!One Dimension Variables/Constants
DOUBLE PRECISION :: G,AU,AbsD,Yr,KE_i,PE_i,PE,KE,E_i,E 

!Arrays for kinematic calculations
DOUBLE PRECISION :: COM(1:3),COV(1:3),ai(1:3,1:9),da(1:3,1:9),m(1:9),D(1:3)

!Position, velocity, and acceleration vectors, with time history
DOUBLE PRECISION :: r(1:3,1:9,-8:1)
DOUBLE PRECISION :: v(1:3,1:9,-8:1)
DOUBLE PRECISION :: a(1:3,1:9,-8:1)

!Variables for Randomising Initial Positions and Velocities
DOUBLE PRECISION :: Abs(1:9),AbsVel(1:9),Val,RN
INTEGER :: Parity(1:2),ParityVel(1:2)

!-------------------------------Set Up------------------------------------!

!Open data files
OPEN(3,file="../Data/DataStable.txt",STATUS = 'REPLACE') !This is for the position vectors against time

!Initialise Arrays, Constants etc..
!Constants
AU = 1.496e11 !AU in metres
G = 6.67e-11  !G in SI
Yr = 3.154e7  !YR in seconds

n=7  !number of bodies in simulation

!1000 seconds timestep
t = 0.       !global time value
dt = 100.    !10 second timesteps
tcount = 0   !Timestep counter

!Initialising Body parameters Sun to Neptune
!Masses in SI

m(1) = 1.99e30  !Sun
m(2) = 5.97e24  !Earth
m(3) = 0.642e24 !Mars
m(4) = 1.898e27 !etc..
m(5) = 568e24
m(6) = 86.8e24
m(7) = 102e24

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

!Absolute Values of Orbital Velocities
absVel(2) = 29.78e3 !Earth
absVel(3) = 24.1e3  !etc...
absVel(4) = 13.07e3
absVel(5) = 9.536e3
absVel(6) = 6.687e3
absVel(7) = 5.372e3

!Producing random positions and velocities in orbit -- Will reintroduce later

!Orbital radii in AU 
r(1,2,-8) = abs(2)*AU
r(1,3,-8) = abs(3)*AU 
r(1,4,-8) = abs(4)*AU
r(1,5,-8) = abs(5)*AU
r(1,6,-8) = abs(6)*AU
r(1,7,-8) = abs(7)*AU


!Orbital velocities - moving parallel to y axis initially
v(2,2,-8) = absVel(2)
v(2,3,-8) = absVel(3) 
v(2,4,-8) = absVel(4)
v(2,5,-8) = absVel(5)
v(2,6,-8) = absVel(6)
v(2,7,-8) = absVel(7)


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
    v(1:3,i,0) = v(1:3,i,-8) - COV(1:3)
end do

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

    write(6,*) tcount

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
write(6,*) r(1:3,7,-8:1)/AU
write(6,*) " " 
write(6,*) v(1:3,7,-8:1) 
write(6,*) " "
write(6,*) a(1:3,7,-8:1) 
write(6,*) " "

!----------------------Predictor-----------------------!
tcountWrite=0
do tcount=0,10000000

    do i=1,n
        r(1:3,i,1) = r(1:3,i,0) + (dt/24)*( -(9*v(1:3,i,-3)) + (37*v(1:3,i,-2)) - (59*v(1:3,i,-1)) + (55*v(1:3,i,0)))
        v(1:3,i,1) = v(1:3,i,0) + (dt/24)*( -(9*a(1:3,i,-3)) + (37*a(1:3,i,-2)) - (59*a(1:3,i,-1)) + (55*a(1:3,i,0)))
    end do
	
	a(1:3,1:n,1)=0.
    do i=1,n !getting acceleration
        do j=1,n
            if (i==j) then
                cycle
            end if 
            D = r(1:3,i,1) - r(1:3,j,1) !Get distance between two objects		
            AbsD = SQRT(SUM(D**2)) !Get absolute distance
	
            a(1:3,i,1) = a(1:3,i,1) - ((G*m(j))/((AbsD)**3))*(D) !Make calculation
				
        end do
    end do
	
	
    if (a(1,1,1) == a(1,1,0)) then
        write (6,*) tcount
        exit
    end if
	
    if (tcountWrite == 10000) then
        write(3,*) (t/Yr), r(1:2,1:n,0)/AU
        tcountWrite = 0
    end if
	
    do i=-8,0,1
        do j=1,n
            r(1:3,j,i) = r(1:3,j,i+1)
            v(1:3,j,i) = v(1:3,j,i+1)
            a(1:3,j,i) = a(1:3,j,i+1)
        end do
    end do

    t = t + dt
    tcountWrite = tcountWrite + 1
end do
	
	
write(6,*) " "
write(6,*) tcount
write(6,*) " "
write(6,*) r(1:3,7,-8:1)/AU
write(6,*) " "
write(6,*) v(1:3,7,-8:1) 
write(6,*) " "
write(6,*) a(1:3,7,-8:1) 
write(6,*) " "

END PROGRAM PCAttempt2