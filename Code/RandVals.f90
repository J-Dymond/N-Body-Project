PROGRAM RandomVal
IMPLICIT NONE

!randomising position in orbit
INTEGER :: Parity(0:1),ParityVel(0:1)
DOUBLE PRECISION :: Val,RN,dis(0:9,0:1),abs(0:9),vel(0:9,0:1),absvel(0:9)
INTEGER :: i,j,n

n = 9
dis = 0

do i=0,9
	abs(i) = i+1
	dis(i,0) = abs(i)
	absvel(i) = 10 - i
	vel(i,1) = absvel(i)
end do

write(6,*) dis
 
!Producing random positions and velocities in orbit
do j=0,n
	!Parity
	do i=0,1
		call random_number(val)
		if (val > 0.5) then
			parity(i) = 1
		else 
			parity(i) = -1
		end if
	end do
	write(6,*) parity

	!Value
	call random_number(RN)
	write(6,*) RN
	
	!Set Values for distance
	dis(j,0) = parity(0)*RN*abs(j)
	dis(j,1) = parity(1)*(SQRT(abs(j)**2 - dis(j,0)**2))
	
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
	
	vel(j,1) = parityVel(1)*RN*absvel(j)
	vel(j,0) = parityvel(0)*(SQRT(absvel(j)**2 - vel(j,1)**2))
	
end do



write(6,*) dis
write(6,*)
write(6,*) vel
do i = 0,n
	write(6,*) sqrt(dis(i,0)**2 + dis(i,1)**2)
end do
do i = 0,n
	write(6,*) sqrt(vel(i,0)**2 + vel(i,1)**2)
end do
END PROGRAM