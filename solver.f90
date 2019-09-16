module solver

implicit none
real               :: xi !dimensionless radius
real               :: dxi !xi step size
integer            :: step !step number
real, dimension(2) :: uOld !Old solution vector
real, dimension(2) :: uNew !New solution vector
real               :: theta !density parameter
real               :: dth !change in theta
real               :: xiMax = 15 !Used to stop Program

contains

subroutine polyInit
	
	implicit none
	
	!Opens up data file to record data
	open(unit=2, file="output.dat", status="unknown")
	write(2, *) "step xi theta dtheta"
	
	!Initializes solution arrays
	uOld = (/ 1.0, 0.0 /)
	uNew = uOld

	!initializes variables
	theta = 1.0
	dth = 0.0
	step = 1
	xi = 0.0
	dxi = 1E-2

end subroutine polyInit

subroutine interrupt

	implicit none

	!Writes out data at each step
	write(2, *) step, xi, theta, dth

end subroutine interrupt

subroutine rk4

	implicit none
	real, dimension(2) :: K1, K2, K3, K4
	real, dimension(2) :: m

	!Calculate 'slopes'
	K1 = diffeq(xi, uOld)
	K2 = diffeq(xi + (dxi/2), uOld + (K1*dxi/2))
	K3 = diffeq(xi + (dxi/2), uOld + (K2*dxi/2))
	K4 = diffeq(xi + dxi, uOld + (K3*dxi))

	!Update solution array
	m = (K1/6) + (K2/3) + (K3/3) + (K4/6)
	uNew = uOld + (m*dxi)

	!update variables
	theta = uNew(1)
	dth = uNew(2)
	xi = xi + dxi
	uOld = uNew
	step = step + 1

end subroutine rk4

subroutine rootCheck(u)

	implicit none
	real, dimension(2), intent(in) :: u
	real :: theta

	theta = u(1)

	!Determine when theta turns negative
	if (theta .le. 0.) then
		
		print*, "Root reached at xi=", xi
		xi = xiMax + 1

	end if

end subroutine rootCheck

!Prevent data leak
subroutine polyEnd

	implicit none

	close(2)

end subroutine polyEnd

function diffeq (xi, u)
!d2th_dxi2 = -(2/xi)dth_dxi - theta^n

	implicit none
	real, dimension(2), intent(in) :: u
	real, dimension(2) :: diffeq
	real :: xi
	real :: n

	n = 0

	diffeq(1) = u(2)
	if (xi .gt. 0) then
		diffeq(2) = -((2./xi)*u(2)) - u(1)**n
	else
		diffeq(2) = -1./3.
	end if

end function diffeq
	 
end module solver
