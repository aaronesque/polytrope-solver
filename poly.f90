program poly

	use solver
	implicit none
	call polyInit
	do while (xi .le. xiMax)

		call interrupt !Record data
		call rk4 !Take next step
		call rootCheck(uNew) !check if zero has beeen crossed

	end do
	call polyEnd

end program poly
