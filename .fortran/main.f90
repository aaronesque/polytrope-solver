program main
  use solver 
  implicit none
  integer :: status = 0
    call ps_init
!    call term
    do while (xi .le. xi_max)
      call ps_interrupt
      call ark4
      call root_check (u_new)

!      if (sho_test == 2) then ! Ends progam if solver reaches root 
!        print*, "Root reached at xi = ", u_xi
!        call exit(status)
!      end if
      call ps_update
    end do
    call ps_exit
end program main


