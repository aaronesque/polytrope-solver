module solver

  ! Uses 

!  use Constants
!  use Parameters
!  use Variables

  ! No implicit typing

  implicit none
  
  integer            :: step          ! step number
  integer            :: data_out  = 2
  real, dimension(2) :: u_old           ! solution vector for satellite (old)
  real, dimension(2) :: u_new           ! solution vector for satellite (new)
  real               :: xi, th, dth
  real               :: dxi             ! time step
  real               :: xi_max           ! max time
  character(len=7)  :: unknown="unknown"
  character(len=3)  :: old="old"
  character(len=6)   :: append="append"
  character(len=*), parameter :: header="step xi theta dtheta"
  integer :: gate = 1 !0

  ! Procedures

contains
   
!***
! ps_init ()
!
! Initializes variables and i/o required for polytrope solution procedures
!***

  subroutine ps_init
  
    implicit none
            
    ! Open output file for data, write header

    open (unit=data_out,file="output.dat",status=unknown)
    write(data_out,*) header

    ! Set initial variables

    u_old   = (/ 1.0,  0.0 /)
    u_new   = u_old 
    
    th      = u_old(1) 
    dth      = u_old(2) 
            
    step   = 1            ! iteration count
    xi      = 0.0            ! total elapsed time
    xi_max   = 10.0
    dxi      = 1E-2

  end subroutine ps_init

!***
! ps_exit ()
!
! Ends program, clears memory
!***

  subroutine ps_exit
         
    implicit none

    close(data_out)
  
  end subroutine ps_exit

!***
! ps_interrupt ()
!
! Outputs calculated info to term and file
!***

  subroutine ps_interrupt 
          
    implicit none
              
    if ( mod(step, 20) .eq. 0 ) then
      write(6,*) step, xi, th, dth
    end if
              
    write(data_out,*) step, xi, th, dth 
           
  end subroutine ps_interrupt

!***
! ps_update ()
!
! New value of last iteration becomes old one for the next. Abscissa step updated.
!***

  subroutine ps_update
    
    implicit none
              
    u_old   = u_new
    step  = step + 1
    xi      = xi + dxi
    !dxi = dxi_scale (xi, dxi, u_old)

  end subroutine ps_update

!***
! rk4 ()
! 
! Takes a step using 4th Order Runge-Kutta integration
!***

  subroutine rk4  ! Takes a step 
    implicit none
    real, dimension(2) :: K1, K2, K3, K4  
    real               :: dxi_1


    K1 = dxi*func(xi, u_old)
    K2 = dxi*func(xi + dxi/2, u_old + K1/2)
    K3 = dxi*func(xi + dxi/2, u_old + K2/2)
    K4 = dxi*func(xi + dxi, u_old + K3)
        
    u_new = u_old + (1./6.)*(K1 + 2*(K2+K3) + K4)

    th = u_new(1)
    dth = u_new(2)

  end subroutine rk4

!***
! ark4 ()
!
! Takes a trial step with rk4() to adapt the step size
!  in order to achieve a desired accuracy
!
!***

  subroutine ark4

    implicit none
    real, dimension(2) :: u_in
    real               :: xi, dxi_test
    real                           :: eps, ftol, frac
    real                           :: u_scale, u_test
    real                           :: dxi_scale
   

   !  
    call rk4
    
    eps = 1E-2
    ftol = 1E-3
    u_scale = eps * 1.
   ! u_test should actually be calling func(), since it's our best way of approximating the error! 
    u_test = abs( u_in(1) ) + abs( dxi_test * u_in(2) )

    frac = u_scale / u_test

    if (frac .gt. 1) then
      ! shrink
      dxi_scale = dxi_test * frac**0.2
    else if (frac .gt. ftol) then
      ! grow
      dxi_scale = dxi_test * frac**0.25
    else
      ! grow more
      dxi_scale = 0.9*eps
    end if
    
    if (dxi_scale .lt. ftol) then
        ! shrunk too much. grow back to ftol
        dxi_scale = ftol
    end if    


  end subroutine ark4


  function dxi_scale (xi, dxi_test, u_in)
    
    implicit none
    real, dimension(2), intent(in) :: u_in
    real, intent(in)               :: xi, dxi_test
    real                           :: eps, ftol, frac
    real                           :: u_scale, u_test
    real                           :: dxi_scale
    
    eps = 1E-2
    ftol = 1E-3
    u_scale = eps * 1.
   ! u_test should actually be calling func(), since it's our best way of approximating the error! 
    u_test = abs( u_in(1) ) + abs( dxi_test * u_in(2) )

    frac = u_scale / u_test

    if (frac .gt. 1) then
      ! shrink
      dxi_scale = dxi_test * frac**0.2
        
    else if (frac .gt. ftol) then
      ! grow
      dxi_scale = dxi_test * frac**0.25

    else
      ! grow more
      dxi_scale = 0.9*eps
    
    end if
      
    if (dxi_scale .lt. ftol) then
        ! shrunk too much. grow back to ftol
        dxi_scale = ftol
    end if    


  end function dxi_scale

!***
! root_check (u_in) 
!
! Checks whether solution function has crossed the horizontal axis
!***

  subroutine root_check (u_in)

    implicit none
    real, dimension(2), intent(in) :: u_in
    real                           :: th

    th = u_in(1)
    
    if (th .le. 0.) then
     
      ! stopping condition
     
      xi = xi_max + 1.
    
    end if

  end subroutine root_check

!***
! func (xi, u_in)
!
! Finds and returns derivatives
!   dy/dxi = dth
!   dz/dxi = f(x, th, dth)
!***   

  function func (xi, u_in) 
    
    implicit none
            ! declaring input and output vars
    real, dimension(2), intent(in)  :: u_in
    real, dimension(2) :: func
    real                            :: y, dy, z, dz
    real                            :: xi

    y = u_in(1)
    z = u_in(2)
            ! local vars
            ! translating
    dy = z

    if (xi .gt. 0.) then
      dz = - y**3. - (2./xi)*z
    else 
      dz = - 1./3.
    end if

            ! computing derivatives
    func(1) = dy
    func(2) = dz
  end function func


end module solver
