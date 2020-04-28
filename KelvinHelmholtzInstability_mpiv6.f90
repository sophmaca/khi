!  iompressible Fluid Dynamics - Finite volume
!  Sophia Macarewich, EECS 587, Fall 2018
!  Kelvin Helmholtz Instability (KHI)
!  compile with gnu (one line):
!  gfortran -o KelvinHelmholtzInstability_mat KelvinHelmholtzInstability_mat.f90
!  mesh.o utils.o types.o
!  compile with: gcc/5.4.0 
!  This is a column major version (decompose into cols)

  program main

  use mpi
  use types, only: dp
  use utils, only: savetxt
  use mesh, only: linspace, meshgrid

  implicit none

! Parameters
  integer, parameter :: Nx = 128, Ny = 128
  real (dp), parameter :: boxSizeX = 1., boxSizeY = 1.
! dx = boxSizeX / Nx ; dy = boxSizeY / Ny ; vol = dx*dy
  real (dp) :: dx, dy, vol
  real (dp), dimension (0:Nx+1,0:Ny+1) :: Y, X
  real (dp), parameter :: pi = 3.1415927
  real (dp) :: t, tEnd, tOut, dt, dt_loc
  logical :: plotThisTurn

  real (dp), parameter :: w0 = 0.1, sigma = 0.05/SQRT(2.), gamma = 5/3.
  real (dp), allocatable, dimension (:,:) :: rho, vx, vy, P

  real (dp), allocatable, dimension (:,:) :: Mass, Momx, Momy, Energy
  integer :: outputCount

  real (dp), allocatable, dimension (:,:) :: rho_gradx, rho_grady, vx_gradx, vx_grady
  real (dp), allocatable, dimension (:,:) :: vy_gradx, vy_grady, P_gradx, P_grady

  real (dp), allocatable, dimension (:,:) :: rho_prime, rho_XL, rho_XR, rho_YL, rho_YR
  real (dp), allocatable, dimension (:,:) :: vx_prime, vx_XL, vx_XR, vx_YL, vx_YR
  real (dp), allocatable, dimension (:,:) :: vy_prime, vy_XL, vy_XR, vy_YL, vy_YR
  real (dp), allocatable, dimension (:,:) :: P_prime, P_XL, P_XR, P_YL, P_YR

  real (dp), allocatable, dimension (:,:) :: rho_Xstar, rho_Ystar, momx_Xstar, momx_Ystar
  real (dp), allocatable, dimension (:,:) :: momy_Xstar, momy_Ystar, en_Xstar, en_Ystar
  real (dp), allocatable, dimension (:,:) :: P_Xstar, P_Ystar

  real (dp), allocatable, dimension (:,:) :: flux_rho_X, flux_rho_Y, flux_momx_X,flux_momx_Y
  real (dp), allocatable, dimension (:,:) :: flux_momy_X, flux_momy_Y, flux_en_X, flux_en_Y,C

  real (dp) :: rho_plot(100,1:Nx,1:Ny)
  integer :: i,j,k
  character(len=70) :: fn

! MPI Parameters

  integer, parameter :: NP = 8                 ! need to change this for # proc
  integer, parameter :: MASTER = 0, FROM_MASTER = 1, FROM_WORKER = 2
  integer :: rank, tag, ierr, source, dest, row, numtasks
  real (dp) :: sizeStrip, finish, start
  integer :: status(MPI_STATUS_SIZE), request

! Cart Parameters

  integer, parameter :: R = -1, L = 1           ! directions for CSHIFT
  integer, parameter :: ndim = 1
  integer :: grid_comm, dims(ndim), coord, s, e
  logical :: periodic(ndim), reorder
  integer :: grid_rank, below, above, grid_id


  call MPI_INIT(ierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, ierr )
!  print *, 'total num tasks= ',numtasks
  print *, 'proc ID= ',rank

  call MPI_Barrier(  MPI_COMM_WORLD, ierr )
  start = MPI_Wtime()

! *********************** EVERY PROC  ***************************
! Redefine processors in a Cartesian topology
    dims(1) = numtasks
    periodic(1) = .TRUE.
    reorder = .TRUE.
    call MPI_Cart_create ( MPI_COMM_WORLD, ndim, dims, periodic, reorder,&
       grid_comm, ierr )
    call MPI_Comm_rank(grid_comm, grid_id, ierr)
    call MPI_Cart_shift(grid_comm, 0, L, below, above, ierr) ! shift row up
!
!  Compute the rows "S" through "E" of data that are to be assigned to
!  this
!  process.
!
  call decomp_band ( Ny, numtasks, grid_id, s, e )
  write (*,*) 'task',grid_id,'has rows',s,'to',e

!
!  Once we have the sizes S and E, we can allocate all variables
!
  allocate ( rho(0:Ny+1,s-1:e+1) )
  allocate ( vx(0:Ny+1,s-1:e+1) )
  allocate ( vy(0:Ny+1,s-1:e+1) )
  allocate ( P(0:Ny+1,s-1:e+1) )
  allocate ( Mass(0:Ny+1,s-1:e+1) )
  allocate ( Momx(0:Ny+1,s-1:e+1) )
  allocate ( Momy(0:Ny+1,s-1:e+1) )
  allocate ( Energy(0:Ny+1,s-1:e+1) )

  allocate ( rho_gradx(0:Ny+1,s-1:e+1) )
  allocate ( rho_grady(0:Ny+1,s-1:e+1) )
  allocate ( vx_gradx(0:Ny+1,s-1:e+1) )
  allocate ( vx_grady(0:Ny+1,s-1:e+1) )
  allocate ( vy_gradx(0:Ny+1,s-1:e+1) )
  allocate ( vy_grady(0:Ny+1,s-1:e+1) )
  allocate ( P_gradx(0:Ny+1,s-1:e+1) )
  allocate ( P_grady(0:Ny+1,s-1:e+1) )

  allocate ( rho_prime(0:Ny+1,s-1:e+1) )
  allocate ( rho_XL(0:Ny+1,s-1:e+1) )
  allocate ( rho_XR(0:Ny+1,s-1:e+1) )
  allocate ( rho_YL(0:Ny+1,s-1:e+1) )
  allocate ( rho_YR(0:Ny+1,s-1:e+1) )
  allocate ( vx_prime(0:Ny+1,s-1:e+1) )
  allocate ( vx_XL(0:Ny+1,s-1:e+1) )
  allocate ( vx_XR(0:Ny+1,s-1:e+1) )
  allocate ( vx_YL(0:Ny+1,s-1:e+1) )
  allocate ( vx_YR(0:Ny+1,s-1:e+1) )
  allocate ( vy_prime(0:Ny+1,s-1:e+1) )
  allocate ( vy_XL(0:Ny+1,s-1:e+1) )
  allocate ( vy_XR(0:Ny+1,s-1:e+1) )
  allocate ( vy_YL(0:Ny+1,s-1:e+1) )
  allocate ( vy_YR(0:Ny+1,s-1:e+1) )
  allocate ( P_prime(0:Ny+1,s-1:e+1) )
  allocate ( P_XL(0:Ny+1,s-1:e+1) )
  allocate ( P_XR(0:Ny+1,s-1:e+1) )
  allocate ( P_YL(0:Ny+1,s-1:e+1) )
  allocate ( P_YR(0:Ny+1,s-1:e+1) )

  allocate ( rho_Xstar(0:Ny+1,s-1:e+1) )
  allocate ( rho_Ystar(0:Ny+1,s-1:e+1) )
  allocate ( momx_Xstar(0:Ny+1,s-1:e+1) )
  allocate ( momx_Ystar(0:Ny+1,s-1:e+1) )
  allocate ( momy_Xstar(0:Ny+1,s-1:e+1) )
  allocate ( momy_Ystar(0:Ny+1,s-1:e+1) )
  allocate ( en_Xstar(0:Ny+1,s-1:e+1) )
  allocate ( en_Ystar(0:Ny+1,s-1:e+1) )
  allocate ( P_Xstar(0:Ny+1,s-1:e+1) )
  allocate ( P_Ystar(0:Ny+1,s-1:e+1) )

  allocate ( flux_rho_X(0:Ny+1,s-1:e+1) )
  allocate ( flux_rho_Y(0:Ny+1,s-1:e+1) )
  allocate ( flux_momx_X(0:Ny+1,s-1:e+1) )
  allocate ( flux_momx_Y(0:Ny+1,s-1:e+1) )
  allocate ( flux_momy_X(0:Ny+1,s-1:e+1) )
  allocate ( flux_momy_Y(0:Ny+1,s-1:e+1) )
  allocate ( flux_en_X(0:Ny+1,s-1:e+1) )
  allocate ( flux_en_Y(0:Ny+1,s-1:e+1) )
  allocate ( C(0:Ny+1,s-1:e+1) )

  t = 0.0D0
  tEnd = 2.0D0
  tOut = 0.01D0
  outputCount = 1
  if (rank .eq. 0) write(*,*) 'tOut=',tOut

! Calc global grid and then select local Y points for task
  dx = 1.0E+00 / real ( Nx + 1 )
  dy = 1.0E+00 / real ( Ny + 1 )
  call meshgrid( linspace(-1.5*dy,boxSizeY-0.5*dy, Ny+2),linspace(0.5*dx,boxSizeX-0.5*dx, Nx+2), Y, X)

! Initalize values in each strip
  call init_band ( Nx, Ny, X, Y, s, e, grid_id, numtasks, &
    dx, dy, rho, vx, vy, P, Mass, Momx, Momy, Energy)

!!  if (rank .eq. 1) write(*,*) 'Rank 1 rho val is',rho

  call MPI_Barrier( grid_comm, ierr )

! Main loop
  do while ( t < tEnd ) 

! Get primitive variables
    call primitiveVar( s, e, Ny, dx, dy, Mass, Momx, Momy, Energy, rho, vx, vy, P)
    if (rank .eq. 1) write(*,*) 'Rank 1 mass min val =',MINVAL(Mass(1:Ny,s:e))

! Get time step of each worker
    if (rank .eq. 1) write(*,*) 'Rank 1 rho min val =',MINVAL(rho(1:Ny,s:e))
    call calc_CFL(s, e, Ny, dx, dy, rho, vx, vy, P, dt_loc)
    write(*,*) 'local dt of rank',grid_id,'is',dt_loc

! MASTER reduces all local dt into global min dt
    call MPI_Allreduce( dt_loc, dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      grid_comm,ierr )
!!    if (rank .eq. 0) write(*,*) 'MASTER calc dt=',dt,'for t=',t

    plotThisTurn = .FALSE.
    ! every proc decides whether or not to plot
    if ( t + dt > outputCount*tOut ) then
      dt = outputCount*tOut - t
      plotThisTurn = .TRUE.
    end if
    if (rank .eq. 0) write(*,*) 'outputCount=',outputCount
    if (rank .eq. 0) write(*,*) 'tOut=',tOut
    if (rank .eq. 0) write(*,*) 'MASTER calc dt=',dt,'for t=',t

! Calculate gradients 
    call calcGradients( dx, dy, rho, vx, vy, P, below, above, &
        s, e, Ny, grid_comm, rho_gradx, rho_grady, vx_gradx, vx_grady, &
        vy_gradx, vy_grady, P_gradx, P_grady )

! Extrapolate to cell faces (in time & space)
    call extrapSpaceTime( dx, dy, rho, vx, vy, P, &
        below, above, dt, s, e, Ny, grid_comm, &
        rho_gradx, rho_grady, vx_gradx, vx_grady, &
        vy_gradx, vy_grady, P_gradx, P_grady, &
        rho_XR, rho_XL, rho_YR, rho_YL, &
        vx_XR, vx_XL, vx_XR, vx_XL, &
        vy_XR, rho_XL, rho_XR, rho_XL, &
        rho_XR, rho_XL, rho_XR, rho_XL )

! Compute star (averaged) states
    rho_Xstar = 0.5*(rho_XL + rho_XR)
    rho_Ystar = 0.5*(rho_YL + rho_YR)
    momx_Xstar = 0.5*(rho_XL * vx_XL + rho_XR * vx_XR)
    momx_Ystar = 0.5*(rho_YL * vx_YL + rho_YR * vx_YR)
    momy_Xstar = 0.5*(rho_XL * vy_XL + rho_XR * vy_XR)
    momy_Ystar = 0.5*(rho_YL * vy_YL + rho_YR * vy_YR)
    en_Xstar = 0.5*( P_XL/(gamma-1)+0.5*rho_XL * (vx_XL**2+vy_XL**2) + P_XR/(gamma-1)+0.5*rho_XR * (vx_XR**2+vy_XR**2))
    en_Ystar = 0.5*( P_YL/(gamma-1)+0.5*rho_YL * (vx_YL**2+vy_YL**2) + P_YR/(gamma-1)+0.5*rho_YR * (vx_YR**2+vy_YR**2))

    P_Xstar = (gamma-1)*(en_Xstar-0.5*(momx_Xstar**2+momy_Xstar**2)/rho_Xstar)
    P_Ystar = (gamma-1)*(en_Ystar-0.5*(momx_Ystar**2+momy_Ystar**2)/rho_Ystar)

! Compute fluxes (local Lax-Friedrichs/Rusanov)
    flux_rho_X = momx_Xstar
    flux_rho_Y = momy_Ystar
    flux_momx_X = momx_Xstar**2/rho_Xstar + P_Xstar
    flux_momx_Y = momy_Ystar * momx_Ystar/rho_Ystar
    flux_momy_X = momx_Xstar * momy_Xstar/rho_Xstar
    flux_momy_Y = momy_Ystar**2/rho_Ystar + P_Ystar
    flux_en_X = (en_Xstar+P_Xstar) * momx_Xstar/rho_Xstar
    flux_en_Y = (en_Ystar+P_Ystar) * momy_Ystar/rho_Ystar

    C = SQRT(gamma*P_XL/rho_XL) + ABS(vx_XL)
    C = MAX( C, SQRT(gamma*P_XR/rho_XR) + ABS(vx_XR))
    C = MAX( C, SQRT(gamma*P_YL/rho_YL) + ABS(vy_YL))
    C = MAX( C, SQRT(gamma*P_YR/rho_YR) + ABS(vy_YR))

    flux_rho_X = flux_rho_X - C * 0.5 * (rho_XL - rho_XR)
    flux_rho_Y = flux_rho_Y - C * 0.5 * (rho_YL - rho_YR)
    flux_momx_X = flux_momx_X - C * 0.5 * (rho_XL * vx_XL - rho_XR * vx_XR)
    flux_momx_Y = flux_momx_Y - C * 0.5 * (rho_YL * vx_YL - rho_YR * vx_YR)
    flux_momy_X = flux_momy_X - C * 0.5 * (rho_XL * vy_XL - rho_XR * vy_XR)
    flux_momy_Y = flux_momy_Y - C * 0.5 * (rho_YL * vy_YL - rho_YR * vy_YR)
    flux_en_X = flux_en_X - C * 0.5 * ( P_XL/(gamma-1)+0.5*rho_XL * (vx_XL**2+vy_XL**2) - &
              (P_XR/(gamma-1)+0.5*rho_XR * (vx_XR**2+vy_XR**2)))
    flux_en_Y = flux_en_Y - C * 0.5 * ( P_YL/(gamma-1)+0.5*rho_YL * (vx_YL**2+vy_YL**2) - & 
              (P_YR/(gamma-1)+0.5*rho_YR * (vx_YR**2+vy_YR**2)))

! Update solution
    call update_soln( Ny, dx, dy, dt, s, e, below, above, grid_comm, &
      flux_rho_X, flux_rho_Y, &
      flux_momx_X, flux_momx_Y, &
      flux_momy_X, flux_momy_Y, &
      flux_en_X, flux_en_Y, &
      Mass, Momx, Momy, Energy )

  if ( plotThisTurn ) then
!! DO NOT DO UNTIL ERROR FIXED***
!!  if ( rank .eq. 0 ) write(*,*) 'gathering all rho rows'
!!    call MPI_Gather( rho(s:e,1:Nx),Nx*(e-s),MPI_DOUBLE_PRECISION,rho_plot,&
!!              Nx*(e-s),MPI_DOUBLE_PRECISION,MASTER,grid_comm,ierr )
    outputCount = outputCount + 1
  if ( rank .eq. 0 ) write(*,*) 't=',t,'complete'
  end if

! Advance time
  t = t + dt

  call MPI_Barrier(  grid_comm, ierr ) 
  end do

  call MPI_Barrier(  grid_comm, ierr )

  if ( rank .eq. 0 ) print *, 'done'

!! DO NOT DO UNTIL ERROR FIXED***
!!  if (rank .eq. 0) then
!!  print *, "MASTER saving rho..."
!!  j = 1
!!  do i=1,100
!!    write(fn,fmt='(i0,a)') j, '.dat'
!!    open (unit=44,file=fn,form='formatted')
!!    write(44, *) rho_plot(i,1:Nx,1:Ny)
!!    close(44)
!!    j = j + 1
!!  end do
!!  end if

! PUT TIMING HERE
  finish = MPI_Wtime()
  if ( rank .eq. 0 ) print *, 'time elapsed=',finish-start

  call MPI_Finalize(ierr)

  end program main

  subroutine decomp_band ( m, n, r, s, e )
  implicit none
  integer :: e, m, n, part, r, s, whole

  whole = m / n
  part = m - whole * n

  if ( r + 1 <= part ) then
    s =    r             * ( whole + 1 ) + 1
    e =  ( r + 1 )       * ( whole + 1 )
  else
    s =   r             * whole + 1 + part
    e = ( r + 1 )       * whole     + part
  end if

  return
  end subroutine decomp_band

  subroutine init_band ( Nx, Ny, X1, Y1, s, e, grid_id, numtasks, &
    dx, dy, rho, vx, vy, P, Mass, Momx, Momy, Energy)
  use mpi
  use types, only: dp

  real (dp), parameter :: pi = 3.1415927
  real (dp), parameter :: w0 = 0.1, sigma = 0.05/SQRT(2.), gamma = 5/3.
  integer :: Nx, Ny, s, e
  real (dp), dimension (0:Ny+1,s-1:e+1) :: x, y
  real (dp), dimension (0:Nx+1,0:Ny+1) :: X1, Y1
  integer :: grid_id, ierr, numtasks
  integer, parameter :: MASTER = 0
  real (dp) :: dx, dy, vol
  integer :: i, j
  real (dp), parameter :: boxSizeX = 1., boxSizeY = 1.
  real (dp), dimension (0:Ny+1,s-1:e+1) :: rho, vx, vy, P
  real (dp), dimension (0:Ny+1,s-1:e+1) :: Mass, Momx, Momy, Energy

  vol = dx*dy

  x = X1(0:Ny+1,s-1:e+1)
  y = Y1(0:Ny+1,s-1:e+1)

  if ( grid_id == MASTER ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INIT_BAND - Master Cartesian process:'
    write ( *, '(a,g14.6)' ) '  The X grid spacing is ', dx
    write ( *, '(a,g14.6)' ) '  The Y grid spacing is ', dy
  end if

! Set initial conditions for KHI
  where ( ABS(y-0.5) < 0.25 )
    rho = 1. + ABS(y-0.5)
    vx = -0.5 + ABS(y-0.5)
  elsewhere
    rho = 1.
    vx = -0.5
  end where

  vy = w0*SIN(4*pi*x) * ( EXP(-(y-0.25)**2/(2 * sigma**2)) + EXP(-(y-0.75)**2/(2*sigma**2)) )
  P = 0*x + 2.5

! Get conserved variables
  Mass = rho * vol
  Momx = rho * vx * vol
  Momy = rho * vy * vol
  Energy = (P/(gamma-1) + 0.5*rho*(vx**2+vy**2))*vol

  return
  end subroutine init_band

  subroutine primitiveVar( s, e, Ny, dx, dy, Mass, Momx, Momy, Energy, rho, vx, vy, P) 
  use types, only: dp
  real (dp) :: dx, dy, vol
  integer :: s, e, Ny
  real (dp), parameter :: gamma = 5/3.
  real (dp), dimension (0:Ny+1,s-1:e+1) :: Mass, Momx, Momy, Energy
  real (dp), dimension (0:Ny+1,s-1:e+1) :: rho, vx, vy, P

  vol = dx*dy

! Get primitive variables
  rho = Mass / vol
  vx = Momx / rho / vol
  vy = Momy / rho / vol
  P = (Energy/vol - 0.5*rho * (vx**2+vy**2)) * (gamma-1)

  return
  end subroutine primitiveVar

  subroutine calc_CFL(s, e, Ny, dx, dy, rho, vx, vy, P, dt_loc)         ! get time step (CFL)
  use types, only: dp
  implicit none
  real (dp) :: dx, dy
  integer :: s, e, Ny
  real (dp), parameter :: courant_fac = 0.4
  real (dp), parameter :: gamma = 5/3.
  real (dp), dimension (0:Ny+1,s-1:e+1) :: rho, vx, vy, P
  real (dp) :: dt_loc

  dt_loc = courant_fac * MINVAL( MIN(dx,dy) / &
        SQRT( gamma*P(1:Ny,s:e)/rho(1:Ny,s:e) ) + SQRT(vx(1:Ny,s:e)**2+vy(1:Ny,s:e)**2) )

!!  ! take care of -Inf issue on proc 1, might lead to issues later ********
!!  if ( dt_loc < 0.01D-2 ) dt_loc = 0.01D-2

  return
  end subroutine calc_CFL

  subroutine calcGradients( dx, dy, rho, vx, vy, P, &
        below, above, s, e, Ny, grid_comm, &
        rho_gradx, rho_grady, &
        vx_gradx, vx_grady, &
        vy_gradx, vy_grady, &
        P_gradx, P_grady )
  use types, only: dp    
  implicit none
  real (dp) :: dx, dy
  integer, parameter :: R = -1, L = 1           ! directions for CSHIFT
  integer :: s, e, Ny
  real (dp), dimension (0:Ny+1,s-1:e+1) :: rho, rho_gradx, rho_grady, rhoUp, rhoDown
  real (dp), dimension (0:Ny+1,s-1:e+1) :: vx , vx_gradx, vx_grady, vxUp, vxDown
  real (dp), dimension (0:Ny+1,s-1:e+1) :: vy, vy_gradx, vy_grady, vyUp, vyDown
  real (dp), dimension (0:Ny+1,s-1:e+1) :: P, P_gradx, P_grady, PUp, PDown
  integer :: below, above, ierr, grid_comm
 
  ! shift band up (up_down = 0) or down (up_down = 1)
  rhoUp = rho
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, rhoUp )

  rhoDown = rho 
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, rhoDown )

  vxUp = vx
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, vxUp )

  vxDown = vx
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, vxDown )

  vyUp = vy
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, vyUp )

  vyDown = vy
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, vyDown )

  PUp = P
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, PUp )

  PDown = P
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, PDown )

  ! split up into rows, so CSHIFT works for dim = 2 (left/right)
  ! calculate gradients
  ! all CSHIFT are along dim = 2 or cols
  rho_gradx = ( rhoDown - rhoUp ) / (2.*dx)
  rho_grady = ( CSHIFT(rho,R) - CSHIFT(rho,L) ) / (2.*dy)
  vx_gradx  = ( vxDown - vxUp ) / (2.*dx)
  vx_grady  = ( CSHIFT(vx,R) - CSHIFT(vx,L) ) / (2.*dy)
  vy_gradx  = ( vyDown - vyUp ) / (2.*dx)
  vy_grady  = ( CSHIFT(vy,R) - CSHIFT(vy,L) ) / (2.*dy)
  P_gradx   = ( PDown - PUp ) / (2.*dx)
  P_grady   = ( CSHIFT(P,R) - CSHIFT(P,L) ) / (2.*dy)

  return
  end subroutine calcGradients

  subroutine extrapSpaceTime( dx, dy, rho, vx, vy, P, &
        below, above, dt, s, e, Ny, grid_comm, &
        rho_gradx, rho_grady, vx_gradx, vx_grady, &
        vy_gradx, vy_grady, P_gradx, P_grady, &
        rho_XR, rho_XL, rho_YR, rho_YL, &
        vx_XR, vx_XL, vx_YR, vx_YL, &
        vy_XR, vy_XL, vy_YR, vy_YL, &
        P_XR, P_XL, P_YR, P_YL )
  use mpi
  use types, only: dp
  implicit none
  integer, parameter :: R = -1, L = 1           ! directions for CSHIFT
  real (dp), parameter :: gamma = 5/3.
  real (dp) :: dt, dx, dy 
  integer :: s, e, Ny
  real (dp), dimension (0:Ny+1,s-1:e+1) :: rho, rho_gradx, rho_grady, rho_prime
  real (dp), dimension (0:Ny+1,s-1:e+1) :: vx, vx_gradx, vx_grady, vx_prime
  real (dp), dimension (0:Ny+1,s-1:e+1) :: vy, vy_gradx, vy_grady, vy_prime 
  real (dp), dimension (0:Ny+1,s-1:e+1) :: P, P_gradx, P_grady, P_prime
  real (dp), dimension (0:Ny+1,s-1:e+1) :: rho_XL, rho_XR, rho_YL, rho_YR
  real (dp), dimension (0:Ny+1,s-1:e+1) :: vx_XL, vx_XR, vx_YL, vx_YR
  real (dp), dimension (0:Ny+1,s-1:e+1) :: vy_XL, vy_XR, vy_YL, vy_YR
  real (dp), dimension (0:Ny+1,s-1:e+1) :: P_XL, P_XR, P_YL, P_YR
  integer :: below, above, ierr, grid_comm

  ! extrapolate to cell faces (in time & space)
  rho_prime = rho - 0.5*dt *( vx * rho_gradx + rho * vx_gradx + vy * rho_grady + rho * vy_grady)
  rho_XL = rho_prime - rho_gradx * dx/2.
  ! shift down row
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, rho_XL )
  rho_XR = rho_prime + rho_gradx * dx/2.
  rho_YL = rho_prime - rho_grady * dy/2.
  rho_YL = CSHIFT(rho_YL,R)
  rho_YR = rho_prime + rho_grady * dy/2.

  vx_prime = vx - 0.5*dt * ( vx * vx_gradx + vy * vx_grady + (1/rho) * P_gradx)
  vx_XL = vx_prime - vx_gradx * dx/2.
  ! shift down row
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, vx_XL )
  vx_XR = vx_prime + vx_gradx * dx/2.
  vx_YL = vx_prime - vx_grady * dy/2.
  vx_YL = CSHIFT(vx_YL,R)
  vx_YR = vx_prime + vx_grady * dy/2.

  vy_prime = vy - 0.5*dt * ( vx * vy_gradx + vy * vy_grady + (1/rho) * P_grady)
  vy_XL = vy_prime - vy_gradx * dx/2.
  ! shift down row
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, vy_XL )
  vy_XR = vy_prime + vy_gradx * dx/2.
  vy_YL = vy_prime - vy_grady * dy/2.
  vy_YL = CSHIFT(vy_YL,R)
  vy_YR = vy_prime + vy_grady * dy/2.

  P_prime = P - 0.5*dt * ( gamma*P * (vx_gradx + vy_grady)  + vx * P_gradx + vy * P_grady )
  P_XL = P_prime - P_gradx * dx/2.
  ! shift down row
  call shift_band ( Ny, s, e, 1, below, above, grid_comm, P_XL )
  P_XR = P_prime + P_gradx * dx/2.
  P_YL = P_prime - P_grady * dy/2.
  P_YL = CSHIFT(P_YL,R)
  P_YR = P_prime + P_grady * dy/2.

  return
  end subroutine 

  subroutine update_soln( Ny, dx, dy, dt, s, e, below, above, grid_comm, &
        flux_rho_X, flux_rho_Y, &
        flux_momx_X, flux_momx_Y, &
        flux_momy_X, flux_momy_Y, &
        flux_en_X, flux_en_Y, &
        Mass, Momx, Momy, Energy )
  use mpi
  use types, only: dp
  implicit none
  real (dp) :: dx, dy, dt
  integer :: s, e, Ny
  integer, parameter :: R = -1, L = 1           ! directions for CSHIFT
  real (dp), dimension (0:Ny+1,s-1:e+1) :: Mass, Momx, Momy, Energy
  real (dp), dimension (0:Ny+1,s-1:e+1) :: flux_rho_X, flux_rho_Y
  real (dp), dimension (0:Ny+1,s-1:e+1) :: flux_momx_X, flux_momx_Y
  real (dp), dimension (0:Ny+1,s-1:e+1) :: flux_momy_X, flux_momy_Y
  real (dp), dimension (0:Ny+1,s-1:e+1) :: flux_en_X, flux_en_Y
  integer :: below, above, ierr, grid_comm
 
  Mass = Mass - dt * dy * flux_rho_X
  ! shift up row
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, flux_rho_X )
  Mass = Mass + dt * dy * flux_rho_X
  Mass = Mass - dt * dx * flux_rho_Y
  Mass = Mass + dt * dx * CSHIFT(flux_rho_Y,L)
  Momx = Momx - dt * dy * flux_momx_X
  ! shift up row
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, flux_momx_X )
  Momx = Momx + dt * dy * flux_momx_X
  Momx = Momx - dt * dx * flux_momx_Y
  Momx = Momx + dt * dx * CSHIFT(flux_momx_Y,L)
  Momy = Momy - dt * dy * flux_momy_X
  ! shift up row
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, flux_momy_X )
  Momy = Momy + dt * dy * flux_momy_X
  Momy = Momy - dt * dx * flux_momy_Y
  Momy = Momy + dt * dx * CSHIFT(flux_momy_Y,L)
  Energy = Energy - dt * dy * flux_en_X
  ! shift up row
  call shift_band ( Ny, s, e, 0, below, above, grid_comm, flux_en_X )
  Energy = Energy + dt * dy * flux_en_X
  Energy = Energy - dt * dx * flux_en_Y
  Energy = Energy + dt * dx * CSHIFT(flux_en_Y,L)

  return
  end subroutine update_soln

  subroutine shift_band ( Ny, s, e, up_down, below, above, grid_comm, a)
  use mpi
  use types, only: dp 

  integer :: Ny, e, s

  real (dp) :: a(0:Ny+1,s-1:e+1)
  integer below, above, tag, ierr, up_down, grid_comm
  integer status(MPI_STATUS_SIZE)

  tag = 0
  call MPI_Send ( a(1,e),   Ny, MPI_DOUBLE_PRECISION, above, tag, grid_comm, ierr )
  call MPI_Recv ( a(1,s-1), Ny, MPI_DOUBLE_PRECISION, below, tag, grid_comm,status, &
     ierr )

  tag = 1
  call MPI_Send ( a(1,s),   Ny, MPI_DOUBLE_PRECISION, below, tag, grid_comm, ierr )
  call MPI_Recv ( a(1,e+1), Ny, MPI_DOUBLE_PRECISION, above, tag, grid_comm,status, &
     ierr )

  if ( up_down == 0) then 
    a = CSHIFT(a,1)
  else if ( up_down == 1) then 
    a = CSHIFT(a,-1) 
  else
    write(*,*) 'no shift in a, check error'
  end if

  return
  end subroutine shift_band

