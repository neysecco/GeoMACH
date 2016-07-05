subroutine computeTip(nD, nu, nv, f0, N, S, inds, Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nu, nv, f0, N, S, inds
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nv) N, S
  !f2py depend(nu,nv) inds
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nu, nv
  double precision, intent(in) ::  f0
  integer, intent(in) ::  N(nv,2,3), S(nv,2,3)
  integer, intent(in) ::  inds(nu,nv,3)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer i, j, k, iD
  double precision den, C(4), u, x_c, f
  real, parameter :: pi = 3.1415927

  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0

  iD = 0

  ! Border: N
  i = 1
  do j=1,nv
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = N(j, 1, k)
     end do
  end do

  ! Border: S
  i = nu
  do j=1,nv
     do k=1,3
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = inds(i, j, k)
        Dj(iD) = S(j, 1, k)
     end do
  end do

  den = 1.0 / (nu-1)
  do i=2,nu-1
     u = den * (i-1)
     do j=1,nv
        x_c = dble(j-1)/dble(nv-1)
        f = f0*(0.5+5*sqrt(dist(x_c)))
        call sparseBezier(u, -f, f, C)
        do k=1,3
           Da(iD+1:iD+4) = C(:)
           Di(iD+1:iD+4) = inds(i, j, k)
           Dj(iD+1) = N(j,1,k)
           Dj(iD+2) = N(j,2,k)
           Dj(iD+3) = S(j,1,k)
           Dj(iD+4) = S(j,2,k)
           iD = iD + 4
        end do
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeTip', iD, nD
  end if

  ! Auxiliary function

  contains

    function dist(x)
      ! This functions gives a number between 0 and 1
      ! that will be applied to the interpolant weights f0
      ! so that the tip becomes rounded when seen from a top view
      
      ! Function inputs
      double precision, intent(in) :: x

      ! Function output
      double precision :: dist

      ! Use naca airfoil polynomial
      dist = 10.0*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x*x*x - 0.1036*x*x*x*x)

    end function dist

end subroutine computeTip
