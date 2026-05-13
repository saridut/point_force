 program main
         
! DESCRIPTION: 
!   This code computes the velocity field due to a point force. 
! 
! METHOD: 
!   A point force pf is provided along with its location r0. The goal is
! to compute the velocity field at different locations due to this point force.
! The domain is given by y = 0, [xmin,xmax] and [zmin,zmax]. At z = 0, there is
! a solid wall. 
! 
! INPUT FILES: 
!   None
        
! OUTPUT FILES: 
!   The output consists of four different files.
!         (i)   X.dat :
!         (ii)  Z.dat :
!         (iii) VX.dat :
!         (iv)  VZ.dat :
 
! HISTORY: 
! Version   Date             Comment 
! -------   --------------   ----------------------------
! 1.0       March 30, 2012   Original code. (Sarit Dutta)
         
! CODE DESCRIPTION: 
!   Language:           Fortran 90. 
         
! DECLARATIONS: 
!   WALL -- Conditional compilation based on whether a wall is present at z = 0
!   or not.
!
! MODULES USED: 
 Use parameters
        
 Implicit none 
         
! Local scalars: 
 real (dbl)     :: xmin, & ! Minimum along x (Supplied by driver script)
                   xmax, & ! Maximum along x (Supplied by driver script)
                   zmin, & ! Minimum along z (Supplied by driver script)
                   zmax, & ! Maximum along z (Supplied by driver script)
                   dx,   & ! Step size along x (Supplied by driver script)
                   dz      ! Step size along z (Supplied by driver script)

 integer        :: nx, & ! # grid points along x
                   nz, & ! # grid points along z
                   i,  & ! loop counter
                   j,  & ! loop counter
                   ierr  ! error flag

 character (20) :: xminstr, & ! String storing xmin from driver script
                   xmaxstr, & ! String storing xmax from driver script
                   zminstr, & ! String storing zmin from driver script
                   zmaxstr, & ! String storing zmax from driver script
                   dxstr,   & ! String storing dx from driver script
                   dzstr,   & ! String storing dz from driver script
                   fmtstrx, & ! Format string 
                   fmtstrz, & ! Format string
                   nxstr,   & ! String containing nx for use in fmtstrx
                   nzstr      ! String containing nz for use in fmtstrz

! Local arrays: 
 real (dbl)  :: pf    (ndim), & ! Point force
                r     (ndim), & ! Any point within the domain
                r0    (ndim), & ! Location of the point force
                div_D (ndim), & ! Divergence of diffusivity tensor
                fcombined (ndim*nbead), & ! Combined force vector
                D(ndim*nbead, ndim*nbead), & ! Diffusivity tensor
                rcombined(ndim, nbead), & ! [r r0]
                A_bdi (nbead*ndim,1) ! The term multiplying dt in Fokker Planck
                                     ! eqn

 real (dbl), allocatable :: X(:), &
                            Z(:), &
                            VX(:,:), &
                            VZ(:,:), &
                            Zmat(:,:)
       
!- End of header --------------------------------------------------------------

#DEFINE WALL

 interface 
         subroutine calc_diffusivity (r, D, ierr)
         use parameters
         real (dbl), dimension(ndim, nbead), intent (in) :: r
         real (dbl), dimension(ndim*nbead, ndim*nbead), intent (out) :: D
         integer, intent (out) :: ierr
         end subroutine

         subroutine calc_div_diffusivity (r, div_D, ierr)
         use parameters
         real (dbl), dimension(ndim), intent (in) :: r
         real (dbl), dimension(ndim), intent (out) :: div_D
         integer, intent (out) :: ierr
         end subroutine
 end interface

 call get_command_argument (1, xminstr)
 call get_command_argument (2, xmaxstr)
 call get_command_argument (3, zminstr)
 call get_command_argument (4, zmaxstr)
 call get_command_argument (5, dxstr)
 call get_command_argument (6, dzstr)

 read (xminstr, *) xmin
 read (xmaxstr, *) xmax
 read (zminstr, *) zmin
 read (zmaxstr, *) zmax
 read (dxstr, *) dx
 read (dzstr, *) dz

 pf = [1._dbl, 0._dbl, 0._dbl]    ! Flip sign of the first element to get pf
                                  ! in opposite direction
 r0 = [-4._dbl, 0._dbl, 5._dbl]*a ! Flip -4 to 4 to get the RHS figure in the
                                  ! paper.

 nx = 1 + int((xmax-xmin)/dx)
 nz = 1 + int((zmax-zmin)/dz)

 allocate (X (nx))
 allocate (Z (nz))
 allocate (VX (nz, nx))
 allocate (VZ (nz, nx))
 allocate (Zmat (nz, nx))

 write (nxstr,*) nx
 write (nzstr,*) nz
 print*, 'nx = ', nx
 print*, 'nz = ', nz

 fmtstrx = '('//trim(adjustl(nxstr))//'(f14.4,1x))'
 fmtstrz = '('//trim(adjustl(nzstr))//'(f14.4,1x))'

 rcombined(:,1) = 0._dbl
 rcombined(:,2) = r0
 fcombined (1:ndim) = 0._dbl
 fcombined (ndim+1:2*ndim) = pf

 open (100, file = 'X1.dat', action = 'write', status = 'unknown')
 open (200, file = 'Z1.dat', action = 'write', status = 'unknown')
 open (300, file = 'VX1.dat', action = 'write', status = 'unknown')
 open (400, file = 'VZ1.dat', action = 'write', status = 'unknown')

 do i = 1, nx
    if (i==nx) then
       X(i) = xmax
    else
       X(i) = xmin + (i-1)*dx
    end if
 end do

 do j = 1, nz
    if (j==nz) then
       Z(j) = zmax
    else
       Z(j) = zmin + (j-1)*dz
    end if
 end do

 do i = 1, nx
    do j = 1, nz
       r = [X(i), 0._dbl, Z(j)]
       rcombined(:,1) = r
       call calc_diffusivity (rcombined, D, ierr)
!      call calc_div_diffusivity (r, div_D, ierr)
       A_bdi = matmul(D, reshape(fcombined,[2*ndim,1]))
!      A_bdi (1:ndim,1) = A_bdi (1:ndim,1)  + div_D
       VX (j,i) = A_bdi(1,1)
       VZ (j,i) = A_bdi(3,1)
!      vx = D(1,4)*pf(1) + D(1,5)*pf(2) + D(1,6)*pf(3)
!      vz = D(3,4)*pf(1) + D(3,5)*pf(2) + D(3,6)*pf(3)
    end do
 end do

 do i = 1, nx
    Zmat(:,i) = Z
 end do

 do j = 1, nz
    write (100, fmtstrx) X
    write (200, fmtstrx) Zmat(j,:)
    write (300,fmtstrx) VX (j,:)
    write (400,fmtstrx) VZ (j,:)
 end do


 close (100)
 close (200)
 close (300)
 close (400)

 end program main

!******************************************************************************
 subroutine calc_diffusivity (r, D, ierr)

 use parameters

 implicit none

 real (dbl), dimension(ndim, nbead), intent (in)            :: r
 real (dbl), dimension(ndim*nbead,ndim*nbead), intent (out) :: D
 integer, intent (out)             :: ierr
 real (dbl), dimension(ndim*nbead) :: omega_w_jj, egn
 real (dbl), dimension(ndim)      :: ri, rj, rij
 real (dbl), dimension(ndim,1)    :: rijmat
 real (dbl), dimension(ndim,ndim) :: rijrij, omega_wij
 real (dbl) :: rijm, irijm, irijm2, C1, C2, consij
 integer    :: i, j, ii, jj, m, n, info

 interface
   function func_eye (n)
     use parameters
     integer, intent (in) :: n
     real (dbl), dimension(n,n) :: func_eye
   end function func_eye
 end interface

 ierr = 0
 D = 0._dbl
 omega_w_jj = 0._dbl

#IFNDEF WALL
!Calculate the RPY tensor (in strictly lower triangular form)
 do j = 1, nbead-1
    rj = r(:,j)
    jj = ndim*(j-1)+1
   do i = j+1, nbead
      ri = r(:,i)
      ii = ndim*(i-1)+1
      rij = ri-rj
      rijm = sqrt(dot_product(rij, rij))

      if (rijm==0._dbl) then
        write(*,*) 'calc_diffusivity: rijm = 0'
        ierr = 7
        return
      endif

      irijm = 1._dbl/rijm
      irijm2 = irijm*irijm
      rijmat = reshape(rij,[ndim,1])
      rijrij = matmul(rijmat,transpose(rijmat))

      if (rijm .ge. 2._dbl*a) then
        C1 =  1._dbl + (2._dbl*a*a/3._dbl)*irijm2
        C2 =  1._dbl - 2._dbl*a*a*irijm2
        consij = 0.75_dbl*a*irijm
      else
        C1 = 1._dbl - 9._dbl*rijm/(32._dbl*a)
        C2 = 3._dbl*rijm/(32._dbl*a)
        consij = 1._dbl
      end if

      D(ii:ii+2,jj:jj+2) = consij*(C1*func_eye(3) + C2*rijrij*irijm2)

   end do
 end do

 D = D + transpose(D)

 do j = 1,ndim*nbead
    D(j,j) = 1._dbl
 end do

#ENDIF

#IFDEF WALL
 do j = 1, nbead-1
    rj = r(:,j)
    jj = ndim*(j-1)+1

    if (rj(3) .ne. 0._dbl) then
       omega_w_jj(jj)   = 0.125_dbl*(a/rj(3))*(a/rj(3))*(a/rj(3)) &
                          - 0.5625_dbl*(a/rj(3))
       omega_w_jj(jj+1)  = 0.125_dbl*(a/rj(3))*(a/rj(3))*(a/rj(3)) &
                          - 0.5625_dbl*(a/rj(3))
       omega_w_jj(jj+2) = 0.5_dbl*(a/rj(3))*(a/rj(3))*(a/rj(3))   &
                          - 1.125_dbl*(a/rj(3))
   end if

   do i = j+1, nbead
      ri = r(:,i)
      ii = ndim*(i-1)+1

      rij = ri-rj
      rijm = sqrt(dot_product(rij, rij))
      if (rijm==0._dbl) then
        write(*, *)  'calc_diffusivity: rijm = 0'
        ierr = 8
        return
      end if
      irijm = 1._dbl/rijm
      irijm2 = irijm*irijm
      rijmat = reshape(rij,[ndim,1])
      rijrij = matmul(rijmat,transpose(rijmat))

      if (rijm .ge. 2._dbl*a) then
        C1 =  1._dbl+(2._dbl*a*a/3._dbl)*irijm2
        C2 =  1._dbl-2._dbl*a*a*irijm2
        consij = 0.75_dbl*a*irijm
      else
        C1 = 1._dbl-9._dbl*rijm/(32._dbl*a)
        C2 = 3._dbl*rijm/(32._dbl*a)
        consij = 1._dbl
      end if

      call calc_omega_w (ri, rj, omega_wij, ierr)

      if (ierr .ne. 0) return

      D(ii:ii+2,jj:jj+2) = consij*(C1*func_eye(3) + C2*rijrij*irijm2) + &
                            omega_wij
   end do
 end do

 rj = r(:,nbead)

 if (rj(3) .ne. 0._dbl) then
    omega_w_jj(ndim*nbead-2) = 0.125_dbl*(a/rj(3))*(a/rj(3))*(a/rj(3)) &
                            - 0.5625_dbl*(a/rj(3))
    omega_w_jj(ndim*nbead-1) = 0.125_dbl*(a/rj(3))*(a/rj(3))*(a/rj(3)) &
                            - 0.5625_dbl*(a/rj(3))
    omega_w_jj(ndim*nbead)   = 0.5_dbl*(a/rj(3))*(a/rj(3))*(a/rj(3))   &
                            - 1.125_dbl*(a/rj(3))
 else
    write(*,*) 'Error: rj(3) <= 0'
    ierr = 11
    return
 end if

 D = D + transpose(D)

 do j = 1, ndim*nbead
    D(j,j) = D(j,j) + 1._dbl + omega_w_jj(j)
 end do

#ENDIF

 end subroutine calc_diffusivity

!******************************************************************************
#IFDEF WALL
!The FPP conditional above ends after subroutine func_delta

 subroutine calc_omega_w (ri, rj, omega_w, ierr)

 use parameters

 implicit none

 real (dbl), dimension(ndim), intent(in) :: ri, rj
 real (dbl), dimension(ndim,ndim), intent(out) :: omega_w
 integer, intent (out)            :: ierr
 real (dbl), dimension(ndim,ndim) :: omega_wc, blake_tensor, RR
 real (dbl), dimension(ndim,1) :: Rmat
 real (dbl), dimension(ndim) :: R
 real (dbl) :: xij, yij, zij, zi, zj, xij2, yij2, zij2, iRij, iRij3, iRij5, &
               iRij7, Rij, Sij, PDij, SDij
 integer  :: i, j

 real (dbl) :: func_delta

 ierr = 0

 xij   = ri(1)-rj(1)
 yij   = ri(2)-rj(2)
 zij   = ri(3)+rj(3)
 zi    = ri(3)
 zj    = rj(3)
 xij2  = xij*xij
 yij2  = yij*yij
 zij2  = zij*zij
 R     = [xij,yij,zij]

 Rij   = sqrt(xij2+yij2+zij2)

 if (Rij==0._dbl) then
   write(*, *)  'calc_omega_w: Rij = 0'
   ierr = 9
   return
 end if

 iRij  = 1._dbl/Rij
 iRij3 = iRij*iRij*iRij
 iRij5 = iRij3*iRij*iRij
 iRij7 = iRij5*iRij*iRij

 omega_wc = 0._dbl

 omega_wc (1,1) = iRij3 - 3._dbl*(xij2+zij2)*iRij5 + 15._dbl*xij2*zij2*iRij7
 omega_wc (2,1) = -3._dbl*xij*yij*iRij5 + 15._dbl*xij*yij*zij2*iRij7
 omega_wc (1,2) = omega_wc (2,1)
 omega_wc (2,2) = iRij3 - 3._dbl*(yij2+zij2)*iRij5 + 15._dbl*yij2*zij2*iRij7
 omega_wc (3,3) = iRij3 + 6._dbl*zij2*iRij5 - 15._dbl*zij2*zij2*iRij7

 omega_wc = 0.75*a*omega_wc

 Rmat = reshape(ri,[ndim,1])
 RR = matmul(Rmat,transpose(Rmat))

 do j = 1,ndim
   do i = 1,ndim
!     Blake tensor from Hoda & Kumar: +/- written as (1-2*delta(j,3)). This
!     is -1 for j = 3.
     Sij = func_delta(i,j)*iRij + R(i)*R(j)*iRij3
     PDij = (1._dbl-2._dbl*func_delta(j,3))*(func_delta(i,j)*iRij3 &
             - 3._dbl*R(i)*R(j)*iRij5)
     SDij = R(3)*PDij + (1._dbl-2._dbl*func_delta(j,3))*(func_delta(j,3)*R(i) &
            -func_delta(i,3)*R(j))*iRij3

     blake_tensor(i,j) = 0.75_dbl*a*(-Sij + 2._dbl*zj*zj*PDij - 2._dbl*zj*SDij)

   end do
 end do

! blake_tensor = 0.75*a*blake_tensor

 omega_w = blake_tensor - (2._dbl*a*a/3._dbl)*omega_wc

 end subroutine calc_omega_w

!******************************************************************************
 function func_delta (x,y)

 use parameters

 implicit none

 integer, intent (in) :: x, y
 real (dbl) :: func_delta

 if (x==y) then
   func_delta = 1._dbl
 else
   func_delta = 0._dbl
 end if

 end function func_delta

!The FPP conditional below starts before subroutine calc_omega_w

#ENDIF

!******************************************************************************
 function func_eye (n)

 use parameters

 implicit none
 integer, intent (in) :: n
 real (dbl), dimension(n,n) :: func_eye
 integer  :: i

 func_eye = 0._dbl
 do i = 1, n
    func_eye(i,i) = 1._dbl
 end do

 end function func_eye

!******************************************************************************
 subroutine calc_div_diffusivity (r, div_D, ierr)

 use parameters

 implicit none

 real (dbl), dimension(ndim), intent (in) :: r
 real (dbl), dimension(ndim), intent (out) :: div_D
 integer, intent (out) :: ierr
 real (dbl) :: izi, izi2
 real (dbl), save :: pf_1, pf_2
 logical, save :: first_call = .true.

 ierr = 0

 if (first_call) then
    pf_1 = 1.125_dbl*a
    pf_2 = 1.5_dbl*a*a*a
    first_call = .false.
 end if

 div_D = 0._dbl

! The following part only if wall exists
#IFDEF WALL
 if (r(3)==0._dbl) then
   write(*, *) 'calc_div_diffusivity = 0'
   ierr = 1
   return
 end if

 izi = 1._dbl/r(3)
 izi2 = izi*izi
 div_D(ndim) = izi2*(pf_1-pf_2*izi2)
#ENDIF

 end subroutine calc_div_diffusivity
!******************************************************************************
