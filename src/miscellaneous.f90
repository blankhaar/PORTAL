subroutine check_neigh(X_phot,node,X_node,neigh_num,neigh_node,X_neigh_node,new_node)
!function to check which node is closest to the photon
!node is the node the photon was closest to in previous step
implicit none
integer node,neigh_num
integer new_node,i

integer, dimension(neigh_num) :: neigh_node

double precision dis,dis_0

double precision, dimension (3) :: x_0,x_n 
double precision, dimension (3) :: X_phot,X_node
double precision, dimension (3,neigh_num) :: X_neigh_node 

!input:
!X_phot       -- photon position
!node         -- current node
!neigh_num    -- number of neighbours
!neigh_node   -- neighbour numbers of current node
!X_neigh_node -- positions of neighbours of current node
!X_node       -- position of current node 

!output:
!new_node     -- closest node to position X_phot 

x_0(:) = X_phot(:) - X_node(:) 

dis_0 = x_0(1)**2.d0 + x_0(2)**2.d0 + x_0(3)**2.d0
new_node = node

do i = 1,neigh_num
  x_n(:) = X_phot(:) - X_neigh_node(:,i)

  dis = x_n(1)**2.d0 + x_n(2)**2.d0 + x_n(3)**2.d0

!find closest neighbour. This is the new_node
!if dis_0 is smallest from the beginning, new_node = node
  if (dis.lt.dis_0) then
    dis_0 = dis
    new_node = neigh_node(i)
  endif
end do 


end subroutine

subroutine find_dir(bvec,gvec,th,ph,dir)
!function to find the direction vector from the magnetic
!field direction and the gauge direction. 
implicit none
integer i

double precision xnorm,th,ph

double precision, dimension(3) :: bvec,gvec,xvec,dir
double precision, dimension(3,3) :: Rth,Rph

!input:
!bvec   -- magnetic field direction
!gvec   -- gauge direction
!th,ph  -- direction angles

!output:
!dir    -- ray direciton vector


!generate x-axis
xvec(1) = bvec(2) * gvec(3) - bvec(3) * gvec(2) 
xvec(2) = bvec(3) * gvec(1) - bvec(1) * gvec(3) 
xvec(3) = bvec(1) * gvec(2) - bvec(2) * gvec(1) 

!normalize
xnorm = dsqrt(xvec(1)**2.d0 + xvec(2)**2.d0 + xvec(3)**2.d0) 
xvec(:) = xvec(:) / xnorm

call rotmat(xvec,th,Rth)
call rotmat(bvec,ph,Rph)

call two_mat_vec_prod(Rph,Rth,bvec,dir)

end subroutine



subroutine find_angles(bvec,gvec,dir,th,ph)
!function to find the angles from the direction vector 
implicit none
integer i

double precision xnorm,th,ph

double precision, dimension(3) :: bvec,gvec,xvec,dir,cbd
double precision, dimension(3,3) :: Rth,Rph

!input:
!bvec   -- magnetic field direction
!gvec   -- gauge direction
!dir    -- ray direciton vector

!output:
!th,ph  -- direction angles


!generate x-axis
xvec(1) = bvec(2) * gvec(3) - bvec(3) * gvec(2)
xvec(2) = bvec(3) * gvec(1) - bvec(1) * gvec(3)
xvec(3) = bvec(1) * gvec(2) - bvec(2) * gvec(1)

!normalize
xnorm = dsqrt(xvec(1)**2.d0 + xvec(2)**2.d0 + xvec(3)**2.d0)
xvec(:) = xvec(:) / xnorm

if (dot_product(bvec,dir).gt.1.d0) then 
  th = dacos(1.d0)
elseif (dot_product(bvec,dir).lt.-1.d0) then
  th = dacos(-1.d0)
else 
  th = dacos(dot_product(bvec,dir))
endif

!cross-product: bvec x dir
cbd(1) = bvec(2)*dir(3) - bvec(3)*dir(2) 
cbd(2) = bvec(3)*dir(1) - bvec(1)*dir(3) 
cbd(3) = bvec(1)*dir(2) - bvec(2)*dir(1) 

ph = datan2(dot_product(xvec,dir),dot_product(xvec,cbd))

!call rotmat(xvec,th,Rth)
!call rotmat(bvec,ph,Rph)
!
!call two_mat_vec_prod(Rph,Rth,bvec,dir)

end subroutine

FUNCTION cross(a, b)
  INTEGER, DIMENSION(3) :: cross
  INTEGER, DIMENSION(3), INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross


subroutine rotmat(nv,th,Rot)
!set up a rotation matrix of rotation vector nv and angle th
implicit none
integer i,j
double precision th,cth,sth,x0,x1,octh
double precision, dimension(3) :: nv
double precision, dimension(3,3) :: Rot

!input:
!nv     -- rotation axis 
!th     -- rotation angle

!output:
!Rot    -- rotation matrix

x0 = 0.d0
x1 = 1.d0

cth = dcos(th)
octh = x1 - cth
sth = dsin(th)

do i = 1,3
  Rot(i,i) = cth + octh*nv(i)*nv(i)
end do

Rot(1,2) = nv(1)*nv(2)*octh - nv(3)*sth
Rot(2,1) = nv(2)*nv(1)*octh + nv(3)*sth

Rot(1,3) = nv(1)*nv(3)*octh + nv(2)*sth
Rot(3,1) = nv(3)*nv(1)*octh - nv(2)*sth

Rot(2,3) = nv(2)*nv(3)*octh - nv(1)*sth
Rot(3,2) = nv(3)*nv(2)*octh + nv(1)*sth

end subroutine

subroutine two_mat_vec_prod(A,B,c,d)
!computes d = A*B*c 
!A,B -- 3x3 matrices
!c,d -- 3   vectors
!A,B,C input
!d output
implicit none
integer i,j,k

double precision, dimension (3) :: c,d
double precision, dimension (3,3) :: A,B

do i = 1,3
  d(i) = 0.d0
  do j = 1,3
    do k = 1,3
      d(i) = d(i) + A(i,j) * B(j,k) * c(k)
    end do
  end do
end do

end subroutine

subroutine setup_count_array(j_lev,nlev,kmax,count_array,i_count)
!subroutine that sets up the count-matrix. The matrix that will keep track of
!what count-number in the stateq-matrix, will correspond to what state (l,j,k)
implicit none
integer nlev,c_dim,i,j,k,i_count,kmax

integer, dimension(nlev) :: j_lev
integer, dimension(nlev) :: count_array 

!input:
!j_lev  -- array with ang. mom. j associated with the levels
!nlev   -- number of levels

!output:
!count_array 	-- count_element where this level starts


count_array(1) = 0
i_count = 0

do i = 1,nlev-1
  j = j_lev(i)
  do k = 0,min(2*j,kmax),2
    i_count = i_count + 1
  end do
  count_array(i+1) = i_count 
end do 

j = j_lev(nlev)
do k = 0,min(2*j,kmax),2
  i_count = i_count + 1
end do

end subroutine 


subroutine matvec_prod(A,b,N,c)
!compute the product A*b = c
!A = NxN matrix
!b = N-vector
!c = N-vector
integer N,I,J
double precision, dimension (N) :: b,c
double precision, dimension (N,N) :: A

do I=1,N
  c(I) = 0.d0
end do

do I = 1,N
  do J = 1,N
    c(I) = c(I) + A(I,J)*b(J)
  end do
end do

end subroutine

subroutine black_body(T,nu0,I_bb)
!black body radiation
implicit none

double precision nu0,I_bb,T,c1,c2

!input:
!T      -- temperature
!nu0    -- frequency

!output:
!I_bb   -- black-body radiation

if (T.eq.0.d0) then
  I_bb = 0.d0
else
  c1 = 1.4744995D-50    !2*h/c^2
  c2 = 4.79924D-11/T    !h/(k_B*T)
  
  I_bb = c1 * nu0**3.d0 * (dexp(c2*nu0) - 1.d0)**(-1.d0)
endif

end subroutine


subroutine make_AC(A,C,lA,lC,nA,nC,nT,nP,nL,E,nH2,Temp,col_temps,n_trans,l_trans,J0,J2,AC)
!subroutine to construct the AC matrix --- the matrix containing all te relevant
!information for the SEE equations
implicit none
integer nA,nC,nT,nL
integer i,j,k,l,n,n_trans,nP

integer, dimension(2,nA) :: lA
integer, dimension(2,nC) :: lC
integer, dimension(2,n_trans) :: l_trans

double precision nH2,Temp,ddT,fac

double precision, dimension(nA) :: A,J0,J2
double precision, dimension(nT,nC,nP) :: C
double precision, dimension(nC) :: C_r,f1,f2,f3,f4
double precision, dimension(nC) :: ff1,ff2 
double precision, dimension(nT) :: col_temps 
double precision, dimension(nL) :: E 
double precision, dimension(5,n_trans) :: AC

!input
!A	    -- einstein A coefficients
!C	    -- collisional rates (s-1 / cm^3)
!lA	    -- radiative transition levels
!lC	    -- collisional transition levels
!nA	    -- number of radiative transitions
!,nC	    -- number of collisional transitions
!nT	    -- number of temperatures in the collision rate-matrix
!nL	    -- number of levels
!E	    -- level energies (in Hz)
!nH2	    -- H2 number density
!Temp	    -- kinetic temperature
!col_temps  -- collision temperatures (for interpolation)
!n_trans    -- total number of unique transitions
!l_trans    -- level numbers of unique transitions
!J0	    -- isotropic radiation intensities
!J2	    -- anisotropic radiation intensities

!output
!AC	    -- transition coefficients (Aji,Cji,nu0,J0,J2)

!interpolate the temperature
if (Temp.le.col_temps(1).and.nP.eq.1) C_r(:) = C(1,:,1)
if (Temp.le.col_temps(1).and.nP.eq.2) C_r(:) = 0.25d0*C(1,:,1) + 0.75d0*C(1,:,2)

do i = 2,nT  
  if (col_temps(i).gt.Temp) then

    ddT = (Temp - col_temps(i-1)) / (col_temps(i) - col_temps(i-1))
    if (nP.eq.1) then
      f1 = C(i-1,:,1)
      f2 = C(i,:,1)
  
      C_r(:) = f1(:) + ddT*(f2(:)-f1(:))
    elseif (nP.eq.2) then !the case for o-H2 and p-H2 collisional partners 
      f1 = C(i-1,:,1)
      f2 = C(i,:,1)

      f3 = C(i-1,:,2)
      f4 = C(i,:,2)

      ff1 = f1(:) + ddT*(f2(:)-f1(:))
      ff2 = f3(:) + ddT*(f4(:)-f3(:))

      C_r(:) = 0.25d0*ff1 + 0.75d0*ff2         !25% para       
                                               !75% ortho 
    endif

    exit

  endif
end do

!we need the collisional rates times the collisional partner 
C_r(:) = C_r(:) * nH2 * 1d-6           !1d-6 factor for conversion nH2 to cm-3
!write(*,*)Temp,nH2,C_r(1)
!write(*,*)C_r(1)
!stop

!now set up the AC matrix
AC(:,:) = 0.d0
do i = 1,n_trans
!transition frequency
  AC(3,i) = E(l_trans(2,i)) - E(l_trans(1,i))
!is this a radiative transition?
  do j = 1,nA
    if (lA(1,j).eq.l_trans(1,i).and.lA(2,j).eq.l_trans(2,i)) then
      AC(1,i) = A(j)
      AC(4,i) = J0(j)
      AC(5,i) = J2(j)
      exit !exit, only one unique transition 
    endif
  end do 
!is this a collisional transition?
  do j = 1,nC
    if (lC(1,j).eq.l_trans(1,i).and.lC(2,j).eq.l_trans(2,i)) then
      AC(2,i) = C_r(j)
      exit !exit, only one unique transition 
    endif
  end do 
end do


end subroutine

subroutine all_transitions(nL,nA,nC,lA,lC,E,l_trans,n_trans)
!find all unique upper-to-lower transitions 
implicit none
integer nL,nA,nC,n_trans
integer i,j,k

integer, dimension(nL,nL) :: T_mat
integer, dimension(2,nA+nC) :: l_trans
integer, dimension(2,nA) :: lA
integer, dimension(2,nC) :: lC

double precision, dimension(nL) :: E

!input:
!nL	-- number of levels
!nC	-- number of collisional transitions
!nA	-- number of radiative transitions
!lA	-- radiative transition levels
!lC	-- collisional transition levels
!E	-- level energies

!output
!l_trans-- all tranistion levels
!n_trans-- total number of transition 

!we construct a T-matrix that is 1 if there is a 
!transition that is lower to upper (energy)
!can be either radiative or collisional

T_mat(:,:) = 0

do i = 1,nA
  j = lA(1,i)
  k = lA(2,i)
  if (E(k).gt.E(j)) then
    T_mat(j,k) = 1 
  else
    write(*,*)'attention'
    T_mat(k,j) = 1
  endif
end do

do i = 1,nC
  j = lC(1,i)
  k = lC(2,i)
  if (E(k).gt.E(j)) then
    T_mat(j,k) = 1 
  else
    write(*,*)'attention'
    T_mat(k,j) = 1
  endif
end do

!now translate the T-matrix to how we want it:
!a transition-level counter 
n_trans = 0
do i = 1,nL
  do j = 1,nL
    if (T_mat(i,j).ne.0) then
      n_trans = n_trans + 1
      l_trans(1,n_trans) = i
      l_trans(2,n_trans) = j
    endif
  end do
end do

end subroutine


subroutine rot_stokes(Q,U,Qnode,Unode,bvec)
!subroutine to rotate the stokes vectors to the appropriate coordinate system
!of the node. Subroutine assumes that the ray-direction is along the z-axis, 
!and the canonical axis is along the x-axis
!tested for some standard cases (test_dipfield.f90) seems to work fine
implicit none

double precision alpha,Q,U,Qnode,Unode

double precision, dimension (3) :: bvec,chi

!input
!Q,U	-- Stokes Q,U parameters
!bvec	-- local magnetic field direction (normalized)

!output
!Qnode 	-- Stokes Q parameter local to node
!Unode	-- Stokes U parameter local to node


!compute the projection of bvec onto the plane of the sky
chi(1) = bvec(1)
chi(2) = bvec(2)
chi(3) = 0.d0 

chi(:) = chi(:) / norm2(chi)

!now rotate the axis system to the local node axis system that is
!gauged on the local magnetic field. Global Q and U are gauged on
!the canonical y-axis.

!chi = |    ---- alpha = 0
!chi = /    ---- alpha = 45
!chi = -    ---- alpha = 90
!chi = \    ---- alpha = -45

alpha = dacos(chi(2))
if (chi(1).lt.0.d0) alpha = -alpha

Qnode = Q * dcos(2.d0*alpha) + U * dsin(2.d0*alpha)
Unode = U * dcos(2.d0*alpha) - Q * dsin(2.d0*alpha)

end subroutine

subroutine rot_stokes_back(Q,U,Qnode,Unode,bvec)
!subroutine to rotate the stokes vectors back to the appropriate coordinate system
!of the node. Subroutine assumes that the ray-direction is along the z-axis, 
!and the canonical axis is along the x-axis
!tested for some standard cases (test_dipfield.f90) seems to work fine
implicit none

double precision alpha,Q,U,Qnode,Unode

double precision, dimension (3) :: bvec,chi

!input
!Q,U	-- Stokes Q,U parameters
!bvec	-- local magnetic field direction (normalized)

!output
!Qnode 	-- Stokes Q parameter local to node
!Unode	-- Stokes U parameter local to node


!compute the projection of bvec onto the plane of the sky
chi(1) = bvec(1)
chi(2) = bvec(2)
chi(3) = 0.d0 

chi(:) = chi(:) / norm2(chi)

!now rotate the axis system to the local node axis system that is
!gauged on the local magnetic field. Global Q and U are gauged on
!the canonical y-axis.

!chi = |    ---- alpha = 0
!chi = /    ---- alpha = 45
!chi = -    ---- alpha = 90
!chi = \    ---- alpha = -45

alpha = dacos(chi(2))
if (chi(1).lt.0) alpha = -alpha

Qnode = Q * dcos(2.d0*alpha) - U * dsin(2.d0*alpha)
Unode = U * dcos(2.d0*alpha) + Q * dsin(2.d0*alpha)

end subroutine

subroutine dipole_magfield(node_num,XYZ,x1,x2,x3,bvec_all)
!subroutine to define the global magnetic field, and assign each node a
!magnetic 
!field direction. 
!This is a dipole magnetic field.  
!dipole magnetic field works as it should. Tested against matlab routine
implicit none
integer i,node_num

double precision x1,x2,x3

double precision, dimension(3) :: bvec,xvec
double precision, dimension(3,node_num) :: bvec_all,XYZ

!input:
!node_num       -- number of nodes
!XYZ            -- node coordinates

!output:
!bvec_all       -- magnetic field orienation for all the nodes

!mvec(1) = 0
!mvec(2) = 0
!mvec(3) = 1.d0

!B(r) ~ 3x(m*x) - m
!x and m are unit vectors

do i = 1,node_num
  xvec(:) = XYZ(:,i)

  xvec(:) = xvec(:)/norm2(xvec(:))

  bvec(1) = 3.d0*xvec(1)*xvec(3)
  bvec(2) = 3.d0*xvec(2)*xvec(3)
  bvec(3) = 3.d0*xvec(3)*xvec(3) - 1.d0

  bvec(:) = bvec(:) / norm2(bvec(:))

  bvec_all(:,i) = bvec(:)
end do

end subroutine

subroutine global_magfield(node_num,XYZ,x1,x2,x3,bvec_all)
!subroutine to define the global magnetic field, and assign each node a
!magnetic field direction. 
implicit none
integer i,node_num

double precision x1,x2,x3

double precision, dimension(3) :: bvec,xvec
double precision, dimension(3) :: bvec_rad,bvec_tor,bvec_pol
double precision, dimension(3,node_num) :: bvec_all,XYZ

!input:
!node_num       -- number of nodes
!XYZ            -- node coordinates
!x1,x2,x3       -- coefficients that determine character of magnetic
!                  field 

!output:
!bvec_all       -- magnetic field orienation for all the nodes

do i = 1,node_num
  xvec(:) = XYZ(:,i)

  bvec_rad(:) = xvec(:)/norm2(xvec(:))
   
  bvec_tor(1) = xvec(2)
  bvec_tor(2) = -xvec(1)
  bvec_tor(3) = 0.d0
    
  bvec_tor(:) = bvec_tor(:) / norm2(bvec_tor(:))

  bvec_pol(1) = xvec(3)*xvec(1)
  bvec_pol(2) = xvec(2)*xvec(3)
  bvec_pol(3) = -(xvec(2)*xvec(2) + xvec(1)*xvec(1))

  bvec_pol(:) = bvec_pol(:) / norm2(bvec_pol(:))

  bvec(:) = x1*bvec_rad(:) + x2*bvec_tor(:) + x3*bvec_pol(:)
  bvec(:) = bvec(:) / norm2(bvec(:)) 
 
  bvec_all(:,i) = bvec(:) 
end do
    
end subroutine

subroutine dust_to_opacities(nu,t_num,dust_table,n_dt,kd_v)
!subroutine to convert the jena_thin_e6 tables to the correct opacities
!for the dust for each transition-frequency 
implicit none
integer i,j,t_num,n_dt

double precision mu_res,mu0,nu_to_mu,kd_m,cgs_to_si

double precision, dimension(t_num) :: nu,kd_v
double precision, dimension(2,n_dt) :: dust_table

!input:
!nu             -- frequencies that are to be investigated
!t_num          -- number of transitions
!dust_table(1,:)-- log10(wavelength in micron )
!dust_table(2,:)-- log10(dust absorption in cm^2/g )
!n_dt           -- number of dust_table inputs

!output:
!kd_v           -- dust opacities at the appropriate frequencies

!we interpolate and extrapolate the log(\lambda) and the log(\kappa) 
!as is common.

nu_to_mu = 2.998E14        !Hz to micron

cgs_to_si = 1E-1           !cm2/g to m2/kg

do i = 1,t_num
  mu0 = nu_to_mu / nu(i)
  mu0 = dlog10(mu0)
!  write(*,*)nu(i),mu0

  do j = 1,n_dt-1
    if (mu0.gt.dust_table(1,j).and.mu0.lt.dust_table(1,j+1)) then
      mu_res = (mu0 - dust_table(1,j)) / (dust_table(1,j+1) - dust_table(1,j))

!      write(*,*)mu0,dust_table(1,j),mu_res      

      kd_m = dust_table(2,j) + mu_res * (dust_table(2,j+1) - dust_table(2,j))
!      write(*,*)kd_m
      exit
    endif
  end do
  if (mu0.gt.dust_table(1,n_dt)) then
    mu_res = (mu0 - dust_table(1,n_dt)) / (dust_table(1,n_dt) - dust_table(1,n_dt-1))
!    write(*,*)mu0,dust_table(1,n_dt),mu_res
 
    kd_m = dust_table(2,n_dt) + mu_res * (dust_table(2,n_dt)-dust_table(2,n_dt-1))
!    write(*,*)kd_m
  endif

  !now, we want to have the bugger in SI
  kd_v(i) = 10**(kd_m) * cgs_to_si
end do

!write(*,*)
end subroutine



  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from flib/hpsort_eps
  !---------------------------------------------------------------------
  subroutine hpsort_eps_epw (n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
  integer, intent(in)   :: n  
  double precision, intent(in)  :: eps
  integer :: ind (n)  
  double precision :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind  
  double precision :: rra  
!
  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
    double precision :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
end subroutine hpsort_eps_epw


subroutine rotate_all_save(node_num,th_inc,ph_az,vel_all,bvec_all,XYZ,XYZ_star,&
&vel_new,bvec_new,XYZ_new,XYZ_star_new)
!rotate axis system to ray-tracing axis system by rotating all the vector
!quantities.
!tested against matlab routines -- rotations are being performed properly
implicit none
integer i,j,k
integer node_num

double precision th_inc,ph_az

double precision, dimension (3,node_num) :: vel_all,XYZ,bvec_all
double precision, dimension (3,node_num) :: vel_new,XYZ_new,bvec_new
double precision, dimension (3,3) :: R_th,R_ph,R_tot
double precision, dimension (3) :: xvec,zvec
double precision, dimension (3) :: XYZ_star,XYZ_star_new

!input
!node_num       -- total number of nodes
!th_inc         -- inclination angle
!ph_az          -- azimuthal angle

!input/output
!XYZ            -- old/new node coordinates
!vel_all        -- old/new 3D node velocities
!bvec_all       -- old/new magnetic field directions


!x' = R_az * R_inc * x
!inclination: rotation around x
xvec(:) = 0.d0
xvec(1) = 1.d0

!azimuth: rotation around z
zvec(:) = 0.d0
zvec(3) = 1.d0

call rotmat(xvec,th_inc,R_th)
call rotmat(zvec,ph_az,R_ph)

R_tot = matmul(R_ph,R_th)

vel_new(:,:) = 0.d0
XYZ_new(:,:) = 0.d0
bvec_new(:,:) = 0.d0

do i = 1,node_num
  do j = 1,3
    do k = 1,3
      vel_new(j,i) = vel_new(j,i) + R_tot(j,k) * vel_all(k,i)
      bvec_new(j,i) = bvec_new(j,i) + R_tot(j,k) * bvec_all(k,i)
      XYZ_new(j,i) = XYZ_new(j,i) + R_tot(j,k) * XYZ(k,i)
    end do
  end do
end do

XYZ_star_new(:) = 0.d0
do j = 1,3
  do k = 1,3
    XYZ_star_new(j) = XYZ_star_new(j) + R_tot(j,k) * XYZ_star(k)
  end do
end do

end subroutine


subroutine adj_list_no_delaunay(in_grid,del_num,del_ne,del_num_ne,max_ne_temp,max_ne,node_num)
!algorithm to construct the neighbour list from the delaunay triangulation
!in a memory efficient manner. 
use SortUnique
implicit none
character(len=100) :: in_grid

integer del_num,node_num,max_ne
integer i,j,max_ne_temp

integer, dimension(3) :: q1,q2,q3,q4 
integer, dimension(4) :: del_tri 
integer, dimension(:), allocatable :: un_ne 
integer, dimension(node_num,max_ne_temp) :: del_ne
integer, dimension(node_num) :: del_num_ne 
integer, dimension(node_num) :: count_ne 

!input
!del_num        -- number of Delaunay Triangles
!node_num       -- number of simulation nodes

!output:
!del_ne         -- list of Delaunay neighbours 
!del_num_ne     -- number of Delaunay neighbours 
!max_ne         -- maximum number of Delauanay neighbours

!we don't need anything from the positions 
open(unit=123,file=in_grid)
read(123,*)
read(123,*)
read(123,*)
read(123,*)
read(123,*)
do i = 1,node_num
  read(123,*)
end do
read(123,*)
read(123,*)

count_ne(:) = 1

q1(1) = 2
q1(2) = 3
q1(3) = 4

q2(1) = 1
q2(2) = 3
q2(3) = 4

q3(1) = 1
q3(2) = 2
q3(3) = 4

q4(1) = 1
q4(2) = 2
q4(3) = 3
!find all the neighbours for each node. Do not distinct if they are already
!there.  
do i = 1,del_num
!for del_tri(1,i) add del_tri(q1,i) to the del_ne(del_tri(1,i),:) list
!  write(*,*)del_tri(:,i)
!  read(123,*)del_tri(1),del_tri(2),del_tri(3),del_tri(4)
  read(123,*)j,del_tri(1),del_tri(2),del_tri(3),del_tri(4)

  del_tri(1) = del_tri(1)+1
  del_tri(2) = del_tri(2)+1
  del_tri(3) = del_tri(3)+1
  del_tri(4) = del_tri(4)+1

  del_ne(del_tri(1),count_ne(del_tri(1)):count_ne(del_tri(1))+2) = del_tri(q1) 
  count_ne(del_tri(1)) = count_ne(del_tri(1))+3 

!and so on for 2,3,4
  del_ne(del_tri(2),count_ne(del_tri(2)):count_ne(del_tri(2))+2) = del_tri(q2) 
  count_ne(del_tri(2)) = count_ne(del_tri(2))+3 

  del_ne(del_tri(3),count_ne(del_tri(3)):count_ne(del_tri(3))+2) = del_tri(q3) 
  count_ne(del_tri(3)) = count_ne(del_tri(3))+3 

  del_ne(del_tri(4),count_ne(del_tri(4)):count_ne(del_tri(4))+2) = del_tri(q4) 
  count_ne(del_tri(4)) = count_ne(del_tri(4))+3 

end do

close(123)
!go through all the nodes, find unique neighbours.
max_ne = 0
do i = 1,node_num

!quick check if max_ne_temp was big enough 
  if (count_ne(i)-1.gt.max_ne_temp) then
    write(*,*)'count_ne=',count_ne(i)-1,'while max_ne_temp=',max_ne_temp
    write(*,*)'adjust max_ne_temp parameter! Make it bigger.'
    stop
  endif


  if (allocated(un_ne)) deallocate (un_ne)
  un_ne = unique(del_ne(i,1:count_ne(i)-1)) 
  del_num_ne(i) = size(un_ne)
  del_ne(i,:) = 0
  del_ne(i,1:del_num_ne(i)) = un_ne 
  if (max_ne.lt.del_num_ne(i)) max_ne = del_num_ne(i)
end do 

end subroutine

subroutine check_neigh_fast(X_phot,dir,X_node,neigh_num,neigh_list,X_neigh,&
&star_mode,R_star,X_star,delta_x,next)
!subroutine to compute the distance to the next neighbour 
implicit none
integer neigh_num,star_mode,i,comp,next

integer, dimension(neigh_num) :: neigh_list

double precision delta_x,new_x,R_star
double precision q1,q2,q4,q5 

double precision, dimension(3) :: X_phot,X_node,dir,X_star
double precision, dimension(3,neigh_num) :: X_neigh

!write(*,*)X_phot
!first we compute the distance until the end of the simulation
call pos_quadratic(norm2(dir)**2.d0,2.d0*dot_product(X_phot,dir),norm2(X_phot)**2.d0-1.d0,delta_x,comp)
!write(*,*)delta_x
next = 0        !if next = 0, then end of ray-path

q4 = norm2(X_node(:)-X_phot)**2.d0
q5 = dot_product(X_node(:)-X_phot(:),dir(:))
 
do i = 1,neigh_num
  q1 = norm2(X_neigh(:,i)-X_phot(:))**2.d0 
  q2 = dot_product(X_neigh(:,i)-X_phot(:),dir(:))
 
!  write(*,*)X_neigh(:,i)

  new_x = 0.5d0 * (q1-q4)/(q2-q5)

  if (new_x.gt.delta_x.or.new_x.le.1E-10) then
    cycle
  else
    delta_x = new_x
    next = neigh_list(i)
  endif
end do 

if (star_mode.eq.1) then
  call pos_quadratic(norm2(dir)**2.d0,2.d0*dot_product(X_phot(:)-X_star(:),dir),&
&norm2(X_phot(:)-X_star(:))**2.d0-R_star**2.0,new_x,comp)
  if (comp.ne.1) then
    if (new_x.lt.delta_x.and.new_x.gt.0.d0) then
!      write(*,*)X_phot
!      write(*,*)dir
!      write(*,*)new_x
!      stop
      delta_x = new_x
      next = 0
    endif
  endif
endif

X_phot = X_phot + delta_x * dir
!write(*,*)X_phot
!
!stop

end subroutine

subroutine pos_quadratic(a,b,c,x1,comp)
!subroutine that returns the minimal and positive solution to the quadratic equation
implicit none
integer comp

double precision x1,x2,a,b,c
double precision dis

dis = b**2.d0 - 4.d0 * a * c

if (dis > 0 ) then
  comp = 0
  x1 = ( -b + dsqrt(dis) ) / (2.d0 * a) 
  x2 = ( -b - dsqrt(dis) ) / (2.d0 * a) 

  if (x2.gt.0) x1 = x2
else
  comp = 1
endif

end subroutine

subroutine legendre_all(nth,nph,n_sphere,xth)
implicit none 
!subroutine to set up the integration of the sphere,
!using a legendre-discretization for theta and using
!a variable number of phi-points (weighed with sin(th)
integer nth,nph,n_sphere
integer i,np,j

double precision pi,ww,w_ph

double precision, dimension (3,nth*nph) :: xth
double precision, dimension (nth) :: x_th
double precision, dimension (nth) :: w_th

pi = dacos(-1.d0)

!xth --- the matrix with the integration parameters  
call cdgqf ( nth, 1, 1.d0, 1.d0, x_th, w_th )

n_sphere = 0
do i = 1,nth

!variable number of phi-points
  np = ceiling(nph*dsin(dacos(x_th(i))))

  w_ph = 2.d0 * pi / np         !maybe also add a random phase-factor for the
                                !series 
  ww = w_ph * w_th(i)
  do j = 1,np
    n_sphere = n_sphere+1

    xth(1,n_sphere) = x_th(i)
    xth(2,n_sphere) = (j-1)*w_ph 
    xth(3,n_sphere) = ww / (4.d0*pi) 
  end do
end do

end subroutine

subroutine star_all(nth,nph,kvec,gvec,bvec,dOm_star,n_star,n_sphere,xth)
implicit none
!subroutine to set up the integration of the sphere,
!using a legendre-discretization for theta and using
!a variable number of phi-points (weighed with sin(th)
integer nth,nph,n_sphere
integer i,np,j,n_star

double precision dth_star,th_star,xx,pi,dOm_star,ww,wph
double precision dth_rest,th_i,th,wth,th_l,th_u,dOm_max

double precision, dimension (3,nth*nph) :: xth
double precision, dimension (3) :: kvec,gvec,bvec,dir 

!R_star -- rotation matrix that gives the rotation-matrix 
!k_star = R_star * z 
!where b_0 is the magnetic-field direction along the z-axis 

pi = dacos(-1.d0)

!first, we'll work from a coordinate system where z-axis is along the 
!star-direction. That means that \int d\phi = 2*pi
dth_star = dOm_star / (2.d0*pi)
th_star  = acos(1-dth_star) 

!how many rays do we want to use for the star integration
dOm_max = 4.d0 * pi / (nth * nph)
n_star = ceiling(dOm_star / dOm_max)

!pick a good distribution of rays on the star-circle
call circle_points(n_star,th_star,dth_star,xth(:,1:n_star))

do i = 1,n_star
  xth(3,i) = xth(3,i) / (4.d0 * pi)
!  write(*,*)xth(3,i)
end do

!now move on to the rest of the simulation

n_sphere = n_star

dth_rest = (pi - th_star)/nth
th_i = th_star + dth_rest/2.d0
do i = 1,nth
  th = (i-1)*dth_rest + th_i 
  if (i.eq.1) then 
    th_u = (i-0.5d0)*dth_rest + th_i 
    wth  = dcos(th_star) - dcos(th_u) 
  elseif (i.eq.nth) then
    th_l = (i-1.5d0)*dth_rest + th_i 
    wth  = dcos(th_l) + 1.d0    !1.d0 = dcos(pi)
  else
    th_u = (i-0.5d0)*dth_rest + th_i
    th_l = (i-1.5d0)*dth_rest + th_i
    wth  = dcos(th_l) - dcos(th_u)
  endif

  !now that we have the weight, evaluate the number of phi-s
  np = ceiling(nph * sin(th))
  wph = 2.d0 * pi / np
  xx = dcos(th)
  ww = wth * wph 
  do j = 1,np
    n_sphere = n_sphere + 1
    if (n_sphere.gt.nph*nth) then
      write(*,*)'error, star integration paramaters do not compute'
      stop
    endif

    xth(1,n_sphere) = xx 
    xth(2,n_sphere) = (j-1)*wph
    xth(3,n_sphere) = wth * wph / (4.d0 * pi)
  end do
end do

!Now, we will express the angles as relative to the magnetic field angles
!by first a conversion to a direction and then convert that direction into
!the magnetic field coordiates  

do i = 1,n_sphere
  call find_dir(kvec,gvec,dacos(xth(1,i)),xth(2,i),dir)
 
  call find_angles(bvec,gvec,dir,xth(1,i),xth(2,i)) 
  xth(1,i) = dcos(xth(1,i))

end do

end subroutine

subroutine circle_points(n,th_s,dth_s,x)
!subroutine that divides a number of points roughly evenly
!on a circle.

implicit none
integer i,n

double precision xx,ga,ww,wph,th_s,dth_s
double precision, parameter :: pi = dacos(-1.d0)

double precision, dimension (3,n) :: x

if (n.eq.1) then
!one point, only centre
  x(1:2,1) = dcos(1.d0)
  x(3,1) = 2.0*pi*dth_s
elseif (n.lt.10) then
!less then ten points, put in circle
  wph = 2.d0 * pi / n
  ww  = wph * dth_s
  xx  = dcos(th_s/2.d0)
  do i = 1,n
    x(1,i) = xx 
    x(2,i) = (i-1)*wph
    x(3,i) = ww
  end do
else
!more than ten points, use Vogel's method
  ga = pi * (3.d0 - dsqrt(5.d0))        !golden angle
  ww = dth_s * 2.d0 * pi / n 
  do i = 1,n
    x(1,i) = dcos(th_s*dsqrt(1.d0*i/n))
    x(2,i) = i * ga
    x(3,i) = ww
  end do
endif

end subroutine

