subroutine ex_elastic(C,lC,nC,nT,col_temps,j_lev,e_levs,e_num,C0)
!subroutine to obtain the elastic collision rates
implicit none
integer nC,nT,nP,e_num
integer dj,d_j,max_dj,min_dj
integer i,j,i1,i2,j1,j2

integer, dimension(2,nC) :: lC
integer, dimension (e_num) :: j_lev 
integer, dimension (nC) :: dj_nums

double precision dE,E1,E2,g1,g2,x,y,Temp,C0

double precision, dimension(nT) :: col_temps
double precision, dimension(nT,nC) :: C
double precision, dimension(e_num) :: e_levs
double precision, dimension(2,2) :: B
double precision, dimension(2) :: a
double precision, allocatable, dimension(:,:) :: const 

!based on Keto & Rybicki 2005
!and de Jong et al. 1970

min_dj = 100
max_dj = 0
do i = 1,nC
  dj = lC(2,i) - lC(1,i)
  if (dj.gt.max_dj) max_dj = dj
  if (dj.lt.min_dj) min_dj = dj
end do 

allocate(const(2,min_dj:max_dj))
!now were gonna find the a(\Delta J) and b(\Delta J) parameters
do dj = min_dj,max_dj
  B(:,:) = 0.d0
  a(:)   = 0.d0

!  write(*,*)dj
  do i = 1,nC
    i1 = lC(1,i)
    i2 = lC(2,i)

    j1 = j_lev(i1)
    j2 = j_lev(i2)

    g1 = 1.d0 * (2.d0 * j1 + 1.d0)
    g2 = 1.d0 * (2.d0 * j2 + 1.d0)

    E1 = e_levs(i1)
    E2 = e_levs(i2)

    d_j = j2 - j1 
    if (dj.eq.d_j) then
      dE = E2 - E1 
      do j = 1,nT
        Temp = col_temps(j) 
        x = 4.779E-11 * dE/Temp   !h/kB * \Delta v / T 
        y = dlog(C(j,i) * ((g2/g1)*(1.d0+x))**(-1.d0))

        B(1,1) = B(1,1) + 1.d0 
        B(1,2) = B(1,2) - dsqrt(x) 
        B(2,2) = B(2,2) + x 
       
        a(1) = a(1) + y
        a(2) = a(2) - dsqrt(x) * y

!        write(*,*)dsqrt(x),y

      end do
    endif 
  end do

  B(2,1) = B(1,2)
 
  !now solve the 2x2 inversion problem
  call invert_22(B) 

  const(1,dj) = B(1,1)*a(1) + B(1,2)*a(2)       !ln(a(\Delta J)) 
  const(2,dj) = B(2,1)*a(1) + B(2,2)*a(2)       !b(\Delta J)
!  write(*,*)
  write(*,*)dj,dexp(const(1,dj)),const(2,dj)
!  write(*,*)
!  write(*,*)
end do

!now let's fit the a-constant 

B(:,:) = 0.d0
do dj = min_dj,max_dj
  B(1,1) = B(1,1) + 1.d0
  B(1,2) = B(1,2) + 1.d0 * dj
  B(2,2) = B(2,2) + 1.d0 * dj*dj

  a(1) = a(1) + 1.d0 * const(1,dj)
  a(2) = a(2) + dj * const(1,dj)
end do

call invert_22(B)

C0 = B(1,1) * a(1) + B(1,2) * a(2)
C0 = dexp(C0)

!write(*,*)C0
end subroutine

subroutine invert_22(A)
!subroutine to invert a 2x2 matrix
implicit none

double precision det,a11

double precision, dimension(2,2) :: A

det = (A(1,1) * A(2,2) - A(2,1) * A(1,2))**(-1.d0)

a11 = A(1,1)

A(1,1) = det * A(2,2)
A(1,2) = -det * A(1,2)  
A(2,1) = -det * A(2,1)  
A(2,2) = det * a11 

end subroutine

subroutine ex_elastic_singleJ(J_state,C,lC,nC,nT,col_temps,j_lev,e_levs,e_num,C0)
!subroutine to obtain the elastic collision rates for a single level J by 
!only considering the temperature variation
implicit none
integer nC,nC_red,nT,nP,e_num
integer dj,d_j,max_dj,min_dj
integer i,j,i1,i2,j1,j2,J_state

integer, dimension(2,nC) :: lC
integer, dimension(2,nC) :: lC_red
integer, dimension(nC) :: red_to_real 
integer, dimension (e_num) :: j_lev 
integer, dimension (nC) :: dj_nums

double precision dE,E1,E2,g1,g2,x,y,Temp,C0

double precision, dimension(nT) :: col_temps
double precision, dimension(nT,nC) :: C
double precision, dimension(e_num) :: e_levs
double precision, dimension(2,2) :: B
double precision, dimension(2) :: a
double precision, allocatable, dimension(:,:) :: const 

!based on Keto & Rybicki 2005
!and de Jong et al. 1970

min_dj = 100
max_dj = 0

nC_red = 0
do i = 1,nC
  if (lC(1,i).eq.J_state) then
    nC_red = nC_red + 1
    dj = lC(2,i) - lC(1,i)
 
    red_to_real(nC_red) = i

    lC_red(1,nC_red) = lC(1,i)
    lC_red(2,nC_red) = lC(2,i)

    if (dj.gt.max_dj) max_dj = dj
    if (dj.lt.min_dj) min_dj = dj
  endif
end do 

allocate(const(2,min_dj:max_dj))
!now were gonna find the a(\Delta J) and b(\Delta J) parameters
do dj = min_dj,max_dj
  B(:,:) = 0.d0
  a(:)   = 0.d0

  write(*,*)dj
  do i = 1,nC_red
    i1 = lC_red(1,i)
    i2 = lC_red(2,i)

    j1 = j_lev(i1)
    j2 = j_lev(i2)

    g1 = 1.d0 * (2.d0 * j1 + 1.d0)
    g2 = 1.d0 * (2.d0 * j2 + 1.d0)

    E1 = e_levs(i1)
    E2 = e_levs(i2)

    d_j = j2 - j1 
    if (dj.eq.d_j) then
      dE = E2 - E1 
      do j = 1,nT
        Temp = col_temps(j) 
        x = 4.779E-11 * dE/Temp   !h/kB * \Delta v / T 
        y = dlog(C(j,red_to_real(i)) * (g2/g1) * ((1.d0+x))**(-1.d0))

        B(1,1) = B(1,1) + 1.d0 
        B(1,2) = B(1,2) - dsqrt(x) 
        B(2,2) = B(2,2) + x 
       
        a(1) = a(1) + y
        a(2) = a(2) - dsqrt(x) * y

        write(*,*)dsqrt(x),y

      end do
    endif 
  end do

  B(2,1) = B(1,2)
 
  !now solve the 2x2 inversion problem
  call invert_22(B) 

  const(1,dj) = B(1,1)*a(1) + B(1,2)*a(2)       !ln(a(\Delta J)) 
  const(2,dj) = B(2,1)*a(1) + B(2,2)*a(2)       !b(\Delta J)
!  write(*,*)
!  write(*,*)const(1,dj),const(2,dj)
!  write(*,*)
!  write(*,*)
end do

!now let's fit the a-constant 

B(:,:) = 0.d0
do dj = min_dj,max_dj
  B(1,1) = B(1,1) + 1.d0
  B(1,2) = B(1,2) + 1.d0 * dj
  B(2,2) = B(2,2) + 1.d0 * dj*dj 

  a(1) = a(1) + 1.d0 * dlog(const(1,dj))
  a(2) = a(2) + dj * dlog(const(1,dj)) 
end do

call invert_22(B)

C0 = B(1,1) * a(1) + B(1,2) * a(2)

end subroutine

subroutine de_polarization(C,lC,nC,nT,e_num,kmax,D0)
!subroutine to obtain the de-polarization rates in the sudden
!approximation 
use fwigxjpf
implicit none
integer nC,nT,e_num
integer i,j,k,l,m,kmax

integer, dimension(2,nC) :: lC
integer, dimension (e_num) :: j_lev 
integer, dimension (nC) :: dj_nums

double precision lj2,fac,w1,w2,w3,fac2

double precision, dimension(nT,nC) :: C
double precision, dimension(nT,e_num,2:kmax) :: D0 
double precision, dimension(nT,e_num) :: C_j0 

!input:
!C      -- collision rates
!lC     -- collisional transitions
!nC     -- number of collisional transitions
!nT     -- number of collision temperatures tabulated
!e_num  -- number of levels
!kmax   -- maximum number of irreducible tensor-rank 

!output:
!D0     -- all depolarization coefficients, temperature dependent       

!based on Morris 1985 -- Eq. (B11)
!same definition of D^(K)(\alpha J) as was 
!used in L&L06, Eq. (7.102) 

do i = 1,nC
  if (lC(1,i).eq.1) C_j0(:,lC(2,i)) = C(:,i)  
end do
!if lC(1,i) is the zero-level, then use these as the sudden-approximation
!paramaters 
!C_j0(:,1) = 0 -- we don't need it

!write(*,*)C_j0(1,2)
!write(*,*)C_j0(1,3)

!only works if j_lev(i) = i-1
D0(:,:,:) = 0.d0
do i = 2,e_num
  j = i-1
  lj2 = (2.d0*j+1.d0)**2.d0
  do k = 2,min(2*j,kmax),2
    do l = 1,nT
      D0(l,i,k) = 0.d0
      !Wigner 3j = 0 if J+J < j
      do m = 1,2*j     
  
        fac = (-1.d0)**(m) * (2.d0*m+1.d0)

        w1 = fwig3jj (2* j , 2* j    , 2* m ,&
        &              2* 0 , 2* 0 , 0    )

        w2 = fwig6jj (2* j , 2* j , 2* m ,&
        &              2* j , 2* j  , 2* 0    )

        w3 = fwig6jj (2* j , 2* j , 2* m ,&
        &              2* j , 2* j  , 2* k    )

        fac2 = w2 - w3

        D0(l,i,k) = D0(l,i,k) + lj2 * fac * w1 * w1 * fac2 * C_j0(l,m)
!        write(*,*)m,k
!        write(*,*)lj2 * fac * w1 * w1 * fac2,C_j0(l,m+1)
      end do
    end do
  end do
end do

end subroutine

subroutine make_depolarization_temp(D0_all,nlev,nT,kmax,col_temps,Temp,nH2,D0)
!interpolate the temperature-dependent de-polarization coefficients
implicit none
integer nlev,nT,kmax,i

double precision Temp,ddT,nH2

double precision, dimension(nT,nlev,2:kmax) :: D0_all
double precision, dimension(nlev,2:kmax) :: D0
double precision, dimension(nlev,2:kmax) :: f1,f2 
double precision, dimension(nT) :: col_temps

!input:
!D0_all         -- all depolarization coefficients, temperature dependent
!nlev           -- number of levels
!nT             -- number of temperatures
!kmax           -- maximum irreducible tensor rank
!col_temps      -- tabulated collision temperatures
!Temp           -- temperature we want coefficients for
!nH2            -- H2 number density

!output
!D0             -- depolarization coefficients for specific temperature


if (Temp.le.col_temps(1)) D0(:,:) = D0_all(1,:,:)

do i = 2,nT
  if (col_temps(i).gt.Temp) then
    ddT = (Temp - col_temps(i-1)) / (col_temps(i) - col_temps(i-1))

    f1(:,:) = D0_all(i-1,:,:)
    f2(:,:) = D0_all(i,:,:)

    D0(:,:) = f1(:,:) + ddT*(f2(:,:) - f1(:,:))
  endif
end do

D0(:,:) = D0(:,:) * nH2 * 1d-6

end subroutine

