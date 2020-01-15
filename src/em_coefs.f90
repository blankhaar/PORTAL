subroutine comp_emcoefs(r1_0,r2_0,r1_2,r2_2,j1,j2,th,nu0,Aji,nsI_tot,nsQ_tot,es_I,es_Q)
use fwigxjpf
!subroutine to compute the emission coefficients
!as of now, only implementation for the ratio nq/ni and eq/ei
!theta angle is required
!tested on a matlab routine

implicit none

integer j1,j2

double precision th,tki,tkq
double precision wjj2_1,wjj2_2
double precision fas,fac1,fac2,w6_1,w6_2
double precision r1_0,r2_0,r1_2,r2_2
double precision g1,g2,nu0,b_to_a,a_to_b
double precision nsI_tot,nsQ_tot
double precision na_I,na_Q,ns_I,ns_Q
double precision es_I,es_Q
double precision pQ,eQ,Aji,Bji,hv_fac 


!input
!r1_0,2		-- lower level rank0,2 population
!r2_0,2		-- upper level rank0,2 population
!j1/j2		-- lower/upper level angular momentum
!th		-- propagation-magnetic field angle
!nu0		-- transition frequency

!output
!oQ		-- \eta_Q/\eta_I -factor
!eQ		-- \epsilon_Q / \epsilon_I - factor 

!compute the w_{jj'}^k-factors
!only required for k=0 and k=2
!Eq.(38) Landi-84

fas = (-1.d0)**(1.d0 + j1 + j2)
fac1 = dsqrt(3.d0 * (2.d0 * j1 + 1.d0))
fac2 = dsqrt(3.d0 * (2.d0 * j2 + 1.d0))

w6_1 = fwig6jj(2*  1 , 2*  1 , 2*  2 ,&
     &         2* j1,  2*  j1, 2*  j2) 

w6_2 = fwig6jj(2*  1 , 2*  1 , 2*  2 ,&
     &         2* j2,  2*  j2, 2*  j1) 

wjj2_1 = fas * fac1 * w6_1
wjj2_2 = fas * fac2 * w6_2

!degeneracy
g1 = 2.d0 * j1 + 1.d0
g2 = 2.d0 * j2 + 1.d0

!
!T^K_Q(i,\Omega) rotation elements from Bommier 1997b
tki = (3.d0 * dcos(th)**2.d0 - 1.d0 )/(2.d0 * dsqrt(2.d0))
tkq = -3.d0 * dsin(th)**2.d0 / (2.d0 * dsqrt(2.d0)) 

!write(*,*)
!write(*,*)tki,tkq
!write(*,*)

!Eqs.(40)-(42) Landi-84
na_I = 1.d0 + wjj2_1 * (r1_2/r1_0) * tki 
na_Q = wjj2_1 * (r1_2/r1_0) * tkq

ns_I = 1.d0 + wjj2_2 * (r2_2/r2_0) * tki 
ns_Q = wjj2_2 * (r2_2/r2_0) * tkq 

a_to_b = 6.781962D49    !best to use L&L relation here
                        !B21 = [c^2 / (2*h * v^3)] * A21
                        !a_to_b = c^2 / (2*h)

!write(*,*)na_I,na_Q
!write(*,*)ns_I,ns_Q


b_to_a = nu0**3.d0 / a_to_b    

Bji = Aji / b_to_a
!write(*,*)
!write(*,*)Aji,Bji

es_I = ns_I 
es_Q = ns_Q 

!now properly substract the absorption and stimulated emission 
!so that the proportionalities stay in-tact

!na ~ B_12 * N1
!ns ~ B_21 * N2 = B_12 * g1 * N2 / g2
!N1 ~ sqrt(g1) * r1_0 
!N2 ~ sqrt(g2) * r2_0 		--- Eq.(43) Landi-84 

nsI_tot = dsqrt(g1)*r1_0 * na_I - dsqrt(g2)*r2_0*ns_I*g1/g2
nsQ_tot = dsqrt(g1)*r1_0 * na_Q - dsqrt(g2)*r2_0*ns_Q*g1/g2

hv_fac = 5.273D-35 * nu0            !(h/4*pi) * nu0

!write(*,*)
!write(*,*)hv_fac

es_I = hv_fac * Aji * es_I * dsqrt(g2) * r2_0 
es_Q = hv_fac * Aji * es_Q * dsqrt(g2) * r2_0 

nsI_tot = hv_fac * (g2/g1)*Bji * nsI_tot        !g2/g1 * B21 = B12
nsQ_tot = hv_fac * (g2/g1)*Bji * nsQ_tot

!write(*,*)nsI_tot,nsQ_tot
!write(*,*)es_I,es_Q
!stop
end subroutine

subroutine get_rhovec(c_dim,count_array,nlev,j_lev,ntrans,trans,AC_trans,t_int,kmax,Temp,D0,rhovec,extra_info,rates)
!subroutine that obtains the rho-vector of a certain node 
!by setting up the SEE and subsequently solving them
implicit none
integer nlev,ntrans,t_int,c_dim
integer i,j,k,jc,ic,INFO,n
integer kmax,extra_info

integer, dimension(nlev) :: j_lev,count_array
integer, dimension(2,ntrans) :: trans
integer, dimension(c_dim) :: IPIV

double precision Temp

double precision, dimension(5,ntrans) :: AC_trans
double precision, dimension(nlev,2:kmax) :: D0 
double precision, dimension(c_dim,c_dim) :: Amat
double precision, dimension(c_dim+1) :: yvec 
double precision, dimension(c_dim+1,c_dim) :: Bmat
double precision, dimension(c_dim+1) :: dvec
double precision, dimension(c_dim) :: WORK,rhovec
double precision, dimension(c_dim,c_dim) :: C_mat
double precision, dimension(c_dim,c_dim) :: B_mat
double precision, dimension(c_dim,c_dim) :: R_mat
double precision, dimension(3,c_dim) :: rates 

!input:
!c_dim          -- dimension of the SEE matrix
!count_array    -- bookkeeping array for the (polarized) levels
!nlev           -- number of levels
!j_lev          -- associated J-angular momentum of the levels
!ntrans         -- number of transitions
!trans          -- associated levels with these transitions
!AC_trans       -- associated rates with these transtions
!t_int          --  

call setup_stateq(c_dim,count_array,nlev,j_lev,trans,ntrans,kmax,&
&AC_trans,D0,Temp,Amat,R_mat,B_mat,C_mat)

!like in Hazel, replace the first row with the level constraint
!that is set at 1d7 --- densities we are interested in are 
!probably lower, so we'll set 1.d3

Amat(1,:) = 0.d0
do i = 1,nlev
  ic = count_array(i) + 1
  jc = j_lev(i)

  Amat(1,ic) = dsqrt(2.d0 * jc + 1.d0)
end do 

!Amat * rhovec = dvec
!dvec(1) = 1.d3 and rest = 0.d0

call DGETRF(c_dim,c_dim,Amat,c_dim,IPIV,INFO)
if (INFO.ne.0) then
  write(*,*)'Density-matrix calculation: failed'
  stop
else
  call DGETRI(c_dim,Amat,c_dim,IPIV,WORK,c_dim,INFO)
  if (INFO.ne.0) then
    write(*,*)'Density-matrix calculation: failed'
    stop
  endif
endif

!see beneath to get all the populations, but we are actually 
!interested only in the populations for our levels
!dvec(1) = 1.d3 and rest = 0.d0
do i = 1,c_dim
  rhovec(i) = Amat(i,1) 
end do

if (extra_info.eq.1) then
  rates(:,:) = 0.d0
  do i = 1,c_dim
    do j = 1,c_dim
      rates(1,i) = rates(1,i) + dabs(R_mat(i,j))*rhovec(j) 
      rates(2,i) = rates(2,i) + dabs(B_mat(i,j))*rhovec(j) 
      rates(3,i) = rates(3,i) + dabs(C_mat(i,j))*rhovec(j) 
!      write(*,*)j,dabs(C_mat(i,j)),rhovec(j)
    end do
!    stop
  end do
endif

end subroutine

