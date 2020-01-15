subroutine eat_lime_dim(in_grid,in_mol,node_num,del_num,e_num,r_num,c_num,c_temp,c_par)
!program to read in the dimensions of the lime-programs. 
!dimensions will be used to allocate the arrays.
!first step in reading in the data.
implicit none

character kasta
character(len=100) :: in_grid,in_mol

integer node_num,del_num
integer i,e_num,r_num,c_num,c_temp,c_par

!output:
!node_num	-- total amount of nodes	
!del_num	-- total amount of delaunay trinagles
!e_num		-- total number of energy levels
!r_num		-- total number of radiative transitions
!c_num		-- total number of collisional transitions
!c_temp		-- number of collisional temperatures

open(unit=11,file=in_grid)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)kasta,node_num

!write(*,*)node_num

do i = 1,node_num
  read(11,*)
end do
read(11,*)
read(11,*)kasta,del_num

!write(*,*)del_num

open(unit=12,file=in_mol)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)e_num
read(12,*)
do i = 1,e_num
  read(12,*)
end do
read(12,*)
read(12,*)r_num
read(12,*)
do i = 1,r_num
  read(12,*)
end do
read(12,*)
read(12,*)c_par
read(12,*)
read(12,*)
read(12,*)
read(12,*)c_num
read(12,*)
read(12,*)c_temp

!write(*,*)e_num,r_num,c_num,c_temp
close(11)
close(12)


end subroutine 

subroutine eat_lime_no_delaunay(in_pop,in_grid,in_mol,in_dust,node_num,del_num,e_num,r_num,c_num,c_temp,&
&XYZ,vel_all,temp_all,lev_p,e_levs,t_Aji,nu,t_Cji,col_temps,c_par,&
&h2_dens,b_all,ab_all,j_lev,r_trans,c_trans,max_ne_temp,maxne,del_num_ne,del_ne,rmax&
&,bulk_node_num,kd_v,dust_ii_gas,dust_p,mass_mol)
!program to read in the data of the lime-programs. 
!dimensions read in before and are used to allocate the arrays.
implicit none

character(len=100) :: in_pop,in_grid,in_mol,in_dust

integer node_num,del_num,e_num,r_num,c_num,c_temp
integer i,j,k,maxne,bulk_node_num,c_par
integer b1,b2,max_ne_temp,c_num_2
integer, parameter :: n_dt=67 

integer, dimension(e_num) :: j_lev 
integer, dimension(2,r_num) :: r_trans 
integer, dimension(2,c_num) :: c_trans 
integer, dimension(node_num) :: del_num_ne 
integer, dimension(node_num,max_ne_temp) :: del_ne 

double precision g1,s_p
double precision a1,a2,a3 
double precision r1,rmax
double precision dust_ii_gas,mass_mol
double precision b_therm

double precision, dimension(3,node_num) :: XYZ
double precision, dimension(3,node_num) :: vel_all
double precision, dimension(2,node_num) :: temp_all
double precision, dimension(e_num,node_num) :: lev_p 
double precision, dimension(node_num) :: dust_p 
double precision, dimension(e_num) :: e_levs 
double precision, dimension(r_num) :: t_Aji,nu,kd_v 
double precision, dimension(2,n_dt) :: dust_table 
double precision, dimension(c_temp,c_num,c_par) :: t_Cji
double precision, dimension(c_temp) :: col_temps 
double precision, dimension(node_num) :: h2_dens
double precision, dimension(node_num) :: b_all
double precision, dimension(node_num) :: ab_all
double precision, dimension(node_num) :: dust_ii_vec 

!input:
!node_num	-- total amount of nodes	
!del_num	-- total amount of delaunay trinagles
!e_num		-- total number of energy levels
!r_num		-- total number of radiative transitions
!c_num		-- total number of collisional transitions
!c_temp		-- number of collisional temperatures

!output:
!

open(unit=11,file=in_grid)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)


do i = 1,node_num
  read(11,*)XYZ(1,i),XYZ(2,i),XYZ(3,i)
end do
read(11,*)
read(11,*)

!do i = 1,del_num
!!  read(11,*)j,del_tri(1,i),del_tri(2,i),del_tri(3,i),del_tri(4,i)
!  read(11,*)del_tri(1,i),del_tri(2,i),del_tri(3,i),del_tri(4,i)
!end do

!0-base to 1-base
!do i = 1,del_num
!  del_tri(1,i) = del_tri(1,i)+1
!  del_tri(2,i) = del_tri(2,i)+1
!  del_tri(3,i) = del_tri(3,i)+1
!  del_tri(4,i) = del_tri(4,i)+1
!end do

open(unit=12,file=in_mol)
read(12,*)
read(12,*)
read(12,*)
read(12,*)mass_mol
read(12,*)
read(12,*)
read(12,*)
do i = 1,e_num
  read(12,*)j,e_levs(i),g1
  j_lev(i) = nint((g1-1.d0)/2.d0)
end do
read(12,*)
read(12,*)
read(12,*)
do i = 1,r_num
  read(12,*)j,r_trans(2,i),r_trans(1,i),t_Aji(i),nu(i)
end do
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)
read(12,*)col_temps(1:c_temp)
read(12,*)
do i = 1,c_num
  read(12,*)j,c_trans(2,i),c_trans(1,i),t_Cji(1:c_temp,i,1)
end do
if (c_par.eq.2) then
  read(12,*)
  read(12,*)
  read(12,*)
  read(12,*)
  read(12,*)
  read(12,*)c_num_2
  read(12,*)
  read(12,*)
  read(12,*)
  do i = 1,c_num_2
    read(12,*)j,b1,b2,t_Cji(1:c_temp,i,2)
  end do
endif

close(11)
close(12)

!write(*,*)t_Cji(1:3,1,1)
!write(*,*)t_Cji(1:3,1,2)
!write(*,*)

dust_ii_vec(:) = 0.d0

open(unit=13,file=in_pop)
!read(13,*)
i = 0

100 CONTINUE
  i = i+1
  read(unit=13,FMT=*,END=200)a1,a2,a3,vel_all(1:3,i),h2_dens(i),temp_all(1,i),temp_all(2,i),b_all(i),&
&ab_all(i),dust_ii_vec(i),g1,lev_p(1:e_num,i)
GO TO 100

200 CONTINUE

close(13)

bulk_node_num = i

!finish with the surface nodes
do j = i+1,node_num
  vel_all(:,j)     = 0.d0
  h2_dens(j)       = 1.d0
  temp_all(:,j)    = 0.d0
  b_all(j)         = 200.d0  
  ab_all(j)        = 1D-3  

  lev_p(1,j) = 1.d0 
  s_p = 1.d0
  do k = 2,e_num
    lev_p(k,j) = 1.d0*dexp(-e_levs(k)*0.527)    !0.527 = 1/cm * c * h / (k*T_CMB)
    s_p = s_p + lev_p(k,j)
  end do
  do k = 1,e_num
    lev_p(k,j) = lev_p(k,j) / s_p
!    write(*,*)lev_p(k,j)
  end do
!  stop
end do

!we might have to add the temperature to the doppler parameter, as it is 
!composed only from the turbulent motion as of yet.

do j = 1,node_num
  b_therm = 2.d0 * 8.314462d0 * temp_all(1,j) / mass_mol
!we use gas-constant because mass_mol is given in Daltons

  b_all(j) = dsqrt(b_all(j)**2.d0 + b_therm)    
!add the thermal broadening to the total broadening and divide
!by two because we want to use the standard-deviation instead
!of the b-parameter
end do

!now the molecular levels times the abundance
do j = 1,node_num
  lev_p(:,j) = lev_p(:,j) * ab_all(j) * h2_dens(j) 
end do
!we're gonna perform the Delaunay triangulation and construct a neighbor list

!first, let us normalize the XYZ coordinates

!find maximum r
rmax = 0.d0
do i = 1,node_num
  r1 = norm2(XYZ(:,i))
  if (r1.gt.rmax) rmax = r1
!  write(*,*)r1
end do  

!now normalize
do i =1,node_num
  XYZ(:,i) = XYZ(:,i) / rmax
end do


call adj_list_no_delaunay(in_grid,del_num,del_ne,del_num_ne,max_ne_temp,maxne,node_num)

!finally, read in the dust_opacity input
open(unit=11,file=in_dust)
do i = 1,n_dt
  read(11,*)dust_table(:,i)
end do

do i = 1,n_dt
  dust_table(1,i) = dlog10(dust_table(1,i))
  dust_table(2,i) = dlog10(dust_table(2,i))
end do

!we can directly convert it to the appropriate dust-opacities at the transition
!frequencies of interest

call dust_to_opacities(nu*1E9,r_num,dust_table,n_dt,kd_v)  


!now we compute the dust density in kg/m3 per node
do i = 1,node_num
  !2.32E-27 = 1.4 Daltons
  !...      = 2.6 Daltons --- prob better for molecules
!  dust_p(i) = h2_dens(i) * 2.32E-27 * dust_ii_vec(i)
  if (dust_ii_vec(i).ne.0.d0) then
    dust_p(i) = h2_dens(i) * 2.32E-27 / dust_ii_vec(i)    !we are reading in gas_ii_dust
  else
    dust_p(i) = 0.d0 
  endif
end do

end subroutine 

