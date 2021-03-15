program main 
use fwigxjpf

!Main file containg "program main" that calls the subroutines to read in LIME
!output, and call the inverse ray-tracing to compute the alignment of all of 
!the energy-levels at all nodes in the simulation

implicit none

character(len=100) :: input_file,output_map
character(len=100) :: c_pol_pops,c_spont_rates,c_rad_rates,c_col_rates,c_anis 
character(len=100) :: in_pop,in_grid,in_mol,in_dust 
character(len=100) :: c_mag

integer node_num,del_num,e_num,r_num,c_num,c_temp
integer i,j,k,i1,maxne,node,node_count,nth,nph,g1,g2,n_trans
integer c_dim,t_int,kmax,im_pix,bulk_node_num
integer verbose_mode,star_mod,k_mu,k_phi 
integer nlevels,nincs,extra_info,c_par,de_pol
integer nvel,n_sphere
integer max_ne_temp,n_depth,n_depth_rt

integer, dimension(23) :: nodes

integer, allocatable, dimension(:) :: j_lev,lev_g
integer, allocatable, dimension(:,:) :: r_trans
integer, allocatable, dimension(:,:) :: c_trans
integer, allocatable, dimension(:,:) :: adj,adj0
integer, allocatable, dimension(:) :: del_num_ne 
integer, allocatable, dimension(:,:) :: del_ne 
integer, allocatable, dimension(:,:) :: del_ne_temp 
integer, allocatable, dimension(:,:) :: all_trans 
integer, allocatable, dimension(:) :: count_array 

double precision x_dip,r1,r2,r3,dx,tx,R_star,T_star
double precision a_to_b,start,finish,I_bb,dust_ii_gas
double precision th,ph,im_size,nu0,d_tot
double precision dOm,th_star,ph_star,mass_mol,x1,x2,x3
double precision, parameter :: pi = dacos(-1.d0)

double precision, dimension(3) :: dir,bvec,gvec
double precision, dimension(3) :: bvec_trans,XYZ_star

double precision, allocatable, dimension(:,:) :: XYZ,vel_all
double precision, allocatable, dimension(:,:) :: temp_all
double precision, allocatable, dimension(:,:) :: lev_p
double precision, allocatable, dimension(:) :: e_levs
double precision, allocatable, dimension(:) :: t_Aji,nu,t_Bij,kd_v
double precision, allocatable, dimension(:,:,:) :: t_Cji
double precision, allocatable, dimension(:) :: dust_p
double precision, allocatable, dimension(:) :: col_temps
double precision, allocatable, dimension(:) :: h2_dens 
double precision, allocatable, dimension(:) :: b_all
double precision, allocatable, dimension(:) :: ab_all
double precision, allocatable, dimension(:,:) :: xth
double precision, allocatable, dimension(:,:) :: bvec_all
double precision, allocatable, dimension(:,:) :: pol_pops 
double precision, allocatable, dimension(:,:) :: JJ_20 
double precision, allocatable, dimension(:,:,:) :: pol_rates

verbose_mode=0
dust_ii_gas = 0.01d0

call getarg(1,input_file)
call getarg(2,output_map)

call getarg(3,in_pop)
call getarg(4,in_grid)
call getarg(5,in_mol)
call getarg(6,in_dust)

open(unit=101,file=input_file)

read(101,*)
read(101,*)star_mod,T_star,R_star,XYZ_star(:)
read(101,*)nth,nph,nvel
read(101,*)x1,x2,x3,x_dip
read(101,*)kmax,extra_info
read(101,*)max_ne_temp,n_depth
close(unit=101)

call fwig_table_init(2*100,9)
call fwig_temp_init(2*100)

!tested and works
call eat_lime_dim(in_grid,in_mol,node_num,del_num,e_num,r_num,c_num,c_temp,c_par)

allocate(j_lev(e_num))     
allocate(lev_g(e_num))     
allocate(r_trans(2,r_num))  
allocate(c_trans(2,c_num)) 

allocate(XYZ(3,node_num))     
allocate(vel_all(3,node_num))     
allocate(temp_all(2,node_num))       
allocate(lev_p(e_num,node_num)) 
allocate(e_levs(e_num))          
allocate(t_Aji(r_num))
allocate(t_Bij(r_num))
allocate(nu(r_num))          
allocate(kd_v(r_num))          
allocate(t_Cji(c_temp,c_num,c_par))   
allocate(col_temps(c_temp))         

allocate(h2_dens(node_num))         
allocate(dust_p(node_num))         
allocate(b_all(node_num))         
allocate(ab_all(node_num))         

allocate(del_num_ne(node_num))
allocate(del_ne_temp(node_num,max_ne_temp))

call eat_lime_no_delaunay(in_pop,in_grid,in_mol,in_dust,node_num,del_num,e_num,r_num,c_num,c_temp,&
&XYZ,vel_all,temp_all,lev_p,e_levs,t_Aji,nu,t_Cji,col_temps,c_par,&
&h2_dens,b_all,ab_all,j_lev,r_trans,c_trans,max_ne_temp,maxne,del_num_ne,del_ne_temp,tx,&
&bulk_node_num,kd_v,dust_ii_gas,dust_p,mass_mol)

!write(*,*)'iii'
nu = 1E9 * nu                  !convert GHz to Hz 
e_levs = 29.9792458E9 * e_levs !convert cm-1 to Hz

R_star = R_star / tx

allocate(del_ne(node_num,maxne))
del_ne(:,:) = 0
do i = 1,node_num
  del_ne(i,1:del_num_ne(i)) = del_ne_temp(i,1:del_num_ne(i))
end do

deallocate(del_ne_temp)
do i = 1,e_num
  lev_g(i) = 2*j_lev(i) + 1
end do

dx = 1E-6

allocate(xth(3,nth*nph))
call legendre_all(nth,nph,n_sphere,xth)
!write(*,*)nth,nph,n_sphere
!write(*,*)
!d_tot = 0.d0
!do i = 1,n_sphere
!  write(*,*)xth(:,i)
!  d_tot = d_tot + xth(3,i)
!end do
!write(*,*)
!write(*,*)d_tot
!write(*,*)
allocate(bvec_all(3,node_num))
!general magnetic field structure
if (x_dip.ne.0.d0)  then
  call dipole_magfield(node_num,XYZ,x1,x2,x3,bvec_all)
else
  call global_magfield(node_num,XYZ,x1,x2,x3,bvec_all)
endif
!write(*,*)bvec_all(:,1)

if (extra_info.eq.1) then
  write(c_mag,fmt="(A)")"/mag_field.dat"
  c_mag = trim(output_map) // trim(c_mag)

  open(unit=30,file=c_mag,recl=8824)
  do node = 1,node_num
    write(30,*)XYZ(:,node),bvec_all(:,node)
  end do
  close(30)
endif

!stop

!canonical gauge direction to measure the angle phi to
call random_number(gvec(1))
call random_number(gvec(2))
call random_number(gvec(3))

gvec(:) = gvec(:) / dsqrt(gvec(1)**2.d0 + gvec(2)**2.d0 + gvec(3)**2.d0)

a_to_b = 1.6179542d57        !c^3/(8*pi*h)

do i = 1,r_num
  g1 = lev_g(r_trans(1,i))
  g2 = lev_g(r_trans(2,i))

  t_Bij(i) = (1.d0*g2 / (1.d0*g1) ) * t_Aji(i) * a_to_b / (nu(i)**3.d0) 
end do

!write(*,*)t_Bij


allocate(all_trans(2,r_num+c_num))

call all_transitions(e_num,r_num,c_num,r_trans,c_trans,e_levs,all_trans,n_trans)


!kmax = 8
allocate(count_array(e_num))

call setup_count_array(j_lev,e_num,kmax,count_array,c_dim)

t_int = 1

allocate(pol_pops(node_num,e_num*2))
allocate(pol_rates(node_num,3,e_num*2))
allocate(JJ_20(node_num,r_num))

pol_pops(:,:)=0.d0
JJ_20(:,:)=0.d0
pol_rates(:,:,:)=0.d0

!$OMP PARALLEL DO
do node = 1,bulk_node_num
  call fwig_temp_init(2*100)
  if (ab_all(node).gt.0.d0) then
  !  call cpu_time(start)
  ! if statement only works if star is in the center
    if (star_mod.eq.0) then
      call trace_node(&
      &node,node_num,&
      &XYZ,bvec_all(:,node),gvec,&
      &del_ne,del_num_ne,maxne,nvel,&
      &nth,nph,xth,n_sphere,&
      &dx,tx,&
      &r_num,r_trans,c_num,c_temp,c_par,c_trans,col_temps,n_trans,all_trans,e_levs,&
      &nu,e_num,lev_p,lev_g,b_all,vel_all,&
      &t_Aji,t_Bij,dust_p,kd_v,t_Cji,h2_dens,temp_all,&
      &j_lev,kmax,n_depth,&
      &c_dim,count_array,&
      &star_mod,R_star,T_star,XYZ_star,&
      &verbose_mode,extra_info,&
      &pol_pops(node,:),pol_rates(node,:,:),JJ_20(node,:))
    else
      if (norm2(XYZ(:,node)).gt.R_star) then
        call trace_node(&
        &node,node_num,&
        &XYZ,bvec_all(:,node),gvec,&
        &del_ne,del_num_ne,maxne,nvel,&
        &nth,nph,xth,n_sphere,&
        &dx,tx,&
        &r_num,r_trans,c_num,c_temp,c_par,c_trans,col_temps,n_trans,all_trans,e_levs,&
        &nu,e_num,lev_p,lev_g,b_all,vel_all,&
        &t_Aji,t_Bij,dust_p,kd_v,t_Cji,h2_dens,temp_all,&
        &j_lev,kmax,n_depth,&
        &c_dim,count_array,&
        &star_mod,R_star,T_star,XYZ_star,&
        &verbose_mode,extra_info,&
        &pol_pops(node,:),pol_rates(node,:,:),JJ_20(node,:))
      endif
    endif
  endif
!  write(*,*)pol_pops(node,1:10)
!  stop

!  call cpu_time(finish)
!  write(*,*)finish-start

end do
!$OMP END PARALLEL DO

write(c_pol_pops,fmt="(A)")"/pol_pops.dat"
write(c_spont_rates,fmt="(A)")"/spont_rates.dat"
write(c_rad_rates,fmt="(A)")"/rad_rates.dat"
write(c_col_rates,fmt="(A)")"/col_rates.dat"
write(c_anis,fmt="(A)")"/rad_anis.dat"

c_pol_pops = trim(output_map) // trim(c_pol_pops)
c_spont_rates = trim(output_map) // trim(c_spont_rates)
c_rad_rates = trim(output_map) // trim(c_rad_rates)
c_col_rates = trim(output_map) // trim(c_col_rates)
c_anis = trim(output_map) // trim(c_anis)

open(unit=25,file=c_pol_pops,recl=8824)
do node = 1,node_num
  write(25,*)pol_pops(node,:)
end do
close(25)

if (extra_info.eq.1) then
  open(unit=26,file=c_spont_rates,recl=8824)
  do node = 1,node_num
    write(26,*)pol_rates(node,1,:)
  end do
  close(26)
  
  open(unit=27,file=c_rad_rates,recl=8824)
  do node = 1,node_num
    write(27,*)pol_rates(node,2,:)
  end do
  close(27)
  
  open(unit=28,file=c_col_rates,recl=8824)
  do node = 1,node_num
    write(28,*)pol_rates(node,3,:)
  end do
  close(28)

  open(unit=29,file=c_anis,recl=8824)
  do node = 1,node_num
    write(29,*)JJ_20(node,:)
  end do
  close(29)

endif

end program


