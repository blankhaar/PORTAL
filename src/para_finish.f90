program main 
use fwigxjpf

!program to do the second part of the PORTAL 
!calculation: using the aligned populations,
!obtain a map of the polarization of an astro
!physical region.

implicit none

character(len=100) :: input_file,output_map
character(len=120) :: file_I,file_Q,file_U,file_tau
character(len=100) :: in_pop,in_grid,in_mol,in_dust
character(len=100) :: in_pol_pop

integer node_num,del_num,e_num,r_num,c_num,c_temp
integer i,j,k,l,i1,maxne,node,node_count,nth,nph,g1,g2,n_trans
integer c_dim,kmax,im_pix,bulk_node_num,pix
integer verbose_mode,star_mod,k_mu,k_phi 
integer nlevels,nincs,nvel,n_rt,c_rt,c_par
integer k1,k2,k3,f1,f2,f3,f4,extra_info,max_ne_temp,n_depth_rt

integer, allocatable, dimension(:,:) :: del_tri
integer, allocatable, dimension(:) :: j_lev,lev_g
integer, allocatable, dimension(:,:) :: r_trans
integer, allocatable, dimension(:,:) :: c_trans
integer, allocatable, dimension(:,:) :: adj,adj0
integer, allocatable, dimension(:) :: del_num_ne 
integer, allocatable, dimension(:,:) :: del_ne 
integer, allocatable, dimension(:,:) :: del_ne_temp 
integer, allocatable, dimension(:,:) :: all_trans 
integer, allocatable, dimension(:) :: count_array 
integer, allocatable, dimension(:) :: l1,l2 
integer, allocatable, dimension(:,:) :: ijk_rt 
integer, allocatable, dimension(:,:,:) :: all_depth_list
integer, allocatable, dimension(:) :: all_node_count

double precision I_to_T,r1,r2,r3,dx,tx,R_star,T_star,x_dip
double precision a_to_b,start,finish,I_bb,dust_ii_gas
double precision th,ph,im_size,nu0,rat,dvel,Q_in,U_in
double precision dOm,th_star,ph_star,mass_mol,x1,x2,x3
double precision alpha,pix0,kasta_a,kasta_b,kasta_c,kasta_d,v0
double precision, parameter :: pi = dacos(-1.d0)

double precision, dimension(3) :: dir,bvec,XYZ_star,XYZ_star_inc

double precision, allocatable, dimension(:,:) :: XYZ,XYZ_inc
double precision, allocatable, dimension(:,:) :: vel_all,vel_inc
double precision, allocatable, dimension(:,:) :: temp_all
double precision, allocatable, dimension(:,:) :: lev_p
double precision, allocatable, dimension(:) :: e_levs
double precision, allocatable, dimension(:) :: t_Aji,nu,kd_v
double precision, allocatable, dimension(:,:,:) :: t_Cji
double precision, allocatable, dimension(:) :: dust_p
double precision, allocatable, dimension(:) :: col_temps
double precision, allocatable, dimension(:) :: h2_dens 
double precision, allocatable, dimension(:) :: b_all
double precision, allocatable, dimension(:) :: ab_all
double precision, allocatable, dimension(:) :: xth,gth,xph,gph 
double precision, allocatable, dimension(:) :: x_in,y_in
double precision, allocatable, dimension(:,:) :: bvec_all,bvec_inc
double precision, allocatable, dimension(:,:) :: all_I_in 
double precision, allocatable, dimension(:,:,:) :: all_pol_rayt
double precision, allocatable, dimension(:) :: Aji_all,kd_all
double precision, allocatable, dimension(:,:) :: pol_rayt 
double precision, allocatable, dimension(:) :: ph_az,th_inc 
double precision, allocatable, dimension(:) :: vel_out 
double precision, allocatable, dimension(:) :: nu0_all 
double precision, allocatable, dimension(:) :: th_all 
double precision, allocatable, dimension(:,:) :: pol_pops 
double precision, allocatable, dimension(:,:,:,:) :: all_I,all_Q,all_U,all_TAU
double precision, allocatable, dimension(:,:,:,:) :: dust_I,dust_Q,dust_U,dust_TAU
double precision, allocatable, dimension(:,:,:) :: array_i,array_q,array_u,array_tau

dust_ii_gas = 0.01d0

pix0 = 0.d0
alpha = 0.15d0

call getarg(1,input_file)
call getarg(2,output_map)

call getarg(3,in_pop)
call getarg(4,in_grid)
call getarg(5,in_mol)
call getarg(6,in_dust)

open(unit=101,file=input_file)
read(101,*) 
read(101,*)nlevels 

allocate(l1(nlevels))
allocate(l2(nlevels))
allocate(nu0_all(nlevels))
allocate(Aji_all(nlevels))
allocate(kd_all(nlevels))


read(101,*)l1
read(101,*)l2
read(101,*)
read(101,*)nincs

allocate(th_inc(nincs))
allocate(ph_az(nincs))

read(101,*)th_inc
read(101,*)ph_az

read(101,*)
read(101,*)nvel,dvel,v0

allocate(vel_out(-nvel:nvel))

read(101,*)
read(101,*)star_mod,T_star,R_star,XYZ_star(:)
read(101,*)x1,x2,x3,x_dip
read(101,*)kmax,extra_info
read(101,*)pix,im_pix,im_size
read(101,*)max_ne_temp,n_depth_rt
close(unit=101)

do i = -nvel,nvel
  vel_out(i) = i * dvel + v0
end do

!write(*,*)'i'

call fwig_table_init(2*100,9)

!tested and works
call eat_lime_dim(in_grid,in_mol,node_num,del_num,e_num,r_num,c_num,c_temp,c_par)

!write(*,*)'ii'

allocate(j_lev(e_num))     
allocate(lev_g(e_num))     
allocate(r_trans(2,r_num))  
allocate(c_trans(2,c_num)) 

allocate(XYZ(3,node_num))     
allocate(vel_all(3,node_num))     
allocate(XYZ_inc(3,node_num))     
allocate(vel_inc(3,node_num))     

allocate(temp_all(2,node_num))       
allocate(lev_p(e_num,node_num)) 
allocate(e_levs(e_num))          
allocate(t_Aji(r_num))
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

do i = 1,nlevels
  nu0_all = e_levs(l2(i)) - e_levs(l1(i))
end do

do i = 1,nlevels
  do j = 1,r_num
    if (l2(i).eq.r_trans(2,j).and.l1(i).eq.r_trans(1,j)) then
      Aji_all(i) = t_Aji(j)
      kd_all(i) = kd_v(j)
      exit
    endif
  end do
end do

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

dx = 1E-4

allocate(bvec_all(3,node_num))
allocate(bvec_inc(3,node_num))
!general magnetic field structure
if (x_dip.ne.0.d0)  then
  call dipole_magfield(node_num,XYZ,x1,x2,x3,bvec_all)
else
  call global_magfield(node_num,XYZ,x1,x2,x3,bvec_all)
endif


allocate(all_trans(2,r_num+c_num))

call all_transitions(e_num,r_num,c_num,r_trans,c_trans,e_levs,all_trans,n_trans)
!write(*,*)'iv'


write (in_pol_pop,fmt='(a)')'/pol_pops.dat' 
in_pol_pop = trim(output_map) // trim(in_pol_pop)
allocate(pol_pops(node_num,e_num*2))
open(unit=25,file=in_pol_pop,recl=8824)
do node = 1,node_num
  read(25,*)pol_pops(node,:)
!now weigh with the abundance and the h2_density
!  pol_pops(node,:) = ab_all(node) * h2_dens(node) * pol_pops(node,:)
!use the lime-results for the populations. But our results for the 
!polarized part.
  pol_pops(node,1) = lev_p(1,node)/dsqrt(1.d0*lev_g(1))      
  do i = 2,e_num
    rat = pol_pops(node,(i-1)*2+2)/pol_pops(node,(i-1)*2+1)
    pol_pops(node,(i-1)*2+2) = rat * lev_p(i,node)/dsqrt(1.d0*lev_g(i)) 
    pol_pops(node,(i-1)*2+1) = lev_p(i,node)/dsqrt(1.d0*lev_g(i))
  end do

end do
close(25)

dir(:) = 0.d0
dir(3) = 1.d0

if (pix.eq.0) then
  allocate(all_I(bulk_node_num,nlevels,nincs,-nvel:nvel))
  allocate(all_Q(bulk_node_num,nlevels,nincs,-nvel:nvel))
  allocate(all_U(bulk_node_num,nlevels,nincs,-nvel:nvel))
  allocate(all_TAU(bulk_node_num,nlevels,nincs,-nvel:nvel))
  
  Q_in = 0.d0
  U_in = 0.d0
  
  
  f1 = 51
  do i = 1,nincs
    call rotate_all_save(node_num,th_inc(i),ph_az(i),vel_all,bvec_all,XYZ,XYZ_star,&
  &vel_inc,bvec_inc,XYZ_inc,XYZ_star_inc) !also rotated star-coordinates required here
  
    write (file_I,fmt='(a,i2.2,a)')'/inc',i,'-xy.dat'
    file_I = trim(output_map) // trim(file_I)
    f1 = f1 + 1
  
    open(unit=f1,file=file_I,recl=8824)
    do node = 1,bulk_node_num
      write(f1,*)tx*XYZ_inc(1,node),tx*XYZ_inc(2,node)
    end do
    close(f1)
  
  
    if (allocated(all_depth_list)) deallocate(all_depth_list)
    if (allocated(all_node_count)) deallocate(all_node_count)
    if (allocated(all_I_in)) deallocate(all_I_in)
    if (allocated(all_pol_rayt)) deallocate(all_pol_rayt)
  
    allocate(all_depth_list(bulk_node_num,2,n_depth_rt))
    allocate(all_node_count(bulk_node_num))
    allocate(all_I_in(nlevels,bulk_node_num))
    allocate(all_pol_rayt(nlevels,4,node_num))
  
    all_pol_rayt(:,:,:) = 0.d0
    all_I_in(:,:) = 0.d0
    all_depth_list(:,:,:) = 0
    all_node_count(:) = 0
  
    !$OMP PARALLEL DO
    do node = 1,bulk_node_num
      call fwig_temp_init(2*100)
  
      call trace_out_raytrace(XYZ_inc(1,node),XYZ_inc(2,node),node_num,XYZ_inc,star_mod,R_star,XYZ_star_inc,T_star,&
  &del_ne,maxne,del_num_ne,dir,dx,n_depth_rt,all_depth_list(node,:,:),all_node_count(node),nlevels,nu0_all,all_I_in(:,node)) 
  
      if (pol_pops(node,1).ne.0.d0) then
        do j = 1,nlevels
          call comp_emcoefs(pol_pops(node,(l1(j)-1)*2+1),pol_pops(node,(l2(j)-1)*2+1),&
  &pol_pops(node,(l1(j)-1)*2+2),pol_pops(node,(l2(j)-1)*2+2),j_lev(l1(j)),j_lev(l2(j))&
  &,dacos(bvec_inc(3,node)),nu0_all(j),t_Aji(l1(j)),all_pol_rayt(j,1,node),all_pol_rayt(j,2,node),&
  &all_pol_rayt(j,3,node),all_pol_rayt(j,4,node))
        end do
      endif 
  
    end do
    !$OMP END PARALLEL DO
  
    !$OMP PARALLEL DO
    do node = 1,bulk_node_num
      do j = 1,nlevels
        do k = -nvel,nvel
          call ray_in_raytrace(vel_out(k),all_node_count(node),&
  &all_pol_rayt(j,:,all_depth_list(node,1,1:all_node_count(node))),b_all(all_depth_list(node,1,1:all_node_count(node)))&
  &,bvec_inc(:,all_depth_list(node,1,1:all_node_count(node))),&
  &dust_p(all_depth_list(node,1,1:all_node_count(node))),alpha,kd_all(j),&
  &temp_all(2,all_depth_list(node,1,1:all_node_count(node))),&
  &vel_inc(3,all_depth_list(node,1,1:all_node_count(node)))&
  &,all_depth_list(node,2,1:all_node_count(node)),tx*dx,nu0_all(j),all_I_in(j,node),all_I(node,j,i,k),&
  &all_Q(node,j,i,k),all_U(node,j,i,k),all_TAU(node,j,i,k),kasta_a,kasta_b,kasta_c,kasta_d)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end do
  
  
  f1 = 11
  do i = 1,nlevels
    I_to_T = 2.045D40 / (nu0_all(i)*nu0_all(i))

    do j = 1,nincs
      do k = -nvel,nvel
  
        write (file_I,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k+nvel+1,'-im.dat'
        file_I = trim(output_map) // trim(file_I)
        f1 = f1 + 1 
  
        open(unit=f1,file=file_I,recl=8824)
        do node = 1,bulk_node_num  
          write(f1,*)I_to_T*all_I(node,i,j,k),I_to_T*all_Q(node,i,j,k),&
&I_to_T*all_U(node,i,j,k),all_TAU(node,i,j,k)
        end do 
        close(f1)
  
      end do
    end do
  end do
else

  n_rt = (im_pix*2 + 1)*(im_pix*2 + 1)
  
  allocate(all_I(n_rt,nlevels,nincs,-nvel:nvel))
  allocate(all_Q(n_rt,nlevels,nincs,-nvel:nvel))
  allocate(all_U(n_rt,nlevels,nincs,-nvel:nvel))
  allocate(all_TAU(n_rt,nlevels,nincs,-nvel:nvel))
  
  allocate(dust_I(n_rt,nlevels,nincs,-nvel:nvel))
  allocate(dust_Q(n_rt,nlevels,nincs,-nvel:nvel))
  allocate(dust_U(n_rt,nlevels,nincs,-nvel:nvel))
  allocate(dust_TAU(n_rt,nlevels,nincs,-nvel:nvel))
  
  allocate(x_in(n_rt))
  allocate(y_in(n_rt))
  
  l = 0
  do i = -im_pix,im_pix
    do j = -im_pix,im_pix
      l = l + 1
      x_in(l) = i * im_size / im_pix
      y_in(l) = j * im_size / im_pix
    end do
  end do
  
  Q_in = 0.d0
  U_in = 0.d0


  do i = 1,nincs
    call rotate_all_save(node_num,th_inc(i),ph_az(i),vel_all,bvec_all,XYZ,XYZ_star,&
  &vel_inc,bvec_inc,XYZ_inc,XYZ_star_inc) !also rotated star-coordinates required here
  
    if (allocated(all_depth_list)) deallocate(all_depth_list)
    if (allocated(all_node_count)) deallocate(all_node_count)
    if (allocated(all_I_in)) deallocate(all_I_in)
    if (allocated(all_pol_rayt)) deallocate(all_pol_rayt)
  
    allocate(all_depth_list(n_rt,2,n_depth_rt))
    allocate(all_node_count(n_rt))
    allocate(all_I_in(nlevels,n_rt))
  
    allocate(all_pol_rayt(nlevels,4,node_num))
  
    all_pol_rayt(:,:,:) = 0.d0
    all_I_in(:,:) = 0.d0
    all_depth_list(:,:,:) = 0
    all_node_count(:) = 0
  
    !$OMP PARALLEL DO
    do node = 1,n_rt
      call trace_out_raytrace(x_in(node),y_in(node),node_num,XYZ_inc,star_mod,R_star,XYZ_star_inc,T_star,&
  &del_ne,maxne,del_num_ne,dir,dx,n_depth_rt,all_depth_list(node,:,:),all_node_count(node),nlevels,nu0_all,all_I_in(:,node))
    end do
    !$OMP END PARALLEL DO
  
    !$OMP PARALLEL DO
    do node = 1,bulk_node_num
      call fwig_temp_init(2*100)
  
      if (pol_pops(node,1).ne.0.d0) then
        do j = 1,nlevels
          call comp_emcoefs(pol_pops(node,(l1(j)-1)*2+1),pol_pops(node,(l2(j)-1)*2+1),&
  &pol_pops(node,(l1(j)-1)*2+2),pol_pops(node,(l2(j)-1)*2+2),j_lev(l1(j)),j_lev(l2(j))&
  &,dacos(bvec_inc(3,node)),nu0_all(j),Aji_all(j),all_pol_rayt(j,1,node),all_pol_rayt(j,2,node),&
  &all_pol_rayt(j,3,node),all_pol_rayt(j,4,node))
        end do
      endif
  
    end do
    !$OMP END PARALLEL DO
  
    !$OMP PARALLEL DO
    do node = 1,n_rt
      do j = 1,nlevels
        do k = -nvel,nvel
          call ray_in_raytrace(vel_out(k),all_node_count(node),&
  &all_pol_rayt(j,:,all_depth_list(node,1,1:all_node_count(node))),b_all(all_depth_list(node,1,1:all_node_count(node)))&
  &,bvec_inc(:,all_depth_list(node,1,1:all_node_count(node))),&
  &dust_p(all_depth_list(node,1,1:all_node_count(node))),alpha,kd_all(j),&
  &temp_all(2,all_depth_list(node,1,1:all_node_count(node))),&
  &vel_inc(3,all_depth_list(node,1,1:all_node_count(node)))&
  &,all_depth_list(node,2,1:all_node_count(node)),tx*dx,nu0_all(j),all_I_in(j,node),all_I(node,j,i,k),&
  &all_Q(node,j,i,k),all_U(node,j,i,k),all_TAU(node,j,i,k),&
  &dust_I(node,j,i,k),dust_Q(node,j,i,k),dust_U(node,j,i,k),dust_TAU(node,j,i,k))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end do

  allocate(array_i(2*im_pix+1,2*im_pix+1,2*nvel+1))
  allocate(array_q(2*im_pix+1,2*im_pix+1,2*nvel+1))
  allocate(array_u(2*im_pix+1,2*im_pix+1,2*nvel+1))
  allocate(array_tau(2*im_pix+1,2*im_pix+1,2*nvel+1))
  
  f1 = 101
  f2 = 201
  f3 = 301
  f4 = 401
  do i = 1,nlevels
  !  freq = nu0_all(i) 
    I_to_T = 2.045D40 / (nu0_all(i)*nu0_all(i))
    do j = 1,nincs
  
      array_i(:,:,:) =0.d0
      array_q(:,:,:) =0.d0
      array_u(:,:,:) =0.d0
      array_tau(:,:,:) =0.d0
  
  !put it in the right array
      l = 0
      do k1 = -im_pix,im_pix
        do k2 = -im_pix,im_pix
          l = l+1
  
          do k3 = -nvel,nvel
            array_i(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = I_to_T*all_I(l,i,j,k3)
            array_q(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = I_to_T*all_Q(l,i,j,k3)
            array_u(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = I_to_T*all_U(l,i,j,k3)
            array_tau(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = all_TAU(l,i,j,k3)
  !          write(*,*)array_i(k1+im_pix+1,k2+im_pix+1,k3+nvel+1)
          end do
        end do
      end do
  !    stop 
  
      write (file_I,fmt='(a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-im_i.fits'
      write (file_Q,fmt='(a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-im_q.fits'
      write (file_U,fmt='(a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-im_u.fits'
      write (file_TAU,fmt='(a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-im_tau.fits'
  
      file_I = trim(output_map) // trim(file_I)
      file_Q = trim(output_map) // trim(file_Q)
      file_U = trim(output_map) // trim(file_U)
      file_TAU = trim(output_map) // trim(file_TAU)

      do k3 = 1,2*nvel+1
  
        write (file_I,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-im_i.dat'
        write (file_Q,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-im_q.dat'
        write (file_U,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-im_u.dat'
        write (file_TAU,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-im_tau.dat'
  
        file_I = trim(output_map) // trim(file_I)
        file_Q = trim(output_map) // trim(file_Q)
        file_U = trim(output_map) // trim(file_U)
        file_TAU = trim(output_map) // trim(file_TAU)
  
        f1 = f1 + 1
        f2 = f2 + 1
        f3 = f3 + 1
        f4 = f4 + 1
  
        open(unit=f1,file=file_I,recl=8824)
        open(unit=f2,file=file_Q,recl=8824)
        open(unit=f3,file=file_U,recl=8824)
        open(unit=f4,file=file_TAU,recl=8824)
  
        do k1 = 1,2*im_pix+1
          write(f1,*)array_i(k1,:,k3)
          write(f2,*)array_q(k1,:,k3)
          write(f3,*)array_u(k1,:,k3)
          write(f4,*)array_tau(k1,:,k3)
        end do
  
        close(f1)
        close(f2)
        close(f3)
        close(f4)
      end do
  
    end do
  end do
  
  deallocate(array_i)
  deallocate(array_q)
  deallocate(array_u)
  deallocate(array_tau)
  
  allocate(array_i(2*im_pix+1,2*im_pix+1,2*nvel+1))
  allocate(array_q(2*im_pix+1,2*im_pix+1,2*nvel+1))
  allocate(array_u(2*im_pix+1,2*im_pix+1,2*nvel+1))
  allocate(array_tau(2*im_pix+1,2*im_pix+1,2*nvel+1))

  f1 = 501
  f2 = 601
  f3 = 701
  f4 = 701
  do i = 1,nlevels
  !  freq = nu0_all(i) 
    I_to_T = 2.045D40 / (nu0_all(i)*nu0_all(i))

    do j = 1,nincs
  
      array_i(:,:,:) =0.d0
      array_q(:,:,:) =0.d0
      array_u(:,:,:) =0.d0
      array_tau(:,:,:) =0.d0
  
  !put it in the right array
      l = 0
      do k1 = -im_pix,im_pix
        do k2 = -im_pix,im_pix
          l = l+1
  
          do k3 = -nvel,nvel
            array_i(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = I_to_T*dust_I(l,i,j,k3)
            array_q(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = I_to_T*dust_Q(l,i,j,k3)
            array_u(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = I_to_T*dust_U(l,i,j,k3)
            array_tau(k1+im_pix+1,k2+im_pix+1,k3+nvel+1) = dust_TAU(l,i,j,k3)
          end do
        end do
      end do
  
      do k3 = 1,2*nvel+1
  
        write (file_I,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-dust_i.dat'
        write (file_Q,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-dust_q.dat'
        write (file_U,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-dust_u.dat'
        write (file_TAU,fmt='(a,i2.2,a,i2.2,a,i2.2,a)')'/inc',j,'-trans',i,'-vel',k3,'-dust_tau.dat'
  
        file_I = trim(output_map) // trim(file_I)
        file_Q = trim(output_map) // trim(file_Q)
        file_U = trim(output_map) // trim(file_U)
        file_TAU = trim(output_map) // trim(file_TAU)
  
        f1 = f1 + 1
        f2 = f2 + 1
        f3 = f3 + 1
        f4 = f4 + 1
  
        open(unit=f1,file=file_I,recl=8824)
        open(unit=f2,file=file_Q,recl=8824)
        open(unit=f3,file=file_U,recl=8824)
        open(unit=f4,file=file_TAU,recl=8824)
  
        do k1 = 1,2*im_pix+1
          write(f1,*)array_i(k1,:,k3)
          write(f2,*)array_q(k1,:,k3)
          write(f3,*)array_u(k1,:,k3)
          write(f4,*)array_tau(k1,:,k3)
        end do
  
        close(f1)
        close(f2)
        close(f3)
        close(f4)
      end do
  
    end do
  end do
endif



end program


