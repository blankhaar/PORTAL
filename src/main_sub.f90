subroutine trace_node(&
&node,node_num,&
&XYZ,bvec,gvec,&
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
&node_pops,node_rates,JJ_20)
!main driver of the node-specific polarized population-calculations
implicit none
integer node_num,e_num,r_num,c_num,c_temp
integer i,maxne,node,nth,nph,n_trans
integer c_dim,t_int,kmax,c_par,n_var
integer verbose_mode,star_mod,k_mu,k_phi
integer extra_info,nvel,n_depth,n_sphere,n_star

integer, dimension(e_num) :: j_lev,lev_g
integer, dimension(2,r_num) :: r_trans
integer, dimension(2,c_num) :: c_trans
integer, dimension(node_num) :: del_num_ne
integer, dimension(node_num,maxne) :: del_ne
integer, dimension(2,r_num+c_num) :: all_trans
integer, dimension(e_num) :: count_array

double precision dx,tx,R_star,T_star,d_tot
double precision dOm,th_star,ph_star,start,finish
double precision, parameter :: pi = dacos(-1.d0)

double precision, dimension(3) :: dir,bvec,gvec,XYZ_star

double precision, dimension(e_num,2:kmax) :: D0
double precision, dimension(3,nth*nph) :: xth,xth_s
double precision, dimension(3,node_num) :: XYZ,vel_all
double precision, dimension(node_num) :: h2_dens,dust_p,b_all
double precision, dimension(2,node_num) :: temp_all
double precision, dimension(e_num,node_num) :: lev_p
double precision, dimension(e_num) :: e_levs
double precision, dimension(r_num) :: t_Aji,nu,t_Bij,kd_v,J0,J2
double precision, dimension(r_num) :: JJ_20
double precision, dimension(c_temp,c_num,c_par) :: t_Cji
double precision, dimension(c_temp) :: col_temps
double precision, dimension(5,n_trans) :: AC
double precision, dimension(c_dim) :: rhovec 
double precision, dimension(3,c_dim) :: rates 
double precision, dimension(2*e_num) :: node_pops 
double precision, dimension(3,2*e_num) :: node_rates

t_int = 1

if (star_mod.eq.1) then
!if node is within star, then computation is meaningless
!  write(*,*)'i'
  dOm = pi * (R_star/norm2(XYZ_star(:)-XYZ(:,node)))**2.d0 
!  write(*,*)dOm

  dir = XYZ_star(:) - XYZ(:,node) 
  dir = dir / norm2(dir)

  xth_s(:,:) = 0.d0
  call star_all(nth,nph,dir,gvec,bvec,dOm,n_star,n_var,xth_s)
!  write(*,*)node,n_star,n_var
!  write(*,*)
!  d_tot = 0.d0
!  do i = 1,n_sphere
!    write(*,*)xth(:,i)
!    d_tot = d_tot + xth(3,i)
!  end do
!  write(*,*)
!  write(*,*)d_tot
  call inverse_rt(del_ne,del_num_ne,node_num,maxne,node,XYZ,nvel,n_var,n_star,bvec,gvec,dx,&
  &tx,r_num,r_trans,nu,e_num,lev_p,lev_g,b_all,vel_all,t_Aji,t_Bij,xth_s(:,1:n_var),&
  &dust_p,kd_v,temp_all(2,:),J0,J2,verbose_mode,star_mod,R_star,T_star,k_mu,k_phi,XYZ_star,n_depth)


else
  call inverse_rt(del_ne,del_num_ne,node_num,maxne,node,XYZ,nvel,n_sphere,n_star,bvec,gvec,dx,&
  &tx,r_num,r_trans,nu,e_num,lev_p,lev_g,b_all,vel_all,t_Aji,t_Bij,xth(:,1:n_sphere),&
  &dust_p,kd_v,temp_all(2,:),J0,J2,verbose_mode,star_mod,R_star,T_star,k_mu,k_phi,XYZ_star,n_depth)
endif

!open(unit=12,file="int_info.dat",recl=8824)
!write(12,*)J0
!write(12,*)J2
!close(12)
!stop

call make_AC(t_Aji,t_Cji,r_trans,c_trans,r_num,c_num,c_temp,c_par,e_num,e_levs,h2_dens(node),&
&temp_all(1,node),col_temps,n_trans,all_trans(:,1:n_trans),J0,J2,AC)

D0(:,:) = 0.d0

rhovec(:) = 0.d0
rates(:,:) = 0.d0

!call cpu_time(start)

call get_rhovec(c_dim,count_array,e_num,j_lev,n_trans,all_trans,AC,t_int,kmax,&
&temp_all(1,node),D0,rhovec,extra_info,rates)

!call cpu_time(finish)

!write(*,*)
!write(*,*)c_dim,kmax,finish-start
!write(*,*)

do i = 1,e_num
  if (j_lev(i).ne.0) then
    node_pops(2*(i-1)+1) = rhovec(count_array(i)+1)
    node_pops(2*(i-1)+2) = rhovec(count_array(i)+2)

    if (extra_info.eq.1) then
      node_rates(:,2*(i-1)+1) = rates(:,count_array(i)+1)
      node_rates(:,2*(i-1)+2) = rates(:,count_array(i)+2)
    endif
  else
    node_pops(2*(i-1)+1) = rhovec(count_array(i)+1)
    node_pops(2*(i-1)+2) = 0.d0 

    if (extra_info.eq.1) then
      node_rates(:,2*(i-1)+1) = rates(:,count_array(i)+1)
      node_rates(:,2*(i-1)+2) = 0.d0 
    endif
  endif

  if (isNaN(node_pops(2*(i-1)+1))) then
    node_pops(2*(i-1)+1)=0.d0
    node_pops(2*(i-1)+2)=0.d0
    node_rates(:,2*(i-1)+2)=0.d0
    node_rates(:,2*(i-1)+2)=0.d0
  endif
end do

if (extra_info.eq.1) JJ_20(:) = J2 / J0

end subroutine


