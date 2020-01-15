subroutine inverse_rt(NN,neigh_list,tot_node,max_ne,node,X_all,nvel,n_sphere,n_star,bvec,gvec,dx,&
&tx,t_num,t_levs,nu,lev_num,lev_p,lev_g,temp_all,vel_all,t_Aji,t_Bij,xth,&
&dust_p,kd_v,bb_temps,J0,J2,verbose_mode,star_mode,R_star,T_star,k_mu,k_phi,X_star,&
n_depth)
!for a certain node construct all the relevant rays, send them out
!and raytrace them to the node. Product is Stokes-I for a range of 
!theta and phi angles. N.B.!: theta and phi are gauged with respect
!to the magnetic field direction of the central node and a canonical
!axis that is use to construct the x-axis (freely choosable because 
!phi is integrated out).    
implicit none
integer tot_node,node,nth,nph,i,j,k,l,verbose_mode
integer t_num,lev_num,max_ne,l1,l2,node_count
integer star_mode,k_mu,k_phi,nvel,n_depth,tn
integer mode,n_sphere,n_star 

integer, dimension (tot_node,max_ne) :: NN
integer, dimension (tot_node) :: neigh_list 
integer, dimension (2,n_depth) :: depth_list 
integer, dimension (2,t_num) :: t_levs 
integer, dimension (lev_num) :: lev_g 

double precision vel_node,pi,mu,ph,dx,tx,I_in,I_out,nu0,c_in,velcorr
double precision R_star,T_star,tau_all

double precision, dimension (3) :: bvec,gvec,dir,X_star
double precision, dimension (3,n_sphere) :: xth
double precision, dimension (3,tot_node) :: X_all
double precision, dimension (tot_node) :: temp_all
double precision, dimension (tot_node) :: bb_temps 
double precision, dimension (tot_node) :: dust_p 
double precision, dimension (3,tot_node) :: vel_all 
double precision, dimension (n_depth) :: vel_ray
double precision, dimension (lev_num,tot_node) :: lev_p 
double precision, dimension (t_num) :: nu,t_Aji,t_Bij,kd_v
double precision, dimension (t_num) :: J0,J2 

!input:
!NN		-- matrix of neighbours for all nodes
!node		-- central node number
!nth		-- number of theta angles
!nph		-- number of phi angles
!bvec		-- magnetic field direction of central node (should be unit)
!gvec		-- gauge vector for definition of phi  (should be unit)
!X_all		-- all coordinates of the nodes 
!dx		-- ray-trace step precision
!tx		-- simulation size 
!max_ne 	-- maximum number of neigbours for all nodes
!tot_node 	-- total number of nodes 
!neigh_list	-- number of neighbours for each node
!t_num		-- number of transitions
!lev_num	-- number of levels
!lev_p		-- density of all levels for all nodes 
!dust_p		-- dust mass density for all nodes 
!kd_v		-- dust opacities at appropriate frequencies  
!nu		-- frequencies of transitions
!temp_all	-- (mass-weighed) temperature of all nodes T = 2*kB*Temp / m
!bb_temps	-- temperature of the black-body temperatures 
!vel_all	-- velocity of all nodes
!t_Aji		-- Einstein A's for all transitions 
!t_Bij		-- Einstein B's for all transitions 
!lev_g		-- degeneracy for all transtions 
!xth,gth	-- cos(th) elements and their associated weights
!xph,gph	-- ph elements and their associated weights

!output:
!J0		-- J_0^0 element of radiation tensor for each transition  
!J2		-- J_0^2 element of radiation tensor for each transition

J0(:) = 0.d0
J2(:) = 0.d0

pi = dacos(-1.d0)

!open(unit=13,file="j_map.dat",recl=8824)
!there is going to be nph*nth rays associated with the central node
do i = 1,n_sphere
  mu = xth(1,i)
  ph = xth(2,i)

!second order legendre polyonmial
  c_in = (3.d0 * mu**2.d0 - 1.d0)/2.d0
!for a certain angle th and ph, find associated ray direction
  call find_dir(bvec,gvec,dacos(mu),ph,dir)

!now trace the photon package to either the outward CMB, source, or an inward source 
  if (star_mode.eq.1.and.i.le.n_star) then
    mode = 1
    call trace_out(node,tot_node,X_all,NN,max_ne,neigh_list,dir,dx,&
&depth_list,node_count,n_depth,mode,X_star,R_star) 
  else
    mode = 0
    call trace_out(node,tot_node,X_all,NN,max_ne,neigh_list,dir,dx,&
&depth_list,node_count,n_depth,mode,X_star,R_star) 
  endif

!we need the velocities projected on the ray-direction
  vel_node = vel_all(1,node)*dir(1) + vel_all(2,node)*dir(2) + vel_all(3,node)*dir(3)

  vel_ray(:) = 0.d0
  do l = 1,node_count
    tn = depth_list(1,l) 
    vel_ray(l) = vel_all(1,tn)*dir(1) + vel_all(2,tn)*dir(2) + vel_all(3,tn)*dir(3)
  end do
!with the path known, we will ray-trace for all the relevant transitions
  do k = 1,t_num   
    nu0 = nu(k)

    l1 = t_levs(1,k) 
    l2 = t_levs(2,k)

!subroutine that determines the intensity of the incoming ray from the point where it starts   
    if (star_mode.eq.1.and.i.le.n_star) then
      call black_body(T_star,nu0,I_in)
    else
      call black_body(2.73d0,nu0,I_in)  
    endif
!subroutine that from the ray-path and the associated node-elements, and input
!intensity
!ray traces towards the central node.   
    call ray_in(node,temp_all(node),vel_node,lev_p(l1,node),lev_p(l2,node),dust_p(node),bb_temps(node),node_count,&
&lev_p(l1,depth_list(1,1:node_count)),lev_p(l2,depth_list(1,1:node_count)),dust_p(depth_list(1,1:node_count)),&
&kd_v(k),bb_temps(depth_list(1,1:node_count)),temp_all(depth_list(1,1:node_count)),vel_ray(1:node_count),&
&depth_list(2,1:node_count),nvel,tx*dx,nu0,t_Aji(k),t_Bij(k),lev_g(l1),lev_g(l2),I_in,I_out,tau_all)

!now we'll add to the integration
!we already have normalized, so no 4*pi factor 
    J0(k) = J0(k) + I_out * xth(3,i)  
    J2(k) = J2(k) + I_out * xth(3,i) * c_in / dsqrt(2.d0)

!    if (k.eq.1) write(13,*)xth(:,i),I_out

  end do
end do 
!close(13)

end subroutine

subroutine trace_out(node,tot_node,X_all,NN,max_ne,neigh_list,dir,dx,depth_list,node_count,n_depth,&
&star_mode,X_star,R_star)
!subroutine that traces from central node to an outward source, listing the
!nodes it passes in depth_list
implicit none
integer tot_node,node,new_node,max_ne,i,star_mode
integer in_sim,neigh_num,k_cur,node_count,node_save,n_depth

integer, dimension (tot_node,max_ne) :: NN
integer, dimension (tot_node) :: neigh_list
integer, dimension (2,n_depth) :: depth_list 

double precision dx,delta_x,R_star
double precision start,finish

double precision, dimension (3) :: X_star,X_phot,dir
double precision, dimension (3,tot_node) :: X_all

!input
!tot_node	-- total node number
!node		-- central node number
!NN     	-- matrix of neighbours for all nodes
!node   	-- central node number
!X_all  	-- all coordinates of the nodes 
!dx     	-- ray-trace step precision
!max_ne		-- maximum number of neighbours
!neigh_list	-- list of number of neighbours for each node
!dir		-- photon direction 
 
!output
!depth_list	-- list with nodes to go through including distance between them
!node_count	-- how many nodes does the ray-trace pass

node_save = node

depth_list(:,:) = 0
depth_list(1,1) = node_save
k_cur = 0
node_count = 0
in_sim = 0

node = node_save

X_phot(:) = X_all(:,node)

!call cpu_time(start)

do while (in_sim.eq.0)
  neigh_num = neigh_list(node)
!  if (star_mode.eq.1) write(*,*)node
  call check_neigh_fast(X_phot,dir,X_all(:,node),neigh_num,NN(node,1:neigh_num),&
&X_all(:,NN(node,1:neigh_num)),star_mode,R_star,X_star,delta_x,new_node)

  if (new_node.eq.0) in_sim=1

  node_count = node_count + 1
!  write(*,*)node_count
  depth_list(2,node_count) = nint(delta_x/dx)
  depth_list(1,node_count+1) = new_node 
  node = new_node
end do

do i = 2,node_count
  depth_list(2,i) = depth_list(2,i-1) + depth_list(2,i)
end do


if (node_count.gt.n_depth) then
  write(*,*)'node_count=',node_count,'while n_depth=',n_depth
  write(*,*)'adjust n_depth parameter! Make it bigger.'
  stop
endif

node = node_save

end subroutine


subroutine ray_in(node,temp_node,vel_node,p1_node,p2_node,dust_node,bb_node,node_count,dens_1,dens_2,dust_dens,kd,&
&bb_list,temp_list,vel_list,length_list,nvel,dx,nu0,Aji,Bij,g1,g2,I_in,I_out,tau_all)
!using the ray-trace list of all the nodes with associated distances, using information about these
!nodes (temp,vel,densities), ray-trace towards the central node 
implicit none
integer node,node_count,g1,g2,i,nv,nvel

integer, dimension(node_count) :: length_list 

double precision Aji,Bij,I_in,I_out,dx,nu0,spi,nbpi
double precision Sv,Iv,kv,ev,p1,p2,tau_v,e_tau,p1_node,p2_node
double precision c1,c2,p1_corr,p2_corr,vel_node,temp_node,kd
double precision alpha_v,j_v,I_bb,dust_node,bb_node
double precision tau_all,corr
double precision cnu,v0,v_node,phi_v0,phi_v,d_v,I_tot

double precision, dimension(node_count) :: dens_1,dens_2,dust_dens,temp_list
double precision, dimension(node_count) :: vel_list,bb_list

!input
!node		--- central node		
!temp_node	--- (mass weighed) temperature of central node
!bb_node	--- temperature of central node
!vel_node	--- velocity of central node
!p1_node	--- lower level density central node
!p2_node	--- upper level density cnetral node
!node_count	--- number of nodes passed in inward tracing 
!dens_1		--- list with lower level densities inward tracing nodes
!dens_2		--- list with upper level densities inward tracing nodes
!dust_dens	--- dust mass densities inward tracing nodes
!kd             --- dust opacity at appropriate transition 
!temp_list	--- list with (mass weighed) temperatures inward tracing nodes
!bb_list	--- list with temperatures inward tracing nodes
!vel_list	--- list with velocities of inward tracing nodes
!length_list	--- list with lengths (in steps) spent in each node
!dx		--- distance weighing of length_list (dx * steps = length)
!nu0		--- central line-frequency
!Aji		--- Einstein Spontaneous emission rate
!Bij		--- Einstein Absorption rate
!g1		--- degeneracy lower level
!g2		--- degeneracy upper level
!I_in		--- incoming intensity

!output:
!I_out		-- line-intensity at the center

c1 = 2.2102191D-42 * nu0   !hv0/c
c2 = 5.272859D-35 * nu0    !hv0/4*pi 

cnu = 299792458.d0 / nu0

I_tot = 0.d0       !\int dv \phi(v) I_v 
d_v  = 2.d0 * temp_node / nvel 


!tau_all = 0.d0
spi = dsqrt(dacos(-1.d0))

do nv = -nvel,nvel

!velocity of this channel and the associated weight-function
!\int dv \phi(v) I_v

  v0 = (2.d0 * temp_node) * nv / nvel
  phi_v0 = (1.d0/(spi*temp_node)) * dexp(-(v0/temp_node)**2.d0)

  v0 = v0 + vel_node 

  Iv = I_in
  do i = node_count,2,-1
  !we start at the node before the end-node and count our way back
    p1 = dens_1(i)
    p2 = dens_2(i) 
  
  !compute the absorption and emission coefficient
  !no line scattering
    kv = c1 * Bij * (p1 - (1.d0*g1)*p2/(1.d0*g2)) 
    ev = c2 * Aji * p2 
  
    call black_body(bb_list(i),nu0,I_bb) 
  
    alpha_v = kd * dust_dens(i) 
    j_v     = alpha_v * I_bb       

    phi_v = cnu*(1.d0/(spi*temp_list(i))) * dexp(-((v0-vel_list(i))/temp_list(i))**2.d0) 
  
    tau_v = (kv*phi_v + alpha_v) * dx * (length_list(i) - length_list(i-1))    !now with dust

!    write(*,'(2x,1000E16.8)')j_v,alpha_v,kv,ev

  !now propagate
    if (tau_v.gt.1D-30) then       !if not, don't even bother -- will result in NaN
    !source function
  !    Sv = ev/kv
      Sv = (ev*phi_v + j_v)/(kv*phi_v + alpha_v)     !now with dust
  
      e_tau = dexp(-tau_v)         
      Iv = Iv * e_tau + Sv * (1.d0 - e_tau) 
    endif
  end do

!  if (isnan(Iv)) stop

  kv = c1 * Bij * (p1_node - (1.d0*g1)*p2_node/(1.d0*g2)) 
  ev = c2 * Aji * p2_node
  
  !dust opacities 
  call black_body(bb_node,nu0,I_bb)
  
  alpha_v = kd * dust_node
  j_v     = alpha_v * I_bb       

  phi_v = phi_v0 * cnu 
 
  !Sv = ev/kv
  Sv = (ev*phi_v + j_v)/(kv*phi_v + alpha_v)
  !tau_v = kv * dx * length_list(1)
  tau_v = (kv*phi_v + alpha_v) * dx * length_list(1)
  
  e_tau = dexp(-tau_v)
  I_out = Iv * e_tau + Sv * (1.d0 - e_tau)

  I_tot = I_tot + phi_v0 * I_out * d_v
end do

I_out = I_tot

end subroutine

