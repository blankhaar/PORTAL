subroutine polarize_prop(eta_I,eta_Q,eps_I,eps_Q,ds,I,Q,U)
!subroutine to do a propagation of polarized radiation over a cell. 
!this method assumes constant conditions over the entire propagation.
!tested against a matlab routine that computes the propagation in a 
!different manner by a finite difference method. 
implicit none

double precision I,Q,U,I0,Q0,ds
double precision eta_I,eta_Q,eps_I,eps_Q
double precision shQ,chQ,e_tau,fac,o11,o12


!input:
!I,Q		-- incoming Stokes I and Q parameters
!eta_I		-- opacity 
!eta_Q		-- polarized opacity 
!eps_I		-- spontaneous emission rate 
!eps_Q		-- polarized emission rate
!ds		-- depth length 

!output:
!I,Q		-- outgoing Stokes parameters

!we use Eq. (8.14) from L&L04
!I(s) = \int_0^s ds' O(s,s') e(s') + O(s,0) I(0)
!where I and e are vectors, and O is a 2x2 matrix
!we use the fact that within a cell e(s') = e. Also
!we use the expression from Eq. (8.28) for O(s,s')
!				  | cosh(eta_Q*(s-s')	-sinh(eta_Q*(s-s')|
!O(s,s') = exp(-eta_I * (s-s')) * |  					  |
!				  |-sinh(eta_Q*(s-s')	 cosh(eta_Q*(s-s')|

!expressions for the integrated elements of O(s,s') are computed using Wolfram
!Alpha. See following links for results:
! https://www.wolframalpha.com/input/?i=%5Cint+exp(-a*(s-x))+*+cosh(b*(s-x))+dx+from+0+to+s 
! https://www.wolframalpha.com/input/?i=%5Cint+exp(-a*(s-x))+*+sinh(b*(s-x))+dx+from+0+to+s 

I0 = I
Q0 = Q

shQ = dsinh(eta_Q*ds)
chQ = dcosh(eta_Q*ds)

e_tau = dexp(-eta_I*ds)

fac = (eta_I**2.d0 - eta_Q**2.d0)**(-1.d0)

o11 = fac*eta_I*(1.d0 - (chQ + (eta_Q/eta_I)*shQ)*e_tau)
o12 =-fac*eta_Q*(1.d0 - (chQ + (eta_I/eta_Q)*shQ)*e_tau)

I = o11*eps_I + o12*eps_Q + e_tau*( chQ*I0 - shQ*Q0)
Q = o12*eps_I + o11*eps_Q + e_tau*(-shQ*I0 + chQ*Q0)
U = U * e_tau


end subroutine 

subroutine trace_out_raytrace(x_co,y_co,tot_node,X_all,star_mod,R_star,XYZ_star,T_star,NN,&
&max_ne,neigh_list,dir,dx,n_depth_rt,depth_list,node_count,nlevels,nu0,I_in)
!subroutine that traces from a starting point on the sphere to the observer on the other side of the
!sphere, listing the nodes it passes in depth_list
implicit none
integer tot_node,node,new_node,max_ne,i,n_depth_rt,star_mod
integer in_sim,neigh_num,k_cur,node_count,node_save
integer nlevels

integer, dimension (tot_node,max_ne) :: NN
integer, dimension (tot_node) :: neigh_list
integer, dimension (2,n_depth_rt) :: depth_list

double precision dx,r,rsave
double precision R_star,T_star
double precision xcorr,ycorr,rcorr
double precision x_co,y_co,z_co

double precision, dimension (3) :: start,XYZ_star
double precision, dimension (3) :: X_phot,dir
double precision, dimension (3,tot_node) :: X_all
double precision, dimension (nlevels) :: nu0 
double precision, dimension (nlevels) :: I_in 

!input
!tot_node       -- total node number
!start          -- starting coordinates 
!X_all          -- all coordinates of the nodes 
!dx             -- ray-trace step precision
!max_ne         -- maximum number of neighbours
!neigh_list     -- list of number of neighbours for each node
!dir            -- photon direction/positive z-direction 

!output
!depth_list     -- list with nodes to go through including distance between them
!node_count     -- how many nodes does the ray-trace pass

!we start by finding the stating point for this calculation.
!

z_co = -1.d0*dsqrt(1.d0 - x_co**2.d0 - y_co**2.d0)
do i = 1,nlevels
  call black_body(2.73d0,nu0(i),I_in(i))
end do

if (star_mod.eq.1) then
  xcorr = x_co - XYZ_star(1)
  ycorr = y_co - XYZ_star(2)

  rcorr = dsqrt(xcorr**2.d0 + ycorr**2.d0)
  if (rcorr.lt.R_star) then
    z_co = XYZ_star(3) + R_star
    do i = 1,nlevels
      call black_body(T_star,nu0(i),I_in(i))
    end do
  endif
endif  

start(1) = x_co
start(2) = y_co
start(3) = z_co

in_sim     = 0
k_cur      = 0
node_count = 0

!find the nearest node for the starting position
!unfortunately, we will have to check every node 
rsave = 10.d0

do i = 1,tot_node
  r = norm2(start(:)-X_all(:,i))
  if (r.lt.rsave) then
    rsave = r
    node = i
  endif
end do

node_save = node

!photon start is the node position
X_phot(:) = start 

do while (in_sim.eq.0)
!evolve photon position
  X_phot(:) = X_phot(:) + dx * dir(:)
  k_cur = k_cur + 1

  if (norm2(X_phot).gt.1.d0) then
    node_count = node_count + 1
    in_sim=1
    depth_list(2,node_count) = k_cur
    depth_list(1,node_count) = node
  endif

  neigh_num = neigh_list(node)

  call check_neigh(X_phot,node,X_all(:,node),neigh_num,NN(node,1:neigh_num),X_all(:,NN(node,1:neigh_num)),new_node)

!  write(*,*)'iii'
  if (new_node.ne.node) then
    node_count = node_count + 1

!how long did the photon stay with the n'th (node_count) node for  
!what is the new node -- save this information for the output
    depth_list(2,node_count) = k_cur
    depth_list(1,node_count) = new_node
    node = new_node
  endif
end do

if (node_count.gt.n_depth_rt) then
  write(*,*)'node_count=',node_count,'while n_depth_rt=',n_depth_rt
  write(*,*)'adjust n_depth_rt parameter! Make it bigger.'
  stop
endif

node = node_save

end subroutine

subroutine ray_in_raytrace(vel_obs,node_count,pol_rayt_list,b_list,bvec_list,&
&dust_dens,alpha,kd,bb_list,vel_list,length_list,dx,nu0,I_in,I_out,Q_out,U_out,tau,&
&I_dust,Q_dust,U_dust,tau_dust)
!using the ray-trace list of all the nodes with associated distances, 
!using information about these nodes (temp,vel) and the earlier computed
!polariztion parameters to ray-trace through the simulation 
implicit none
integer node,node_count,g1,g2,i

integer, dimension(node_count) :: length_list

double precision I_in,I_out,Q_in,Q_out,U_in,U_out,dx,nu0
double precision v_corr,vel_obs,delta_x,I_corr,kd
double precision eta_I,eta_Q,eps_I,eps_Q,tau,cnu
double precision I_dust,Q_dust_in,Q_dust_out,U_dust_in,U_dust_out
double precision sig_dust_0,sig_dust_2,mu,mu2,sI,sQ,I_bb,e_dust
double precision tau_dust,tau_dust_tot,alpha,Q_dust,U_dust
double precision phi,q_phi_rot,u_phi_rot 

double precision, dimension(node_count) :: bb_list,vel_list,b_list,dust_dens
double precision, dimension(3,node_count) :: bvec_list
double precision, dimension(4,node_count) :: pol_rayt_list 

!input
!vel_obs        --- velocity of observation 
!node_count     --- number of nodes passed in inward tracing 
!pol_rayt_list  --- list of polarization coefficients 
!b_list         --- list with (mass weighed) temperatures inward tracing nodes
!bvec_list      --- list of magnetic field vectors 
!vel_list       --- list with velocities of inward tracing nodes
!length_list    --- list with lengths (in steps) spent in each node
!dx             --- distance weighing of length_list (dx * steps = length)
!nu0            --- central line-frequency
!I_in           --- incoming Stokes I 
!Q_in           --- incoming Stokes Q 
!U_in           --- incoming Stokes U 

!output:
!I_out           --- outgoing Stokes I 
!Q_out           --- outgoing Stokes Q 
!U_out           --- outgoing Stokes U 

cnu = 2.99792E8 / nu0

I_out = I_in
Q_in = 0.d0
U_in = 0.d0
tau = 0.d0
tau_dust_tot = 0.d0

I_dust = I_in
sig_dust_0 = 0.d0
sig_dust_2 = 0.d0
Q_dust_in = 0.d0
U_dust_in = 0.d0


do i = 1,node_count-1 
!correct for the velocity-shift --- Maxwell-Boltzmann distribution from the Temperature
!temperature is in fact a temperature that is weighed already over particle size 
  v_corr = (1.d0/(b_list(i)*dsqrt(dacos(-1.d0))))*dexp((-(vel_obs - vel_list(i))**2.d0)/(b_list(i)**2.d0))    

!rotate Q and U parameters from the global axis system 
!to the node axis system 
  call rot_stokes(Q_in,U_in,Q_out,U_out,bvec_list(:,i))
!  call rot_stokes(Q_dust_in,U_dust_in,Q_dust_out,U_dust_out,bvec_list(:,i))

  mu = bvec_list(3,i)
  mu2 = 1.d0 - mu * mu          !gamma = angle of B with plane of the sky  

!we have already computed the propagation parameters
  eta_I = cnu * v_corr * pol_rayt_list(1,i) 
  eta_Q = cnu * v_corr * pol_rayt_list(2,i) 
  eps_I = cnu * v_corr * pol_rayt_list(3,i) 
  eps_Q = cnu * v_corr * pol_rayt_list(4,i) 

!add the dust part

!compute the optical depth for this cell-path
  delta_x = dx * (length_list(i+1) - length_list(i))
!now propagate
  if (delta_x*eta_I.gt.1D-30) then       !if not, don't even bother -- will result in NaN
  !source function
 
    call polarize_prop(eta_I,eta_Q,eps_I,eps_Q,delta_x,I_out,Q_out,U_out) 
    tau = tau + eta_I * delta_x

  endif

  tau_dust = kd * dust_dens(i) * delta_x
  if (tau_dust.gt.1D-30) then
    call black_body(bb_list(i),nu0,I_bb) 

    e_dust = dexp(-tau_dust)
    
    sI = 0.5d0 * alpha * (mu2 - 2.d0 / 3.d0 )
    sQ = alpha * mu2
     
    I_dust = I_dust*e_dust + I_bb*(1-e_dust)
    sig_dust_0 = sig_dust_0*e_dust + I_bb*(1-e_dust)
    sig_dust_2 = sig_dust_2*e_dust + sI*I_bb*(1-e_dust)

    phi = dacos(bvec_list(2,i)/dsqrt(bvec_list(1,i)**2.d0 + bvec_list(2,i)**2.d0))
    if (bvec_list(1,i).lt.0.d0) phi = -phi
  
    q_phi_rot = dcos(2.d0 * phi)
    u_phi_rot = dsin(2.d0 * phi)
 
    Q_dust_in = Q_dust_in*e_dust + sQ*I_bb*(1-e_dust) * q_phi_rot
    U_dust_in = U_dust_in*e_dust + sQ*I_bb*(1-e_dust) * u_phi_rot
    tau_dust_tot = tau_dust_tot + tau_dust
  endif

!  write(*,*)i,eta_I,eta_I*delta_x
!rotate the Q and U parameters from the node axis system
!to the global axis system
  call rot_stokes_back(Q_out,U_out,Q_in,U_in,bvec_list(:,i))
!  call rot_stokes_back(Q_dust_out,U_dust_out,Q_dust_in,U_dust_in,bvec_list(:,i))


end do

!correct for the black-body background radiation.
call black_body(2.73d0,nu0,I_corr)

I_out = I_out - I_corr
Q_out = Q_in
U_out = U_in

I_dust = I_dust - I_corr
if (I_dust.eq.0.d0) then
  Q_dust = 0.d0
  U_dust = 0.d0
  tau_dust = 0.d0
else 
  Q_dust = I_dust * Q_dust_in/(sig_dust_0 - sig_dust_2)
  U_dust = I_dust * U_dust_in/(sig_dust_0 - sig_dust_2)
  tau_dust = tau_dust_tot
endif
!convert to temperature units?

!write(*,*)I_out,Q_out,U_out
!write(*,*)I_dust,Q_dust,U_dust

end subroutine


 
