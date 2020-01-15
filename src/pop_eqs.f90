subroutine t_abs_const(j,k,jl,kl,t0,t2)
!subroutine to compute the (t^0)_{jk}^{j_l k_l} and (t^2)_{jk}^{j_l k_l} 
!constants for absorption t_A from L&L04 Eqs.(7.20a). 
!tested against matlab routines
use fwigxjpf

implicit none
integer j,k,jl,kl
integer dj,dk
integer m,ml,q

double precision t0,t2
double precision w1,w2,w3,ww,w4_0,w4_2
double precision ljl,fac1,fas,fac

!input:
!j	-- level 1 angular momentum
!k	-- level 1 rank
!jl	-- level 2 (lower) angular momentum
!kl	-- level 2 rank

!output:
!t0	-- (t^0)_{jk}^{j_l k_l} factor
!t0	-- (t^2)_{jk}^{j_l k_l} factor

dj = jl - j
dk = kl - k

t0 = 0.d0
t2 = 0.d0

if (abs(dj).le.1) then
!compute some factors that are required and summarize them in fac
  ljl  = 2.d0 * jl + 1.d0
  fac1 = dsqrt(3.d0 * (2.d0*k+1.d0) * (2.d0*kl+1.d0) )
  fas  = (-1.d0)**(1.d0 - dj)

  fac = ljl * fac1 * fas

  if (dk.eq.0.or.abs(dk).eq.2) then
    !sum for the t0 and t2 constants
    do m = -j,j
      do q = -1,1
        ml = q + m 
!these are the Wigner-3j symbols that are computed using the fwigxjpf
!programs -- times two because it also treats half-integers and has
!integers as input
        w1 = fwig3jj (2* j , 2* j    , 2* k ,&
        &              2* m , 2* (-m) , 0    )
 
        w2 = fwig3jj (2* jl , 2* jl    , 2* kl ,&
        &              2* ml , 2* (-ml) , 0    )

        w3 = fwig3jj (2* j , 2* jl    , 2* 1 ,&
        &              2* (-m) , 2* ml , 2* (-q)  )

        w4_0 = fwig3jj (2* 1 , 2* 1    , 2* 0 ,&
        &              2* q , 2* (-q) , 0    )

        w4_2 = fwig3jj (2* 1 , 2* 1    , 2* 2 ,&
        &              2* q , 2* (-q) , 0    )

        ww = w1*w2*w3*w3

!only if k=kl do we have a contribution from the J_0^0 element   
        if (dk.eq.0) t0 = t0 + ww*w4_0      
        t2 = t2 + ww*w4_2
      end do
    end do

!extra sqrt(5) factor from the 2*kr+1
    t0 = t0 * fac
    t2 = t2 * dsqrt(5.d0) * fac

  endif
endif

end subroutine

subroutine t_stim_const(j,k,ju,ku,t0,t2)
!subroutine to compute the (t^0)_{jk}^{j_l k_l} and (t^2)_{jk}^{j_l k_l} 
!constants for stimulated emission from L&L04 Eqs.(7.20c). 
!tested against matlab routines
use fwigxjpf

implicit none
integer j,k,ju,ku
integer dj,dk
integer m,mu,q

double precision t0,t2
double precision w1,w2,w3,ww,w4_0,w4_2
double precision lju,fac1,fas,fac

!input:
!j	-- level 1 angular momentum
!k	-- level 1 rank
!ju	-- level 2 (upper) angular momentum
!ku	-- level 2 rank

!output:
!t0	-- (t^0)_{jk}^{j_l k_l} factor
!t0	-- (t^2)_{jk}^{j_l k_l} factor

dj = ju - j
dk = ku - k

t0 = 0.d0
t2 = 0.d0

if (abs(dj).le.1) then
!compute some factors that are required and summarize them in fac
  lju  = 2.d0 * ju + 1.d0
  fac1 = dsqrt(3.d0 * (2.d0*k+1.d0) * (2.d0*ku+1.d0) )
  fas  = (-1.d0)**(1.d0 - dj)

  fac = lju * fac1 * fas

  if (dk.eq.0.or.abs(dk).eq.2) then
    !sum for the t0 and t2 constants
    do m = -j,j
      do q = -1,1
        mu = m - q 
!these are the Wigner-3j symbols that are computed using the fwigxjpf
!programs -- times two because it also treats half-integers and has
!integers as input 
        w1 = fwig3jj (2* j , 2* j    , 2* k ,&
        &              2* m , 2* (-m) , 0    )
 
        w2 = fwig3jj (2* ju , 2* ju    , 2* ku ,&
        &              2* mu , 2* (-mu) , 0    )

        w3 = fwig3jj (2* ju , 2* j    , 2* 1 ,&
        &              2* (-mu) , 2* m , 2* (-q)  )

        w4_0 = fwig3jj (2* 1 , 2* 1    , 2* 0 ,&
        &              2* q , 2* (-q) , 0    )

        w4_2 = fwig3jj (2* 1 , 2* 1    , 2* 2 ,&
        &              2* q , 2* (-q) , 0    )

        ww = w1*w2*w3*w3

!only if k=ku do we have a contribution from the J_0^0 element   
        if (dk.eq.0) t0 = t0 + ww*w4_0      
        t2 = t2 + ww*w4_2
      end do
    end do

    t0 = t0 * fac
    t2 = t2 * dsqrt(5.d0) * fac

  endif
endif

end subroutine

subroutine t_spont_const(j,k,ju,ku,t0)
!subroutine to compute the (t^0)_{jk}^{j_l k_l} and (t^2)_{jk}^{j_l k_l} 
!constants for spontaneous emission from L&L04 Eqs.(7.20b). 
!tested against matlab routines
use fwigxjpf

implicit none
integer j,k,ju,ku
integer dj,dk
integer m,mu,q

double precision t0
double precision w1,lju,fas

!input:
!j      -- level 1 angular momentum
!k      -- level 1 rank
!ju     -- level 2 (upper) angular momentum
!ku     -- level 2 rank

!output:
!t0     -- (t^0)_{jk}^{j_l k_l} factor


dj = ju - j
dk = ku - k

t0 = 0.d0

!delta_kku and selection rule
if (abs(dj).le.1.and.dk.eq.0) then
!required factors
  lju  = 2.d0 * ju + 1.d0
  fas  = (-1.d0)**(1.d0 + 1.d0*j + 1.d0*ju + 1.d0*k)

!Wigner 6j-symbol
  w1 = fwig6jj (2* ju , 2* ju , 2* k ,&
  &              2* j , 2* j  , 2* 1    )
 
  t0 = lju * fas * w1 

endif

end subroutine


subroutine r_abs_const(j,k,ju,kp,r0,r2)
!subroutine to compute the (r^0)_{jk}^{j_u k_p} and (r^2)_{jk}^{j_u k_p} 
!constants for absorption r_A from L&L04 Eqs.(7.20d). 
!tested against matlab routines
use fwigxjpf

implicit none
integer j,k,ju,kp
integer dj,dk
integer m,mu,q

double precision r0,r2
double precision w1,w2,w3,ww,w4_0,w4_2
double precision lj,fac1,fas,fac

!input:
!j	-- level 1 angular momentum
!k	-- level 1 rank
!ju	-- level 2 (upper) angular momentum
!kp	-- level 2 rank (k-prime)

!output:
!r0	-- (r^0)_{jk}^{j_l k_l} factor
!r0	-- (r^2)_{jk}^{j_l k_l} factor

dj = ju - j
dk = kp - k

r0 = 0.d0
r2 = 0.d0

if (abs(dj).le.1) then
!compute some factors that are required and summarize them in fac
  lj  = 2.d0 * j + 1.d0
  fac1 = dsqrt(3.d0 * (2.d0*k+1.d0) * (2.d0*kp+1.d0) )

  if (dk.eq.0.or.abs(dk).eq.2) then
    !sum for the r0 and r2 constants
    do m = -j,j
      do q = -1,1
        mu = m-q 

        fas = (-1.d0)**(1.d0 + q)
        fac = lj * fac1 * fas

!these are the Wigner-3j symbols that are computed using the fwigxjpf
!programs -- times two because it also treats half-integers and has
!integers as input
        w1 = fwig3jj (2* j , 2* j    , 2* k ,&
        &              2* m , 2* (-m) , 0    )
 
        w2 = fwig3jj (2* j , 2* j    , 2* kp ,&
        &              2* m , 2* (-m) , 0    )

        w3 = fwig3jj (2* ju , 2* j    , 2* 1 ,&
        &              2* (-mu) , 2* m , 2* (-q)  )

        w4_0 = fwig3jj (2* 1 , 2* 1    , 2* 0 ,&
        &              2* q , 2* (-q) , 0    )

        w4_2 = fwig3jj (2* 1 , 2* 1    , 2* 2 ,&
        &              2* q , 2* (-q) , 0    )

        ww = w1*w2*w3*w3

!only if k=kl do we have a contribution from the J_0^0 element   
        if (dk.eq.0) r0 = r0 + fac*ww*w4_0      
        r2 = r2 + fac*ww*w4_2
      end do
    end do

!extra sqrt(5) factor from the 2*kr+1
    r2 = r2 * dsqrt(5.d0) 

  endif
endif

end subroutine

subroutine r_stim_const(j,k,jl,kp,r0,r2)
!subroutine to compute the (r^0)_{jk}^{j_l k_l} and (r^2)_{jk}^{j_l k_l} 
!constants for stimulated emission from L&L04 Eqs.(7.20f). 
!tested against matlab routines
use fwigxjpf

implicit none
integer j,k,jl,kp
integer dj,dk
integer m,ml,q

double precision r0,r2
double precision w1,w2,w3,ww,w4_0,w4_2
double precision lj,fac1,fas,fac

!input:
!j	-- level 1 angular momentum
!k	-- level 1 rank
!jl	-- level 2 (lower) angular momentum
!kp	-- level 2 rank (k-prime)

!output:
!t0	-- (r^0)_{jk}^{j_u k_p} factor
!t0	-- (r^2)_{jk}^{j_u k_p} factor

dj = jl - j
dk = kp - k

r0 = 0.d0
r2 = 0.d0

if (abs(dj).le.1) then
!compute some factors that are required and summarize them in fac
  lj  = 2.d0 * j + 1.d0
  fac1 = dsqrt(3.d0 * (2.d0*k+1.d0) * (2.d0*kp+1.d0) )

  if (dk.eq.0.or.abs(dk).eq.2) then
    !sum for the r0 and r2 constants
    do m = -j,j
      do q = -1,1
        ml = m + q 

        fas = (-1.d0)**(1.d0 - q)
        fac = lj * fac1 * fas

!these are the Wigner-3j symbols that are computed using the fwigxjpf
!programs -- times two because it also treats half-integers and has
!integers as input 
        w1 = fwig3jj (2* j , 2* j    , 2* k ,&
        &              2* m , 2* (-m) , 0    )
 
        w2 = fwig3jj (2* j , 2* j    , 2* kp ,&
        &              2* m , 2* (-m) , 0    )

        w3 = fwig3jj (2* j , 2* jl    , 2* 1 ,&
        &              2* (-m) , 2* ml , 2* (-q)  )

        w4_0 = fwig3jj (2* 1 , 2* 1    , 2* 0 ,&
        &              2* q , 2* (-q) , 0    )

        w4_2 = fwig3jj (2* 1 , 2* 1    , 2* 2 ,&
        &              2* q , 2* (-q) , 0    )

        ww = w1*w2*w3*w3

!only if k=kp do we have a contribution from the J_0^0 element   
        if (dk.eq.0) r0 = r0 + fac*ww*w4_0      
        r2 = r2 + fac*ww*w4_2
      end do
    end do

    r2 = r2 * dsqrt(5.d0)

  endif
endif

end subroutine

subroutine setup_stateq(c_dim,count_array,nlev,j_lev,trans,ntrans,kmax,&
&AC_trans,D0,Temp,Amat,A_mat,B_mat,C_mat)
!subroutine to setup up the statistical equilibrium equation matrix
!subroutine still has to be polished and tested
!calls to the constant-calculations can be replaced by pre-calculated
!matrices to save computational effort 
!no de-polarization term D^(K) built in as of yet

!The Hazel program has a built-in kmax=6 -- we have kmax=8

!The Hazel program substitutes the first row of the Amat matrix by
!an equation that says \sum_i sqrt(2*j_i+1) * rho^0_i = 1 
!maybe we can do it in a different way

implicit none
integer i,c_dim,nlev,ntrans
integer i1,j1,k1
integer i2,j2,k2
integer i3,j3,q
integer c1,c2,kmax
integer c_up_trans,c_down_trans,t_up,t_down 
integer, dimension(2,ntrans) :: trans 
integer, dimension(ntrans) :: up_trans 
integer, dimension(ntrans) :: down_trans 
integer, dimension(3,c_dim) :: C_all
integer, dimension(nlev) :: j_lev,count_array 

double precision A21,B21,B12,C21,C12,a_to_b
double precision ts0,ts2,te0,ta0,ta2,r0,r2,JR0,JR2 
double precision t_rad,t_col,nu0,col_abs,col_stim
double precision rad_abs,rad_stim,g1,g2,g3,Temp
double precision col_dep

double precision, dimension(nlev,2:kmax) :: D0
double precision, dimension(5,ntrans) :: AC_trans 
double precision, dimension(5,ntrans) :: AC_up_trans 
double precision, dimension(5,ntrans) :: AC_down_trans 
double precision, dimension(c_dim,c_dim) :: Amat 
double precision, dimension(c_dim,c_dim) :: A_mat 
double precision, dimension(c_dim,c_dim) :: B_mat 
double precision, dimension(c_dim,c_dim) :: C_mat 

!input
!c_dim		-- dimension of Amat
!count_array	-- counting array that finds initial counting 
!		   number associated with certain level
!nlev		-- number of levels
!trans		-- list of transitions (only lower --> upper)
!ntrans		-- total number of transitions. N.B. these are 
!		   collisional and radiative transitions combined
!AC_trans	-- transition coefficients (Aij,Cij,nu0,J0,J2)
!Temp		-- temperature required to obtain C12 from 
!                  Einstein-Milne relation
!D0             -- Collisional de-polarization rates per level

!output
!Amat		-- statistical equillibrium matrix

a_to_b = 6.781962D49    !best to use L&L relation here
                        !B21 = [c^2 / (2*h * v^3)] * A21
                        !a_to_b = c^2 / (2*h)

Amat(:,:) = 0.d0
A_mat(:,:) = 0.d0
B_mat(:,:) = 0.d0
C_mat(:,:) = 0.d0

do i1 = 1,nlev
!find relevant transitions
!divide relevant transitions up into excitations
!and de-excitations
  c_up_trans = 0
  c_down_trans = 0
  do i2 = 1,nlev
    do q = 1,ntrans
      if (trans(1,q).eq.i1.and.trans(2,q).eq.i2) then
        c_up_trans = c_up_trans + 1
        up_trans(c_up_trans) = i2  
        AC_up_trans(1,c_up_trans) = AC_trans(1,q) !A21 
        AC_up_trans(2,c_up_trans) = AC_trans(2,q) !C21
        AC_up_trans(3,c_up_trans) = AC_trans(3,q) !nu0
        AC_up_trans(4,c_up_trans) = AC_trans(4,q) !J0
        AC_up_trans(5,c_up_trans) = AC_trans(5,q) !J2
      endif

      if (trans(1,q).eq.i2.and.trans(2,q).eq.i1) then
        c_down_trans = c_down_trans + 1
        down_trans(c_down_trans) = i2  
        AC_down_trans(1,c_down_trans) = AC_trans(1,q) 
        AC_down_trans(2,c_down_trans) = AC_trans(2,q)
        AC_down_trans(3,c_down_trans) = AC_trans(3,q)
        AC_down_trans(4,c_down_trans) = AC_trans(4,q)
        AC_down_trans(5,c_down_trans) = AC_trans(5,q)
      endif
    end do
  end do 

!  write(*,*)i1,c_up_trans,c_down_trans
!  write(*,*)up_trans(1:c_up_trans)
!  write(*,*)down_trans(1:c_down_trans)
!  write(*,*)
  
  j1 = j_lev(i1)
  g1 = 2.d0 * j1 + 1.d0 

!now we have c_up_trans relevant upper transitions
  do t_up = 1,c_up_trans
!   j1 -- lower 
!   j2 -- upper

    i2 = up_trans(t_up)
    j2 = j_lev(i2)          
 
    g2 = 2.d0 * j2 + 1.d0 

    A21 = AC_up_trans(1,t_up) 
    C21 = AC_up_trans(2,t_up)
    nu0 = AC_up_trans(3,t_up)
    JR0  = AC_up_trans(4,t_up)
    JR2  = AC_up_trans(5,t_up)

    B21 = a_to_b * A21 / nu0**3.d0  

!    write(*,*)i1,j1,i2,j2
!    write(*,*)A21,C21,B21,nu0,JR0,JR2
!    stop

    c1 = count_array(i1)
    do k1 = 0,min(2*j1,kmax),2
      c1 = c1 + 1
      c2 = count_array(i2) 
      do k2 = 0,min(2*j2,kmax),2
        c2 = c2 + 1  
        !upper transition factors
        call t_stim_const(j1,k1,j2,k2,ts0,ts2)
        call t_spont_const(j1,k1,j2,k2,te0)

        !radiative coefficent
        t_rad = B21 * (ts0 * JR0 + ts2 * JR2) + A21 * te0 

        !collisional coefficient only if k1=k2
        t_col = 0.d0
        if (k1.eq.k2) then
!currently (as in GK), only the rank-0 contribution of the collisional
!operator to the rate equations
          if (k1.eq.0) t_col = dsqrt(g2/g1)*C21          !no degeneracy terms? 
        endif

        Amat(c1,c2) = t_rad + t_col         !Amat is only defined once
        A_mat(c1,c2) = A_mat(c1,c2) + A21 * te0
        B_mat(c1,c2) = B_mat(c1,c2) + B21 * (ts0 * JR0 + ts2 * JR2) 
        C_mat(c1,c2) = C_mat(c1,c2) + t_col 
      end do
    end do
  end do

!now we have c_down_trans relevant down transitions
  do t_down = 1,c_down_trans
!   j1 -- upper
!   j2 -- lower

    i2 = down_trans(t_down)
    j2 = j_lev(i2)          
    g2 = 2.d0 * j2 + 1.d0

    A21  = AC_down_trans(1,t_down) 
    C21  = AC_down_trans(2,t_down)
    nu0  = AC_down_trans(3,t_down)
    JR0  = AC_down_trans(4,t_down)
    JR2  = AC_down_trans(5,t_down)

    B21 = a_to_b * A21 / nu0**3.d0 
    B12 = B21 * g1 / g2 

    C12 = C21 * g1 * dexp(-4.79924d-11 * nu0 / Temp) / g2 !constant is h/kB 

    c1 = count_array(i1)
    do k1 = 0,min(2*j1,kmax),2
      c1 = c1 + 1
      c2 = count_array(i2) 
      do k2 = 0,min(2*j2,kmax),2
        c2 = c2 + 1  
        !upper transition factors
        call t_abs_const(j1,k1,j2,k2,ta0,ta2)

        !radiative coefficent
        t_rad = B12 * (ta0 * JR0 + ta2 * JR2) 

        !collisional coefficient only if k1=k2
        t_col = 0.d0
        if (k1.eq.k2) then
!currently (as in GK), only the rank-0 contribution of the collisional
!operator to the rate equations
          if (k1.eq.0) t_col = dsqrt(g2/g1)*C12          
        endif

        Amat(c1,c2) = t_rad + t_col         !Amat is only defined once
        B_mat(c1,c2) = B_mat(c1,c2) + t_rad 
        C_mat(c1,c2) = C_mat(c1,c2) + t_col 
      end do
    end do
  end do

!inner level depedancies
  i2 = i1
  j2 = j1

  c1 = count_array(i1)
  do k1 = 0,min(2*j1,kmax),2
    c1 = c1 + 1
    c2 = count_array(i2)
    do k2 = 0,min(2*j2,kmax),2
      c2 = c2 + 1

      rad_abs = 0.d0
      col_abs = 0.d0
      !all excitations away to upper levels 
      do t_up = 1,c_up_trans
        !j1 -- lower
        !j3 -- upper

        i3 = up_trans(t_up)
        j3 = j_lev(i3)
        g3 = 2.d0 * j3 + 1.d0

        A21  = AC_up_trans(1,t_up)
        C21  = AC_up_trans(2,t_up)
        nu0  = AC_up_trans(3,t_up)
        JR0  = AC_up_trans(4,t_up)
        JR2  = AC_up_trans(5,t_up)
    
        B21 = a_to_b * A21 / nu0**3.d0 
        B12 = B21 * g3 / g1
        C12 = C21 * g3 * dexp(-4.79924d-11 * nu0 / Temp) / g1 !constant is h/kB (K/Hz) 

        call r_abs_const(j1,k1,j3,k2,r0,r2)
        rad_abs = rad_abs + B12 * (JR0 * r0 + JR2 * r2)

        B_mat(c1,c2) = B_mat(c1,c2) - B12 * (JR0 * r0 + JR2 * r2)

        if (k1.eq.k2) then
          col_abs = col_abs + C12

          C_mat(c1,c2) = C_mat(c1,c2) - C12

!          if (k1.eq.0) write(*,*)j1,j3,B12*(JR0*r0+JR2*r2),r0

        endif
      end do

      !all de-excitations to lower levels 
      rad_stim = 0.d0      !also includes spontaneous emission
      col_stim = 0.d0 
      do t_down = 1,c_down_trans
        i3 = down_trans(t_down)
        j3 = j_lev(i3)

        A21 = AC_down_trans(1,t_down)
        C21 = AC_down_trans(2,t_down)
        nu0 = AC_down_trans(3,t_down)
        JR0  = AC_down_trans(4,t_down)
        JR2  = AC_down_trans(5,t_down)
    
        B21 = a_to_b * A21 / nu0**3.d0 

        call r_stim_const(j1,k1,j3,k2,r0,r2)

        rad_stim = rad_stim + B21 * (JR0 * r0 + JR2 * r2)
        B_mat(c1,c2) = B_mat(c1,c2) - B21 * (JR0 * r0 + JR2 * r2) 

        if (k1.eq.k2) then
          rad_stim = rad_stim + A21  
          col_stim = col_stim + C21

          A_mat(c1,c2) = A_mat(c1,c2) - A21
          C_mat(c1,c2) = C_mat(c1,c2) - C21 

!          if (k1.eq.0) write(*,*)j1,j3,B21*(JR0*r0+JR2*r2),r0
        endif
      end do

      if (k1.eq.k2.and.k2.ge.2) then
        !collisional de-polarization rates
        col_dep = D0(i1,k1) 
        C_mat(c1,c2) = C_mat(c1,c2) - col_dep
      endif
!note the minus sign -- Amat is only defined once 
      Amat(c1,c2) = -1.d0*(rad_abs + rad_stim + col_abs + col_stim + col_dep) 
    end do
  end do
end do

end subroutine

