      subroutine lagrange_1(q,n,g,calc_grad_g,eps,iter_max)
      implicit none
c-----solve the problem
c     J(q^) -> min
c     g(q)=0
c------------------------------
c     where J(q^)=sum{j=1,n}(q_i**2)
c----------------------------------------------------------
c     method
c
c     Lagendre function:
c     L(q,lambda)=J(q^)+lambda*g(q)
c
c     q^_k+1=q^_k-alpha_k *grad( L(q^k,lamda_k) )
c
c     lamda_k+1=lambda_k*alpha_k*d_L/d_labmda(q_k,lambda_k)
c
c     In our case:
c     J(q^)=(q^)**2
c     grad(L)= 2*q^+lambda*grad(g(q^))
c
c     d_L/d_labmda = g(q^)
c
c     q^k+1 = q^_k - alpha_k*(2*q^_k+lambda*grad(g(q^)))
c     lambda_k+1 = lambda_k + alpha_k*g(q^_k)
c
c     The condition to find alpha_k:
c     Let introduce the function 
c     f_k(alpha) = J(q^_k - alpha * grad(J(q^_k)))

c     f_k(alpha) = J(q^_k - alpha * (2*q^_k+lambda*grad(g(q^))))

c     alpha will determine from the condition
c     f_k(alpha_k) = min[f_k(alpha)] at alpha_k >0
c
c------------------------------------------------------------------   
c
c     The equation to get monimal falue of f_k(alpha) is 
c     d_f_k(alpha)/d_alpha=0
c
c     In our case d_f_k(alpha)/d_alpha=0 at 
c
c     alpha_k=sum{i=1,n}(2*q_i_k**2+lambda_k*q_i_k*d_g(u)/d_q_i)/
c          [ sum{i=1,}(2*q_i_k+lambda_k*d_g(q^)/d_q_i)**2)]
c----------------------------------------------------------------

c-----input
      integer n,  ! dimension of u vector
     &iter_max    !maximalnumber of iterations
 
      real*8 eps ! accuracy for condition to finish calculations
                  ! at |g|<eps
c-----ouptut
      real*8 q(n) 

c-----extenals    
      real*8 g    ! function of condition g(u)=0
      external g,calc_grad_g

c-----local
      real*8 lambda_k,lambda_k_p,    ! Legendre multiplier
     &alpha_k,sum1,sum2,p,g_k,g_k_old

      real*8, dimension(n) :: q_k,q_k_p, grad_g_k

      integer k,        ! iteration number
     &iter,i 

c----------------------------------------------------
c     zero iteration values
c----------------------------------------------------- 
      do i=1,n
        q_k(i)=0.d0
        q_k_p(i)=0.d0
      enddo
      g_k=g(q_k,n)
      g_k_old=g_k

      lambda_k=0.d0
      iter=0
c-----------------------------------------------------
 10   continue
      iter=iter+1

      g_k=g(q_k_p,n)

      if( abs(g_k).gt.abs(g_k_old)) then
        do i=1,n
          q(i)=q_k(i)
        enddo
        goto 20
      endif
      
      if( abs(g_k).lt.eps) then
        do i=1,n
          q(i)=q_k_p(i)
        enddo
        goto 20
      endif

      do i=1,n
         q_k(i)=q_k_p(i)
      enddo
  
      if (iter.gt.iter_max) then
        write(*,*)'in lagrange_1 number of iterations>iter_max'
        q(i)=q_k(i)
        goto 20
      endif

      lambda_k=lambda_k_p
c--------------------------------------------------   
      call calc_grad_g(q,grad_g_k,n)

      sum1=0.d0
      sum2=0.d0

      do i=1,n
         p=2*q_k(i)+lambda_k*grad_g_k(i)
         sum1=sum1 + q_k(i)*p
         sum2=sum2 + p*p
      enddo 

      alpha_k=sum1/sum2

      do i=1,n
         p=2*q_k(i)+lambda_k*grad_g_k(i)
         q_k_p(i)=q_k(i)-alpha_k*p
      enddo 
 
      lambda_k_p=lambda_k-alpha_k*g_k
            
      goto 10

 20   continue

      return
      end

      subroutine correction_lagrange_1(z,r,phi,cnz,cnr,cm,iter_max)
c---------------------------------------------------------
c     ray corrections using lagrange_1
c---------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
c-----input - output  
      real*8 z,r,phi,cnz,cnr,cm
      integer iter_max !maximal number of iteration

c-----local    
      integer n    
      real*8  eps         ! accuracy
      real*8, dimension(6) :: u
      real*8, dimension(2) :: q  !at n=2
      
      real*8 z_l,r_l,phi_l,cnz_l,cnr_l,cm_l
      common /cor_lagrange_1/ z_l,r_l,phi_l,cnz_l,cnr_l,cm_l

c-----externals
      real*8 b,g
      external b,g,grad_g
      
      
c-----put ray coordinates to common block  /cor_lagrange_1/
      z_l=z
      r_l=r
      phi_l=phi
      cnz_l=cnz
      cnr_l=cnr
      cm_l=cm

      u(1)=z
      u(2)=r
      u(3)=phi
      u(4)=cnz
      u(5)=cnr
      u(6)=cm


      bmod=b(z,r,phi)
      n=2
      eps=prmt4 

      call lagrange_1(q,n,g,grad_g,eps,iter_max)


      return
      end

      return
      end


      real*8 function g(q,n)
c-----calculate condition function 
c     g=hamiltonian(z,r,phi,cnz+q(1),cnr+q(2),cm)
      implicit none 
      include 'param.i'
      include 'one.i'
c-----input      
      integer n
      real*8 q(n)

      real*8 z_l,r_l,phi_l,cnz_l,cnr_l,cm_l
      common /cor_lagrange_1/ z_l,r_l,phi_l,cnz_l,cnr_l,cm_l

c-----local
      real*8 gam,cnz_cor,cnr_cor
c-----external
      real*8 gamma1,hamilt1

      cnz_cor=cnz_l + q(1)
      cnr_cor=cnr_cor + q(2)

      gam=gamma1(z_l,r_l,phi_l,cnz_cor,cnr_cor,cm_l)
      
      g=hamilt1(z_l,r_l,phi_l,cnz_cor,cnr_cor,cm_l)

      return
      end

      subroutine calc_grad_g(q,grad_g,n)
c-----calculate derivatives from g=hamiltonian
c     grad_g(1)=d_hamiltonian/d_z
c     grad_g(2)=d_hamiltonian/d_r
c     
c     at z_l,r_l,phi_l,cm_l
c     and cnz_cor=cnz_l + q(1),  cnr_cor=cnr_l + q(2)
c
c     z_l,r_l,phi_l,cnz_l,cnr_l,cm_l are given in common//cor_lagrange_1/
c
      implicit none

      include 'param.i'
      include 'one.i'

c-----input      
      integer n
      real*8 q(n)

      real*8 z_l,r_l,phi_l,cnz_l,cnr_l,cm_l
      common /cor_lagrange_1/ z_l,r_l,phi_l,cnz_l,cnr_l,cm_l

c-----output
      real*8 grad_g(n)
      
c-----local
      real*8 cnz_cor,cnr_cor,t
      real*8, dimension(6) :: u,deru

c-----external dddrz1(t,u,deru)
      
      cnz_cor=cnz_l + q(1)
      cnr_cor=cnr_l + q(2)

      u(1)=z_l
      u(2)=r_l
      u(3)=phi_l 
      u(4)=cnz_cor
      u(5)=cnr_cor
      u(6)=cm_l

      call dddrz1(t,u,deru)

      grad_g(1)=deru(1)
      grad_g(2)=deru(2)
      
      
      return
      end
