c-----------------------------------------------------------------------
c                                                                       
c   purpose             - zeros of an analytic complex function         
c                           using the muller method with deflation      
c                                                                       
c   usage               - call muller (f,eps,nsig,kn,nguess,n,x,itmax,  
c                           infer,ier)                                  
c                                                                       
c   arguments    f      - a complex function subprogram, f(z), written  
c                           by the user specifying the equation whose   
c                           roots are to be found.  f must appear in    
c                           an external statement in the calling pro-   
c                           gram.                                       
c                eps    - 1st stopping criterion.  let fp(z)=f(z)/p     
c                           where p = (z-z(1))*(z-z(2))*,,,*(z-z(k-1))  
c                           and z(1),...,z(k-1) are previously found    
c                           roots.  if ((cdabs(f(z)).le.eps) .and.       
c                           (cdabs(fp(z)).le.eps)), then z is accepted   
c                           as a root. (input)                          
c                nsig   - 2nd stopping criterion.  a root is accepted   
c                           if two successive approximations to a given 
c                           root agree in the first nsig digits. (input)
c                             note. if either or both of the stopping   
c                             criteria are fulfilled, the root is       
c                             accepted.                                 
c                kn     - the number of known roots which must be stored
c                           in x(1),...,x(kn), prior to entry to muller 
c                nguess - the number of initial guesses provided. these 
c                           guesses must be stored in x(kn+1),...,      
c                           x(kn+nguess).  nguess must be set equal     
c                           to zero if no guesses are provided. (input) 
c                n      - the number of new roots to be found by        
c                           muller (input)                              
c                x      - a complex vector of length kn+n.  x(1),...,   
c                           x(kn) on input must contain any known       
c                           roots.  x(kn+1),..., x(kn+n) on input may,  
c                           on user option, contain initial guesses for 
c                           the n new roots which are to be computed.   
c                           if the user does not provide an initial     
c                           guess, zero is used.                        
c                           on output, x(kn+1),...,x(kn+n) contain the  
c                           approximate roots found by muller.          
c                itmax  - the maximum allowable number of iterations    
c                           per root (input)                            
c                infer  - an integer vector of length kn+n.  on         
c                           output infer(j) contains the number of      
c                           iterations used in finding the j-th root    
c                           when convergence was achieved.  if          
c                           convergence was not obtained in itmax       
c                           iterations, infer(j) will be greater than   
c                           itmax (output).                             
c                ier    - error parameter (output)                      
c                         warning error                                 
c                           ier = 33 indicates failure to converge with-
c                             in itmax iterations for at least one of   
c                             the (n) new roots.                        
c                                                                       
c                                                                       
c   remarks      muller always returns the last approximation for root j
c                in x(j). if the convergence criterion is satisfied,    
c                then infer(j) is less than or equal to itmax. if the   
c                convergence criterion is not satisified, then infer(j) 
c                is set to either itmax+1 or itmax+k, with k greater    
c                than 1. infer(j) = itmax+1 indicates that muller did   
c                not obtain convergence in the allowed number of iter-  
c                ations. in this case, the user may wish to set itmax   
c                to a larger value. infer(j) = itmax+k means that con-  
c                vergence was obtained (on iteration k) for the defla-  
c                ted function                                           
c                              fp(z) = f(z)/((z-z(1)...(z-z(j-1)))      
c                                                                       
c                but failed for f(z). in this case, better initial      
c                guesses might help or, it might be necessary to relax  
c                the convergence criterion.                             
c                                                                       
C ENM -- made all variables double precision that were real
C ENM 3 jan 2006 -- For use with prep3d.f iabsorp=12 option,
C   need to modify so that the function f takes two input variables
C   but only one is used, so no need to rewrite all the dispfun code
C   in eric_disp.f. (uses dum, which is the variable for eps in dispfun)

c-----------------------------------------------------------------------
c                                                                       
      subroutine muller (f,eps,nsig,kn,nguess,n,x,itmax,infer,ier)
      implicit double precision (a-h,o-z)                                  
c      dimension           x(1),infer(1)  
      dimension           x(*),infer(*)  
      double precision    rzero,rten,rhun,rp01,ax,eps1,qz,eps,tpq       
      complex*16          x,d,dd,den,fprt,frt,h,rt,t1,t2,t3,            
     1                    tem,z0,z1,z2,bi,xx,xl,y0,y1,y2,x0,          
     2                    zero,p1,one,four,p5,f,dum(1:3,1:3)
      external f
      data                zero/(0.0d0,0.0d0)/,p1/(0.1d0,0.0d0)/,                
     1                    one/(1.0d0,0.0d0)/,four/(4.0d0,0.0d0)/,               
     2                    p5/(0.5d0,0.0d0)/,                                
     3                    rzero/0.0d0/,rten/10.0d0/,rhun/100.0d0/,            
     4                    ax/0.1d0/,ickmax/3/,rp01/0.01d0/                  
c                                  first executable statement           
      ier = 0      

c      write(*,*)'Muller eps,nsig,kn,nguess,n,itmax,ier',              
c     &                  eps,nsig,kn,nguess,n,itmax,ier
        
c      print *,'Muller: zero: ',zero,' p1: ',p1,' one: ',one
c      print*,'n',n 
      if (n .lt. 1) go to 9005                                          
      eps1 = rten **(-nsig)                                             
      eps1 = min(eps1,rp01)
c      print *,'eps1'                                            
c                                  set number of iterations             
      knp1 = kn+1                                                       
      knpn = kn+n                                                       
      knpng = kn+nguess

c      print *,'knpn',knpn
                                                 
      do 5 i=1,knpn                                                     
         infer(i) = 0                                                   
         if (i .gt. knpng) x(i) = zero                                  
    5 continue                       
c      write(*,*)'infer',infer
c      write(*,*)'x',x
                                   
      l= knp1                                                           
   10 jk = 0                                                            
      ick = 0                                                           
      xl = x(l)                                                         
   15 ic = 0                                                            
      h = ax                                                            
      h = p1*h                                                          
      if (cdabs(xl) .gt. ax) h = p1*xl                                   
c                                  first three points are               
c                                    xl+h,  xl-h,  xl                   
      rt = xl+h                                                         
      assign 20 to nn                                                   
      go to 50                                                          
   20 z0 = fprt                                                         
      y0 = frt                                                          
      x0 = rt                                                           
      rt = xl-h                                                         
      assign 25 to nn                                                   
      go to 50                                                          
   25 z1 = fprt                                                         
      y1 = frt                                                          
      h = xl-rt                                                         
      d = h/(rt-x0)                                                     
      rt = xl                                                           
      assign 30 to nn                                                   
      go to 50                                                          
   30 z2 = fprt                                                         
      y2 = frt                                                          
c                                  begin main algorithm                 
   35 dd = one + d                                                      
      t1 = z0*d*d                                                       
      t2 = z1*dd*dd                                                     
      xx = z2*dd                                                        
      t3 = z2*d                                                         
      bi = t1-t2+xx+t3                                                  
      den = bi*bi-four*(xx*t1-t3*(t2-xx))                               
c                                  use denominator of maximum amplitude 
      t1 = cdsqrt(den)                                                   
      qz = rhun*max(cdabs(bi),cdabs(t1))                               
      t2 = bi + t1                                                      
      tpq = cdabs(t2)+qz                                                 
      if (tpq .eq. qz) t2 = zero                                        
      t3 = bi - t1                                                      
      tpq = cdabs(t3) + qz                                               
      if (tpq .eq. qz) t3 = zero                                        
      den = t2                                                          
      qz = cdabs(t3)-cdabs(t2)                                            
      if (qz .gt. rzero) den = t3                                       
c                                  test for zero denominator            
      assign 30 to nn                                                   
      if (cdabs(den) .eq. rzero) go to 65                               
      d = -xx/den                                                       
      d = d+d                                                           
      h = d*h                                                           
      rt = rt + h                                                       
c                                  check convergence of the first kind  
c                                                                       
      if (cdabs(h) .le. eps1*max(cdabs(rt),ax)) go to 70               
      if (ic .ne. 0) go to 15                                           
      assign 40 to nn                                                  
      go to 50                                                          
   40 qz = cdabs(fprt)-cdabs(z2)*rten                                    
      if (qz .ge. rzero) go to 45                                       
      z0 = z1                                                           
      z1 = z2                                                           
      z2 = fprt                                                         
      y0 = y1                                                           
      y1 = y2                                                           
      y2 = frt                                                          
      go to 35                                                          
c                                  take remedial action to induce       
c                                    convergence                        
   45 continue                                                          
      d = d*p5                                                          
      h = h*p5                                                          
      rt = rt-h                                                         
   50 jk = jk+1      
c      write(*,*)'jk,itmax',jk,itmax                                                   
      if (jk .gt. itmax) go to 75  
c      write(*,*)'rt',rt
c      write(*,*)'dum',dum  
cSAP110323                                   
c      frt = f(rt,dum)
      frt = f(rt)
c      write(*,*)'frt',frt 
      fprt = frt         
                                               
c                                  test to see if first root is being   
c                                     determined                        
      if (l .eq. 1) go to 60                                            
c                                  compute deflated function            
      lm1 = l-1                                                         
      do 55 i=1,lm1                                                     
         tem = rt - x(i)                                                
         if (cdabs(tem) .eq. rzero) go to 65                            
   55 fprt = fprt/tem                                                   
   60 continue                                                          
c                                  check convergence of the second kind 
c                                                                       
      if (cdabs(fprt) .le. eps .and. cdabs(frt) .le. eps) go to 80        
      go to nn,(20,25,30,40)                                            
   65 continue                                                          
      if (ic .ne. 0) go to 15                                           
      tem = rten*eps1                                                   
      if (cdabs(rt) .gt. ax) tem = tem*rt                                
      rt = rt+tem                                                       
      d = (h+tem)*d/h                                                   
      h = h+tem                                                         
      go to 50                                                          
c                                  check solution                       
   70 continue                                                          
      if (ic .ne. 0) go to 80                                           
      ic = 1                                                            
      z0 = y1                                                           
      z1 = y2                                                           
cSAP110322
c      z2 = f(rt,dum)
      z2 = f(rt)

      xl = rt                                                           
      ick = ick+1                                                       
      if (ick .le. ickmax) go to 35                                     
c                                  warning error, itmax = maximum       
      jk = itmax + jk                                                   
   75 ier = 33                                                          
c                                  a root has been found                
   80 x(l) = rt                                                         
C********* ROOT FOUND, SHOW MESSAGE
c      print *,'++++++MULLER ROOT: ',rt
      infer(l) = jk                                                     
      l = l+1                                                           
      if (l .le. knpn) go to 10                                         
      if (ier .eq. 0) go to 9005                                        
 9000 continue                                                          
 9005 return                                                            
      end                                                               
