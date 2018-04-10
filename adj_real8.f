CBH Notes: Sept, 2002
C          Search for "e-16" to find points where precision
C          was changed to "e-8", for purposes of getting code
C          running at 32-bit precision.
C          WASN'T ABLE TO GET THIS WORKING (GOT some NaN).
C          Should revert to e-16 for double precision.

CBH   Used replace real*8  real*8 < adj_cft_bh.f > adj_real8_bh.f
CBH   replaced 1.e-8 with 1.d-16 by hand.



C FTANGLE v1.53, created with UNIX on "Saturday, September 23, 1995 at 12:43." 
C  COMMAND LINE: "ftangle -mCRAY -=n=#.cft -=c=/dev/null adj"
C  RUN TIME: "Monday, March 4, 1996 at 11:45."
C  WEB FILE:    "adj.web"
C  CHANGE FILE: (none)
C* 14: * 
*line 477 "adj.web"
      program adjoint
      implicit none
      integer nmax,imax,kmax,tmax,lmax,tskip
      real*8 t,c2,ze,zi,eps,psivar,kappa,ipola,umax,dt,alpha,rho
      real*8 cpu,cpu0
cSm020909
c      external operinit,second,runa
	external operinit,runa
C* :14 * 
C* 14: * 
*line 490 "adj.web"
      call operinit
      write(6,*)'Enter '//'psivar,kappa,ipola'
C* :14 * 
C* 14: * 
*line 501 "adj.web"
      read(5,*,end=90001)psivar,kappa,ipola
      if(int(psivar+0.5d0).EQ.2)then
      ipola=8*atan(1.0d0)*29.1d0
      end if
90000 continue
      write(6,*)'Enter t,c2,ze,zi,eps,nmax,umax,imax,kmax,dt,tmax,alpha,
     &rho,lmax'
C* :14 * 
C* 14: * 
*line 518 "adj.web"
      read(5,*,end=90001)t,c2,ze,zi,eps,nmax,umax,imax,kmax,dt,tmax,alph
     &a,rho,lmax
      if(t.LE.0.OR.c2.LE.0)goto 90001
      write(6,*)t,c2,ze,zi,eps,nmax,umax,imax,kmax,dt,tmax,alpha,rho,lma
     &x
      write(6,*)psivar,kappa,ipola
      tskip=10
cSm020909
c      call second(cpu0)
      call runa(t,c2,ze,zi,eps,psivar,kappa,ipola,nmax,umax,imax,kmax,dt
     &,tmax,alpha,rho,lmax,tskip)
cSm020909
c      call second(cpu)
c      write(6,'(1x,a,f10.2,a)')'CPU time: ',cpu-cpu0,' secs'
      goto 90000
90001 continue
      stop
      end
C* 15: * 
*line 553 "adj.web"
C* :15 * 
C* 15: * 
*line 568 "adj.web"
C* :15 * 
C* 16: * 
*line 571 "adj.web"
C* :16 * 
C* 16: * 
*line 573 "adj.web"
      subroutine operinit
      implicit none
      open(unit=20,file='adjout',form='unformatted',status='unknown')
      return
      end
C* :16 * 
C* 16: * 
*line 582 "adj.web"
C* :16 * 
C* 17: * 
*line 585 "adj.web"
C* :17 * 
C* 17: * 
*line 611 "adj.web"
C* :17 * 
C* 18: * 
*line 614 "adj.web"
C* :18 * 
C* 18: * 
*line 640 "adj.web"
C* :18 * 
C* 24: * 
*line 756 "adj.web"
      subroutine runa(t,c2,ze,zi,eps,psivar,kappa,ipola,nmax,umax,imax,k
     &max,dt,tmax,alpha,rho,lmax,tskip)
      implicit none
      integer nmax,imax,kmax,tmax,lmax,tskip
      real*8 t,c2,ze,zi,eps,psivar,kappa,ipola,umax,dt,alpha,rho
      integer nmax1,imax1,lmaxa,lmaxa1
C* 25: * 
*line 813 "adj.web"
      integer memsize,worksize,pwork
C* 30: * 
*line 934 "adj.web"
      integer pua,pug,pth0a,pth0g,pcsth0a,pchi
C* :30 * 
C* 34: * 
*line 989 "adj.web"
      integer plegjy
C* :34 * 
C* 44: * 
*line 1150 "adj.web"
      integer pchil,plegp
C* :44 * 
C* 51: * 
*line 1264 "adj.web"
      integer pduug,pfu,pdtt,pfm,ppsid
C* :51 * 
C* 72: * 
*line 1843 "adj.web"
      integer plam1a,plam2g,pharr
C* :72 * 
C* 92: * 
*line 2261 "adj.web"
      integer puinv,pth0inv
C* :92 * 
C* 100: * 
*line 2424 "adj.web"
      integer presid,presida
C* :100 * 
*line 814 "adj.web"
      
C* :25 * 
C* 25: * 
*line 820 "adj.web"
C* :25 * 
C* 25: * 
*line 824 "adj.web"
C* :25 * 
C* 25: * 
*line 831 "adj.web"
C* :25 * 
C* 25: * 
*line 832 "adj.web"
cSm020909
      integer memory_max
      parameter (memory_max=5000000)
      real*8 memory(0:memory_max)
 
cSm020909
c      real*8 memory
c      dimension memory(0:*)
c      integer malloc 

      integer ibaseaddr
      real*8 baseaddr1
cSm020909
      real*8 baseaddr
cSm020909 Warning: because of COMMON the alignment of object is 
c inconsistent with its type   [BASEADDR1] real*8 baseaddr1 
      common /dptr/ baseaddr,baseaddr1
cSm020909
c      pointer(baseaddr,memory)
      integer status
c_cray      external hpalloc,hpdeallc
c_pc        external malloc,free
c      external malloc,free
C* :25 * 
C* 25: * 
*line 838 "adj.web"
C* :25 * 
C* 25: * 
*line 845 "adj.web"
C* :25 * 
C* 25: * 
*line 846 "adj.web"
C* :25 * 
*line 764 "adj.web"
      
      external runb
      
      nmax1=nmax
      imax1=imax
C* :24 * 
C* 24: * 
*line 776 "adj.web"
      if(mod(nmax1,2).EQ.1)nmax1=nmax1+1
      if(mod(imax1,2).EQ.0)imax1=imax1+1
C* :24 * 
C* 24: * 
*line 780 "adj.web"
      if(mod(lmax,2).EQ.0)lmax=lmax-1
      lmaxa=(lmax+1)/2
      lmaxa1=lmaxa
      if(lmaxa1.EQ.0)lmaxa1=1
      memsize=0
      worksize=0
C* 31: * 
*line 938 "adj.web"
      pua=memsize
      memsize=memsize+nmax
      
      pug=memsize
      memsize=memsize+nmax+1
      
      pth0a=memsize
      memsize=memsize+imax
      
      pth0g=memsize
      memsize=memsize+imax+1
      
      pcsth0a=memsize
      memsize=memsize+imax
      
      pchi=memsize
      memsize=memsize+imax1*nmax
      
C* :31 * 
C* 35: * 
*line 993 "adj.web"
      plegjy=memsize
      memsize=memsize+nmax1*7*4*(lmaxa1+1)
      
C* :35 * 
C* 38: * 
*line 1033 "adj.web"
      worksize=max(worksize,2*(max(2*lmaxa-1,0)+2+1)*7*2)
C* :38 * 
C* 41: * 
*line 1102 "adj.web"
      worksize=max(worksize,nmax+(nmax1+1)*7*2)
C* :41 * 
C* 45: * 
*line 1154 "adj.web"
      pchil=memsize
      memsize=memsize+nmax1*lmaxa1
      
      plegp=memsize
      memsize=memsize+imax1*lmaxa1
      
C* :45 * 
C* 48: * 
*line 1198 "adj.web"
      worksize=max(worksize,imax)
C* :48 * 
C* 52: * 
*line 1268 "adj.web"
      pduug=memsize
      memsize=memsize+nmax+1
      
      pfu=memsize
      memsize=memsize+nmax
      
      pdtt=memsize
      memsize=memsize+nmax
      
      pfm=memsize
      memsize=memsize+nmax
      
      ppsid=memsize
      memsize=memsize+nmax1*5*2
      
C* :52 * 
C* 58: * 
*line 1459 "adj.web"
      worksize=max(worksize,nmax*2+nmax+(nmax1+1)*7*2)
C* :58 * 
C* 73: * 
*line 1847 "adj.web"
      plam1a=memsize
      memsize=memsize+imax
      
      plam2g=memsize
      memsize=memsize+imax+1
      
      pharr=memsize
      memsize=memsize+(((lmaxa1-1)*lmaxa1*(2*lmaxa1-1))/6+(lmaxa1-1)*lma
     &xa1+lmaxa1)
      
C* :73 * 
C* 76: * 
*line 1905 "adj.web"
      worksize=max(worksize,2*kmax+((((lmaxa1+3)-1)*(lmaxa1+3))/2+(lmaxa
     &1+1)))
C* :76 * 
C* 90: * 
*line 2250 "adj.web"
      worksize=max(worksize,3*imax)
C* :90 * 
C* 93: * 
*line 2265 "adj.web"
      puinv=memsize
      memsize=memsize+nmax1*3
      
      pth0inv=memsize
      memsize=memsize+imax1*3
      
C* :93 * 
C* 101: * 
*line 2428 "adj.web"
      presid=memsize
      memsize=memsize+imax1*nmax
      
      presida=memsize
      memsize=memsize+imax1*nmax
      
C* :101 * 
C* 101: * 
*line 2432 "adj.web"
      worksize=max(worksize,imax1*nmax)
C* :101 * 
C* 101: * 
*line 2434 "adj.web"
C* :101 * 
C* 109: * 
*line 2574 "adj.web"
      worksize=max(worksize,imax)
C* :109 * 
*line 793 "adj.web"
      
      pwork=memsize
      memsize=memsize+worksize
      
C* 26: * 
*line 872 "adj.web"
C* :26 * 
C* 26: * 
*line 879 "adj.web"
c_cray      call hpalloc(baseaddr,1*memsize,status,0)
c_cray      if(status.NE.0)then
c_cray      write(6,*)'Couldn''t allocate memory: ',memsize
c_cray      return
c_cray      end if
c_pc     Using unix library libU77.a
c_pc     baseaddr=malloc(8*memsize)
c_pc     if (baseaddr.eq.0) stop 'failed malloc in runa'

cBH  FOLLOWING LINE IS FOR 8 BYTE WORDS.  PROBABLY 4* OK HERE,
cBH  (OR, 1*, see above) UNLESS RESET real*8 ==> real*8*8
cSm020909
c      baseaddr=malloc(8*memsize)
cSm020909
      baseaddr=8*memsize !it can be redused
      
      baseaddr1=baseaddr
      ibaseaddr=baseaddr
      write(*,*)'runa: baseaddr,ibaseaddr',baseaddr,ibaseaddr
      if (baseaddr.eq.0) stop 'failed malloc in runa'
C* :26 * 
C* 26: * 
*line 885 "adj.web"
C* :26 * 
*line 796 "adj.web"
      write(*,*)'in runa memsize,memory_max',memsize,memory_max
      call runb(t,c2,ze,zi,eps,psivar,kappa,ipola,nmax,umax,imax,kmax,dt
     &,tmax,alpha,rho,lmax,tskip,nmax1,imax1,lmaxa,lmaxa1,worksize,memor
     &y(pua),memory(pug),memory(pth0a),memory(pth0g),memory(pcsth0a),mem
     &ory(pchi),memory(plegjy),memory(pchil),memory(plegp),memory(pduug)
     &,memory(pfu),memory(pdtt),memory(pfm),memory(ppsid),memory(plam1a)
     &,memory(plam2g),memory(pharr),memory(puinv),memory(pth0inv),memory
     &(presid),memory(presida),memory(pwork))
C* 27: * 
*line 887 "adj.web"
C* :27 * 
C* 27: * 
*line 889 "adj.web"
c_cray      call hpdeallc(baseaddr,status,0)
c_pc        call free(baseaddr)

cSm020909
c      call free(baseaddr)

C* :27 * 
C* 27: * 
*line 891 "adj.web"
C* :27 * 
*line 806 "adj.web"
      
      return
      end
C* :24 * 
C* 28: * 
*line 896 "adj.web"
      subroutine runb(t,c2,ze,zi,eps,psivar,kappa,ipola,nmax,umax,imax,k
     &max,dt,tmax,alpha,rho,lmax,tskip,nmax1,imax1,lmaxa,lmaxa1,worksize
     &,ua,ug,th0a,th0g,csth0a,chi,legjy,chil,legp,duug,fu,dtt,fm,psid,la
     &m1a,lam2g,harr,uinv,th0inv,resid,resida,work)
      implicit none
      integer worksize
      real*8 work
      dimension work(0:worksize-1)
C* 29: * 
*line 925 "adj.web"
      integer nmax,imax,nmax1,imax1,n,i
      real*8 du,umax,ua,ug,dth0,th0max,th0a,th0g,csth0a,chi
      dimension ua(0:nmax-1),ug(0:nmax),th0a(0:imax-1),th0g(0:imax),csth
     &0a(0:imax-1),chi(0:imax1-1,0:nmax-1)
      real*8 pi
C* :29 * 
C* 33: * 
*line 981 "adj.web"
      real*8 c,c2
      real*8 legjy
      integer lmax,lmaxa,lmaxa1
      dimension legjy(0:nmax1-1,0:6,0:3,0:lmaxa1)
      external initjy
C* :33 * 
C* 43: * 
*line 1144 "adj.web"
      real*8 chil,legp
      dimension chil(0:nmax1-1,1:lmaxa1),legp(0:imax1-1,1:lmaxa1)
      external initp
C* :43 * 
C* 50: * 
*line 1249 "adj.web"
      real*8 t,ze,zi
      real*8 duug,fu,dtt
      dimension duug(0:nmax),fu(0:nmax-1),dtt(0:nmax-1)
      real*8 fm
      dimension fm(0:nmax-1)
      real*8 psid
      dimension psid(0:nmax1-1,0:4,0:1)
      external maxwel,isotrop
C* :50 * 
C* 60: * 
*line 1492 "adj.web"
      real*8 psivar,kappa,ipola
      real*8 eps
C* :60 * 
C* 71: * 
*line 1837 "adj.web"
      integer kmax
      real*8 lam1a,lam2g,harr
      dimension lam1a(0:imax-1),lam2g(0:imax),harr((((lmaxa1-1)*lmaxa1*(
     &2*lmaxa1-1))/6+(lmaxa1-1)*lmaxa1+lmaxa1))
C* :71 * 
C* 91: * 
*line 2255 "adj.web"
      real*8 uinv,th0inv,rho
      dimension uinv(0:nmax1-1,-1:1),th0inv(0:imax1-1,-1:1)
      external matinit
C* :91 * 
C* 99: * 
*line 2408 "adj.web"
      integer tmax,time,tskip
      real*8 resid,abserr,relerr,maxerr,alpha,dt,resida
      dimension resid(0:imax1-1,0:nmax-1),resida(0:imax1-1,0:nmax-1)
C* :99 * 
C* 99: * 
*line 2413 "adj.web"
C$$$      external hdfout
C* :99 * 
C* 99: * 
*line 2415 "adj.web"
      external decomp,reactall,residue,error
C* :99 * 
C* 99: * 
*line 2417 "adj.web"
      integer ifail
      external md03uaf
C* :99 * 
C* 99: * 
*line 2422 "adj.web"
C* :99 * 
C* 106: * 
*line 2521 "adj.web"
      real*8 curv,current,cond,e
      external current
C* :106 * 
*line 909 "adj.web"
      
C* 32: * 
*line 948 "adj.web"
      du=umax/nmax
      do n=0,nmax-1
      ua(n)=(n+0.5d0)*du
      ug(n)=n*du
      end do
      ug(nmax)=umax
      pi=4*atan(1.0d0)
cBH      th0max=atan2(1.0d0,sqrt(2*eps/(1-eps)))
      th0max=atan2(1.0d0,sqrt(2*eps/(1-eps)))
      dth0=th0max/imax
      do i=0,imax-1
      th0a(i)=(i+0.5d0)*dth0
      csth0a(i)=cos(th0a(i))
      th0g(i)=i*dth0
      end do
      th0g(imax)=th0max
      do n=0,nmax-1
      do i=0,imax-1
      chi(i,n)=0
      end do
      end do
C* :32 * 
C* 36: * 
*line 997 "adj.web"
      c=sqrt(c2)
      lmaxa=(lmax+1)/2
      call initjy(ua,c,nmax,nmax1,lmaxa,lmaxa1,legjy,work(0),work(2*(max
     &(2*lmaxa-1,0)+2+1)*7),max(2*lmaxa-1,0)+2)
C* :36 * 
C* 46: * 
*line 1160 "adj.web"
      call initp(csth0a,imax,imax1,lmaxa,lmaxa1,legp,work)
C* :46 * 
C* 53: * 
*line 1277 "adj.web"
      call maxwel(ua,du,c,nmax,fm,t)
      do n=0,nmax-1
      chil(n,1)=1
      end do
      call isotrop(ua,du,c,fm,chil(0,1),nmax,nmax1,legjy(0,0,0,0),psid,z
     &e,zi,duug,fu,dtt,work)
C* :53 * 
C* 74: * 
*line 1854 "adj.web"
      call maginit(eps,psivar,kappa,ipola,th0a,th0g,lam1a,lam2g,harr,ima
     &x,lmaxa1,kmax,work(0),work(kmax),work(2*kmax))
C* :74 * 
C* 94: * 
*line 2271 "adj.web"
      call matinit(du,ua,ug,nmax,nmax1,dth0,th0a,th0g,imax,imax1,duug,fu
     &,dtt,lam1a,lam2g,uinv,th0inv,rho)
C* :94 * 
*line 911 "adj.web"
      
      do time=0,tmax
C* 102: * 
*line 2437 "adj.web"
      call decomp(chi,nmax,nmax1,imax,imax1,lmaxa,lmaxa1,th0a,dth0,legp,
     &chil)
      call reactall(ua,du,c,fm,t,chil,nmax,nmax1,imax,imax1,lam1a,legjy(
     &0,0,0,1),psid,legp,harr,ze,lmaxa,lmaxa1,resid,work(0),work(nmax),w
     &ork(2*nmax))
C* :102 * 
C* 103: * 
*line 2445 "adj.web"
C* :103 * 
C* 103: * 
*line 2451 "adj.web"
      call residue(ua,nmax,nmax1,csth0a,imax,imax1,chi,dtt,c,lam1a,uinv,
     &th0inv,resid)
C* :103 * 
C* 103: * 
*line 2453 "adj.web"
      if(mod(time,tskip).EQ.0)then
      call error(du,ua,umax,nmax,dth0,th0a,th0max,imax,imax1,resid,chi,a
     &bserr,relerr,maxerr,work(0),work(imax),work(2*imax))
C* :103 * 
C* 103: * 
*line 2457 "adj.web"
C$$$      call hdfout(resid,chi,du,umax,nmax,dth0,imax,imax1,time)
C* :103 * 
C* 103: * 
*line 2459 "adj.web"
      end if
C* :103 * 
C* 104: * 
*line 2463 "adj.web"
C* :104 * 
C* 104: * 
*line 2465 "adj.web"
      ifail=0
      call md03uaf(imax,nmax,imax1,th0inv(0,-1),th0inv(0,0),th0inv(0,1),
     &uinv(0,-1),uinv(0,0),uinv(0,1),-dt,alpha,time,resid,resida,work,if
     &ail)
      if(ifail.NE.0)then
      write(6,*)'md03uaf failed: ',ifail,time
      end if
C* :104 * 
C* 104: * 
*line 2478 "adj.web"
      do n=0,nmax-1
      do i=0,imax-1
      chi(i,n)=chi(i,n)+dt*resid(i,n)
      end do
      end do
C* :104 * 
C* 107: * 
*line 2531 "adj.web"
      if(mod(time,tskip).EQ.0)then
      curv=current(du,ua,c,nmax,dth0,th0a,imax,imax1,fm,chi,work)
      e=t/harr((((1-1)*1*(2*1-1))/6+(1-1)*1+1))
      cond=curv*zi/(t*sqrt(t)*e)
      write(6,*)time,cond,abserr,relerr,maxerr
      end if
C* :107 * 
*line 914 "adj.web"
      
      end do
C* 111: * 
*line 2618 "adj.web"
      write(20)nmax,imax,kmax,lmax
      write(20)real(t),real(c2),real(ze),real(zi),real(eps)
      write(20)real(umax),real(th0max),real(du),real(dth0),real(rho)
      write(20)real(dt),tmax,real(alpha)
      write(20)real(cond),real(abserr),real(relerr)
      write(20)(real(ua(n)),n=0,nmax-1)
      write(20)(real(ug(n)),n=0,nmax)
      write(20)(real(th0a(n)),n=0,imax-1)
      write(20)(real(th0g(n)),n=0,imax)
      write(20)(real(fm(n)),n=0,nmax-1)
      write(20)((real(chi(i,n)),i=0,imax-1),n=0,nmax-1)
C* :111 * 
*line 917 "adj.web"
      
      return
      end
C* :28 * 
C* 37: * 
*line 1004 "adj.web"
      subroutine initjy(ua,c,nmax,nmax1,lmaxa,lmaxa1,legjy,jarray,djarra
     &y,lmax1)
      implicit none
      integer nmax,nmax1,lmaxa,lmaxa1,n,mult,la,l,a,lmax1
      real*8 ua,legjy,c,jarray,djarray
      dimension ua(0:nmax-1),legjy(0:nmax1-1,0:6,0:3,0:lmaxa1)
      dimension jarray(-lmax1-1:lmax1,0:6),djarray(-lmax1-1:lmax1,0:6)
      external legjn
      do n=0,nmax-1
      call legjn(ua(n),c,max(2*lmaxa-1,0),jarray,djarray,lmax1)
      do la=0,lmaxa
      l=max(max(2*la-1,0),0)
      mult=(-1)**(l+1)
      do a=0,6
      legjy(n,a,0,la)=jarray(l,a)
      legjy(n,a,1,la)=mult*jarray(-l-1,a)
      legjy(n,a,2,la)=djarray(l,a)
      legjy(n,a,3,la)=mult*djarray(-l-1,a)
      end do
      end do
      end do
      return
      end
C* :37 * 
C* 40: * 
*line 1056 "adj.web"
      subroutine pot(ua,du,c,fm,chil,nmax,nmax1,legjy,psid,fla,ijy)
      implicit none
      integer p0,p02,p022,p1,p11
      parameter(p0=0,p02=1,p022=2,p1=3,p11=4)
      integer j0,j1,j2,j02,j11,j22,j022
      parameter(j0=0,j1=1,j2=2,j02=3,j11=4,j22=5,j022=6)
      integer nmax,nmax1,n,a
      real*8 ua,du,c,fm,chil,legjy
      dimension ua(0:nmax-1),fm(0:nmax-1),chil(0:nmax-1),legjy(0:nmax1-1
     &,0:6,0:3)
      real*8 psid
      dimension psid(0:nmax1-1,0:4,0:1)
      real*8 fla,ijy
      dimension fla(0:nmax-1),ijy(0:nmax1,0:6,0:1)
C* 42: * 
*line 1109 "adj.web"
      do n=0,nmax-1
      fla(n)=fm(n)*chil(n)*ua(n)**2/sqrt(1+(ua(n)/c)**2)*du/2
      end do
      do a=0,6
      do n=0,nmax-1
      ijy(n+1,a,0)=fla(n)*legjy(n,a,0)
      ijy(n,a,1)=fla(n)*legjy(n,a,1)
      end do
      ijy(0,a,0)=0
      ijy(nmax,a,1)=0
      end do
      do n=1,nmax
      do a=0,6
      ijy(n,a,0)=ijy(n-1,a,0)+ijy(n,a,0)
      end do
      end do
      do n=nmax-1,0,-1
      do a=0,6
      ijy(n,a,1)=ijy(n+1,a,1)+ijy(n,a,1)
      end do
      end do
      do a=0,6
      do n=0,nmax-1
      ijy(n,a,0)=ijy(n,a,0)+ijy(n+1,a,0)
      ijy(n,a,1)=ijy(n,a,1)+ijy(n+1,a,1)
      end do
      end do
C* :42 * 
*line 1073 "adj.web"
      
      do n=0,nmax-1
      psid(n,p0,0)=legjy(n,j0,1)*ijy(n,j0,0)+ijy(n,j0,1)*legjy(n,j0,0)
      psid(n,p02,0)=legjy(n,j0,1)*ijy(n,j02,0)+legjy(n,j02,1)*ijy(n,j2,0
     &)+ijy(n,j0,1)*legjy(n,j02,0)+ijy(n,j02,1)*legjy(n,j2,0)
      psid(n,p022,0)=legjy(n,j0,1)*ijy(n,j022,0)+legjy(n,j02,1)*ijy(n,j2
     &2,0)+legjy(n,j022,1)*ijy(n,j2,0)+ijy(n,j0,1)*legjy(n,j022,0)+ijy(n
     &,j02,1)*legjy(n,j22,0)+ijy(n,j022,1)*legjy(n,j2,0)
      psid(n,p1,0)=legjy(n,j1,1)*ijy(n,j1,0)+ijy(n,j1,1)*legjy(n,j1,0)
      psid(n,p11,0)=legjy(n,j1,1)*ijy(n,j11,0)+legjy(n,j11,1)*ijy(n,j1,0
     &)+ijy(n,j1,1)*legjy(n,j11,0)+ijy(n,j11,1)*legjy(n,j1,0)
      psid(n,p0,1)=legjy(n,j0,3)*ijy(n,j0,0)+ijy(n,j0,1)*legjy(n,j0,2)
      psid(n,p02,1)=legjy(n,j0,3)*ijy(n,j02,0)+legjy(n,j02,3)*ijy(n,j2,0
     &)+ijy(n,j0,1)*legjy(n,j02,2)+ijy(n,j02,1)*legjy(n,j2,2)
      psid(n,p022,1)=legjy(n,j0,3)*ijy(n,j022,0)+legjy(n,j02,3)*ijy(n,j2
     &2,0)+legjy(n,j022,3)*ijy(n,j2,0)+ijy(n,j0,1)*legjy(n,j022,2)+ijy(n
     &,j02,1)*legjy(n,j22,2)+ijy(n,j022,1)*legjy(n,j2,2)
      psid(n,p1,1)=legjy(n,j1,3)*ijy(n,j1,0)+ijy(n,j1,1)*legjy(n,j1,2)
      psid(n,p11,1)=legjy(n,j1,3)*ijy(n,j11,0)+legjy(n,j11,3)*ijy(n,j1,0
     &)+ijy(n,j1,1)*legjy(n,j11,2)+ijy(n,j11,1)*legjy(n,j1,2)
      end do
      return
      end
C* :40 * 
C* 47: * 
*line 1166 "adj.web"
      subroutine initp(x,imax,imax1,lmaxa,lmaxa1,legp,leg0)
      implicit none
      integer imax,imax1,lmaxa,lmaxa1,i,l,la
      real*8 x,legp
      dimension x(0:imax-1),legp(0:imax1-1,1:lmaxa1)
      real*8 leg0
      dimension leg0(0:imax-1)
      if(lmaxa.LT.1)return
      la=1
      do i=0,imax-1
      leg0(i)=1
      legp(i,la)=x(i)
      end do
      do la=2,lmaxa
      l=max(2*la-1,0)
      do i=0,imax-1
      leg0(i)=((2*l-3)*x(i)*legp(i,la-1)-(l-2)*leg0(i))/(l-1)
      legp(i,la)=((2*l-1)*x(i)*leg0(i)-(l-1)*legp(i,la-1))/l
      end do
      end do
      return
      end
C* :47 * 
C* 49: * 
*line 1213 "adj.web"
      subroutine decomp(chi,nmax,nmax1,imax,imax1,lmaxa,lmaxa1,th0a,dth0
     &,legp,chil)
      implicit none
      integer nmax,nmax1,imax,imax1,lmaxa,lmaxa1,la,l,n,i
      real*8 chi,th0a,dth0,legp,chil,sn
      dimension chi(0:imax1-1,0:nmax-1),th0a(0:imax-1),legp(0:imax1-1,1:
     &lmaxa1),chil(0:nmax1-1,1:lmaxa1)
      do la=1,lmaxa
      l=max(2*la-1,0)
      do n=0,nmax-1
      chil(n,la)=0
      end do
      do i=0,imax-1
      sn=sin(th0a(i))*legp(i,la)
      do n=0,nmax-1
      chil(n,la)=chil(n,la)+chi(i,n)*sn
      end do
      end do
      do n=0,nmax-1
      chil(n,la)=(2*l+1)*dth0*chil(n,la)
      end do
      end do
      return
      end
C* :49 * 
C* 54: * 
*line 1290 "adj.web"
      subroutine isotrop(ua,du,c,fm,chil,nmax,nmax1,legjy,psid,ze,zi,duu
     &g,fu,dtt,work)
      implicit none
      integer p0,p02,p022,p1,p11
      parameter(p0=0,p02=1,p022=2,p1=3,p11=4)
      integer nmax,nmax1,n
      real*8 ua,c,psid,duug,fu,dtt,pi,g,u,ze,zi,du,fm,chil,legjy
      dimension ua(0:nmax-1),fm(0:nmax-1),chil(0:nmax-1),legjy(0:nmax1-1
     &,0:6,0:3),psid(0:nmax1-1,0:4,0:1),duug(0:nmax),fu(0:nmax-1),dtt(0:
     &nmax-1)
      real*8 work
      dimension work(0:nmax+(nmax1+1)*7*2-1)
      external pot
      call pot(ua,du,c,fm,chil,nmax,nmax1,legjy,psid,work(0),work(nmax))
      pi=4*atan(1.0d0)
      do n=0,nmax-1
      u=ua(n)
      g=sqrt(1+(u/c)**2)
      duug(n)=ze*4*pi*g/u*(2*g**2*psid(n,p02,1)-u*psid(n,p0,0)-8*g**2/c*
     &*2*psid(n,p022,1)+8*u/c**4*psid(n,p022,0))
      dtt(n)=ze*4*pi/(g*u)*(-g**2*psid(n,p02,1)-u/c**2*psid(n,p02,0)+4*g
     &**2/c**2*psid(n,p022,1)-4*u/c**4*psid(n,p022,0))+zi*g/(2*u)
      fu(n)=ze*4*pi*g*(-psid(n,p1,1)+2/c**2*psid(n,p11,1))
      end do
      
      duug(nmax)=(3*duug(nmax-1)-duug(nmax-2))/2
      do n=nmax-1,1,-1
      duug(n)=(duug(n)+duug(n-1))/2
      end do
      duug(0)=duug(0)
      return
      end
C* :54 * 
C* 55: * 
*line 1338 "adj.web"
      subroutine maxwel(ua,du,c,nmax,fm,t)
      implicit none
      integer nmax,n
      real*8 ua,du,c,fm,t,sum,pi
      dimension ua(0:nmax-1),fm(0:nmax-1)
      sum=0
      do n=0,nmax-1
      fm(n)=exp(-ua(n)**2/(sqrt(1+(ua(n)/c)**2)+1)/t)
      sum=sum+ua(n)**2*fm(n)
      end do
      pi=4*atan(1.0d0)
      sum=4*pi*du*sum
      do n=0,nmax-1
      fm(n)=fm(n)/sum
      end do
      return
      end
C* :55 * 
C* 56: * 
*line 1369 "adj.web"
      subroutine reaction(ua,du,c,fm,t,chil,nmax,nmax1,legjy,psid,ze,rea
     &ctl,l,work)
      implicit none
      integer p0,p02,p022,p1,p11
      parameter(p0=0,p02=1,p022=2,p1=3,p11=4)
      integer nmax,nmax1,n,l
      real*8 ua,du,c,legjy,psid,ze,reactl,fm,chil,t,pi,u,g
      dimension ua(0:nmax-1),legjy(0:nmax1-1,0:6,0:3),psid(0:nmax1-1,0:4
     &,0:1),reactl(0:nmax-1),fm(0:nmax-1),chil(0:nmax-1)
      real*8 work
      dimension work(0:nmax+(nmax1+1)*7*2-1)
      external pot
      call pot(ua,du,c,fm,chil,nmax,nmax1,legjy,psid,work(0),work(nmax))
      pi=4*atan(1.0d0)
      do n=0,nmax-1
      u=ua(n)
      g=sqrt(1+(u/c)**2)
      reactl(n)=reactl(n)+ze*4*pi*(fm(n)*chil(n)/g-u/t*psid(n,p1,1)-2/(c
     &**2*g)*psid(n,p1,0)+2*u/(c**2*t)*psid(n,p11,1)+u/t*psid(n,p0,1)-(u
     &**2/(g*t**2)-1/t)*psid(n,p0,0)+(2*g*u/t**2-2*u/(c**2*t))*psid(n,p0
     &2,1)-(l*(l+1)/(g*t**2)-2/(c**2*t))*psid(n,p02,0)-8*g*u/(c**2*t**2)
     &*psid(n,p022,1)+4*((l+2)*(l-1)/g+2*g)/(c**2*t**2)*psid(n,p022,0))
      end do
      return
      end
C* :56 * 
C* 57: * 
*line 1417 "adj.web"
      subroutine reactall(ua,du,c,fm,t,chil,nmax,nmax1,imax,imax1,lam1a,
     &legjy,psid,legp,harr,ze,lmaxa,lmaxa1,react,reactl,chila,work)
      implicit none
      integer nmax,nmax1,imax,imax1,lmaxa,lmaxa1
      real*8 ua,du,c,fm,t,chil,lam1a,legjy,psid,legp,harr,ze,react
      dimension ua(0:nmax-1),fm(0:nmax-1),chil(0:nmax1-1,1:lmaxa1),lam1a
     &(0:imax-1),legjy(0:nmax1-1,0:6,0:3,1:lmaxa1),psid(0:nmax1-1,0:4,0:
     &1),legp(0:imax1-1,1:lmaxa1),harr((((lmaxa1-1)*lmaxa1*(2*lmaxa1-1))
     &/6+(lmaxa1-1)*lmaxa1+lmaxa1)),react(0:imax1-1,0:nmax-1)
      integer la,la1,la2,n,i
      real*8 reactl,chila,work
      dimension reactl(0:nmax-1),chila(0:nmax-1),work(0:nmax+(nmax1+1)*7
     &*2-1)
      external reaction
      do n=0,nmax-1
      do i=0,imax-1
      react(i,n)=0
      end do
      end do
      do la1=1,lmaxa
C* 59: * 
*line 1464 "adj.web"
      do n=0,nmax-1
      reactl(n)=0
      end do
      do la=la1,lmaxa
      do n=0,nmax-1
      chila(n)=0
      end do
      do la2=1,la
      do n=0,nmax-1
      chila(n)=chila(n)+harr((((la-1)*la*(2*la-1))/6+(la1-1)*la+la2))*ch
     &il(n,la2)
      end do
      end do
      call reaction(ua,du,c,fm,t,chila,nmax,nmax1,legjy(0,0,0,la),psid,z
     &e,reactl,max(2*la-1,0),work)
      end do
C* :59 * 
*line 1439 "adj.web"
      
      do n=0,nmax-1
      do i=0,imax-1
      react(i,n)=react(i,n)+reactl(n)*legp(i,la1)
      end do
      end do
      end do
      do n=0,nmax-1
      do i=0,imax-1
      react(i,n)=react(i,n)/lam1a(i)
      end do
      end do
      return
      end
C* :57 * 
C* 61: * 
*line 1499 "adj.web"
      subroutine psipolf(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipo
     &lr)
      implicit none
      real*8 dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipolr
      external psipola,psipolb,psipolc
      if(int(psivar+0.5d0).EQ.0)then
      call psipola(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipolr)
      else if(int(psivar+0.5d0).EQ.1)then
      call psipolb(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipolr)
      else
      call psipolc(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipolr)
      end if
      return
      end
C* :61 * 
C* 62: * 
*line 1518 "adj.web"
      subroutine psipola(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipo
     &lr)
      implicit none
      real*8 dr,z,psivar,kappa,ipola,r0,psipol1,psipol,psipolz,psipolr
      psipol1=1
      psipol=-psipol1*(dr**2+(z/kappa)**2)
      psipolz=-psipol1*2*z/kappa**2
      psipolr=-psipol1*2*dr
      return
      end
C* :62 * 
C* 63: * 
*line 1534 "adj.web"
      subroutine psipolb(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipo
     &lr)
      implicit none
cSm0209009
c      real*8 dr,z,psivar,kappa,ipola,r0,a,alpha,r,psipol1,psipol,psipolz,p
c     &sipolr
	real*8 dr,z,psivar,kappa,ipola,r0,a,alpha,r,psipol1,psipol,psipolz
     &,psipolr

      psipol1=1
      a=0
      alpha=kappa**2*(a**2-r0**2)/(4*r0**2)
      r=r0+dr
      psipol=-psipol1*((a**2-r**2)*z**2+alpha*(dr*(r0+r))**2)
      psipolz=-psipol1*2*(a**2-r**2)*z
      psipolr=-psipol1*2*(-r*z**2+2*alpha*r*dr*(r0+r))
      return
      end
C* :63 * 
C* 64: * 
*line 1554 "adj.web"
      subroutine psipolc(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipo
     &lr)
      implicit none
      real*8 dr,z,psivar,kappa,ipola,r0,psipol,psipolr,psipolz
      real*8 rv,zv,numer,numerr,numerz,denom,denomr,denomz,temp,tempr
      integer k,m,l
      integer nnum,nden
      parameter(nnum=3,nden=2)
C* 65: * 
*line 1607 "adj.web"
      real*8 ac
      dimension ac(((nnum+1)*(nnum+2))/2+((nden+1)*(nden+2))/2)
      save ac
      data ac/0.000510d0,-0.003820d0,-0.029876d0,-0.003174d0,0.142820d0,
     &5.484722d0,-0.000972d0,0.099224d0,0.000000d0,0.000000d0,0.001580d0
     &,0.000021d0,-0.015488d0,0.000110d0,-0.009744d0,1.000000d0/
C* :65 * 
*line 1561 "adj.web"
      
      rv=dr*(dr+2*r0)
      zv=z**2
      numer=0
      numerr=0
      numerz=0
      k=0
      do m=nnum,0,-1
      temp=0
      tempr=0
      do l=nnum-m,0,-1
      k=k+1
      temp=temp*rv+ac(k)
      if(l.GT.0)tempr=tempr*rv+l*ac(k)
      end do
      numer=numer*zv+temp
      numerr=numerr*zv+tempr
      if(m.GT.0)numerz=numerz*zv+m*temp
      end do
      denom=0
      denomr=0
      denomz=0
      do m=nden,0,-1
      temp=0
      tempr=0
      do l=nden-m,0,-1
      k=k+1
      temp=temp*rv+ac(k)
      if(l.GT.0)tempr=tempr*rv+l*ac(k)
      end do
      denom=denom*zv+temp
      denomr=denomr*zv+tempr
      if(m.GT.0)denomz=denomz*zv+m*temp
      end do
      psipol=numer/denom
      psipolr=2*(dr+r0)*(numerr-numer*denomr/denom)/denom
      psipolz=2*z*(numerz-numer*denomz/denom)/denom
      return
      end
C* :64 * 
C* 66: * 
*line 1620 "adj.web"
      function ipol(psivar,kappa,ipola,psipol)
      real*8 ipol,psivar,kappa,ipola,psipol
      ipol=ipola
      return
      end
C* :66 * 
C* 67: * 
*line 1631 "adj.web"
      function magfield(dr,z,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipo
     &l)
      implicit none
cSm020909
c      real*8 magfield,ipol,dr,z,psivar,kappa,ipola,r0,bz,br,bzeta,dbzeta,b
c     &0,psipol,psipolz,psipolr,pi,l
	real*8 magfield,ipol,dr,z,psivar,kappa,ipola,r0,bz,br,bzeta,dbzeta
     &,b0,psipol,psipolz,psipolr,pi,l
      external psipolf,ipol
      call psipolf(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipolr)
      pi=4*atan(1.0d0)
      l=2*pi*(r0+dr)
      bz=psipolr/l
      br=-psipolz/l
      b0=ipol(psivar,kappa,ipola,psipol)/(2*pi*r0)
      dbzeta=-dr/(r0+dr)*b0
      bzeta=b0+dbzeta
      magfield=(bz**2+br**2+dbzeta*(b0+bzeta))/(sqrt(bz**2+br**2+bzeta**
     &2)+b0)
      return
      end
C* :67 * 
C* 68: * 
*line 1656 "adj.web"
      function newt(thpol,psipol0,rpol1,psivar,kappa,ipola,r0,drc,zc)
      implicit none
cSm020909
c      real*8 newt,thpol,psipol0,rpol1,psivar,kappa,ipola,r0,drc,zc,drpol,p
c     &sipol,psipolz,psipolr,dr,z,rpol,rteps
	real*8 newt,thpol,psipol0,rpol1,psivar,kappa,ipola,r0,drc,zc,drpol
     &,psipol,psipolz,psipolr,dr,z,rpol,rteps
      integer i,jj
      logical conv
      external psipolf
      if(rpol1.EQ.0)then
      newt=0
      return
      end if
cBH      rteps=sqrt(1.0e-16)
cBH      rteps=sqrt(1.0e-8)
      rteps=sqrt(1.0d-16)
      conv=.FALSE.
      jj=0
      rpol=rpol1
      do i=1,10
      dr=drc+rpol*cos(thpol)
      z=zc+rpol*sin(thpol)
      call psipolf(dr,z,psivar,kappa,ipola,r0,psipol,psipolz,psipolr)
      drpol=-(psipol-psipol0)/(psipolr*cos(thpol)+psipolz*sin(thpol))
      rpol=rpol+drpol
      if(rpol.LT.0)then
      write(6,*)'newt: no convergence: drpol, rpol',drpol,rpol
      newt=0
      end if
      if(conv)then
      jj=jj+1
      if(jj.EQ.3)then
      newt=rpol
      return
      end if
      else
      conv=abs(drpol).LT.rteps*rpol
      end if
      end do
      write(6,*)'newt: no convergence: drpol, rpol',drpol,rpol
      newt=rpol
      return
      end
C* :68 * 
C* 69: * 
*line 1711 "adj.web"
      subroutine psiinit(eps,psivar,kappa,ipola,r0,psipol0,dra,za,drc,zc
     &,b0,dbmin,dbmax,epsg)
      implicit none
cSm020909
c      real*8 eps,psivar,kappa,ipola,r0,psipol0,dra,za,drc,zc,b0,dbmin,dbma
c     &x,bz,br,dbzeta,magfield,pi,dr1,dr2,drb,eps1,eps2,newt,z0,epsg,dumm
c     &y
      real*8 eps,psivar,kappa,ipola,r0,psipol0,dra,za,drc,zc,b0,dbmin,
     &dbmax,bz,br,dbzeta,magfield,pi,dr1,dr2,drb,eps1,eps2,newt,z0,epsg
     &,dummy
      integer i
      external magfield,newt
      pi=4*atan(1.0d0)
      z0=0.0d0
      if(eps.EQ.0)then
      dra=0
      za=z0
      dbmin=magfield(dra,za,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      dbmax=dbmin
      drc=dra
      zc=za
      epsg=0
      return
      end if
      dr1=r0*1.2d0*eps
      dbmin=magfield(dr1,z0,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      drb=-newt(pi,psipol0,dr1,psivar,kappa,ipola,r0,0.0d0,0.0d0)
      dbmax=magfield(drb,z0,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      eps1=(dbmax-dbmin)/(2*b0+dbmax+dbmin)-eps
      dr2=r0*0.8d0*eps
      do i=1,50
      dbmin=magfield(dr2,z0,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      drb=-newt(pi,psipol0,-drb,psivar,kappa,ipola,r0,0.0d0,0.0d0)
      dbmax=magfield(drb,z0,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      eps2=(dbmax-dbmin)/(2*b0+dbmax+dbmin)-eps
      if(abs(eps2).LT.abs(eps1))then
      dummy=eps2
      eps2=eps1
      eps1=dummy
      dummy=dr2
      dr2=dr1
      dr1=dummy
      end if
      write(6,*)'eps1,dr1=',eps1,dr1
cBH      if(abs(eps1).LE.100*1.0e-16*eps)then
cBH      if(abs(eps1).LE.100*1.0e-8*eps)then
      if(abs(eps1).LE.100*1.0d-16*eps)then
      dbmin=magfield(dr1,z0,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      drb=-newt(pi,psipol0,-drb,psivar,kappa,ipola,r0,0.0d0,0.0d0)
      dbmax=magfield(drb,z0,psivar,kappa,ipola,r0,bz,br,dbzeta,b0,psipol
     &0)
      dra=dr1
      za=z0
      drc=(dra+drb)/2
      zc=z0
      epsg=(dra-drb)/(2*r0+dra+drb)
      return
      end if
      dr2=dr1-eps1*(dr2-dr1)/(eps2-eps1)
      end do
      write(6,*)'psiinit: no convergence'
      return
      end
C* :69 * 
C* 70: * 
*line 1779 "adj.web"
      subroutine magfldar(kmax,deltab,ds,len,qs,psivar,kappa,ipola,r0,ps
     &ipol0,drc,zc,b0,dbmin,dbmax,dra,za)
      implicit none
      real*8 deltab,ds,len,qs
      integer kmax,k
      dimension deltab(0:kmax-1),ds(0:kmax-1)
      external newt,magfield
cSm020909
c      real*8 newt,magfield,psipol0,pi,rpol,rpol1,thpol0,thpol,dthpol,b0,db
c     &min,dbmax,dbtot,bz,br,dbzeta,bthpol,psivar,kappa,ipola,r0,drc,zc,d
c     &ra,za,dummy,dummya
	real*8 newt,magfield,psipol0,pi,rpol,rpol1,thpol0,thpol,dthpol,b0,
     &dbmin,dbmax,dbtot,bz,br,dbzeta,bthpol,psivar,kappa,ipola,r0,drc,zc
     &,dra,za,dummy,dummya
      pi=4*atan(1.0d0)
      dthpol=2*pi/kmax
      rpol1=sqrt((dra-drc)**2+(za-zc)**2)
      if(rpol1.GT.0)then
      thpol0=atan2(za-zc,dra-drc)
      else
      thpol0=0
      end if
      len=0
      qs=0
      do k=0,kmax-1
      thpol=thpol0+(k+0.25d0)*dthpol
      rpol=newt(thpol,psipol0,rpol1,psivar,kappa,ipola,r0,drc,zc)
      dbtot=magfield(drc+rpol*cos(thpol),zc+rpol*sin(thpol),psivar,kappa
     &,ipola,r0,bz,br,dbzeta,dummy,dummya)
      if(dbtot.LT.dbmin.OR.dbtot.GT.dbmax)then
      write(6,*)'Sanity check on magnetic field failed:',dbtot,dbmin,dbm
     &ax
      end if
      deltab(k)=(dbtot-dbmin)/(b0+dbmin)
      if(rpol.EQ.0)then
cSm020909
c      ds(k)=1/real*8(kmax)
	 ds(k)=1/dfloat(kmax)
      else
      bthpol=-br*sin(thpol)+bz*cos(thpol)
      ds(k)=rpol*(b0+dbtot)/bthpol*dthpol
      end if
      len=len+ds(k)
      qs=qs+((b0+dbzeta)/(b0+dbtot))*ds(k)/(r0+drc+rpol*cos(thpol))
      rpol1=rpol
      end do
      qs=qs/(2*pi)
      do k=0,kmax-1
      ds(k)=ds(k)/len
      end do
      return
      end
C* :70 * 
C* 75: * 
*line 1861 "adj.web"
      subroutine maginit(eps,psivar,kappa,ipola,th0a,th0g,lam1a,lam2g,ha
     &rr,imax,lmaxa1,kmax,deltab,ds,garr)
      implicit none
      integer kmax,imax,lmaxa1,i
cSm020909
c      real*8 eps,psivar,kappa,ipola,r0,psipol0,dra,za,drc,zc,b0,dbmin,dbma
c     &x,epsg,len,qs,th0a,th0g,lam1a,lam2g,harr,deltab,ds,garr
      real*8 eps,psivar,kappa,ipola,r0,psipol0,dra,za,drc,zc,b0,dbmin,
     &dbmax,epsg,len,qs,th0a,th0g,lam1a,lam2g,harr,deltab,ds,garr

      dimension th0a(0:imax-1),th0g(0:imax),lam1a(0:imax-1),lam2g(0:imax
     &),harr((((lmaxa1-1)*lmaxa1*(2*lmaxa1-1))/6+(lmaxa1-1)*lmaxa1+lmaxa
     &1))
      dimension deltab(0:kmax-1),ds(0:kmax-1),garr(((((lmaxa1+3)-1)*(lma
     &xa1+3))/2+(lmaxa1+1)))
      real*8 lam1,lam2
      external lam1,lam2,hinit
      if(int(psivar+0.5d0).EQ.0)then
      r0=1
      else if(int(psivar+0.5d0).EQ.1)then
      r0=1
      else
      r0=6.35
      end if
      call psiinit(eps,psivar,kappa,ipola,r0,psipol0,dra,za,drc,zc,b0,db
     &min,dbmax,epsg)
      call magfldar(kmax,deltab,ds,len,qs,psivar,kappa,ipola,r0,psipol0,
     &drc,zc,b0,dbmin,dbmax,dra,za)
      write(6,*)'eps,epserr,epsg,len,qs,psipol',eps,(dbmax-dbmin)/(2*b0+
     &dbmax+dbmin)-eps,epsg,len,qs,psipol0
      do i=0,imax-1
      lam1a(i)=lam1(deltab,ds,th0a(i),kmax)
      end do
      do i=0,imax
      lam2g(i)=lam2(deltab,ds,th0g(i),kmax)
      end do
      call hinit(deltab,ds,kmax,lmaxa1,harr,garr)
      return
      end
C* :75 * 
C* 78: * 
*line 1923 "adj.web"
      function lam1(deltab,ds,th0,kmax)
      implicit none
      integer kmax,k
      real*8 lam1,deltab,ds,th0,sum,costh
      dimension deltab(0:kmax-1),ds(0:kmax-1)
      sum=0
      do k=0,kmax-1
      costh=sqrt(1-min(1.0d0,deltab(k)*tan(th0)**2))
      sum=sum+ds(k)/costh
      end do
      lam1=sum
      return
      end
C* :78 * 
C* 79: * 
*line 1943 "adj.web"
      function lam2(deltab,ds,th0,kmax)
      implicit none
      integer kmax,k
      real*8 lam2,deltab,ds,th0,sum,costh
      dimension deltab(0:kmax-1),ds(0:kmax-1)
      sum=0
      do k=0,kmax-1
      costh=sqrt(1-min(1.0d0,deltab(k)*tan(th0)**2))
      sum=sum+ds(k)*costh/(1+deltab(k))
      end do
      lam2=sum
      return
      end
C* :79 * 
C* 80: * 
*line 1974 "adj.web"
      subroutine hinit(deltab,ds,kmax,lmaxa1,harr,garr)
      implicit none
      integer kmax,k,lmaxa1,la,la1,la2,l,l1
      real*8 harr,deltab,ds,b,garr
      dimension harr((((lmaxa1-1)*lmaxa1*(2*lmaxa1-1))/6+(lmaxa1-1)*lmax
     &a1+lmaxa1))
      dimension garr(((((lmaxa1+3)-1)*(lmaxa1+3))/2+(lmaxa1+1)))
      dimension deltab(0:kmax-1),ds(0:kmax-1)
C* 81: * 
*line 2012 "adj.web"
      do la=1,lmaxa1
      do la1=1,la
      do la2=1,la
      harr((((la-1)*la*(2*la-1))/6+(la1-1)*la+la2))=0
      end do
      end do
      end do
      do la=-2,lmaxa1-1
      garr(((((la+3)-1)*(la+3))/2+(0+1)))=0
      end do
      do la=-1,lmaxa1-1
      garr(((((la+3)-1)*(la+3))/2+(la+1+1)))=0
      garr(((((la+3)-1)*(la+3))/2+(la+2+1)))=0
      end do
      garr(((((1+3)-1)*(1+3))/2+(1+1)))=1
C* :81 * 
*line 1983 "adj.web"
      
      do k=0,kmax-1
      b=1+deltab(k)
C* 82: * 
*line 2032 "adj.web"
      do la=2,lmaxa1
      l=max(2*la-1,0)
      do la1=1,la
      l1=max(2*la1-1,0)
      garr(((((la+3)-1)*(la+3))/2+(la1+1)))=(2*(2*l-3)*(l**2-3*l+1)*garr
     &(((((la-1+3)-1)*(la-1+3))/2+(la1+1)))-(l-3)*(l-2)*(2*l-1)*garr((((
     &(la-2+3)-1)*(la-2+3))/2+(la1+1))))/(l*(l-1)*(2*l-5))+(2*l-3)*(2*l-
     &1)*b*((l1+1)*(l1+2)*garr(((((la-1+3)-1)*(la-1+3))/2+(la1+1+1)))/((
     &2*l1+3)*(2*l1+5))-2*(l1**2+l1-1)*garr(((((la-1+3)-1)*(la-1+3))/2+(
     &la1+1)))/((2*l1-1)*(2*l1+3))+l1*(l1-1)*garr(((((la-1+3)-1)*(la-1+3
     &))/2+(la1-1+1)))/((2*l1-3)*(2*l1-1)))/(l*(l-1))
      end do
      end do
C* :82 * 
*line 1987 "adj.web"
      
      do la=1,lmaxa1
      do la1=1,la
      do la2=1,la
      harr((((la-1)*la*(2*la-1))/6+(la1-1)*la+la2))=harr((((la-1)*la*(2*
     &la-1))/6+(la1-1)*la+la2))+ds(k)*b*garr(((((la+3)-1)*(la+3))/2+(la1
     &+1)))*garr(((((la+3)-1)*(la+3))/2+(la2+1)))
      end do
      end do
      end do
      end do
      do la=1,lmaxa1
      do la1=1,la
      do la2=1,la
cSm020909
c      harr((((la-1)*la*(2*la-1))/6+(la1-1)*la+la2))=(2*max(2*la-1,0)+1)/
c     &real*8(2*max(2*la2-1,0)+1)*harr((((la-1)*la*(2*la-1))/6+(la1-1)*la+l
c     &a2))
      harr((((la-1)*la*(2*la-1))/6+(la1-1)*la+la2))=(2*max(2*la-1,0)+1)/
     &dfloat(2*max(2*la2-1,0)+1)*harr((((la-1)*la*(2*la-1))/6+(la1-1)*la
     &+la2))
      end do
      end do
      end do
      return
      end
C* :80 * 
C* 83: * 
*line 2054 "adj.web"
C* :83 * 
C* 83: * 
*line 2087 "adj.web"
C* :83 * 
C* 84: * 
*line 2090 "adj.web"
C* :84 * 
C* 84: * 
*line 2092 "adj.web"
C* :84 * 
C* 84: * 
*line 2095 "adj.web"
      subroutine residue(ua,nmax,nmax1,csth0a,imax,imax1,chi,dtt,c,lam1a
     &,uinv,th0inv,resid)
      implicit none
      integer nmax,nmax1,imax,imax1,n,i
      real*8 ua,csth0a,chi,dtt,c,lam1a,resid,uinv,th0inv
      dimension ua(0:nmax-1),chi(0:imax1-1,0:nmax-1),csth0a(0:imax-1),dt
     &t(0:nmax-1),lam1a(0:imax-1),resid(0:imax1-1,0:nmax-1),uinv(0:nmax1
     &-1,-1:1),th0inv(0:imax1-1,-1:1)
C* 88: * 
*line 2196 "adj.web"
      do n=0,nmax-1
      do i=0,imax-1
      resid(i,n)=resid(i,n)+csth0a(i)/lam1a(i)*ua(n)/sqrt(1+(ua(n)/c)**2
     &)
      end do
      end do
C* :88 * 
*line 2105 "adj.web"
      
      
      do n=0,nmax-1
      do i=0,imax-1
      resid(i,n)=resid(i,n)*ua(n)**2/dtt(n)+chi(i,n-1)*uinv(n,-1)+chi(i,
     &n)*uinv(n,0)+chi(i,n+1)*uinv(n,1)+chi(i-1,n)*th0inv(i,-1)+chi(i,n)
     &*th0inv(i,0)+chi(i+1,n)*th0inv(i,1)
      end do
      end do
      return
      end
C* :84 * 
C* 84: * 
*line 2124 "adj.web"
C* :84 * 
C* 89: * 
*line 2206 "adj.web"
      subroutine error(du,ua,umax,nmax,dth0,th0a,th0max,imax,imax1,resid
     &,chi,abserr,relerr,maxerr,aerr,rerr,merr)
      implicit none
      integer nmax,imax,imax1,i,n
      real*8 du,ua,umax,dth0,th0a,th0max,resid,chi,abserr,relerr,maxerr
      dimension ua(0:nmax-1),th0a(0:imax-1),resid(0:imax1-1,0:nmax-1),ch
     &i(0:imax1-1,0:nmax-1)
      real*8 aerr,rerr,merr
      dimension aerr(0:imax-1),rerr(0:imax-1),merr(0:imax-1)
      do i=0,imax-1
      aerr(i)=0
      rerr(i)=0
      merr(i)=0
      end do
      do n=0,nmax-1
      do i=0,imax-1
      aerr(i)=aerr(i)+resid(i,n)**2*ua(n)**2
cBH      rerr(i)=rerr(i)+(resid(i,n)/max(1.0e-16,abs(chi(i,n))))**2*ua(n)**
cBH     &2
cBH      rerr(i)=rerr(i)+(resid(i,n)/max(1.0e-8,abs(chi(i,n))))**2*ua(n)**2
      rerr(i)=rerr(i)+(resid(i,n)/max(1.0e-16,abs(chi(i,n))))**2*ua(n)**
     &2

      merr(i)=max(merr(i),abs(resid(i,n)))
      end do
      end do
      abserr=0
      relerr=0
      maxerr=0
      do i=0,imax-1
      abserr=abserr+sin(th0a(i))*aerr(i)
      relerr=relerr+sin(th0a(i))*rerr(i)
      maxerr=max(maxerr,merr(i))
      end do
      
      abserr=sqrt(abserr*du*dth0/(umax**3/3*(1-cos(th0max))))
      relerr=sqrt(relerr*du*dth0/(umax**3/3*(1-cos(th0max))))
      return
      end
C* :89 * 
C* 95: * 
*line 2278 "adj.web"
      subroutine matinit(du,ua,ug,nmax,nmax1,dth0,th0a,th0g,imax,imax1,d
     &uug,fu,dtt,lam1a,lam2g,uinv,th0inv,rho)
      implicit none
      integer nmax,nmax1,imax,imax1,n,i
cSm020909
c      real*8 du,ua,ug,dth0,th0a,th0g,duug,fu,dtt,lam1a,lam2g,uinv,th0inv,r
c     &ho
	real*8 du,ua,ug,dth0,th0a,th0g,duug,fu,dtt,lam1a,lam2g,uinv,th0inv
     &,rho
      dimension ua(0:nmax-1),ug(0:nmax),th0a(0:imax-1),th0g(0:imax),duug
     &(0:nmax),fu(0:nmax-1),dtt(0:nmax-1),lam1a(0:imax-1),lam2g(0:imax),
     &uinv(0:nmax1-1,-1:1),th0inv(0:imax1-1,-1:1)
C* 96: * 
*line 2301 "adj.web"
      do n=0,nmax-1
      uinv(n,-1)=ug(n)**2*duug(n)/(ua(n)**2*du**2)
      uinv(n,1)=ug(n+1)**2*duug(n+1)/(ua(n)**2*du**2)
      uinv(n,0)=-uinv(n,-1)-uinv(n,1)
      uinv(n,-1)=uinv(n,-1)-fu(n)/(2*du)
      uinv(n,1)=uinv(n,1)+fu(n)/(2*du)
      end do
      n=0
      uinv(n,0)=uinv(n,0)-uinv(n,-1)
      uinv(n,-1)=0
      n=nmax-1
      uinv(n,-1)=uinv(n,-1)-rho*uinv(n,1)
      uinv(n,0)=uinv(n,0)+(1+rho)*uinv(n,1)
      uinv(n,1)=0
      do n=0,nmax-1
      uinv(n,-1)=uinv(n,-1)*ua(n)**2/dtt(n)
      uinv(n,0)=uinv(n,0)*ua(n)**2/dtt(n)
      uinv(n,1)=uinv(n,1)*ua(n)**2/dtt(n)
      end do
C* :96 * 
*line 2290 "adj.web"
      
C* 97: * 
*line 2327 "adj.web"
      do i=0,imax-1
      th0inv(i,-1)=sin(th0g(i))*lam2g(i)/(lam1a(i)*sin(th0a(i))*dth0**2)
      th0inv(i,1)=sin(th0g(i+1))*lam2g(i+1)/(lam1a(i)*sin(th0a(i))*dth0*
     &*2)
      th0inv(i,0)=-th0inv(i,-1)-th0inv(i,1)
      end do
      i=0
      i=imax-1
      th0inv(i,0)=th0inv(i,0)-th0inv(i,1)
      th0inv(i,1)=0
C* :97 * 
*line 2292 "adj.web"
      
      return
      end
C* :95 * 
C* 98: * 
*line 2355 "adj.web"
C* :98 * 
C* 98: * 
*line 2405 "adj.web"
C* :98 * 
C* 105: * 
*line 2490 "adj.web"
C* :105 * 
C* 105: * 
*line 2518 "adj.web"
C* :105 * 
C* 108: * 
*line 2541 "adj.web"
      function current(du,ua,c,nmax,dth0,th0a,imax,imax1,fm,chi,cur)
      implicit none
      integer nmax,imax,imax1,n,i
      real*8 current,du,ua,c,dth0,th0a,fm,chi
      dimension ua(0:nmax-1),th0a(0:imax-1),fm(0:nmax-1),chi(0:imax1-1,0
     &:nmax-1)
      real*8 cur,pi
      dimension cur(0:imax-1)
      pi=4*atan(1.0d0)
      do i=0,imax-1
      cur(i)=0
      end do
      do n=0,nmax-1
      do i=0,imax-1
      cur(i)=cur(i)+chi(i,n)*fm(n)*ua(n)**3/sqrt(1+(ua(n)/c)**2)
      end do
      end do
      current=0
      do i=0,imax-1
      current=current+sin(th0a(i))*cos(th0a(i))*cur(i)
      end do
      current=4*pi*du*dth0*current
      return
      end
C* :108 * 
C* 110: * 
*line 2580 "adj.web"
C* :110 * 
C* 110: * 
*line 2582 "adj.web"
c$$$      subroutine hdfout(resid,chi,du,umax,nmax,dth0,imax,imax1,time)
c$$$      implicit none
c$$$      integer nmax,imax,imax1,time,minchar,maxchar,intervals,pixsize,n,i
c$$$     &,k,d8pimg,d8aimg,d8spal,dpgpal,ret
c$$$      external d8pimg,d8aimg,d8spal,dpgpal
c$$$      real*8 resid,du,umax,dth0,chi
c$$$      dimension resid(0:imax1-1,0:nmax-1),chi(0:imax1-1,0:nmax-1)
c$$$      parameter(minchar=2,maxchar=255,intervals=20,pixsize=40000)
c$$$      integer DFTAG_RLE,DFTAG_IMC
c$$$      parameter(DFTAG_RLE=11,DFTAG_IMC=12)
c$$$      character*1 ch(0:pixsize-1)
c$$$      character*1 palette(3,256)
c$$$      k=0
c$$$      do n=nmax-1,0,-1
c$$$      do i=imax-1,0,-1
c$$$      ch(k)=char(max(minchar,min(maxchar,int(log10(max(1.0e-16,abs(resid
c$$$     &(i,n))))*intervals+maxchar))))
c$$$      k=k+1
c$$$      end do
c$$$      end do
c$$$      if(time.EQ.0)then
c$$$      ret=dpgpal('rainbow.pal',palette)
c$$$      ret=d8spal(palette)
c$$$      ret=d8pimg('adj.hdf',ch,imax,nmax,DFTAG_RLE)
c$$$      else
c$$$      ret=d8aimg('adj.hdf',ch,imax,nmax,DFTAG_RLE)
c$$$      end if
c$$$      return
c$$$      end
C* :110 * 
C* 110: * 
*line 2615 "adj.web"
C* :110 * 
C* 113: * 
*line 2668 "adj.web"
      subroutine legjs(z,lmax,amax,jarray,lmax1)
      implicit none
      integer l,a,lmax,amax,lmax1,lstart
      real*8 z,jarray,z2,g,ja0,ja1,ja2,ratio,sigma
      dimension jarray(-lmax1-1:lmax1,0:amax)
      if(amax.LT.0.OR.lmax.LT.0)return
      if(lmax.GT.lmax1)return
      z2=z**2
      g=sqrt(1+z2)
C* 114: * 
*line 2696 "adj.web"
      
      do a=1,amax
      ja2=1
      ja1=g
      do l=a-3,lmax-1,-1
      ja0=g*ja1-z2*(l-a+2)*(l+a+2)*ja2/((2*l+3)*(2*l+5))
      ja2=ja1
      ja1=ja0
      end do
      l=min(lmax,a-1)
      jarray(l,a)=ja2
      l=l-1
      jarray(l,a)=ja1
      do l=min(lmax,a-1)-2,max(-a,-lmax-1),-1
      jarray(l,a)=g*jarray(l+1,a)-z2*(l-a+2)*(l+a+2)*jarray(l+2,a)/((2*l
     &+3)*(2*l+5))
      end do
      end do
C* :114 * 
*line 2681 "adj.web"
      
C* 115: * 
*line 2722 "adj.web"
      
      do a=0,min(amax,lmax)
      l=-a-1
      jarray(l,a)=1
      l=l-1
      if(l.GE.-lmax-1)then
      jarray(l,a)=g
      end if
      do l=-a-3,-lmax-1,-1
      jarray(l,a)=g*jarray(l+1,a)-z2*(l-a+2)*(l+a+2)*jarray(l+2,a)/((2*l
     &+3)*(2*l+5))
      end do
      end do
C* :115 * 
*line 2683 "adj.web"
      
C* 116: * 
*line 2742 "adj.web"
      a=0
      if(abs(z).LE.lmax+1)then
cBH      lstart=lmax+log(1.0e-16)/(2*log((1.0e-16+abs(z))/(g+1)))+10
cBH      lstart=lmax+log(1.0e-8)/(2*log((1.0e-8+abs(z))/(g+1)))+10
      lstart=lmax+log(1.0d-16)/(2*log((1.0d-16+abs(z))/(g+1)))+10
      ratio=2/(g+1)
      do l=lstart,lmax-1,-1
      ratio=1/(g-z2*(l-a+2)*(l+a+2)*ratio/((2*l+3)*(2*l+5)))
      end do
      jarray(lmax,a)=ratio
      ja0=1
      do l=lmax-2,-1,-1
      jarray(l+1,a)=ja0
      ja0=g*jarray(l+1,a)-z2*(l-a+2)*(l+a+2)*jarray(l+2,a)/((2*l+3)*(2*l
     &+5))
      end do
      ja0=1/ja0
      do l=0,lmax
      jarray(l,a)=ja0*jarray(l,a)
      end do
      else
      sigma=log(g+abs(z))
      l=0
      jarray(l,a)=sigma/abs(z)
      do l=1,lmax
      jarray(l,a)=(g*jarray(l-1,a)-jarray(l-2,a))*(2*l-1)*(2*l+1)/(z2*(l
     &-a)*(l+a))
      end do
      end if
C* :116 * 
*line 2685 "adj.web"
      
C* 117: * 
*line 2779 "adj.web"
      do a=1,min(lmax,amax)
      l=lmax
      jarray(l,a)=((2*l+1)*jarray(l-1,a-1)-(l+1-a)*g*jarray(l,a-1))/(l+a
     &)
      l=lmax-1
      if(l.EQ.0)then
      ja0=((2*l+1)-(l+1-a)*g*jarray(l,a-1))/(l+a)
      else
      ja0=((2*l+1)*jarray(l-1,a-1)-(l+1-a)*g*jarray(l,a-1))/(l+a)
      end if
      do l=lmax-2,a-1,-1
      jarray(l+1,a)=ja0
      ja0=g*jarray(l+1,a)-z2*(l-a+2)*(l+a+2)*jarray(l+2,a)/((2*l+3)*(2*l
     &+5))
      end do
      ja0=1/ja0
      do l=a,lmax
      jarray(l,a)=ja0*jarray(l,a)
      end do
      end do
C* :117 * 
*line 2687 "adj.web"
      
      return
      end
C* :113 * 
C* 118: * 
*line 2807 "adj.web"
C* :118 * 
C* 118: * 
*line 2838 "adj.web"
C* :118 * 
C* 120: * 
*line 2872 "adj.web"
      subroutine legjn(u,c,lmax,jarray,djarray,lmax1)
      implicit none
      integer j0,j1,j2,j02,j11,j22,j022
      parameter(j0=0,j1=1,j2=2,j02=3,j11=4,j22=5,j022=6)
      integer l,a,lmax,amax,lmax1
      real*8 u,c,jarray,djarray,scale,g
      dimension jarray(-lmax1-1:lmax1,0:6),djarray(-lmax1-1:lmax1,0:6)
      external legjs
      amax=2
      call legjs(u/c,lmax+2,amax,jarray,lmax1)
      scale=1
      do l=1,lmax+2
      scale=scale*u/(2*l+1)
      do a=0,amax
      jarray(l,a)=scale*jarray(l,a)
      end do
      end do
      scale=1
      do l=-1,-lmax-2,-1
      scale=scale/u*(2*l+3)
      do a=0,amax
      jarray(l,a)=scale*jarray(l,a)
      end do
      end do
      do l=-lmax-2,lmax
      jarray(l,j02)=u/2*jarray(l+1,1)
      jarray(l,j11)=u/2*jarray(l+1,0)
      jarray(l,j22)=u/2*jarray(l+1,1)+(u/c)**2/2*jarray(l+2,0)
      jarray(l,j022)=u**2/8*jarray(l+2,0)
      end do
      g=sqrt(1+(u/c)**2)
      do a=0,6
      do l=-lmax-1,lmax
      djarray(l,a)=jarray(l-1,a)/g-(l+1)/u*jarray(l,a)
      end do
      end do
      return
      end
C* :120 * 
*line 540 "adj.web"
      
C* :14 * 
C* 14: * 
*line 547 "adj.web"
      
C* :14 * 
      
