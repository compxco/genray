c*****************************************************************
c*                   Spline programs                             *
c*****************************************************************
! YuP[2020-08-20] Renamed b into bsp, to avoid conflict with function b()
      subroutine iac1r(x,nx,nx4,fx,lx,mxa,fxa,mxb,fxb,tx,cx)
cSm030221
         implicit none !double precision(a-h,o-z)

         integer lx,mxa,mxb,nx,nx4
         double precision cx,fx,fxa,fxb,tx,x
cSm011006
c         dimension cx(nx4),fx(nx),tx(nx4),x(nx)
         dimension cx(*),fx(*),tx(*),x(*)

         external iac1r1,iac1r2,iac1r3,iac1r4
         integer i,j,k,k1,k2,k3
         double precision a,bsp,iac1r1
         dimension a(3),bsp(3)
         data k1/1/,k2/2/,k3/3/
       call iac1r2(x,nx,nx4,mxa,tx)
       if (mxa .lt. k3) go to 10
          i = nx - k2
          cx(k1) = fx(i)
          i = i + k1
          cx(k2) = fx(i)
          i = nx + k2 + k1
          cx(i) = fx(k2)
          i = i + k1
          cx(i) = fx(k3)
          go to 20
   10  continue
       call iac1r3(tx,nx,nx4,cx,lx,mxa,mxb,a,bsp)
   20  continue
       do 30 i = k1, nx
          j = i + k2
          cx(j) = fx(i)
   30  continue
       if (mxa .gt. k2) go to 40
          call iac1r4(tx,nx,nx4,cx,lx,mxa,fxa,mxb,fxb,a,bsp)
   40  continue
       i = nx + k2
       do 50 j = k1, i
          k = j + k1
          cx(j) = iac1r1(tx,cx,nx4,lx,k)
   50  continue
       return
      end subroutine iac1r
      
      
      
      double precision function iac1r1(tx,cx,nx4,lx,jx)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer jx,lx,nx4
         double precision cx,tx
cSm011006
c         dimension cx(*),tx(*)
         dimension cx(nx4),tx(nx4)

         integer i,j,k1
         double precision a,bsp,c,c1,c112,c13,c240,c3,c438,d
         data c1/1.0d0/,c112/112.0d0/,c13/13.0d0/,c240/240.0d0/,
     1        c3/3.0d0/,
     1        c438/438.0d0/,k1/1/
       a = cx(jx)
       if (lx .lt. k1) go to 20
          i = jx + k1
          j = jx - k1
          if (lx .gt. k1) go to 10
             bsp = tx(i) - tx(jx)
             c = tx(jx) - tx(j)
             d = c*c/bsp/(bsp + c)
             c = bsp*bsp/c/(bsp + c)
             bsp = c1/c/d
             a = (bsp*a - c*cx(j) - d*cx(i))/c3
             go to 20
   10  continue
       bsp = cx(j) + cx(i)
       i = i + k1
       j = j - k1
       c = cx(j) + cx(i)
       a = (c438*a - c112*bsp + c13*c)/c240
   20  continue
       iac1r1 = a
       return
      end function iac1r1
      
      
      subroutine iac1r2(x,nx,nx4,mx,tx)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer mx,nx,nx4
         double precision tx,x
cSm011006
c         dimension tx(nx4),x(nx4)
         dimension tx(*),x(*)
         integer i,j,k1,k2,k3
         double precision a
         data k1/1/,k2/2/,k3/3/
       do 10 i = k1, nx
          j = i + k2
          tx(j) = x(i)
   10  continue
       if (mx .lt. k3) go to 20
          a = x(nx) - x(k1)
          j = nx - k2
          tx(k1) = x(j) - a
          j = j + k1
          tx(k2) = x(j) - a
          j = nx + k3
          tx(j) = x(k2) + a
          j = j + k1
          tx(j) = x(k3) + a
          go to 30
   20  continue
       a = x(k2) - x(k1)
       tx(k1) = x(k1) - k2*a
       tx(k2) = x(k1) - a
       i = nx - k1
       a = x(nx) - x(i)
       i = nx + k3
       tx(i) = x(nx) + a
       i = i + k1
       tx(i) = x(nx) + k2*a
   30  continue
       return
      end subroutine iac1r2
      
      
      subroutine iac1r3(tx,nx,nx4,cx,lx,mxa,mxb,ac,bc)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer lx,mxa,mxb,nx,nx4
         double precision ac,bc,cx,tx
cSm011006
c        dimension ac(3),bc(3),cx(nx4),tx(nx4)
         dimension ac(3),bc(3),cx(*),tx(*)
         external iac1r1,ias1r2
         integer i,j,k,k0,k1,k2,k3,k4,l
         double precision a,bsp,c,c0,c1,d,e,f,g,iac1r1,ias1r2
         data c0/0.0d0/,c1/1.0d0/,k0/0/,k1/1/,k2/2/,k3/3/,k4/4/
       i = k4 + k2
       do 10 j = k1, i
          cx(j) = c0
          k = nx + i - k1 - j
          cx(k) = c0
   10  continue
       i = nx + k1
       j = i + k1
       k = j + k1
       l = k + k1
       bsp = c1
       cx(k1) = c1
       a = ias1r2(tx,cx,nx4,k4,mxa,tx(k3),k3,k1)
       if (lx .eq. k0) go to 30
          bsp = iac1r1(tx,cx,nx4,lx,k2)
          cx(k1) = c0
          cx(k2) = c1
          c = iac1r1(tx,cx,nx4,lx,k3)
          if (mxa .gt. k0) go to 20
             d = ias1r2(tx,cx,nx4,k4,k0,tx(k4),k4,k1)
             cx(k2) = 0
             ac(k2) = c1/c/d
             go to 30
   20  continue
       e = ias1r2(tx,cx,nx4,k4,k0,tx(k3),k3,k1)
       f = ias1r2(tx,cx,nx4,k4,mxa,tx(k3),k3,k1)
       cx(k2) = c0
       cx(k1) = c1
       g = ias1r2(tx,cx,nx4,k4,k0,tx(k3),k3,k1)
       d = g*f - e*a
       ac(k3) = g/d/c
       ac(k2) = a/d/c
       a = g
   30  continue
       cx(k1) = c0
       ac(k1) = c1/a/bsp
       cx(j) = c1
       a = ias1r2(tx,cx,nx4,k4,mxb,tx(j),i,k1)
       if (lx .eq. k0) go to 50
          cx(j) = c0
          cx(l) = c1
          bsp = iac1r1(tx,cx,nx4,lx,k)
          cx(l) = c0
          cx(k) = c1
          c = iac1r1(tx,cx,nx4,lx,j)
          cx(k) = c0
          cx(i) = c1
          if (mxb .gt. k0) go to 40
             d = ias1r2(tx,cx,nx4,k4,mxb,tx(i),i,k1)
             cx(i) = c0
             bc(k2) = c1/c/d
             go to 50
   40  continue
       e = ias1r2(tx,cx,nx4,k4,k0,tx(j),i,k1)
       f = ias1r2(tx,cx,nx4,k4,mxb,tx(j),i,k1)
       cx(i) = c0
       cx(j) = c1
       g = ias1r2(tx,cx,nx4,k4,k0,tx(j),i,k1)
       d = g*f - e*a
       bc(k3) = g/d/c
       bc(k2) = a/d/c
       a = g
   50  continue
       cx(j) = c0
       bc(k1) = c1/a/bsp
       return
      end subroutine iac1r3
      
      
      subroutine iac1r4(tx,nx,nx4,cx,lx,mxa,fxa,mxb,fxb,ac,bc)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer lx,mxa,mxb,nx,nx4
         double precision ac,bc,cx,fxa,fxb,tx
cSm011006
c         dimension ac(3),bc(3),cx(nx4),tx(nx4)
         dimension ac(3),bc(3),cx(*),tx(*)

         external iac1r1,ias1r2
         integer i,j,k,k0,k1,k2,k3,k4,l,m
         double precision a,bsp,c,c0,iac1r1,ias1r2
         dimension a(4)
         data c0/0.0d0/,k0/0/,k1/1/,k2/2/,k3/3/,k4/4/
       cx(k1) = c0
       cx(k2) = c0
       k = nx + k3
       l = k + k1
       cx(k) = c0
       cx(l) = c0
       m = mxa
       do 10 i = k1, k4
          j = i + k1
          a(i) = iac1r1(tx,cx,nx4,lx,j)
   10  continue
       if (m .gt. k0 .or. lx .eq. k0) go to 30
          c = ias1r2(tx,a,nx4,k4,k0,tx(k4),k3,k0)
          cx(k2) = (cx(k4) - c)*ac(k2)
   20  continue
       a(k1) = iac1r1(tx,cx,nx4,lx,k2)
       a(k2) = iac1r1(tx,cx,nx4,lx,k3)
   30  continue
       c = ias1r2(tx,a,nx4,k4,m,tx(k3),k3,k0)
       if (m .gt. k0) go to 40
          j = k2 - lx
          cx(j) = (cx(k3) - c)*ac(k1)
          go to 60
   40  continue
       if (lx .eq. k1) go to 50
          cx(k2) = (fxa - c)*ac(k1)
          go to 60
   50  continue
       bsp = ias1r2(tx,a,nx4,k4,k0,tx(k3),k3,k0)
       cx(k2) = ac(k3)*(fxa - c) - ac(k2)*(cx(k3) - bsp)
       m = k0
       go to 20
   60  continue
       m = mxb
       k = k - k1
       l = k - k1
       do 70 i = k1, k4
          j = k - k3 + i
          a(i) = iac1r1(tx,cx,nx4,lx,j)
   70  continue
       if (m .gt. k0 .or. lx .eq. k0) go to 90
          c = ias1r2(tx,a,nx4,k4,k0,tx(l),l,k0)
          cx(j) = (cx(l) - c)*bc(k2)
   80  continue
       a(k3) = iac1r1(tx,cx,nx4,lx,k)
       a(k4) = iac1r1(tx,cx,nx4,lx,j)
   90  continue
       c = ias1r2(tx,a,nx4,k4,m,tx(k),l,k0)
       if (m .gt. k0) go to 100
          j = j + lx
          cx(j) = (cx(k) - c)*bc(k1)
          go to 120
  100  continue
       if (lx .eq. k1) go to 110
          cx(j) = (fxb - c)*bc(k1)
          go to 120
  110  continue
       bsp = ias1r2(tx,a,nx4,k4,k0,tx(k),l,k0)
       cx(j) = bc(k3)*(fxb - c) - bc(k2)*(cx(k) - bsp)
       m = k0
       go to 80
  120  continue
       return
      end subroutine iac1r4

      double precision function ias1r(tx,nx,nx4,cx,idx,xx)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer idx,nx,nx4
         double precision cx,tx,xx
cSm011006
c         dimension cx(nx4),tx(nx4)
         dimension cx(*),tx(*)
         external ias1r1,ias1r2
         integer i,ias1r1,k1,k4
         double precision ias1r2
         data k1/1/,k4/4/
cc	 write(*,*)'in ias1r nx=',nx,'nx4=',nx4,'idx=',idx,'xx=',xx
c	 write(*,*)(tx(i),i=1,nx4)
c	 write(*,*)(cx(i),i=1,nx)
cc	 write(*,*)'tx(i)'
cc	 write(*,*)tx
cc	 write(*,*)'cx(i)'
cc	 write(*,*)cx
       i = ias1r1(tx,nx,nx4,xx)
       ias1r = ias1r2(tx,cx,nx4,k4,idx,xx,i,k1)
       return
      end function ias1r
      
      
      integer function ias1r1(tx,nx,nx4,xx)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer nx,nx4
         double precision tx,xx
cSm011006
c         dimension tx(nx4)
         dimension tx(*)

         integer i,j,k,k1,k2
         data k1/1/,k2/2/
       i = k2*k2
       j = nx + k1
   10  continue
       k = (i + j)/k2
       if (xx .ge. tx(k)) go to 20
          j = k
          if (j .gt. i) go to 10
             k = k - k1
             go to 30
   20  continue
       i = k + k1
       if (j .ge. i) go to 10
   30  continue
       ias1r1 = k
       return
      end function ias1r1
      
      
      double precision function ias1r2(tx,cx,nx4,kx,idx,xx,ix,jx)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer idx,ix,jx,kx,nx4
         double precision cx,tx,xx
cSm011006
c         dimension cx(nx4),tx(nx4)
         dimension cx(*),tx(*)

         integer i,j,k,k1,l,m,n,n1
         double precision a,bsp,c
         dimension a(8),bsp(7),c(7)
         data k1/1/
       i = jx*(ix - kx + k1)
       do 10 j = k1, kx
          k = i + j
          a(j) = cx(k)
   10  continue
       if (idx .lt. k1) go to 40
          do 30 i = k1, idx
             j = kx - i
             do 20 k = k1, j
                l = ix + k
                m = l - j
                n = k + k1
                a(k) = (a(n) - a(k))/(tx(l) - tx(m))*j
   20        continue
   30     continue
   40  continue
       i = kx - idx - k1
       if (i .lt. k1) go to 80
          do 50 j = k1, i
             k = ix + j
             bsp(j) = tx(k) - xx
             k = ix + k1 - j
             c(j) = xx - tx(k)
   50     continue
          n1 = kx - k1
          i = idx + k1
          do 70 j = i, n1
             k = kx - j
             do 60 l = k1, k
                m = k + k1 - l
                n = l + k1
                a(l) = (a(n)*c(m) + a(l)*bsp(l))/(c(m) + bsp(l))
   60        continue
   70     continue
   80  continue
       ias1r2 = a(k1)
       return
      end function ias1r2
      
      
      subroutine iac2r(x,nx,y,ny,fxy,nfx,nfy,lx,ly,mxa,fxa,mxb,fxb,
     1                 mya,fya,myb,fyb,tx,ty,cxy,ncx,ncy,nry,cy)
	implicit none
         integer lx,ly,mxa,mxb,mya,myb,ncx,ncy,nfx,nfy,nx,ny,nry
         double precision cxy,cy,fxa,fxb,fxy,fya,fyb,tx,ty,x,y

cSm030221
c         dimension cxy(ncx,ncy),cy(nry),fxa(ny),fxb(ny),fxy(nfx,nfy),
c     1             fya(nx),fyb(nx),tx(ncx),ty(ncy),x(nx),y(ny)
         dimension cxy(ncx,*),cy(*),fxa(*),fxb(*),fxy(nfx,*),
     1             fya(*),fyb(*),tx(*),ty(*),x(*),y(*)
         external iac1r1,iac1r2,iac1r3,iac1r4
         integer i,j,k,k0,k1,k2,k3,l,m,n
         double precision a,bsp,fa,fb,iac1r1
         dimension a(3),bsp(3)
         data k0/0/,k1/1/,k2/2/,k3/3/
       call iac1r2(x,nx,ncx,mxa,tx)
       call iac1r2(y,ny,ncy,mya,ty)
       do 20 i = k1, nx
          j = k2 + i
          do 10 k = k1, ny
             l = k2 + k
             cxy(j,l) = fxy(i,k)
   10     continue
   20  continue
       m = nx + k2
       n = ny + k2
       if (mxa .gt. k2) go to 50
          call iac1r3(tx,nx,ncx,cy,lx,mxa,mxb,a,bsp)
          do 40 i = k3, n
             j = i - k2
             if (mxa .ge. k1) fa = fxa(j)
             if (mxb .ge. k1) fb = fxb(j)
             do 30 j = k3, m
                cy(j) = cxy(j,i)
   30        continue
             call iac1r4(tx,nx,ncx,cy,lx,mxa,fa,mxb,fb,a,bsp)
             cxy(k1,i) = cy(k1)
             cxy(k2,i) = cy(k2)
             k = m + k1
             cxy(k,i) = cy(k)
             k = k + k1
             cxy(k,i) = cy(k)
   40     continue
   50  continue
       if (mya .gt. k2) go to 110
          call iac1r3(ty,ny,ncy,cy,ly,mya,myb,a,bsp)
          do 70 i = k3, m
             j = i - k2
             if (mya .ge. k1) fa = fya(j)
             if (myb .ge. k1) fb = fyb(j)
             do 60 j = k3, n
                cy(j) = cxy(i,j)
   60        continue
             call iac1r4(ty,ny,ncy,cy,ly,mya,fa,myb,fb,a,bsp)
             cxy(i,k1) = cy(k1)
             cxy(i,k2) = cy(k2)
             k = n + k1
             cxy(i,k) = cy(k)
             k = k + k1
             cxy(i,k) = cy(k)
   70     continue
          if (mxa .gt. k2) go to 130
             call iac1r3(ty,ny,ncy,cy,ly,k0,k0,a,bsp)
             do 100 i = k1, k2
                k = k3 - i
                do 80 j = k3, n
                   cy(j) = cxy(k,j)
   80           continue
                call iac1r4(ty,ny,ncy,cy,ly,k0,fa,k0,fb,a,bsp)
                cxy(k,k1) = cy(k1)
                cxy(k,k2) = cy(k2)
                l = n + k1
                cxy(k,l) = cy(l)
                l = l + k1
                cxy(k,l) = cy(l)
                k = m + i
                do 90 j = k3, n
                   cy(j) = cxy(k,j)
   90           continue
                call iac1r4(ty,ny,ncy,cy,ly,k0,fa,k0,fb,a,bsp)
                cxy(k,k1) = cy(k1)
                cxy(k,k2) = cy(k2)
                l = n + k1
                cxy(k,l) = cy(l)
                l = l + k1
                cxy(k,l) = cy(l)
  100        continue
             m = m + k2
             n = n + k2
             go to 150
  110  continue
       k = k3
       if (mxa .lt. k3) k = k - k2
       if (mxa .lt. k3) m = m + k2
       do 120 i = k, m
          j = n - k1
          cxy(i,k1) = cxy(i,j)
          j = j - k1
          cxy(i,k2) = cxy(i,j)
          j = k3 + k1
          l = n + k1
          cxy(i,l) = cxy(i,j)
          j = j + k1
          l = l + k1
          cxy(i,l) = cxy(i,j)
  120  continue
  130  continue
       n = n + k2
       if (mxa .gt. k2) go to 150
          do 140 i = k1, n
             j = m - k1
             cxy(k1,i) = cxy(j,i)
             j = j - k1
             cxy(k2,i) = cxy(j,i)
             l = k3 + k1
             j = m + k1
             cxy(j,i) = cxy(l,i)
             l = l + k1
             j = j + k1
             cxy(j,i) = cxy(l,i)
  140     continue
          m = m + k2
  150  continue
       l = n - k2
       do 180 i = k1, m
          do 160 j = k1, n
             cy(j) = cxy(i,j)
  160     continue
          do 170 j = k1, l
             k = k1 + j
             cxy(i,j) = iac1r1(ty,cy,ncy,ly,k)
  170     continue
  180  continue
       n = m - k2
       do 210 i = k1, l
          do 190 j = k1, m
             cy(j) = cxy(j,i)
  190     continue
          do 200 j = k1, n
             k = k1 + j
             cxy(j,i) = iac1r1(tx,cy,ncx,lx,k)
  200     continue
  210  continue
       return
      end subroutine iac2r
      
      
      double precision function ias2r(tx,nx,ty,ny,cxy,
     1 ncx,ncy,idx,idy,xx,yy)
cSm030221
         implicit none !double precision(a-h,o-z)
         integer idx,idy,ncx,ncy,nx,ny
         double precision cxy,tx,ty,xx,yy
cSm020321
c         dimension cxy(ncx,ncy),tx(ncx),ty(ncy)
          dimension cxy(ncx,*),tx(*),ty(*)

         external ias1r1,ias1r2
         integer i,j,k,k0,k1,k4,l,m,n,ias1r1
         double precision a,bsp,ias1r2
         dimension a(4),bsp(4)
         data k0/0/,k1/1/,k4/4/
       i = ias1r1(tx,nx,ncx,xx)
       j = ias1r1(ty,ny,ncy,yy)
       do 20 k = k1, k4
          l = i - k4 + k + k1
          do 10 m = k1, k4
             n = j - k4 + m + k1
             a(m) = cxy(l,n)
   10     continue
          bsp(k) = ias1r2(ty,a,ncy,k4,idy,yy,j,k0)
   20  continue
       ias2r = ias1r2(tx,bsp,ncx,k4,idx,xx,i,k0)
       return
      end function ias2r
      
      
c*****************************************************************
c*              End of Spline programs                           *
c*****************************************************************

c      modified spline function to use the arrasys with maximal
c      dimensions (nxa,nya) >(nx,ny)
      subroutine iac2r_Sm(x,nx,y,ny,fxy,nfx,nfy,lx,ly,mxa,fxa,mxb,fxb,
     1                 mya,fya,myb,fyb,tx,ty,cxy,ncx,ncy,nry,cy,
     1                 nxa,nx4a)
c-----------------------------------------------------------
c      The additional input parameters nxa,nx4a
c      are to use this subroutine when the arrays dimensions
c      are bigger then the used dimensions.
c      nxa is a maximal value nx: nxa.ge.nx
c      nx4x is a maximal value of nx4: nx4a.ge.nx4
c----------------------------------------------------------               
	implicit none
         integer nxa,nx4a
         integer lx,ly,mxa,mxb,mya,myb,ncx,ncy,nfx,nfy,nx,ny,nry
         double precision cxy,cy,fxa,fxb,fxy,fya,fyb,tx,ty,x,y
cSm030218
c         dimension cxy(ncx,ncy),cy(nry),fxa(ny),fxb(ny),fxy(nfx,nfy),
c     1             fya(nx),fyb(nx),tx(ncx),ty(ncy),x(nx),y(ny)
         dimension cxy(nx4a,*),cy(*),fxa(*),fxb(*),fxy(nxa,*),
     1             fya(*),fyb(*),tx(*),ty(*),x(*),y(*)

         external iac1r1,iac1r2,iac1r3,iac1r4
         integer i,j,k,k0,k1,k2,k3,l,m,n
          double precision a,bsp,fa,fb,iac1r1
         dimension a(3),bsp(3)
         data k0/0/,k1/1/,k2/2/,k3/3/
       call iac1r2(x,nx,ncx,mxa,tx)
       call iac1r2(y,ny,ncy,mya,ty)
       do 20 i = k1, nx
          j = k2 + i
          do 10 k = k1, ny
             l = k2 + k
             cxy(j,l) = fxy(i,k)
   10     continue
   20  continue
       m = nx + k2
       n = ny + k2
       if (mxa .gt. k2) go to 50
          call iac1r3(tx,nx,ncx,cy,lx,mxa,mxb,a,bsp)
          do 40 i = k3, n
             j = i - k2
             if (mxa .ge. k1) fa = fxa(j)
             if (mxb .ge. k1) fb = fxb(j)
             do 30 j = k3, m
                cy(j) = cxy(j,i)
   30        continue
             call iac1r4(tx,nx,ncx,cy,lx,mxa,fa,mxb,fb,a,bsp)
             cxy(k1,i) = cy(k1)
             cxy(k2,i) = cy(k2)
             k = m + k1
             cxy(k,i) = cy(k)
             k = k + k1
             cxy(k,i) = cy(k)
   40     continue
   50  continue
       if (mya .gt. k2) go to 110
          call iac1r3(ty,ny,ncy,cy,ly,mya,myb,a,bsp)
          do 70 i = k3, m
             j = i - k2
             if (mya .ge. k1) fa = fya(j)
             if (myb .ge. k1) fb = fyb(j)
             do 60 j = k3, n
                cy(j) = cxy(i,j)
   60        continue
             call iac1r4(ty,ny,ncy,cy,ly,mya,fa,myb,fb,a,bsp)
             cxy(i,k1) = cy(k1)
             cxy(i,k2) = cy(k2)
             k = n + k1
             cxy(i,k) = cy(k)
             k = k + k1
             cxy(i,k) = cy(k)
   70     continue
          if (mxa .gt. k2) go to 130
             call iac1r3(ty,ny,ncy,cy,ly,k0,k0,a,bsp)
             do 100 i = k1, k2
                k = k3 - i
                do 80 j = k3, n
                   cy(j) = cxy(k,j)
   80           continue
                call iac1r4(ty,ny,ncy,cy,ly,k0,fa,k0,fb,a,bsp)
                cxy(k,k1) = cy(k1)
                cxy(k,k2) = cy(k2)
                l = n + k1
                cxy(k,l) = cy(l)
                l = l + k1
                cxy(k,l) = cy(l)
                k = m + i
                do 90 j = k3, n
                   cy(j) = cxy(k,j)
   90           continue
                call iac1r4(ty,ny,ncy,cy,ly,k0,fa,k0,fb,a,bsp)
                cxy(k,k1) = cy(k1)
                cxy(k,k2) = cy(k2)
                l = n + k1
                cxy(k,l) = cy(l)
                l = l + k1
                cxy(k,l) = cy(l)
  100        continue
             m = m + k2
             n = n + k2
             go to 150
  110  continue
       k = k3
       if (mxa .lt. k3) k = k - k2
       if (mxa .lt. k3) m = m + k2
       do 120 i = k, m
          j = n - k1
          cxy(i,k1) = cxy(i,j)
          j = j - k1
          cxy(i,k2) = cxy(i,j)
          j = k3 + k1
          l = n + k1
          cxy(i,l) = cxy(i,j)
          j = j + k1
          l = l + k1
          cxy(i,l) = cxy(i,j)
  120  continue
  130  continue
       n = n + k2
       if (mxa .gt. k2) go to 150
          do 140 i = k1, n
             j = m - k1
             cxy(k1,i) = cxy(j,i)
             j = j - k1
             cxy(k2,i) = cxy(j,i)
             l = k3 + k1
             j = m + k1
             cxy(j,i) = cxy(l,i)
             l = l + k1
             j = j + k1
             cxy(j,i) = cxy(l,i)
  140     continue
          m = m + k2
  150  continue
       l = n - k2
       do 180 i = k1, m
          do 160 j = k1, n
             cy(j) = cxy(i,j)
  160     continue
          do 170 j = k1, l
             k = k1 + j
             cxy(i,j) = iac1r1(ty,cy,ncy,ly,k)
  170     continue
  180  continue
       n = m - k2
       do 210 i = k1, l
          do 190 j = k1, m
             cy(j) = cxy(j,i)
  190     continue
          do 200 j = k1, n
             k = k1 + j
             cxy(j,i) = iac1r1(tx,cy,ncx,lx,k)
  200     continue
  210  continue
      
       return
      end subroutine iac2r_Sm



      double precision function ias2r_Sm(tx,nx,ty,ny,cxy,
     1 ncx,ncy,idx,idy,xx,yy,nx4a)
c-----------------------------------------------------------
c      The additional input parameters nx4a
c      are to use this subroutine when the arrays dimensions
c      are bigger then the used dimensions:(nxa,nya).ge.(nx,ny)
c      nxa is a maximal value nx: nxa.ge.nx
c      nx4x is a maximal value of nx4: nx4a.ge.nx4
c----------------------------------------------------------
      implicit none !double precision(a-h,o-z) 
         integer idx,idy,ncx,ncy,nx,ny, nx4a
         double precision cxy,tx,ty,xx,yy
cSm030218
c         dimension cxy(ncx,ncy),tx(ncx),ty(ncy)
c         dimension cxy(ncx,*),tx(*),ty(*)
         dimension cxy(nx4a,*),tx(*),ty(*)
        

         external ias1r1,ias1r2
         integer i,j,k,k0,k1,k4,l,m,n,ias1r1
         double precision a,bsp,ias1r2
         dimension a(4),bsp(4)
         data k0/0/,k1/1/,k4/4/

c       write(*,*)'ncx,ncy,nx4a',ncx,ncy,nx4a
c       do i=1,nx4a
c         do j=1,nx4a
c         write(*,*)'i,j,cxy(i,j)',i,j,cxy(i,j)
c         enddo
c       enddo


       i = ias1r1(tx,nx,ncx,xx)
       j = ias1r1(ty,ny,ncy,yy)
       do 20 k = k1, k4
          l = i - k4 + k + k1
          do 10 m = k1, k4
             n = j - k4 + m + k1
             a(m) = cxy(l,n)
   10     continue
          bsp(k) = ias1r2(ty,a,ncy,k4,idy,yy,j,k0)
   20  continue
       ias2r_Sm = ias1r2(tx,bsp,ncx,k4,idx,xx,i,k0)
       return
      end function ias2r_Sm
