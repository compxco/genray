
c..................................................................
c
c   Include file for dskin.f, 
c     which reads diskf distribution function file from CQL3D.
c
c   parameters iya,jxa,lrza,ngena need to be set 
c     in accord with iy,jx,lrz,ngen in diskf file.
c
c   Note: You might need to know that storage is slightly 
c         different here than in CQL3D.
c
c   It uses pointers
c   Maximal dimentions were iya,jxa,lrza,ngena
c   Now jxa,lrza,ngena are not used.
c   and dimentions are iya,jx,lrz,ngen
c   Now jxa1 ia not used, the code uses jx1=jx+1  
c..................................................................

      real*8, pointer ::
     &x(:),                   !(jxa)
     &y(:,:),                 !(iya,lrza),
     &rovera(:),              !(lrza), 
     &rya(:),                 !(lrza), =rovera?
     &elecfld(:),             !(lrza),
     &bthr(:),                !(lrza),
     &btoru(:),               !(lrza),
     &bthr0(:),               !(lrza),
     &btor0(:),               !(lrza),
     &reden(:,:),             !(lrza,ngena),
     &temp(:,:),              !(lrza,ngena)
     &bnumb(:),               !(ngena),
     &fmass(:),               !(ngena)
     &f(:,:,:,:),             !(iya,jxa,lrza,ngena)  
     &df_dpitch(:,:,:,:),     !(iya1,jxa,lrza,ngena),iya1=iya+1
     &df_dx(:,:,:,:)          !(iya,jxa1,lrza,ngena),jxa1=jxa+1

      integer*4, pointer ::
     &iy_(:),                  !(lrza),
     &itl(:),                 !(lrza),
     &itu(:)                  !(lrza)

      real*8 radmin,vnorm,vmaxdvt,eovedd
      integer iy !max iy_ 
c------------------------------------------------
c    double precision x(jxa),y(iya,lrza),rovera(lrza),elecfld(lrza),
c     +       bthr(lrza),btoru(lrza),bthr0(lrza),btor0(lrza),
c     +       reden(lrza,ngena),temp(lrza,ngena)
c      real*8 radmin,vnorm,vmaxdvt,eovedd
c      double precision bnumb(ngena),fmass(ngena)
c      integer*4 iy_(lrza),itl(lrza),itu(lrza)
c      double precision f(iya,jxa,lrza,ngena)

c-----derivatives from the 3D distribution function
c     by the pitch angle df_dpich and the momentum  u=p/m df_dx
c     This derivatives are calculated as a numerical derivative in the
c     inter-medium points i+0.5
c     df_dx(i)=df_dx(x_(i-0.5))=f(i)-f(i-1))/(x(i)-x(i-1)) i=2,jx
c     df_dx(1)=df_dx(x(1))=0 is the derivative in the first point of x mesh
c     df_dx(jx+1)=df(x(jx))=0 is the derivative in the last point of x mesh

c      double precision  df_dpitch(iya1,jxa,lrza,ngena)
c      double precision  df_dx(iya,jxa1,lrza,ngena)
 
      COMMON /dskincomm/
     &x,                   !(jxa)
     &y,                   !(iya,lrza),
     &rovera,              !(lrza),
     &rya,                 !(lrza),
     &elecfld,             !(lrza),
     &bthr,                !(lrza),
     &btoru,               !(lrza),
     &bthr0,               !(lrza),
     &btor0,               !(lrza),
     &reden,               !(lrza,ngena),
     &temp,                !(lrza,ngena)
     &bnumb,               !(ngena),
     &fmass,               !(ngena)
     &f,                   !(iya,jxa,lrza,ngena)  
     &df_dpitch,           !(iya1,jxa,lrza,ngena)
     &df_dx,               !(iya,jxa1,lrza,ngena)
     &radmin,vnorm,vmaxdvt,eovedd,
     &iy_,                 !(lrza)
     &itl,                 !(lrza),
     &itu,                 !(lrza)
     &iy
c      COMMON /dskincomm/
c     & x,
c     & y,
c     & rovera,
c     & elecfld,bthr,btoru,bthr0,btor0,
c     & reden,temp,radmin,vnorm,vmaxdvt,eovedd,
c     & bnumb,fmass,itl,itu,
c     + f,df_dpitch,df_dx
c     + ,iy_ ! for test






  




