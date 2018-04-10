c The data for the parameters of the plasma species.
c They are used in hot non-relativistic dispersion.
c Here nbulka is the max number of the plasma species. 
c      integer nbulka !it is in param.i
c      parameter(nbulka=5)
      double precision nllc,massc,tec,xc,yc,tpopc,vflowc
      integer nbulkc  
      common/nperpcom/nllc,massc(nbulka),tec(nbulka),
     .xc(nbulka),yc(nbulka),tpopc(nbulka),vflowc(nbulka),
     .nbulkc
