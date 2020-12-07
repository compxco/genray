*****************************************************
c----------------subroutine MK_GRAP----------------*
c Prepares for file drawgenr.in for xdraw              *
c INPUT: from common blocks 'gr' and 'write'        *
c        iray -number of the ray at antenna         *
c*****************************************************
      SUBROUTINE MK_GRAP(iray)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      include 'param.i'
      INTEGER I,J
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'

      if(myrank.ne.0) return

      write(*,*)'mk_grap iray=',iray
      if(iray.eq.1) then
c         CALL ASSIGN("assign -F f77 -N ieee u:11",ier)
         OPEN(11,file='genray.bin',form='unformatted')
      end if
      write(*,*)'mk_grap nrayelt=',nrayelt
      if (nrayelt.eq.0)then
         goto 70
      endif
      if (iray.gt.1) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
      DO 10 I=1,nrayelt
       WRITE(11) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
10    CONTINUE
      WRITE(11)

      DO 30 J=1,3
       DO 20 I=1,NP+1
        WRITE(11) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
20     CONTINUE
       WRITE(11)
30    CONTINUE

      DO 50 J=4,NL
       DO 40 I=1,NP+1
        WRITE(11) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
40     CONTINUE
       WRITE(11)
50    CONTINUE
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
c      p=0.0
      DO 61 I=1,nrayelt
       WRITE(11) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(p),REAL(p),REAL(p),REAL(p),
c     2  REAL(p),REAL(p),REAL(p)
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
61    CONTINUE
      WRITE(11)
c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      CLOSE(11)
      RETURN
      END


c----------------subroutine MK_GRAPH----------------
c Prepares for file drawgenr.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c        iray -number of the ray at antenna        
c        n_wall is a number of points in arrays
c                which set the wall form
c*****************************************************
cSAP090313
c      SUBROUTINE MK_GRAPH(iray,nray,ifreq,nfreq,nbulk)
cSAP110211
c      SUBROUTINE MK_GRAPH(iray,nray,ifreq)  
      SUBROUTINE MK_GRAPH_old_110111(iray,nray,ifreq)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J

c-----input
      integer iray,nray,ifreq

      integer isave

cSm061209
      integer isave_freq

      save isave

cSm061209
      save isave_freq

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'
cSAP090313
      include 'one.i'
      include 'fourb.i'
      include 'three.i'

      data isave /0/
      data isave_freq/0/

      if(myrank.ne.0) return

      if (isave.lt.1) then
        isave=iray
      endif

cSm061229
      if(isave_freq.lt.1) then
        isave_freq=ifreq
      endif

      write(*,*)'!!!!!!in mk_graph iray,nrayelt,isave,isave_freq',
     + iray,nrayelt,isave,isave_freq

       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:82",ier)
         write(*,*)'mk_grph before open 82 and 81'
         OPEN(82,file='genray.bin',form='unformatted')
cBH111031         open(81,file='genray.doc')
         open(85,file='absorp.bin',form='unformatted')
         open(86,file='absorp.doc')
         open(75,file='efield.bin',form='unformatted')
         open(76,file='eps_r.bin',form='unformatted')
         open(77,file='eps_i.bin',form='unformatted')
         open(78,file='em_res.bin',form='unformatted')
         open(79,file='cd.bin',form='unformatted')
      end if
cSm030508
c      if (nrayelt.eq.0)then
      write(*,*)'nrayelt,iray,isave,ifreq',nrayelt,iray,isave,ifreq
      if (nrayelt.lt.2)then
        if((iray.gt.isave).or.(ifreq.gt.1)) then
          write(*,*)'before goto 70'
          goto 70  
        endif
        
      endif
     
c      write(*,*)'iray,isave,ifreq,isave_freq',
c     &            iray,isave,ifreq,isave_freq

      if (((iray.gt.isave).or.(ifreq.gt.1)).and.
     &    (ifreq.gt.isave_freq))then
c        data for noncentral antenna rays      
c         write(*,*)'before goto 60' 
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(15(1pe10.3))
 7    format(7(1pe10.3))
 2    format(/)
      write(*,*)'data for central ray nrayelt',nrayelt
      
      if(nrayelt.eq.0) nrayelt=1  ! YuP[07-2017] added: 
         ! could happen if the root Nper at the launch was not found.
         ! Set nrayelt to 1, to avoid out-of-bounds problem.
         
      DO 10 I=1,nrayelt
cSm070720
        if (nrayelt.eq.1) goto 11

c        write(*,*)'2before 82 I',I
c        write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),wye(I),wyi(I)'
c        write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
c     3  REAL(wyi(i))

c        write(*,*)'I,wal_emis(I)',I,wal_emis(I)
        WRITE(82) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),REAL(wye(I)),
     3  REAL(wyi(I)),real(wal_emis(I)),real(wj_emis(I))

cBH111031        WRITE(81,1)wr(I),wz(I),xarr(I),yarr(I),
cBH111031     1  ws(I),delpwr(I),rez(I),spsi(I),
cBH111031     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I),
cBH111031     3  wal_emis(I),wj_emis(I)
         
c       for ioin ion resonance for nbulk >2

c        write(*,*)'nbulk,I,wyi(I),wxi2(I),wyi2(I),wxi(I)',
c     &             nbulk,I,wyi(I),wxi2(I),wyi2(I),wxi(I)

c        if (nbulk.gt.2) then
        if (nbulk.gt.3) then
           p12=wyi(I)*wxi2(I)/(wyi2(I)*wxi(I))

           write(*,*)'p12',p12

           del=1.d0/(1.d0+p12)
           yii=wyi(I)*wyi2(I)*((1.d0-del)*wyi(I)+del*(wyi2(I)))/
     .          ((1.d0-del)*wyi2(I)+del*(wyi(I)))   
        else
           p12=0.d0
           del=1.d0
           yii=0.d0
        endif

c        write(*,*)'mk_graph bef 85'
                       
        WRITE(85)REAL(ws(I)),REAL(delpwr(I)),REAL(spsi(I)),REAL(wye(I)),
     3  REAL(wyi(I)),REAL(wyi2(I)),REAL(dsqrt(yii))

c        write(*,*)'mk_graph bef 86'

        WRITE(86,7)ws(I),delpwr(I),spsi(I),wye(I),
     3  wyi(I),wyi2(I),dsqrt(dabs(yii))

c        write(*,*)'mk_graph bef 75'

        write(75)real(ws(I)),
     +   real(dreal(cwexde(I))),real(dimag(cwexde(I))),
     +   real(dreal(cweyde(I))),real(dimag(cweyde(I))),
     +   real(dreal(cwezde(I))),real(dimag(cwezde(I))),
     +   real(fluxn(I)),REAL(spsi(I))

c        write(*,*)'mk_graph bef 76'

        write(76)real(ws(I)),
     &  real(dreal(w_ceps(1,1,I))),real(dreal(w_ceps(1,2,I))),
     &  real(dreal(w_ceps(1,3,I))),
     &  real(dreal(w_ceps(2,1,I))),real(dreal(w_ceps(2,2,I))),
     &  real(dreal(w_ceps(2,3,I))),
     &  real(dreal(w_ceps(3,1,I))),real(dreal(w_ceps(3,2,I))),
     &  real(dreal(w_ceps(3,3,I)))


c        write(*,*)'mk_graph bef 77'

        write(77)real(ws(I)),
     &  real(dimag(w_ceps(1,1,I))),real(dimag(w_ceps(1,2,I))),
     &  real(dimag(w_ceps(1,3,I))),
     &  real(dimag(w_ceps(2,1,I))),real(dimag(w_ceps(2,2,I))),
     &  real(dimag(w_ceps(2,3,I))),
     &  real(dimag(w_ceps(3,1,I))),real(dimag(w_ceps(3,2,I))),
     &  real(dimag(w_ceps(3,3,I)))

c        write(*,*)'mk_graph bef 78'

        WRITE(78) REAL(wr(I)),REAL(wye(I)),
     &  real(wp_perpmax_dmvt(I,2)),
     &  real(wp_parmin_dmvt(I,2)),real(wp_parmax_dmvt(I,2)),
     &  real(wp_0_dmvt(I,2)),real(wdel_p_dmvt(I,2)),
     &  REAL(wnpar(I))

c      write(*,*)'mk_graph i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)',
c     &i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)

c       write(*,*)'mk_graph bef 79'

       WRITE(79) REAL(ws(I)),REAL(wr(I)),REAL(spsi(I)),REAL(wye(I)),
     &  real(eff(I)),REAL(wnpar(I)),real(delpow_e_ar(I)),
     &  real(delcur_par_ar(I)),real(wtheta_pol(I))
     

10    CONTINUE
cyup      write(*,*)'after 10'
      WRITE(82)
      write(85)
      write(75)
      write(78)
      write(79)
cBH111031      WRITE(81,2)

cSm070720
 11   continue
      DO 30 J=1,3 
       DO 20 I=1,NP+1
         WRITE(82) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

cBH111031        WRITE(81,1)AR(J,I),AZ(J,I),XT(J,I),
cBH111031     1   YT(J,I),
cBH111031     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
cBH111031     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
cBH111031     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt),
cBH111031     5   wal_emis(nrayelt),wj_emis(nrayelt)
20     CONTINUE
       WRITE(82)
cBH111031       WRITE(81,2)
30    CONTINUE

      DO 50 J=4,NL
       DO 40 I=1,NP+1
         WRITE(82) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
    
cBH111031         WRITE(81,1)AR(J,I),AZ(J,I),
cBH111031     1   XT(3,NP+1),YT(3,NP+1),
cBH111031     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
cBH111031     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
cBH111031     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt),
cBH111031     5   wal_emis(nrayelt),wj_emis(nrayelt)
40     CONTINUE
       WRITE(82)
cBH111031       WRITE(81,2)
50    CONTINUE
     
cSAP090313 write wall coordinates

      if (n_wall.gt.0) then
         do i=1,n_wall  
         WRITE(82) 
     1    REAL(r_wall(i)*100.0),REAL(z_wall(i)*100.0),
     &    REAL(xarr(nrayelt)),REAL(yarr(nrayelt)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3    REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5    real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         enddo

         WRITE(82)

      endif

      if (n_wall.gt.1) then        
         do m=0,max_limiters          
            do i=1,n_wall_add(m)           
               WRITE(82) 
     &         REAL(r_wall_add(i,m)*100.0),
     &         REAL(z_wall_add(i,m)*100.0),
     &         REAL(xarr(nrayelt)),REAL(yarr(nrayelt)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5         real(wal_emis(nrayelt)),real(wj_emis(nrayelt))  
           enddo ! i=1,n_wall_add(m)       
 
           WRITE(82)
         enddo !m
      endif
         
cSAP090413     
      p=pi/180.d0
      do m=1,max_limiters
         r_min=1.d0
         do i=1,n_limiter(m)
            if(r_min.gt.r_limiter(i,m)) r_min=r_limiter(i,m)
         enddo
         r_max=rr(nxeqd)*100.d0
         r_min=r_min*100.d0
         write(*,*)'degree phi_limiter(1,m),phi_limiter(2,m)',
     &              phi_limiter(1,m),phi_limiter(2,m)
         write(*,*)'radian phi_limiter(1,m)*p,phi_limiter(2,m)*p',
     &              phi_limiter(1,m)*p,phi_limiter(2,m)*p

         x_lim_min=r_min*dcos(phi_limiter(1,m)*p)
         x_lim_max=r_max*dcos(phi_limiter(1,m)*p)
         y_lim_min=r_min*dsin(phi_limiter(1,m)*p)
         y_lim_max=r_max*dsin(phi_limiter(1,m)*p)
 
         write(*,*)'mk_graph r_min,r_max',r_min,r_max
         write(*,*)'x_lim_min,y_lim_min',x_lim_min,y_lim_min
         write(*,*)'x_lim_max,y_lim_max',x_lim_max,y_lim_max

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_min),REAL(y_lim_min),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_max),REAL(y_lim_max),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         WRITE(82)

         x_lim_min=r_min*dcos(phi_limiter(2,m)*p)
         x_lim_max=r_max*dcos(phi_limiter(2,m)*p)
         y_lim_min=r_min*dsin(phi_limiter(2,m)*p)
         y_lim_max=r_max*dsin(phi_limiter(2,m)*p)

         write(*,*)'x_lim_min,y_lim_min',x_lim_min,y_lim_min
         write(*,*)'x_lim_max,y_lim_max',x_lim_max,y_lim_max

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_min),REAL(y_lim_min),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_max),REAL(y_lim_max),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         WRITE(82)
      enddo !m
c-------------------------------------------------------------------          
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cSm070729
      if (nrayelt.eq.0) nrayelt=1  ! YuP[07-2017] added: 
         ! could happen if the root Nper at the launch was not found.
         ! Set nrayelt to 1, to avoid out-of-bounds problem.
         
      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt

c       write(*,*)'befor 82 i',i
c       write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I))'
c       write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I))


       WRITE(82) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
     3  REAL(wyi(I)),real(wal_emis(I)),real(wj_emis(I))

cBH111031       WRITE(81,1)wr(I),wz(I),xarr(I),yarr(I),
cBH111031     1  ws(I),delpwr(I),rez(I),spsi(I),
cBH111031     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I),
cBH111031     3  wal_emis(I),wj_emis(I)

        write(75)real(ws(I)),
     +   real(dreal(cwexde(I))),real(dimag(cwexde(I))),
     +   real(dreal(cweyde(I))),real(dimag(cweyde(I))),
     +   real(dreal(cwezde(I))),real(dimag(cwezde(I))),
     +   real(fluxn(I)),REAL(spsi(I))

        WRITE(78) REAL(wr(I)),REAL(wye(I)),
     &  real(wp_perpmax_dmvt(I,2)),
     &  real(wp_parmin_dmvt(I,2)),real(wp_parmax_dmvt(I,2)),
     &  real(wp_0_dmvt(I,2)),real(wdel_p_dmvt(I,2)),
     &  REAL(wnpar(I))

c      write(*,*)'mk_graph i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)',
c     &i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)

       WRITE(79) REAL(ws(I)),REAL(wr(I)),REAL(spsi(I)),REAL(wye(I)),
     &  real(eff(I)),REAL(wnpar(I)),real(delpow_e_ar(I)),
     &  real(delcur_par_ar(I)),real(wtheta_pol(I))
61    CONTINUE
c      write(*,*)'after 61'
      WRITE(82)
cBH111031      WRITE(81,2)
      write(75)
      WRITE(78)
      write(79) 

cSm070720 
 62   CONTINUE

c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      if (iray.eq.nray) then
       
c      write(*,*)'mk_graph iray,nraym,ifreq,nfreq',
c     & iray,nray,ifreq,nfreq
      if ((iray.eq.nray).and.(ifreq.eq.nfreq)) then
c       write(*,*)'mk_graph before close 82'
       CLOSE(82)
cBH111031       CLOSE(81)
       close(85)
       close(86)
       close(75)
       close(76)
       close(77)
       close(78)
       close(79)
c       write(*,*)'mk_graph after close 82 and 81'
      endif
      RETURN
      END

				    

c----------------subroutine MK_GRAPc----------------
c Prepares for file drawgenc.in for xdraw              
c  +  plots of contours 1/Yc_i=n
c INPUT: from common blocks 'gr' and 'write'        
c         iray -number of the ray at antenna         
c    	  iwopen and iwj are in common one.i         
c iwopen =1 mk_grapc will calculate open contours wb_c=n (using contrb1)
c         2 mk_grapc will calculate close contours wb_c=n (using contrb2)
c iwj    =mk_grapc will calculate contours wb_cj=n.
c         Here j gives the plasma  component, 
c        must be.le.nbulk, j=1 for the electron gyrofrequency
c*****************************************************
      SUBROUTINE MK_GRAPc(iray,iwopen,iwj)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'
      
      if(myrank.ne.0) return

c      write(*,*)'!!!!!!in mk_graphc iray,nrayelt',iray,nrayelt
      if(iray.eq.1) then
c         CALL ASSIGN("assign -F f77 -N ieee u:84",ier)
         OPEN(84,file='genrac.bin',form='unformatted')
         open(83,file='genrac.doc')

      end if
      if (nrayelt.eq.0)then
         goto 70
      endif
      if (iray.gt.1) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(11(1pe10.3))
 2    format(/)
      
      do i=1,nrayelt
       WRITE(83,1)wr(I),wz(I),xarr(I),yarr(I),
     1 ws(I),delpwr(I),rez(I),spsi(I),
     2 wnpar(I),wnper(I),salphal(I)
      enddo
      DO 10 I=1,nrayelt
cSm070720
        if (nrayelt.eq.1) goto 11

       WRITE(84) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
10    CONTINUE
      WRITE(84)
      WRITE(83,2)
cc      write(*,*)'in mk_graph before 30'
cc       read(*,*)

cSm070720
 11     continue

      DO 30 J=1,3
       DO 20 I=1,NP+1
        WRITE(84) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
        WRITE(83,1)AR(J,I),AZ(J,I),XT(J,I),
     1   YT(J,I),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt)
20     CONTINUE
       WRITE(84)
       WRITE(83,2)
30    CONTINUE
cc      write(*,*)'in mk_grapc after 30'
cc      read(*,*)
cc      write(*,*)'in mk_grapc before 50'
      DO 50 J=4,NL
       DO 40 I=1,NP+1
        WRITE(84) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt))
        WRITE(83,1)AR(J,I),AZ(J,I),
     1   XT(3,NP+1),YT(3,NP+1),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt)
40     CONTINUE
       WRITE(84)
       WRITE(83,2)
50    CONTINUE
c      write(*,*)'in mk_graphc after 50'
c	write(*,*)'in mk_graphc before call contrb1 iwopen,iwj',
c     *  iwopen,iwj 
        if (iwopen.eq.1) then
c	  open contours R=R(z), iwj are in common one.i
          call contrb1
	else 
          if(iwopen.eq.2) then
c           close contours rho=rho(theta), iwj are in common one
cSm030514
            write(*,*)'in mk_graph iwopen=',iwopen
            call contrb2
            write(*,*)'after contrb2'
          endif
 	endif
c	write(*,*)'in mk_graphc after call contrb1'
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
cc      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
c      p=0.0
cc      write(*,*)'in mk_graph before 61'
cSm070720
        if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt
       WRITE(84) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(p),REAL(p),REAL(p),REAL(p),
c     2  REAL(p),REAL(p),REAL(p)
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I))
cc       write(*,*)'in mk_graph 61 i',i
cc       WRITE(*,1)wr(I),wz(I),xarr(I),yarr(I),
cc     1  ws(I),delpwr(I),rez(I),spsi(I),
cc     2  wnpar(I),wnper(I),salphal(I)
       WRITE(83,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I)
61    CONTINUE
cc      write(*,*)'in mk_graph after 61 30'
      WRITE(84)
      WRITE(83,2)

cSm070720 
 62   CONTINUE


c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      CLOSE(84)
c      CLOSE(83)
      RETURN
      END


c----------------subroutine mk_gronetwo-------------*
c Prepares for file drawonet.in for xdraw            *
c INPUT: from common block  'onetwo.i'              *
c****************************************************
      SUBROUTINE mk_gronetwo
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      include 'onetwo.i'
cSAP090404
      include 'rho.i'
      include 'one.i'

      if(myrank.ne.0) return
  
c         CALL ASSIGN("assign -F f77 -N ieee u:11",ier)
      OPEN(11,file='onetwo.bin',form='unformatted')
cBH110313      OPEN(21,file='onetwo.dat')
20    format(8e12.4)
cBH110313      WRITE(21,*)'rho spower powden powden_e
cBH110313     1powden_i powden_cl'
      hrho=1.d0/(NR-1)

cSAP090404
      p=1.d0/(hrho*dsqrt(1.d4*areatot/pi)) !1/[cm]

      DO 10 I=1,NR
         rho=hrho*(I-1)
         WRITE(11) REAL(rho),REAL(spower(I)),REAL(powden(I)),
     1   REAL(powden_e(I)),REAL(powden_i(I)),REAL(powden_cl(I)),
cSAP090404
     &   real(spower(I)*1.d-7*p),real(rho*dsqrt(1.d4*areatot/pi))   
cBH110313         WRITE(21,20)rho,spower(I),powden(I),
cBH110313     1   powden_e(I),powden_i(I),powden_cl(I)
   
10    CONTINUE
      WRITE(11)
     
      CLOSE(11)
cBH110313      CLOSE(21)

      write(*,*)'mk_graph.f: mk_gronetwo created onewto.bin onetwo.dat'
    
      RETURN
      END


c----------------subroutine MK_GRAPT----------------*
c Prepares for file drawgenr.in for xdraw toray_data  *
c INPUT: from common blocks 'gr' and 'write'        *
c        iray -number of the ray at antenna         *
c****************************************************
      SUBROUTINE MK_GRAPT (iray,nray)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J

      integer isave
      save isave

      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank

      include 'gr.i'
      include 'write.i'

      if(myrank.ne.0) return

      if (isave.lt.1) then
        isave=iray
      endif

c      write(*,*)'!!!!!!in mk_graph iray,nrayelt,isave',
c     + iray,nrayelt,isave

c      if(iray.eq.1) then
       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:52",ier)
c         write(*,*)'mk_grph before open 52 and 51'
         OPEN(52,file='genrayt.bin',form='unformatted')
         open(51,file='genrayt.doc')
         open(55,file='absorpt.bin',form='unformatted')
         open(56,file='absorpt.doc')
      end if
      if (nrayelt.eq.0)then
         goto 70
      endif
c      if (iray.gt.1) then
      if (iray.gt.isave) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(11(1pe10.3))
 7    format(7(1pe10.3))
 2    format(/)
c      write(*,*)'in mk_graph nrayelt=',nrayelt,'iray',iray

      DO 10 I=1,nrayelt

c        write(*,*)'before 82 I',I
c        write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),wye(I),wyi(I)'
c        write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
c     3  REAL(wyi(i))

cSm070720
        if (nrayelt.eq.1) goto 11

        WRITE(52) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),REAL(wye(I)),
     3  REAL(wyi(I))
        WRITE(51,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I)
         
c       for ion resonance for nbulk >2
        p12=wyi(I)*wxi2(I)/(wyi2(I)*wxi(I))
        del=1.d0/(1.d0+p12)
        yii=wyi(I)*wyi2(I)*((1.d0-del)*wyi(I)+del*(wyi2(I)))/
     .   ((1.d0-del)*wyi2(I)+del*(wyi(I)))   
                   
        WRITE(55)REAL(ws(I)),REAL(delpwr(I)),REAL(spsi(I)),REAL(wye(I)),
     3  REAL(wyi(I)),REAL(wyi2(I)),REAL(dsqrt(yii))
        WRITE(56,7)ws(I),delpwr(I),spsi(I),wye(I),
     3  wyi(I),wyi2(I),dsqrt(dabs(yii))
10    CONTINUE
      WRITE(52)
      write(55)
      WRITE(51,2)
c      write(*,*)'in mk_graph before 30'
cc      read(*,*)

cSm070720
 11   continue

      DO 30 J=1,3
       DO 20 I=1,NP+1
        WRITE(52) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt))
        WRITE(51,1)AR(J,I),AZ(J,I),XT(J,I),
     1   YT(J,I),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt)
20     CONTINUE
       WRITE(52)
       WRITE(51,2)
30    CONTINUE
c      write(*,*)'in mk_graph before 50'
      DO 50 J=4,NL
       DO 40 I=1,NP+1
        WRITE(52) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt))
        WRITE(51,1)AR(J,I),AZ(J,I),
     1   XT(3,NP+1),YT(3,NP+1),
     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt)
40     CONTINUE
       WRITE(52)
       WRITE(51,2)
50    CONTINUE
c      write(*,*)'in mk_graph after 50'
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
cc      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cc      write(*,*)'in mk_graph before 61'

cSm070720
      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt

       WRITE(52) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
     3  REAL(wyi(I))
       WRITE(51,1)wr(I),wz(I),xarr(I),yarr(I),
     1  ws(I),delpwr(I),rez(I),spsi(I),
     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I)
61    CONTINUE
      WRITE(52)
      WRITE(51,2)
cSm070720 
 62   continue
c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
      if (iray.eq.nray) then
       CLOSE(52)
       CLOSE(51)
       close(55)
       close(56)


c       write(*,*)'mk_grapht  after close 52 and 51'
      endif
      RETURN
      END



c----------------subroutine mkgrtool-------------*
c Prepares for file tool.bin for xdraw              *
c INPUT: from common blocks: onetwo, one, 'five' *
c************************************************
      SUBROUTINE mkgrtool
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
            
      include 'onetwo.i'
      include 'one.i'
      include 'five.i'

      if(myrank.ne.0) return

c     CALL ASSIGN("assign -F f77 -N ieee u:11",ier)
      write(*,*)'in mkgrtool'
      OPEN(11,file='tool.bin',form='unformatted')
cBH110313      OPEN(21,file='tool.dat')
20    format(8e12.4)
      
      hrho=1.d0/(NR-1)
cBH110313      WRITE(21,*) 'rho,dens_e,temp_e,bmod,ye,
cBH110313     1 xe,ur,cutoffp,cutoffm,rx'
      thetax=0.d0

      write(*,*)' in mkgrtool before do 30'

c-----write data to tool.bin and tool.dat files

      DO 30 I=1,NR
         rho=hrho*(I-1)

         psix=psi_rho(rho)
         call zr_psith(psix,thetax,zx,rx)
         bmod=b(zx,rx,0.d0)
         ye=y(zx,rx,0.d0,1)
         xe=x(zx,rx,0.d0,1)
         ur=xe+ye**2
         cnpar=0.d0
         cutoffp=xe/(1.d0+ye)+cnpar**2
         if (ye.ne.0.d0) then
          cutoffm=xe/(1.d0-ye)+cnpar**2
         else
          cutoffm=0.d0
         endif
         
         WRITE(11) REAL(rho),REAL(densrho(rho,1)),
     1   REAL(temperho(rho,1)),REAL(bmod),REAL(ye),
     1   REAL(xe),REAL(ur),REAL(cutoffp),REAL(cutoffm),
     1   REAL(rx)
    
cBH110313         WRITE(21,20)rho,densrho(rho,1),
cBH110313     1   temperho(rho,1),bmod,ye,
cBH110313     1   xe,ur,cutoffp,cutoffm,
cBH110313     1   rx
 30   CONTINUE  

      write(*,*)' in mkgrtool after do 30'

      WRITE(11)
     
      CLOSE(11)
cBH110313      CLOSE(21)
     
      OPEN(11,file='tool1.bin',form='unformatted')
cBH110313      OPEN(21,file='tool1.dat')

c-----data for x-mode optimal parameters
      OPEN(31,file='tool2.bin',form='unformatted')
      OPEN(41,file='tool2.dat')
c---------------------------------------


      N_r=50
c      write(*,*)'rmax,rmin',rmax,rmin
      hr=(rmax-rmin)/(N_r-1)
cBH110313      WRITE(21,*) 'r,dens_e,temp_e bmod,y_e,x_e,ur,cutoofp,
cBH110313     1cutoffm,omega_p/omega_b'
 
      nmax=0
      mmax=0

      write(*,*)' in mkgrtool before do 10'

      DO 10 I=1,N_r
         r=rmin+hr*(I-1)
c         write(*,*)'i,r',i,r
         bmod=b(0.d0,r,0.d0)
c         write(*,*)'bmod',bmod
         ye=y(0.d0,r,0.d0,1)
         xe=x(0.d0,r,0.d0,1)
         ur=xe+ye**2
c         write(*,*)'ye,xe,ur',ye,xe,ur

         cnpar=0.d0
         cutoffp=xe/(1.d0+ye)+cnpar**2
         if (ye.ne.0.d0) then
          cutoffm=xe/(1.d0-ye)+cnpar**2
         else
          cutoffm=0.d0
         endif

c--------sqrt(x(r))=omega_p(r)/omega=eta o-mode cutoff for eta=1. 
c        y(r)=omega_b(r)/omega=1/n  'n' harmonics EC resonance condition
         eta=0.95d0  


         WRITE(11) REAL(r),REAL(densrho(rho,1)),
     1   REAL(temperho(rho,1)),REAL(bmod),REAL(ye),
     1   REAL(xe),REAL(ur),REAL(cutoffp),REAL(cutoffm),
     1   REAL(dsqrt(xe)/(eta*ye))
        
cBH110313         WRITE(21,20)r,densrho(rho,1),
cBH110313     1   temperho(rho,1),bmod,ye,xe,ur,cutoffp,cutoffm,
cBH110313     1   dsqrt(xe)/(eta*ye)
         
         n=dsqrt(xe)/(eta*ye)
         if(n.gt.nmax) then
           nmax=n
           r_l=r
         endif        
         
         WRITE(31) REAL(r),REAL(densrho(rho,1)),
     1   REAL(temperho(rho,1)),REAL(bmod),REAL(ye),
     1   REAL(xe),REAL(ur),REAL(cutoffp-1.d0),REAL(cutoffm-1.d0),
     1   REAL(1.d0/(xe/(1.d0-cnpar**2)-1.d0))
        
         WRITE(41,20)r,densrho(rho,1),
     1   temperho(rho,1),bmod,ye,xe,ur,cutoffp-1.,cutoffm-1.,
     1   (1.d0/(xe/(1.d0-cnpar**2)-1.d0))

         m=1.d0/(xe/(1.d0-cnpar**2)-1.d0)
         
         if(iabs(m).gt.mmax) then
           mmax=iabs(m)

           if(m.lt.0)then
             i_mmax=-1
           else
             i_mmax=1
           endif

           r_lm=r
c           write(*,*)'mmax,r_lm',mmax,r_lm
         endif
c         write(*,*)'0 r_lm',r_lm
     
10    CONTINUE

      write(*,*)' in mkgrtool after do 10'

c      write(*,*)'1 r_lm',r_lm     
      WRITE(11)
      WRITE(31)

      CLOSE(11)
cBH110313      CLOSE(21)

      CLOSE(31)
      CLOSE(41)
c      write(*,*) 'mk_graph.f mkgrtool nmax',nmax
c      write(*,*) 'mk_graph.f mkgrtool mmax,i_mmax',mmax,i_mmax
      z=0.d0
c      write(*,*)'2 r_lm',r_lm    

c      write(*,*)' in mkgrtool before do i=1,nmax nmax=',nmax

c      write(*,*)'eta',eta
 
      goto 100

      do i=1,nmax
       write(*,*)'i',i
       call solvropt(i,z,eta,r_l,rmax,r_opt_o)
       bmod=b(z,r_opt_o,0.d0)
c       write(*,*)'mk_graph.f mkgrtool i,r_opt_o,rho',i,r_opt_o,rho
       
       den_opt=dense(z,r_otp_o,0.d0,1) ! the optimal density
c------omega_optimal=omega_pe/eta
       friq_o=dsqrt(806.2*den_opt/eta**2)
c       write(*,*)'den_opt,friq_o',den_opt,friq_o

      enddo

      write(*,*)' in mkgrtool after do i=1,nmax'


c      write(*,*)'3 r_lm',r_lm     
      mmax=mmax*i_mmax
      k1=i_mmax
      k2=i_mmax
      cnpar=0.d0

      write(*,*)' in mkgrtool before do i=k1,nmax,k2'

      do i=k1,mmax,k2
        write(*,*)'i',i
c       write(*,*)'mk_graph before solvropx i=',i
c       write(*,*)'r_lm,rmax',r_lm,rmax
       call solvropx(i,z,r_lm,rmax,cnpar,r_opt_x)
       bmod=b(z,r_opt_x,0.d0)
c       write(*,*)'mk_graph.f mkgrtool i,r_opt_x,rho',i,r_opt_x,rho
       
       den_opt=dense(z,r_otp_x,0.d0,1)
c------omega_optimal=omega_pe/eta
       friq_x=dsqrt(806.2*den_opt/((1-cnpar**2)*(1.d0+1.d0/i)))
c       write(*,*)'den_opt,friq_x',den_opt,friq_x

      enddo

c      write(*,*)' in mkgrtool after do i=k1,nmax,k2'
 100  continue
      write(*,*) 'end of mkgrtool'

      RETURN
      END


      subroutine solvropt(n,z,eta,r_left,r_right,r_opt_o)
c     calculates the root r_opt_o from the equation( for o-mode)
c     sqrt(x_e(r,z)/(eta*y_e(r,z))=n
c     using the binary method
c     on (r_left<r<r_right) interval
c     input:
c     n is a number of harmonics,
c     z is vericle variable
c     eta is a parameter ~<1,
c     r_left,r_right are the boundaries of the interval
c
c     output:
c     r_opt_o is the root of the equation 
      
      IMPLICIT double precision (a-h,o-z)
      
      include 'param.i'
      include 'pgconst.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
cSm030514      
c      common /cf_r_o/z_loc,n_loc,eta_loc
      common /cf_r_o/z_loc,eta_loc,n_loc

      external f_r_o 

      z_loc=z
      n_loc=n
      eta_loc=eta

      racc=1.d-6
      
      r_opt_o=rtbis(f_r_o,r_left,r_right,racc)
           
      return
      end 

 
      double precision function f_r_o(r)
      implicit double precision (a-h,o-z)      
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
cSm030514
c      common /cf_r_o/z_loc,n_loc,eta_loc
      common /cf_r_o/z_loc,eta_loc,n_loc

c-----input
c     r is a major radius
c     z_loc is a vertical distance: from common /cf_r_o/
c     n_loc is the EC hatmonic number: from common /cf_r_o/
c     eta_loc ia a parameter xe=eta**2: from common /cf_r_o/

c      write(*,*)'f_r_o ,z_loc,n_loc,eta_loc,r',z_loc,n_loc,eta_loc,r
      z_loc1=z_loc
      bmod=b(z_loc1,r,0.d0)      
      ye=y(z_loc1,r,0.d0,1)
      xe=x(z_loc1,r,0.d0,1)
      
      n_loc1=n_loc
      eta_loc1=eta_loc
      if (ye.ne.0.d0) then
        f_r_o=dsqrt(xe)/(eta_loc1*ye)-n_loc1
      else
        f_r_o=1000.d0
      endif

      return
      end

      subroutine solvropx(n,z,r_left,r_right,cnpar,r_opt_x)
c     calculates the root r_opt_o from the equation( for x-mode)
c     1/(x/(1-N_par**2)-1)=n
c     using the binary method
c     on (r_left<r<r_right) interval
c     input:
c     n is a number og harmonics,
c     z is vericle variable
c     r_left,r_right are the boundaries of the interval
c     cnpar is N_par
c     output:
c     r_opt_o is the root of the equation 
      
      IMPLICIT double precision (a-h,o-z)
      
      include 'param.i'
      include 'pgconst.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      
c      common /cf_r_x/z_loc,n_loc,cnpar_l
      common /cf_r_x/z_loc,cnpar_l,n_loc


      external f_r_x 

      z_loc=z
      n_loc=n
      cnpar_l=cnpar 

      racc=1.d-6

      r_opt_x=rtbis(f_r_x,r_left,r_right,racc)

           
      return
      end

      double precision function f_r_x(r)
      implicit double precision (a-h,o-z)      
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
c      common /cf_r_x/z_loc,n_loc,cnpar_l
      common /cf_r_x/z_loc,cnpar_l,n_loc


c-----input
c     r is a major radius
c     z_loc is avertical distance: from common /cf_r_x/
c     n_loc is the EC hatmonic number: from common /cf_r_x/
c     cnpar_l is N_parrael from common/cf_r_x/

c      write(*,*)'f_r_x ,z_loc,n_loc,cnpar_l,r',
c     *z_loc,n_loc,cnpar_l,r
      z_loc1=z_loc
      bmod=b(z_loc1,r,0.d0)      
      ye=y(z_loc1,r,0.d0,1)
      xe=x(z_loc1,r,0.d0,1)
      
      n_loc1=n_loc
      cnpar_l1=cnpar_l

      if((1.d0-cnpar_l1**2.ne.0.d0).and.
     *((xe/(1.d0-cnpar_l1**2)-1).ne.0.d0))then
        f_r_x=1/(xe/(1.d0-cnpar_l1**2)-1)-n_loc1
      else
        f_r_x=1000.d0
      endif

      return
      end




c----------------subroutine MK_GR3----------------*
c Prepares for file drawgr3d.in for xdraw              *
c INPUT: from common blocks 'gr' and 'write'        *
c         iray -number of the ray at antenna         *
c*****************************************************
      SUBROUTINE MK_GR3d(iray,nray)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J

      integer isave
      save isave

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'

      if(myrank.ne.0) return

      if (isave.lt.1) then
        isave=iray
      endif

c      write(*,*)'!!!!!!in mk_gr3d iray,nrayelt,isave',
c     + iray,nrayelt,isave

       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:82",ier)
c         write(*,*)'mk_gr3d before open'
         OPEN(71,file='gr3d_1.bin',form='unformatted')
         open(73,file='gr3d_2.bin',form='unformatted')         
cBH111031         open(72,file='gr3d_1,doc')
cBH111031         open(74,file='gr3d_2,doc')
      end if
 11   format(10(1pe10.3))
 9    format(9(1pe10.3))
 1    format(/)

      if (nrayelt.eq.0)then
         goto 70
      endif

      write(*,*)'mk_gr3d iray,myrank,nrayelt',iray,myrank,nrayelt
      if (iray.gt.isave) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
cSm070720
      if (nrayelt.eq.1) goto 12

      DO 10 I=1,nrayelt

        WRITE(71) real(ws(I)),REAL(seikon(I)),real(spsi(I)),
     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))
  
cBH111031        write(72,11)real(ws(I)),REAL(seikon(I)),real(spsi(I)),
cBH111031     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
cBH111031     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))

        write(73)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
     +  REAL(salphac(I)),REAL(salphal(I))

cBH111031        write(74,9)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
cBH111031     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
cBH111031     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
cBH111031     +  REAL(salphac(I)),REAL(salphal(I))

10    CONTINUE
      WRITE(71)
      write(73)
cBH111031      write(72,1)
cBH111031      write(74,1)

cSm070720
 12   continue

      goto 70
c  end data for the central ray
c--------------------------------------------------------------
60    continue
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cSm070720
      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt
        WRITE(71) real(ws(I)),REAL(seikon(I)),real(spsi(I)),
     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))
  
cBH111031        write(72,11)real(ws(I)),REAL(seikon(I)),real(spsi(I)),
cBH111031     +  REAL(wr(I)),real(wphi(I)),REAL(wz(I)),
cBH111031     2  REAL(wnpar(I)),REAL(wnper(I)),real(delpwr(I)),real(sdpwr(I))

        write(73)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
     +  REAL(salphac(I)),REAL(salphal(I))

cBH111031        write(74,9)real(ws(I)),real(wdnpar(I)),real(dreal(cwexde(I))),
cBH111031     +  real(dreal(cweyde(I))),real(dreal(cwezde(I))),
cBH111031     +  real(fluxn(I)),real(sbtot(I)),real(sene(I)),
cBH111031     +  REAL(salphac(I)),REAL(salphal(I))




61    CONTINUE
      
      WRITE(71)
      write(73)
cBH111031      write(72,1)
cBH111031      write(74,1)

cSm070720
 62   continue

c  end data for the noncentral rays
c---------------------------------------------------------------------
70    continue
      if (iray.eq.nray) then
       CLOSE(71)
cBH111031       CLOSE(72)
       CLOSE(73)
cBH111031       CLOSE(74)      
      endif

      write(*,*)'end MK_GR3d' 

      RETURN
      END

	


c----------------subroutine MK_gremis----------------
c Prepares for file drawem.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c         iray -number of the ray at antenna         
c*****************************************************
      subroutine mk_gremis(iray,nray,ifreq,nfreq)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J

      integer isave
cSm070108
      integer isave_freq

      save isave  
      save isave_freq

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'

      data isave /0/
      data isave_freq/0/

      if(myrank.ne.0) return

      if (isave.lt.1) then
        isave=iray
      endif

cSm061229
      if(isave_freq.lt.1) then
        isave_freq=ifreq
      endif

c      write(*,*)'!!!!!!in mk_gremis iray,nrayelt_emis,isave',
c     + iray,nrayelt,isave

       if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:90",ier)
c         write(*,*)'mk_gremis before open 90 and 91'
         OPEN(90,file='emis.bin',form='unformatted')
         open(91,file='emis.doc') 
      end if

c      write(*,*)'nrayelt_emis,iray,isave,ifreq',
c     &           nrayelt_emis,iray,isave,ifreq
cSm070198
c      if (nrayelt_emis.eq.0)then
c      if (nrayelt_emis.lt.2)then
      if (nrayelt_emis.lt.3)then
        if((iray.gt.isave).or.(ifreq.gt.1)) then
c         write(*,*)'before goto 70'
         goto 70
        endif
      endif
 
c      if (iray.gt.1) then
c      if (iray.gt.isave) then
cSm070108
c      if ((iray.gt.isave).or.(ifreq.gt.1)) then 

c        write(*,*)'iray,isave,ifreq,isave_freq',
c     &             iray,isave,ifreq,isave_freq
      

       if (((iray.gt.isave).or.(ifreq.gt.1)).and.
     &    (ifreq.gt.isave_freq))then
c        data for noncentral antenna rays
c         write(*,*)'before goto 60' 
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(13(1pe10.3))
 2    format(/)
    
c      write(*,*)'data for central ray nrayelt_emis',nrayelt_emis
cSm070720
      if (nrayelt_emis.lt.3) goto 11

      DO 10 I=1,nrayelt_emis-1
c       DO 10 I=1,nrayelt_emis

         ds=wsn(I+1)-wsn(I)

        if(ds.gt.1.d-18) then
          win_sn_ds =win_sn(I)/ds
          win_0_ds=win_0(I)/ds
        else
          ds_old=wsn(I)-wsn(I-1)
          win_sn_ds =win_sn(I-1)/ds_old
          win_0_ds=win_0(I-1)/ds_old
        endif

c        write(*,*)'I,wsn(I+1),wsn(I),ds',I,wsn(I+1),wsn(I),ds
c        write(*,*)'wr_em(I),wz_em(I),wphi_em(I)',
c     &  wr_em(I),wz_em(I),wphi_em(I)
c        write(*,*)'wsn(I),wal_emis(I),wj_emis(I),wnray(I)',
c     &  wsn(I),wal_emis(I),wj_emis(I),wnray(I)
cSm050913
c        write(*,*)'win_0(I)',win_0(I)
c        write(*,*)'win_0(I)/ds',win_0(I)/ds 
c        write(*,*)'win_sn_ds,win_0_ds',win_sn_ds,win_0_ds
c        write(*,*)'wtemp_em(I),wtemp_rad_em(I)',
c     &             wtemp_em(I),wtemp_rad_em(I)  
c        write(*,*)'wi_0sn(I),wtaun_em(I),win_0(I)',
c     &             wi_0sn(I),wtaun_em(I),win_0(I)

c        write(*,*)'wtemp_em(I),wtemp_rad_em(I),wi_0sn(I),wtaun_em(I)',  
c     &  wtemp_em(I),wtemp_rad_em(I),wi_0sn(I),wtaun_em(I)
  
        WRITE(90) REAL(wr_em(I)),REAL(wz_em(I)),real(wphi_em(I)),
     +  REAL(wsn(I)),
     +  REAL(wal_emis(I)),REAL(wj_emis(I)),REAL(wnray(I)),
cSm050913
c     +  real(win_sn(I)/ds),real(win_0(I)/ds),
     +  real(win_sn_ds),real(win_0_ds),
     +  REAL(wtemp_em(I)),real(wtemp_rad_em(I)),
     +  real(wi_0sn(I)),real(wtaun_em(I)),
     +  real(win_0(I)),
cSm070120
     +  real(wj_emis(I)*dexp(-wtaun_em(I)))

        WRITE(91,1) wr_em(I),wz_em(I),wphi_em(I),
     +  wsn(I),
     +  wal_emis(I),wj_emis(I),wnray(I),
cSm050913
c     +  win_sn(I)/ds,win_0(I)/ds,
     +  win_sn_ds,win_0_ds,
     +  wtemp_em(I),wtemp_rad_em(I),
     +  wi_0sn(I),wtaun_em(I)
              
10    CONTINUE
c      write(*,*)'after 10'
      WRITE(90)
      WRITE(91,2)

cSm070720
 11   continue

      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue

c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
c      write(*,*)'!!!!mk_gremis after 60,nrayelt_emis=',nrayelt_emis

cSm070720
      if (nrayelt_emis.lt.3) goto 62
 
      DO 61 I=1,nrayelt_emis-1
c      DO 61 I=1,nrayelt_emiscSm070720
      if (nrayelt_emis.lt.3) goto 62
 

        ds=wsn(I+1)-wsn(I)

        if(ds.gt.1.d-18) then
          win_sn_ds =win_sn(I)/ds
          win_0_ds=win_0(I)/ds
        else
          ds_old=wsn(I)-wsn(I-1)
          win_sn_ds =win_sn(I-1)/ds_old
          win_0_ds=win_0(I-1)/ds_old
        endif

c        write(*,*)'I,wsn(I+1),wsn(I),ds',I,wsn(I+1),wsn(I),ds
c        write(*,*)'wr_em(I),wz_em(I),wphi_em(I)',
c     &  wr_em(I),wz_em(I),wphi_em(I)
c        write(*,*)'wsn(I),wal_emis(I),wj_emis(I),wnray(I)',
c     &  wsn(I),wal_emis(I),wj_emis(I),wnray(I)
c        write(*,*)'win_sn_ds,win_0_ds',
c     &  win_sn_ds,win_0_ds
c        write(*,*)'wtemp_em(I),wtemp_rad_em(I),wi_0sn(I),wtaun_em(I)',  
c     &  wtemp_em(I),wtemp_rad_em(I),wi_0sn(I),wtaun_em(I)
  
        WRITE(90) REAL(wr_em(I)),REAL(wz_em(I)),real(wphi_em(I)),
     +  REAL(wsn(I)),
     +  REAL(wal_emis(I)),REAL(wj_emis(I)),REAL(wnray(I)),
cSm050913
c     +  real(win_sn(I)/ds),real(win_0(I)/ds),
     +  real(win_sn_ds),real(win_0_ds),
     +  REAL(wtemp_em(I)),real(wtemp_rad_em(I)),
     +  real(wi_0sn(I)),real(wtaun_em(I))
     +  ,real(win_0(I)),
cSm070120
     +  real(wj_emis(I)*dexp(-wtaun_em(I)))

        WRITE(91,1) wr_em(I),wz_em(I),wphi_em(I),
     +  wsn(I),
     +  wal_emis(I),wj_emis(I),wnray(I),
     +  win_sn(I)/ds,win_0(I)/ds,
     +  wtemp_em(I),wtemp_rad_em(I),
     +  wi_0sn(I),wtaun_em(I)
  
61    CONTINUE
c      write(*,*)'after 61'
      WRITE(90)
      WRITE(91,2)

cSm070720
 62   continue

c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      if (iray.eq.nray) then
c      write(*,*)'mk_graph iray,nraym,ifreq,nfreq',
c     & iray,nray,ifreq,nfreq
      if ((iray.eq.nray).and.(ifreq.eq.nfreq)) then
c       write(*,*)'mk_graph before close 90'
       CLOSE(90)
       CLOSE(91)     

c       write(*,*)'mk_gremis after close 90 and 91'
      endif
      RETURN
      END


c----------------subroutine mk_gremfr----------------
c Prepares for file drawemfr.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c         iray -number of the ray at antenna         
c*****************************************************
      subroutine mk_gremfr(iray,nray)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J

      integer isave
      save isave

      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'write.i'

      if(myrank.ne.0) return

      if (isave.lt.1) then
        isave=iray
      endif

c      write(*,*)'!!!!!!in mk_gremfr iray,isave',
c     + iray,isave

      if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:94",ier)
c         write(*,*)'mk_gremfr before open 93 and 94'
         OPEN(95,file='emfr.bin',form='unformatted')
         open(93,file='emfr.doc') 
      end if
    
      if (iray.gt.isave) then
c        data for noncentral antenna rays
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(10(1pe10.3))
 2    format(/)
      write(*,*)'in mk_gremfr iray',iray,'nfreq',nfreq
      write(*,*)'mk_gremf freqncy0',freqncy0


      DO 10 I=1,nfreq

c        write(*,*)'in mk_gremfr I,wfreq(I),wtemp_rad_fr(I),
c     +  wi_0(iray,I),wtemp_rad_fr_wall(I)',
c     +  I,wfreq(I),wtemp_rad_fr(I),wi_0(iray,I),
c     +  wtemp_rad_fr_wall(I)

c        write(*,*)'in mk_gremfr wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
c     +  wrho0_em(I)',
c     +  wtemp_pl_fr(I),
c     +  wr0_em(I),wz0_em(I),wrho0_em(I)

        WRITE(95) REAL(wfreq(I)/freqncy0),Real(wtemp_rad_fr(I)),
     +  REAL(wi_0(iray,I)),real(wtemp_rad_fr_wall(I)),
     +  real(wtemp_pl_fr(I)),real(wr0_em(I)),real(wz0_em(I)),
     +  real(wrho0_em(I)),real(wr_2nd_harm(I)),
     +  real(wtemp_2nd_harm(I)),real(wtau_em(iray,I))      

        WRITE(93,1) wfreq(I)/freqncy0,wtemp_rad_fr(I),wi_0(iray,I),
     +  wtemp_rad_fr_wall(I),wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
     +  wrho0_em(I),wr_2nd_harm(I),wtemp_2nd_harm(I),
     +  wtau_em(iray,I)    


10    CONTINUE
      WRITE(95)
      WRITE(93,2)
     
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue

c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
      
      DO 61 I=1,nfreq
        
        WRITE(95) REAL(wfreq(I)/freqncy0),REAL(wtemp_rad_fr(I)),
     +  REAL(wi_0(iray,I)),real(wtemp_rad_fr_wall(I)),
     +  real(wtemp_pl_fr(I)),real(wr0_em(I)),real(wz0_em(I)),
     +  real(wrho0_em(I)),real(wr_2nd_harm(I)),
     +  real(wtemp_2nd_harm(I)),real(wtau_em(iray,I))    
        WRITE(93,1) wfreq(I),wtemp_rad_fr(I),wi_0(iray,I),
     +  wtemp_rad_fr_wall(I),wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
     +  wrho0_em(I),wr_2nd_harm(I),wtemp_2nd_harm(I),
     +  wtau_em(iray,I)
    
61    CONTINUE
      WRITE(95)
      WRITE(93,2)
c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue

      if (iray.eq.nray) then
         CLOSE(95)
         CLOSE(93)
         write(*,*)'mk_gremfr after close 95 and 93'
      endif

      RETURN
      END


      subroutine map_b(n_r,n_z)
c-----creates 2D array B(r,z) for plotting B contours at the mesh
c     r_min<r_i<r_max
c     input:
c     n_r n_z are the number of r and z mesh points
c
c     from five.i
c     rmin ,rmax are the min and max values of the major radius

c                at the limitter  
c     zmax,zmin are the min and max values of vertical coordinate z
c                at the limitter
      !implicit none
      implicit double precision (a-h,o-z)
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'five.i'
      include 'three.i'     
c-----input
      integer n_r,n_z
c      double precision rmin,rmax,zmax,zmin
c-----external
c     zpzmlim
      real*8  b
c-----local
      real*8 step_r,step_z,r,z,zp,zm,phi
      integer i,j

c     for pgplot
      parameter(n_r_a=100,n_z_a=100) !are the max number of mesh points
      parameter(n_contour=10) !is the number of contours
c      parameter(n_rho_a=10)! is the number of points in small radius
c      parameter(n_rho_a=201)! is the number of points in small radius
      parameter(n_rho_a=ndensa)! is the number of points in small radius
      parameter(n_contour_wide=30) !is the number of contours

      REAL*4 f2d_b(n_r_a,n_z_a), ! magnetic field
     . f2d_bp(n_r_a,n_z_a), ! poloidal magnetic field
     . f2d_bt(n_r_a,n_z_a),  ! toroidal magnetic field
     . f2d_psi(n_r_a,n_z_a),  !poloidal flux 
     . f2d_xe(n_r_a,n_z_a), f2d_ye(n_r_a,n_z_a) !xe,ye
      REAL*4 x_axis(n_r_a),y_axis(n_z_a) ! mesh
      
      REAL*4 contour_b(n_contour),contpur_bt(n_contour),
     .contour_bp(n_contour) !the values of contours 
      REAL*4 contour_wide(n_contour_wide)
      real*8 f1d_bz(n_r_a),f1d_bt(n_r_a),dx_axis(n_r_a),
     .     f1d_xe(n_r_a),f1d_ye(n_r_a)
      real*8 f1d_dens_e(n_rho_a),f1d_temp_e(n_rho_a),
     .f1d_rho(n_rho_a) 
      real*8 sign_dp
      
      if(myrank.ne.0) return
              
c      write(*,*)'map_b rmax,rmin,zmax,zmin',rmax,rmin,zmax,zmin
      step_r=(rmax-rmin)/(n_r-1)
      step_z=(zmax-zmin)/(n_z-1)
      phi=0.d0
      
      if (n_r.gt.n_r_a) then
         write(*,*)'in file mk_graph.f in map_b n_r.gt.n_r_a'
         write(*,*)'n_r=',n_r,'n_r_a=',n_r_a
         write(*,*)'change parameter n_r_a in map_b'
         stop
      endif

      if (n_z.gt.n_z_a) then
         write(*,*)'in file mk_graph.f in map_b n_z.gt.n_z_a'
         write(*,*)'n_z=',n_z,'n_z_a=',n_z_a
         write(*,*)'change parameter n_z_a in map_b'
         stop
      endif
      
      do i=1,n_r
         x_axis(i)=rmin+step_r*(i-1)
c         dx_axis(i)=x_axis(i)
c         write(*,*)'i,x_axis(i),dx_axis(i)',i,x_axis(i),dx_axis(i)
      enddo

      write(*,*)'map_b n_r',n_r
      write(*,*)'x_axis',x_axis
 

      do j=1,n_z
         y_axis(j)=zmin+step_z*(j-1)
      enddo

      write(*,*)'map_b n_z',n_z
      write(*,*)'y_axis',y_axis
              
      do i=1,n_r
c         r=x_axis(i)
         r=rmin+step_r*(i-1)
c--------zp zm are the values of the top and bottom points
c        of the limitter at the given r 

         write(*,*)'map_b i,r',i,r
         call zpzmlim(r,zp,zm)

         write(*,*)'zp,zm',zp,zm

         bmod=b(0.d0,r,phi)

         write(*,*)'bmod,bz,br',bmod,bz,br

         f1d_bz(i)=bz
         f1d_bt(i)=bphi
         f1d_xe(i)=x(0.d0,r,phi,1)
         f1d_ye(i)=y(0.d0,r,phi,1)

         write(*,*)'rho,f1d_xe(i),f1d_ye(i)',rho,f1d_xe(i),f1d_ye(i)                

         do j=1,n_z
c         do j=n_z/2,n_z/2
            z=y_axis(j)
            z=zmin+step_z*(j-1)
c            write(*,*)'j,z',j,z
            bmod=b(z,r,phi)
   
            if ((z.ge.zm).and.(z.le.zp)) then
c              the point is in the plasma
               sign_dp=(dpdrd*bz-dpdzd*br)
               if (sign_dp.ge.0.d0) then
                sign_dp=1.d0
               else
                sign_dp=-1.d0
               endif
               f2d_b(i,j)=real(bmod)
               f2d_bt(i,j)=real(bphi)
               f2d_psi(i,j)=real(psid)
               f2d_bp(i,j)=sqrt(real(bz**2+br**2))*sign_dp
c               write(*,*) 'in map_b before f2d_xe=real(x) '
               f2d_xe(i,j)=real(x(z,r,phi,1))
               f2d_ye(i,j)=real(y(z,r,phi,1))
            else
c              the point is outside the plasma
c               f2d_b(i,j)=0.0
c               f2d_bt(i,j)=0.0
c               f2d_bp(i,j)=0.0      
c               write(*,*)'mk_graph.f in map_b the point outside plasma'
               !YuP[07-2017] sign_dp was not defined. Not sure about this def:  
               sign_dp=(dpdrd*bz-dpdzd*br)
               if (sign_dp.ge.0.d0) then
                sign_dp=-1.d0 !YuP not sure
               else
                sign_dp=+1.d0 !YuP not sure
               endif
               f2d_b(i,j)=real(bmod)
               f2d_bt(i,j)=real(bphi)
               f2d_bp(i,j)=sqrt(real(bz**2+br**2))*sign_dp
               f2d_psi(i,j)=real(psid)
c               write(*,*) 'in map_b before f2d_xe=real(x) '
               f2d_xe(i,j)=real(x(z,r,phi,1))
               f2d_ye(i,j)=real(y(z,r,phi,1))
            endif
         enddo
      enddo

      write(*,*)'map_b abefore contour2d' 
cSm030519 (comment contour2d only for test) 
      call contour2d(n_r,n_z,f2d_b,x_axis,y_axis,
     .'B total [Tl]', 'R [m]','Z [m]',
     .n_contour,contour,
     .'name_param',param,0)

      call contour2d(n_r,n_z,f2d_bt,x_axis,y_axis,
     .'B toroidal [Tl]', 'R [m]','Z [m]',
     .n_contour,contour,
     .'name_param',param,0)

      call contour2d(n_r,n_z,f2d_bp,x_axis,y_axis,
     .'B poloidal [Tl]', 'R [m]','Z [m]',
     .n_contour,contour,
     .'name_param',param,0)

      call contour2d(n_r,n_z,f2d_psi,x_axis,y_axis,
     .'poloidal psi', 'R [m]','Z [m]',
     .n_contour,contour,
     .'name_param',param,0)

      call contour2d(n_r,n_z,f2d_xe,x_axis,y_axis,
     .'Xe', 'R [m]','Z [z]',
     .n_contour,contour,
     .'name_param',param,0)

      call contour2d(n_r,n_z,f2d_ye,x_axis,y_axis,
     .'Ye', 'R [m]','Z [m]',
     .n_contour,contour,
     .'name_param',param,0)


      do i=1,n_r
         dx_axis(i)=x_axis(i)
      enddo     
      ymin=0.d0
      ymax=0.d0
 
cSm030519 (comment plot1dt only for test)
      call plot1dt(dx_axis,f1d_bz,0,0,n_r,1,'linlin',0.d0,0.d0,
     .'R [m] at z=0','B_z [Tl]')

 10   format(5(' ',1pe10.3))
      open(7,file='bxy_r.dat')
      write(7,*)'r   bz   bt   ye   xe'
      do i=1,n_r
         write(7,10)x_axis(i),f1d_bz(i),f1d_bt(i),f1d_ye(i),f1d_xe(i)
      enddo
      close(7)

cSm030519 (comment plot1dt only for test)
      call plot1dt(dx_axis,f1d_bt,0,0,n_r,1,'linlin',0.d0,0.d0,
     .'R [m] at z=0','B_toroidal [Tl]')
              
      call plot1dt(dx_axis,f1d_xe,0,0,n_r,1,'linlin',0.d0,0.d0,
     .'R at z=0','Xe')
      
      call plot1dt(dx_axis,f1d_ye,0,0,n_r,1,'linlin',0.d0,0.d0,
     .'R [m] at z=0','Ye')

      step_rho=1.d0/dfloat(n_rho_a-1)
      do i=1,n_rho_a
         f1d_rho(i)=step_rho*(i-1)
         f1d_dens_e(i)=densrho(f1d_rho(i),1)
         f1d_temp_e(i)=temperho(f1d_rho(i),1)
      enddo
cSm030519 (comment plot1dt only for test)
      call plot1dt(f1d_rho,f1d_dens_e,0,0,n_rho_a,1,'linlin',0.d0,0.d0,
     .'small radius','electron density [10**13/cm**3]')


      call plot1dt(f1d_rho,f1d_temp_e,0,0,n_rho_a,1,'linlin',0.d0,0.d0,
     .'small radius','electron temperature [KeV]')

c--------------------------------------------------------------
c      create poloidal flux countours at full eqdsk box
c--------------------------------------------------------------
      step=xdimeqd/(n_r-1)
      do i=1,n_r
         x_axis(i)=redeqd+step*(i-1)
      enddo
      step=ydimeqd/(n_z-1)
      do i=1,n_z
         y_axis(i)=ymideqd+step*(i-1)-ydimeqd*0.5d0
      enddo

      do i=1,n_r
         r=x_axis(i)
         do j=1,n_z
            z=y_axis(j)
            f2d_psi(i,j)=fpsi(r,z)
         enddo
      enddo
      call contour2d(n_r,n_z,f2d_psi,x_axis,y_axis,
     .'poloidal flux', 'R [m]','Z [m]',
     .n_contour_wide,contour_wide,'name_param',param,1)

c      stop 'mk_graph.f in map_b'
      return
      end
     
      
      subroutine read_emfr_bin(iray)
c-----reads emfr.bin file        
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'write.i'
      REAL*4 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10

      OPEN(95,file='emfr.bin',form='unformatted',status='old')
      write(*,*)' in read_emfr_bin nfreq=',nfreq
      DO I=1,nfreq
      write(*,*)'I',I
      READ(95) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
c      write(*,*)' a1,a2,a3,a4,a5,a6,a7,a8,a9,a10',
c     +a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
c     +   wfreq(I),wtemp_rad_fr(I),
c     +  wi_0(iray,I),wtemp_rad_fr_wall(I),
c     +  wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
c     +  wrho0_em(I),wr_2nd_harm(I),
c     +  wtemp_2nd_harm(I)

c      write(*,*)'wfreq(I),wtemp_rad_fr(I),
c     +  wi_0(iray,I),wtemp_rad_fr_wall(I),
c     +  wtemp_pl_fr(I)),wr0_em(I),wz0_em(I),
c     +  wrho0_em(I)),wr_2nd_harm(I),
c     +  wtemp_2nd_harm(I)',
c     +  wfreq(I),wtemp_rad_fr(I),
c     +  wi_0(iray,I),wtemp_rad_fr_wall(I),
c     +  wtemp_pl_fr(I),wr0_em(I),wz0_em(I),
c     +  wrho0_em(I),wr_2nd_harm(I),
c     +  wtemp_2nd_harm(I)
c      READ(95)

      enddo  
      close(95)

      RETURN
      END
       

      subroutine map_d_cold(z,r,phi,npar,n_nperp,nperp_min,nperp_max,
     &name_param,param,n_param)
c-----creates----------------------------------------------------------
c     1D array D(nperp) for plotting of the cold plasma dispersion function
c     D(nperp) versus perpendicular refractive index N_perpendicular=nperp
c     at the given point (r,z,phi) and given N_parallel=npar
c
c     Plots D_cold(ReN_perp) to plot.ps file
c----------------------------------------------------------------------------
      implicit none
      include 'pgconst.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input:
      integer n_nperp                         !number of points in nperp mesh
      double precision nperp_min,nperp_max,   !minimal and maximal nperp values
     &z,r,phi, ! space point coordinates
     &npar     ! parallel refractive index
      integer n_param                    !number of input parameters
      character*(*) name_param(*)        !names  of the input parameters
      REAL*4 param(*)                      !values of the input parameters    
c-----locals
      integer n_nperp_a !maximal number of n_nperp
c      parameter(n_nperp_a=10000)
      parameter(n_nperp_a=100000)
      integer i

      double precision d_ar(n_nperp_a), !array of dispersion function values
     &nperp_ar(n_nperp_a)   !mesh of nperp

c      real d_ar(n_nperp_a), !array of dispersion function values
c     &nperp_ar(n_nperp_a)   !mesh of nperp

      double precision step,nperp
      double complex disp_func

      if(myrank.ne.0) return
      
      write(*,*)'in subroutine map_d_cold z,r,phi',z,r,phi
      write(*,*)'npar,n_nperp,nperp_min,nperp_max',
     &npar,n_nperp,nperp_min,nperp_max


      if(n_nperp.gt.n_nperp_a)then
        write(*,*)'in map_d_cold n_nperp>n_nperp_a'
        write(*,*)'n_nperp,n_nperp_a',n_nperp,n_nperp_a
        write(*,*)'Please increase n_nperp_a or decrease n_nperp'
        write(*,*)'and recompile the code'
        stop 'in map_d_cold'
      endif

      step=(nperp_max-nperp_min)/dfloat(n_nperp-1)

c      write(*,*)'mk_graph.f step',step
   
      do i=1,n_nperp
        nperp=nperp_min+(i-1)*step           
c        write(*,*)'i,nperp',i,nperp    
        call d_cold(z,r,phi,npar,nperp,disp_func)
        nperp_ar(i)=nperp
        d_ar(i)=dreal(disp_func)
c        write(*,*)'nperp_ar(i),d_ar(i)',nperp_ar(i),d_ar(i)
      enddo

c      call plot1dt(nperp_ar,d_ar,0,0,n_nperp,1,'linlin',0.d0,0.d0,
c     &'cold plasma dispersion function',
c     &'n_perpendicular','dispersion function')

      call plot1dt_param(nperp_ar,d_ar,0,0,n_nperp,1,'linlin',0.d0,0.d0,
     &'cold plasma dispersion function',
     &'n_perpendicular','cold dispersion function',n_param,name_param,
     &param)

      return
      end

      subroutine d_cold(z,r,phi,cnpar,cnper,disp_func)
c-----calculates cold plasma dispersion function disp_func
c     in space point (z,r,phi)
c      with refractive index components: cnpar,cnper
      implicit none
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'eps.i'
      double complex aK(3,3)
      double precision z,r,phi,cnpar,cnper
      double complex disp_func
c-----external
      double complex dcold_rlt
      integer i,j

      if(myrank.ne.0) return
      
      do i=1,3
        do j=1,3
cRobtAndre140412           aK(3,3)=dcmplx(0.d0,0.d0)
           aK(i,j)=dcmplx(0.d0,0.d0)
        enddo
      enddo 
      call tensrcld(z,r,phi) 
      disp_func=dcold_rlt(reps,aK,cnpar,cnper)
        
      return
      end

      subroutine wave_normal_surface(z,r,phi,cnpar_ray,cnper_ray,
     &t,n_gam,ib)
c-----calculates and creates plot of the wave normal surface 
c     for cold plasma dispersion
c     in z,r,phi point
c     
      implicit none 
      include 'pgconst.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input
      real*8 z,r,phi, ! space coordinates 
     &cnpar_ray,cnper_ray, ! parallel and perpendicular refractive index
                           ! components along the ray
     &t               ! poloidal distance of the ray point
      integer n_gam,  ! the number of angle points
     &ib ! 
c-----local
      integer n_gam_a !maximal value of n_gam
      parameter (n_gam_a=500)
      real*8 gam,     !the angle between B_field and k_vector
     & gam_ar(n_gam_a),! array of angles
     & step,           ! angle step      
     & nperp_1,nperp_2, ! roots of the dispersion relation
     & ns_1_ar(n_gam_a),! first root array
     & ns_2_ar(n_gam_a),! second root array
     & x_ar_1(n_gam_a),z_ar_1(n_gam_a),!projections of the 
     & x_ar_2(n_gam_a),z_ar_2(n_gam_a),!the wave normal vector
     & x_ar_3(n_gam_a,2),z_ar_3(n_gam_a,2),!the wave normal vector
     & pi,
     & ad,bd,cd,         ! dispersion relation coefficients
     & d4,d2,d0,det,cn1,ds,dc,ds2,dc2,ds4,
     & xe,xi,ye,yb,dele,delib,dn

      integer n_param                    !number of parameters at figure
      character*(20) name_param(4)       !names  of parameters at figure
      REAL*4 param(4)                      !values of parameters at figure 

c-----for Drawing Markers:
      integer n_marker        ! number of points to be marked,
      parameter  (n_marker=1)   
      REAL*4  x_marker(n_marker),y_marker(n_marker)
      integer n_symbol_marker
      real*8 gam_ray
c-----externals
      real*8 x,y

      integer i,i_sign   
      integer ibmx !local 
     
      if(myrank.ne.0) return
      
      i_sign=1

      if (n_gam.gt.n_gam_a) then
         write(*,*)'in wave_normal_surface'
         write(*,*)'n_gam > n_gam_a'
         write(*,*)'n_gam=',n_gam,'n_gam_a=',n_gam_a
         write(*,*)'Please change parameter n_gam_a'
         write(*,*)' in wave_normal_surface and recompile the code'
         stop 'in wave_normal_surface'
      endif


      do i=1,n_gam
         ns_1_ar(i)=0.d0
         ns_2_ar(i)=0.d0
         x_ar_1(i)=0.d0
         x_ar_2(i)=0.d0
         z_ar_1(i)=0.d0
         z_ar_2(i)=0.d0
      enddo     

      pi=4.d0*datan(1.d0)

      step=2.d0*pi/dfloat(n_gam-1)
      
      do i=1,n_gam
        gam=step*(i-1)
        gam_ar(i)=gam
        ds=dsin(gam)
        dc=dcos(gam) 
        ds2=ds*ds
        dc2=dc*dc
        ds4=ds2*ds2
       
        call abc(z,r,phi,ds2,dc2,ad,bd,cd) 
c-------dispersion relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
        d4=ad
        d2=bd
        d0=cd
        det=d2*d2-4.d0*d4*d0
        cn1=dsqrt(d2*d2-4.d0*d4*d0)

        if (ib.eq.1) then
c----------ib.eq.1 electron resonance condition may be
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
           write(*,*)'wave_normal_surface xe,ye',xe,ye	 
           dele=1.d0-ye
           if (dele.lt.0.d0) i_sign=-1 ! not used
c---------------------------------------------------------------------
c          calculate roots of the dispersion
c          relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
c           write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
c           write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
           ns_1_ar(i)=(-d2+cn1)/(2.d0*d4)
           ns_2_ar(i)=(-d2-cn1)/(2.d0*d4)         
	  
        endif ! ib.eq.1

        if (ib.gt.1) then
c----------ib.gt.1 iones resonance condition may be
           ibmx=ib !min(ib,nbulk) ! safety check: not to exceed nbulk
           yb=y(z,r,phi,ibmx)
           delib=1.d0-yb     
           if (delib.lt.0.d0) i_sign=-1 ! not used
c--------------------------------------------------------------------
c          calculate roots of the dispersion
c          relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c--------------------------------------------------------------------
           cn1=dsqrt(d2*d2-4.d0*d4*d0)
c           write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
c           write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)       
	   ns_1_ar(i)=(-d2+cn1)/(2.d0*d4)
           ns_2_ar(i)=(-d2-cn1)/(2.d0*d4)     
        endif ! ib.gt.1

c        write(*,*)'i, ns_1_ar(i)',i, ns_1_ar(i)
        if (ns_1_ar(i).gt.0.d0) then
           dn=1.d0/dsqrt(ns_1_ar(i))
           z_ar_1 (i)=dc*dn
           x_ar_1 (i)=ds*dn
        endif 
c        write(*,*)'i,x_ar_1(i),z_ar_1(i)',i,x_ar_1(i),z_ar_1(i)

c        write(*,*)'i, ns_2_ar(i)',i, ns_2_ar(i)
        if(ns_2_ar(i).gt.0.d0) then
           dn=1.d0/dsqrt(ns_2_ar(i))
           z_ar_2 (i)=dc*dn
           x_ar_2 (i)=ds*dn
        endif
c        write(*,*)'i,x_ar_2(i),z_ar_2(i)',i,x_ar_2(i),z_ar_2(i)

      enddo !i   
c-----test hamilt=0?
c     hamtest=d4*cn2*cn2+d2*cn2+d0
  
c      do i=1,n_gam
c        write(*,*)'i,gam_ar(i),ns_1_ar(i)',i,gam_ar(i),ns_1_ar(i)
c      enddo
c      do i=1,n_gam
c        write(*,*)'i,gam_ar(i),ns_2_ar(i)',i,gam_ar(i),ns_2_ar(i)
c      enddo

c      write(*,*)'mk_graph.f in subroutine wave_normal_surface'

c     do i=1,n_gam
c       write(*,*)'i,gam_ar(i),x_ar_1(i),z_ar_1(i)',
c     &i,gam_ar(i),x_ar_1(i),z_ar_1(i)
c      enddo

c      do i=1,n_gam
c       write(*,*)'i,gam_ar(i),x_ar_2(i),z_ar_2(i)',
c     &i,gam_ar(i),x_ar_2(i),z_ar_2(i)
c      enddo

      n_param=4
      name_param(1)='r [m]'
      name_param(2)='pol. dist [m]'
c----------------------------------------------------------------------------
c     n**2 at the dipersion curve VS gam (theta=gam)
c----------------------------------------------------------------------------
      name_param(3)='n**2 at ray'
      name_param(4)='theta  at ray'
      param(1)=r
      param(2)=t
      param(3)=(cnpar_ray**2+cnper_ray**2)
      gam_ray=dacos(cnpar_ray/dsqrt(cnpar_ray**2+cnper_ray**2))
      param(4)=gam_ray

c-----with marker

      n_symbol_marker=9
      gam_ray=dacos(cnpar_ray/dsqrt(cnpar_ray**2+cnper_ray**2))
      x_marker(1)=gam_ray
      y_marker(1)=param(3)
        
      call plot1dt_marker_param(gam_ar,ns_1_ar,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'n**2 1=(-b+sqrt(b**2-4*a*c))/(2*a)',
     &'theta','n**2 1',
     &n_param,name_param,param)

      call plot1dt_marker_param(gam_ar,ns_2_ar,0,0,n_gam,2,
     & n_marker,x_marker,y_marker,n_symbol_marker,    
     &'linlin',0.d0,0.d0,
     &'n**2 2=(-b-sqrt(b**2-4*a*c))/(2*a) ',
     &'theta','n**2 2',
     &n_param,name_param,param)

c--------------------------------------------------------------------
c     wave normals
c--------------------------------------------------------------------
      name_param(3)='n_perp/n**2 at ray'
      name_param(4)='n_par/n**2  at ray'
      param(3)=cnper_ray/(cnpar_ray**2+cnper_ray**2)
      param(4)=cnpar_ray/(cnpar_ray**2+cnper_ray**2)
                              
c-----with marker

      n_symbol_marker=9
      x_marker(1)=param(3)
      y_marker(1)=param(4)

      if(dabs(x_ar_1(n_gam/2)-x_ar_1(1)).gt.1.d-13) then
        call plot1dt_marker_param(x_ar_1,z_ar_1,0,0,n_gam,2,
     &  n_marker,x_marker,y_marker,n_symbol_marker,    
     &  'linlin',0.d0,0.d0,
     &'wave normal surface cold n**2 1 root','n_perp/n**2','n_par/n**2'
     &  ,n_param,name_param,param)
      endif

      if(dabs(x_ar_2(n_gam/2)-x_ar_2(1)).gt.1.d-13) then
        call plot1dt_marker_param(x_ar_2,z_ar_2,0,0,n_gam,2,
     &  n_marker,x_marker,y_marker,n_symbol_marker,    
     &  'linlin',0.d0,0.d0,
     &'wave normal surface cold n**2 2 root','n_perp/n**2','n_par/n**2'
     &  ,n_param,name_param,param)
      endif

c      stop 'subroutine wave_normal_surface'
      return
      end

      subroutine wave_ray_normal_surface(z,r,phi,cnpar_ray,cnper_ray,
     &t,n_gam,ib)
c-----calculates and creates plot of the wave normal surface 
c     and the reciprocal ray refractive index (1/N_r) surface
c     for cold plasma dispersion
c     in z,r,phi point
c     
      implicit none 
      include 'pgconst.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c-----input
      real*8 z,r,phi, ! space coordinates 
     &cnpar_ray,cnper_ray, ! parallel and perpendicular refractive index
                           ! components along the ray
     &t               ! poloidal distance of the ray point
      integer n_gam,  ! the number of angle points
     &ib ! 
c-----local
      integer n_gam_a !maximal value of n_gam
      parameter (n_gam_a=500)
      real*8 gam,     !the angle between B_field and k_vector
     & gam_ar(n_gam_a),! array of angles
     & step,           ! angle step      
     & nperp_1,nperp_2, ! roots of the dispersion relation
     & ns_1_ar(n_gam_a),! first root array
     & ns_2_ar(n_gam_a),! second root array
     & x_ar_1(n_gam_a),z_ar_1(n_gam_a),!projections of the 
     & x_ar_2(n_gam_a),z_ar_2(n_gam_a),!the wave normal vector
     & x_ar_3(n_gam_a,2),z_ar_3(n_gam_a,2),!the wave normal vector
     & x_Nray_ar_1(n_gam_a),z_Nray_ar_1(n_gam_a),!projections of the
                                                 !the wave ray normal vector
     & x_Nray_ar_2(n_gam_a),z_Nray_ar_2(n_gam_a),!projections of the
                                                 !the wave ray normal vector
     & pi,
     & ad,bd,cd,         ! dispersion relation coefficients
     & d4,d2,d0,det,cn1,ds,dc,ds2,dc2,ds4,
     & xe,xi,ye,yb,dele,delib,dn,
     & omegpedce,omegdce,cnpar_loc,cnper_loc,
     & N_ray, !ray refractive index
     & cnpar_Nray,cnper_Nray,gam_Nray,
     & p     

      integer n_param                    !number of parameters at figure
      character*(22) name_param(4)       !names  of parameters at figure
      REAL*4 param(4)                      !values of parameters at figure 

c-----for Drawing Markers:
      integer n_marker        ! number of points to be marked,
      parameter  (n_marker=1)   
      REAL*4  x_marker(n_marker),y_marker(n_marker)
      integer n_symbol_marker
      real*8 gam_ray
c-----externals
      real*8 x,y,rrind    
 
      integer i,i_sign    
      integer ibmx !local

      if(myrank.ne.0) return
      
      i_sign=1
      
      if (n_gam.gt.n_gam_a) then
         write(*,*)'in wave_normal_surface'
         write(*,*)'n_gam > n_gam_a'
         write(*,*)'n_gam=',n_gam,'n_gam_a=',n_gam_a
         write(*,*)'Please change parameter n_gam_a'
         write(*,*)' in wave_normal_surface and recompile the code'
         stop 'in wave_normal_surface'
      endif

      do i=1,n_gam
         ns_1_ar(i)=0.d0
         ns_2_ar(i)=0.d0
         x_ar_1(i)=0.d0
         x_ar_2(i)=0.d0
         z_ar_1(i)=0.d0
         z_ar_2(i)=0.d0
         x_Nray_ar_1(i)=0.d0
         x_Nray_ar_2(i)=0.d0
         z_Nray_ar_1(i)=0.d0
         z_Nray_ar_2(i)=0.d0
      enddo     

      xe=x(z,r,phi,1)
      ye=y(z,r,phi,1)
      omegpedce=dsqrt(xe)/dabs(ye)      !omega_pe/omega_ce
      omegdce=1.d0/dabs(ye)             !omega/omega_ce

      pi=4.d0*datan(1.d0)

      step=2.d0*pi/dfloat(n_gam-1)

      do i=1,n_gam
        gam=step*(i-1)
        gam_ar(i)=gam
        ds=dsin(gam)
        dc=dcos(gam) 
        ds2=ds*ds
        dc2=dc*dc
        ds4=ds2*ds2

        call abc(z,r,phi,ds2,dc2,ad,bd,cd) 

c-------dispersion relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
        d4=ad
        d2=bd
        d0=cd
        det=d2*d2-4.d0*d4*d0
        cn1=dsqrt(d2*d2-4.d0*d4*d0)

        if (ib.eq.1) then
c----------ib.eq.1 electron resonance condition may be
           xe=x(z,r,phi,1)
           ye=y(z,r,phi,1)
c           write(*,*)'wave_normal_surface xe,ye',xe,ye	 
           dele=1.d0-ye
           if (dele.lt.0.d0) i_sign=-1 ! not used
c---------------------------------------------------------------------
c          calculate roots of the dispersion
c          relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c---------------------------------------------------------------------
c           write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
c           write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)
           ns_1_ar(i)=(-d2+cn1)/(2.d0*d4)
           ns_2_ar(i)=(-d2-cn1)/(2.d0*d4)         
	  
        endif ! ib.eq.1

        if (ib.gt.1) then
c----------ib.gt.1 iones resonance condition may be
           ibmx=ib !min(ib,nbulk) ! safety check: not to exceed nbulk
           yb=y(z,r,phi,ibmx)
           delib=1.d0-yb     
           if (delib.lt.0.d0) i_sign=-1 ! not used
c--------------------------------------------------------------------
c          calculate roots of the dispersion
c          relation  n**2=(-b+ioxm*sqrt(b**2-4*a*c))/(2*a)
c--------------------------------------------------------------------
           cn1=dsqrt(d2*d2-4.d0*d4*d0)
c           write(*,*)'(-d2+cn1)/(2.d0*d4)',(-d2+cn1)/(2.d0*d4)
c           write(*,*)'(-d2-cn1)/(2.d0*d4)',(-d2-cn1)/(2.d0*d4)       
	   ns_1_ar(i)=(-d2+cn1)/(2.d0*d4)
           ns_2_ar(i)=(-d2-cn1)/(2.d0*d4)     
        endif ! ib.gt.1

c        write(*,*)'i, ns_1_ar(i)',i, ns_1_ar(i)
        if (ns_1_ar(i).gt.0.d0) then
           dn=1.d0/dsqrt(ns_1_ar(i))
           z_ar_1 (i)=dc*dn
           x_ar_1 (i)=ds*dn
c------------------------------------------------------------------------
c          calculate ray refractive index: cnray
c--------------------------------------------------------------
           cnpar_loc=dc/dn
           cnper_loc=ds/dn

           N_ray=rrind(cnpar_loc,cnper_loc,omegdce,omegpedce)

           if (i.gt.1) then
            z_Nray_ar_1 (i-1)=dc/N_ray
            x_Nray_ar_1 (i-1)=ds/N_ray
           endif
        endif 

        if(ns_2_ar(i).gt.0.d0) then
           dn=1.d0/dsqrt(ns_2_ar(i))
           z_ar_2 (i)=dc*dn
           x_ar_2 (i)=ds*dn
c------------------------------------------------------------------------
c          calculate ray refractive index: cnray
c--------------------------------------------------------------
           cnpar_loc=dc/dn
           cnper_loc=ds/dn
           
           N_ray=rrind(cnpar_loc,cnper_loc,omegdce,omegpedce)
c           write(*,*)'N_ray',N_ray
           if (i.gt.1) then
             z_Nray_ar_2 (i-1)=dc/N_ray
             x_Nray_ar_2 (i-1)=ds/N_ray
           endif
        endif
c        write(*,*)'i,x_ar_2(i),z_ar_2(i)',i,x_ar_2(i),z_ar_2(i)
       
      enddo !i=1,n_gam

c      stop 'wave_ray_normal_surface'

      n_param=4
      name_param(1)='r [m]'
      name_param(2)='pol. dist [m]'
c----------------------------------------------------------------------------
c     n**2 is at the dipersion curve VS (theta=gam)
c----------------------------------------------------------------------------
      param(1)=r
      param(2)=t

c--------------------------------------------------------------------
c     wave ray normal surface
c--------------------------------------------------------------------
c-----with marker
      name_param(3)='sin(theta)/Nray at ray'
      name_param(4)='cos(theta)/Nray at ray'
      N_ray=rrind(cnpar_ray,cnper_ray,omegdce,omegpedce)
      param(3)=cnper_ray/dsqrt(cnpar_ray**2+cnper_ray**2)/N_ray
      param(4)=cnpar_ray/dsqrt(cnpar_ray**2+cnper_ray**2)/N_ray
      
      n_symbol_marker=9

      gam_ray=dacos(cnpar_ray/dsqrt(cnpar_ray**2+cnper_ray**2))     
      x_marker(1)=param(3)
      y_marker(1)=param(4)

      write(*,*)'x_marker(1),y_marker(1) ',x_marker(1),y_marker(1)

c      do i=1,n_gam-1
c        write(*,*)'i,x_Nray_ar_1(i),z_Nray_ar_1(i)',   
c     &             i,x_Nray_ar_1(i),z_Nray_ar_1(i)
c      enddo

      if(dabs(x_Nray_ar_1(n_gam/2)-x_Nray_ar_1(1)).gt.1.d-13) then
        call plot1dt_marker_param(x_Nray_ar_1,z_Nray_ar_1,0,0,n_gam-1,
     &  2,n_marker,x_marker,y_marker,n_symbol_marker,    
     &  'linlin',0.d0,0.d0,
c     &  'wave normal surface 1/N_ray cold n**2 1 root',
c     &  'recip ray ref surface 1/N_r, cold n**2 1 root',
     &  '1/(ray refr index, N_r), cold n**2 1st root',
     &  'sin(theta)/N_ray',
     &  'cos(theta)/N_ray',
     &  n_param,name_param,param)
      endif 

c      do i=1,n_gam-1
c        write(*,*)'i,x_Nray_ar_2(i),z_Nray_ar_1(i)',   
c     &             i,x_Nray_ar_2(i),z_Nray_ar_1(i)
c      enddo
c      write(*,*)'n_gam/2=',n_gam/2
c      write(*,*)'x_Nray_ar_2(n_gam/2),x_Nray_ar_2(1)',
c     &x_Nray_ar_2(n_gam/2),x_Nray_ar_2(1)

c      write(*,*)'x_Nray_ar_2(n_gam/2)-x_Nray_ar_2(1)',
c     &x_Nray_ar_2(n_gam/2)-x_Nray_ar_2(1)
    
      if(dabs(x_Nray_ar_2(n_gam/2)-x_Nray_ar_2(1)).gt.1.d-13) then
c        write(*,*)'mk_graph.f before  plot1dt_marker_param N 2'
        call plot1dt_marker_param(x_Nray_ar_2,z_Nray_ar_2,0,0,n_gam-1,
     &  2,n_marker,x_marker,y_marker,n_symbol_marker,    
     & 'linlin',0.d0,0.d0,
c     & 'wave normal surface 1/N_ray cold n**2 2 root',
c     &  'recip ray ref surface 1/N_r, cold n**2 2 root',
     &  '1/(ray refr index, N_r), cold n**2 2nd root',
     & 'sin(theta)/N_ray',
     & 'cos(theta)/N_ray',
     &  n_param,name_param,param)
c        write(*,*)'mk_graph.f after plot1dt_marker_param N 2'
      endif 
 
      return
      end


      subroutine map_d_hot(z,r,phi,npar,n_nperp,nperp_min,nperp_max,
     &name_param,param,n_param)
c-----creates----------------------------------------------------------
c     1D array D(nperp) for plotting hot plasma dispersion function
c     Re(D(nperp)) on perpendicular refractive index N_perpendicular=nperp
c     in given the point (r,z,phi) and given N_parallel=npar
c
c     Plots ReD_hot(ReN_perp) to plot.ps file
c----------------------------------------------------------------------------
      implicit none
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'ions.i'
cfor test
      include 'eps.i'
c-----input:
      integer n_nperp                         !number of points in nperp mesh
      double precision nperp_min,nperp_max,   !minimal and maximal nperp values
     &z,r,phi, ! space point coordinates
     &npar     ! parallel refractive index
      integer n_param                         !number of input parameters
      character*(*) name_param(n_param)        !names  of the input parameters
      REAL*4 param(n_param)                      !values of the input parameters    
c-----locals
      integer n_nperp_a !maximal number of n_nperp
c      parameter(n_nperp_a=10000)
      parameter(n_nperp_a=1000000)
      integer i

      double precision d_ar(n_nperp_a), !array of dispersion function values
     &nperp_ar(n_nperp_a)   !mesh of nperp

c      real d_ar(n_nperp_a), !array of dispersion function values
c     &nperp_ar(n_nperp_a)   !mesh of nperp

      double precision step,nperp
      double complex disp_func

      double precision x_ar(nbulka),y_ar(nbulka),
     &t_av_ar(nbulka),tpop_ar(nbulka),vflow_ar(nbulka),te

      double complex K_sum(3,3),d_complex

      integer iherm_loc,k1,k2


c-----externals
      double precision b,x,y,tempe,tpoprho,vflowrho
      double complex dhot_sum

      if(myrank.ne.0) return
      
      write(*,*)'in subroutine map_d_hot z,r,phi',z,r,phi
      write(*,*)'npar,n_nperp,nperp_min,nperp_max',
     &npar,n_nperp,nperp_min,nperp_max


      if(n_nperp.gt.n_nperp_a)then
        write(*,*)'in map_d_hot n_nperp>n_nperp_a'
        write(*,*)'n_nperp,n_nperp_a',n_nperp,n_nperp_a
        write(*,*)'Please increase n_nperp_a or decrease n_nperp'
        write(*,*)'and recompile the code'
        stop 'in map_d_hot'
      endif

      bmod=b(z,r,phi)

      do i=1,nbulk
        x_ar(i)=x(z,r,phi,i)
	y_ar(i)=y(z,r,phi,i)
        if(i.eq.1) y_ar(1)=-y_ar(1)
	te=tempe(z,r,phi,i) ! kev
	t_av_ar(i)=te*1000.d0      ! ev 
        tpop_ar(i)=tpoprho(rho,i)
        vflow_ar(i)=vflowrho(rho,i)
      enddo

c      write(*,*)'dmas',dmas
ctest      
      call tensrcld(z,r,phi) !tensor reps is in one.i  
cendtest
      step=(nperp_max-nperp_min)/dfloat(n_nperp-1)

      write(*,*)'n_nperp, step',n_nperp,step
      iherm_loc=1

      write(*,*)' nbulk,npar',nbulk,npar
      do i=1,nbulk
        write(*,*)'=i,dmas(i),x_ar(i),y_ar(i)',i,dmas(i),x_ar(i),y_ar(i)
        write(*,*)'t_av_ar(i),tpop_ar(i),vflow_ar(i)',
     &             t_av_ar(i),tpop_ar(i),vflow_ar(i)
      enddo
      do i=1,n_nperp
        nperp=nperp_min+(i-1)*step           
c        call d_cold(z,r,phi,npar,nperp,disp_func)
        nperp_ar(i)=nperp
        d_complex=dhot_sum(nbulk,dmas,x_ar,y_ar,t_av_ar,tpop_ar,
     &   vflow_ar,npar,nperp,iherm_loc,K_sum)
        d_ar(i)=dreal(d_complex)

        write(*,*)'i,nperp_ar(i),d_ar(i)',i,nperp_ar(i),d_ar(i)

ctest
c        do k1=1,3
c           do k2=1,3
c              write(*,*)'k1,k2,K_sum(k1,k2),reps(k1,k2)',
c     &                   k1,k2,K_sum(k1,k2),reps(k1,k2)
c           enddo
c        enddo
cendtest

      enddo

c      write(*,*)'in map_d_hot'
c      write(*,*)'n_param',n_param

c      do i=1,n_param
c         write(*,*)'i, name_param(i) ',i,name_param(i)
c         write(*,*)'param(i) ',param(i)
c      enddo

c      write(*,*)'map_d_hot before  plot1dt_param'

c      call plot1dt(nperp_ar,d_ar,0,0,n_nperp,1,'linlin',0.d0,0.d0,
c     &'n_perpendicular','Re(hot dispersion function)')

      call plot1dt_param(nperp_ar,d_ar,0,0,n_nperp,1,'linlin',0.d0,
     &0.d0,'hot plasma dispersion function','n_perpendicular',
     &'Re(Dhot)',n_param,name_param,param)

c      write(*,*)'map_d_hot after plot1dt_param'

      return

      end

     
      double precision function freq_p(z,r,phi,i)
c--------------------------------------------------------
c     calculates plasma frequency [GHZ]
c-----------------------------------------------------
      implicit none
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'ions.i'
c--------------------------------------------
c     input
      double precision z,r,phi
      integer i !number of plasma component
c-----externals
      double precision dense

c-----locals
      double precision den,mass_e,charge_electron

      den=dense(z,r,phi,i) !10**13 cm**-3
      if(den.lt.0.d0)then
         den=0.d0
      endif

      mass_e=9.1094d-28         !electron rest mass (g)
      charge_electron=4.8032d-10 !electron charge (statcoulomb)  

c     greq_pe=sqrt(4*pi*n_e*e**2/m_e)/2pi      HZ
c     freq_pi=sqrt(4*pi*n_i*Z**2*e**2/m_i)/2pi      

      freq_p=dsqrt(4.d0*pi*den*1.d13*(charge_electron*charge(i))**2/
     &             (mass_e*dmas(i)))/(2.d0*pi)*1.d-9          !GHZ
     
    
      return
      end

       
      double precision function freq_c(z,r,phi,i)
c--------------------------------------------------------
c     calculates gyrofrequency [GHZ]
c-----------------------------------------------------
      implicit none
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
      include 'ions.i'
c--------------------------------------------
c     input
      double precision z,r,phi
      integer i !number of plasma component 

c-----locals
      double precision den,mass_e,charge_electron,clight

      mass_e=9.1094d-28          !electron rest mass (g)
      charge_electron=4.8032d-10 !electron charge (statcoulomb) 
      clight=2.99792458d10       !light speed (cm/sec)

c     freq_ce=(e*B/m_e*clight)/2pi      HZ
c     freq_ci=(Z*e*B/m_i*clihgt)/2pi      

      freq_c=charge_electron*charge(i)*bmod*1.d4/
     &       (mass_e*dmas(i)*clight)/(2.d0*pi)*1.d-9   

      return
      end
 
    

      subroutine plot_fcefuh(z_freq,r_freq,alpha_freq,
     &beta_freq,dist_freq,
     &nsteps_freq,n_ec_harmonics_freq,npar_freq,
     &max_plot_freq)
 
C     straightpath_plot_fcefuh to plot f_ce*j, f_uh, f_pe
C     with a path through the plasma (midplane for example)  

      !implicit none
      implicit double precision (a-h,o-z)
      include 'pgconst.i'

      include 'param.i' ! specifies code parameters 
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
c      include 'commons.i'
      include 'one.i'
      include 'three.i'
      include 'ions.i'
c-----input
      double precision
     & r_freq,z_freq, !straight line edge point [m]
     & dist_freq,     !straight line length [m]
     & alpha_freq,    !toroidal angle [degree] of straight line
     & beta_freq,     !angle between Z direction and straigth
                      !line direction [degree]
     &npar_freq,      !N_parallel to plot X mode cutoff
     &max_plot_freq   ! maximal frquency at plot GHZ

      integer nsteps_nfreq, !the number of points along the line
     &n_ec_harmonics_nfreq  !number of EC harmonic
c-----external
      double precision b,freq_c,freq_p
     
      integer maxnstep
      parameter (maxnstep=1000) !max value of nsteps_freq

      double precision raxis,zaxis,
     &fpi,fci,flh,fx1,fx2
      REAL*4 ss(maxnstep),fuhs(maxnstep)
      REAL*4 fces(maxnstep),rs(maxnstep)
      REAL*4 fpes(maxnstep),rhos(maxnstep),den_e_s(maxnstep)
      REAL*4 xmin,xmax,ymin,ymax,newx(maxnstep),newy(maxnstep)
     
      integer xtitlelen,titlelen,PGOPEN,nchoice,j,k
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if(myrank.ne.0) return

      if (nsteps_freq.gt.maxnstep) then
         write(*,*)'mk_graph.f in plot_fcefuh'
         write(*,*)'nsteps_freq.gt.maxnstep'
         write(*,*)'it should be nsteps_freq.le.maxnstep'
         write(*,*)'nsteps_freq,maxnstep',nsteps_freq,maxnstep
         write(*,*)'Please reduse nsteps_freq in gebray.dat'
         stop 'in plot_fcefuh'
      endif

c-----open PGplot
c     call plotinit

C-------- find magnetic axis (psi=0, theta=0)
      call zr_psith(psimag,0.0,zaxis,raxis)
C**** FOR ARIESST: Magnetic axis at: R= 4.678357129999998  Z= -4.316799816857955E-18

c-----plot magnetic field contours to plot.ps file

C****** convert alpha, beta into radians from degrees
      alpha=alpha_freq*PI/180.0
      beta=beta_freq*PI/180.0
c      print *,'PI= ',PI
C****** starting coordinates in x,y,z
      xc=r_freq
      yc=0
      zc=z_freq
      sc=0
      ds=dist_freq/nsteps_freq
      dx=ds*sin(beta)*cos(alpha)
      dy=ds*sin(beta)*sin(alpha)
      dz=ds*cos(beta)
    
      do j=1,nsteps_freq 
C------- keep track of coordinates x,y,z, then convert to r,z,phi for
C------- calling bmod, dense, tempe
         r=sqrt(xc**2+yc**2)
         z=zc
         phi=atan2(yc,xc)
         bmod=b(z,r,phi)         
c         print *,'---r=',r,' z=',z,' phi=',phi,' rho=',rho,'--'
C----------- now calculate cold plasma stuff
         fces(j)=freq_c(z,r,phi,1)
         fpes(j)=freq_p(z,r,phi,1)
         fuhs(j)=sqrt(fces(j)**2+fpes(j)**2)
         ss(j)=sc
         rs(j)=r
         rhos(j)=rho
         den_e_s(j)=densrho(rho,1)
c         print *,'sc=',sc,'fces(j)=',fces(j),'fpes(j)=',fpes(j)
C**** step forward
         sc=sc+ds
         xc=xc+dx
         yc=yc+dy
         zc=zc+dz
      enddo

c      print *,'psimag=',psimag
c      print *,'Magnetic axis at: R=',raxis,' Z=',zaxis
    
c------------------------------------------------------------------
c     calulate minimal and maximal frequencies: ymin,ymax [GHZ]
c------------------------------------------------------------------
      ymin=fces(1) 
      ymax=fces(1)

      do j=1,nsteps_freq 
         if (fpes(j).lt.ymin) ymin=fpes(j)
         if (fces(j).lt.ymin) ymin=fces(j)
         if (fuhs(j).lt.ymin) ymin=fuhs(j)

         if (fpes(j).gt.ymax) ymax=fpes(j)
         if (fces(j)*n_ec_harmonics_freq.gt.ymax) then
             ymax=fces(j)*n_ec_harmonics_freq
         endif
         if (fuhs(j).gt.ymax) ymax=fuhs(j)
      enddo



C**** find limits for plot of real parts
C**** was ss before
        xmin=rs(nsteps_freq)
        xmax=rs(1)
C******* (reverse)
c        print *,'xmin=',xmin,' xmax=',xmax
C**** machine axis is at reqd (not magnetic axis with shafranov shift)
c        print *,'reqd: ',reqd
C***** Plot freq0, and dimensionless scale, and eta
        ymin=0.0      ! GHz for y axis
        ymax=max_plot_freq

        CALL PGSLW(4)
        CALL PGSCH(R41P5)
        CALL PGENV(xmin,xmax,ymin,ymax,0,1)
c        print *,'Past PGENV '
C*** Sets limit for plot.
        xtitle='R along midplane (m)'
        xtitlelen=LEN_TRIM(xtitle)
        titlestr=''
        titlelen=LEN_TRIM(titlestr)
c        print *,'xtitlelen: ',xtitlelen,' titlelen: ',titlelen
        CALL PGLAB(xtitle(1:xtitlelen),'Frequency (GHz)'
     +  ,titlestr(1:titlelen))
        
        print *,'fces: ',fces(1),fces(nsteps_freq)
        print *,'fuhs: ',fuhs(1),fuhs(nsteps_freq)
        print *,'fpes: ',fpes(1),fpes(nsteps_freq)
        print *,'rs: ',rs(1),rs(nsteps_freq)

C*** note on screen, BLACK=0 means background color
C*** WHITE=1 means default foreground color
        CALL PGSLS(1)
        CALL PGSCI(1)
        CALL PGSLW(6) 
        CALL PGLINE(nsteps_freq,rs,fces)
C*** plot fundamental and up to 5th harmonic
        DO j=2,n_ec_harmonics_freq
           DO k=1,nsteps_freq
              newy(k)=fces(k)*j
           ENDDO
           CALL PGLINE(nsteps_freq,rs,newy)
           print *,'newy: ',newy(1),newy(nsteps_freq)
        ENDDO
        print *,'WHITE: fce and harmonics'

        CALL PGSLS(3)
        CALL PGSCI(RED)
        CALL PGLINE(nsteps_freq,rs,fpes)
        print *,'RED: fpes'

        CALL PGSLS(2)
        CALL PGSCI(BLUE)
        CALL PGLINE(nsteps_freq,rs,fuhs)
        print *,'BLUE: fuhs'

        CALL PGSLS(4)
        CALL PGSCI(GREEN)
        newx(1)=raxis
        newx(2)=raxis
        newy(1)=ymin
        newy(2)=ymax
        CALL PGLINE(2,newx,newy)

C***** magnetic axis
C** Now add labels   :BEED TO ADJUST FOR .not.ARIES (BH)
        CALL PGSCI(WHITE)
        CALL PGTEXT(R44,R431,'f\dce\u')
        CALL PGTEXT(R43P8,R470,'2 f\dce\u')
        CALL PGTEXT(R43P8,R4110,'3 f\dce\u')
        CALL PGSCI(RED)
        CALL PGTEXT(R41P4,R445,'f\dpe\u')
        CALL PGSCI(BLUE)
        CALL PGTEXT(R42,R4165,'f\duh\u')
        CALL PGSCI(GREEN)
        CALL PGTEXT(R44P75,R410,'R\dm\u')
c        call PGEND

c--------------------------------------------------------
c     write frequencies to frequency.bin file
c-----------------------------------------------
c$$$      OPEN(10,file='frequency.bin',form='unformatted') 
c$$$
c$$$      do j=1,nsteps_freq 
c$$$        write(10) real(rs(j)),real(fpes(j)),real(fces(j)),
c$$$     &   real(fuhs(j)),real(2*fces(j)),
c$$$     &   real(3*fces(j)),real(4*fces(j)),real(5*fces(j)),
c$$$     &   real(6*fces(j))
c$$$      enddo
c$$$
c$$$      CLOSE(10)


      OPEN(10,file='freqelec.bin',form='unformatted') 

      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(fpes(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(fuhs(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(2*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
       do j=1,nsteps_freq 
        write(10) real(rs(j)),real(3*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(4*fces(j)),real(rhos(j)),  
     &            real(den_e_s(j))
      enddo   
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(5*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo
      WRITE(10)
      do j=1,nsteps_freq 
        write(10) real(rs(j)),real(6*fces(j)),real(rhos(j)),
     &            real(den_e_s(j))
      enddo    
      WRITE(10)
      do j=1,nsteps_freq
        fx1=0.5d0*(-fces(j)+dsqrt(fces(j)**2+4.d0*fpes(j)**2/
     &                           (1.d0-npar_freq**2))) 
        write(10) real(rs(j)),real(fx1),real(rhos(j)),
     &            real(den_e_s(j))
c        write(*,*)'xe-(1+ye)*(1-npar_freq**2)'
c      write(*,*)'j,rs(j),fx1',j,rs(j),fx1
c      write(*,*)'fpes(j)/fx1)**2-(1.d0+fces(j)/fx1)*(1.d0-npar_freq**2)'
c      write(*,*)(fpes(j)/fx1)**2-(1.d0+fces(j)/fx1)*(1.d0-npar_freq**2)
      enddo 
      WRITE(10)
      do j=1,nsteps_freq
        fx2=0.5d0*(fces(j)+dsqrt(fces(j)**2+4.d0*fpes(j)**2/
     &                           (1.d0-npar_freq**2))) 
        write(10) real(rs(j)),real(fx2),real(rhos(j)),
     &            real(den_e_s(j))
c      write(*,*)'j,rs(j),fx2',j,rs(j),fx2
c      write(*,*)'fpes(j)/fx2)**2-(1.d0-fces(j)/fx2)*(1.d0-npar_freq**2)'
c      write(*,*)(fpes(j)/fx2)**2-(1.d0-fces(j)/fx2)*(1.d0-npar_freq**2)
      enddo 

      CLOSE(10)


      OPEN(10,file='freqion.bin',form='unformatted') 
c-----fpi for the first ion component
      do j=1,nsteps_freq 
        fpi=charge(2)*fpes(j)/dsqrt(dmas(2))
        write(10) real(rs(j)),real(fpi),real(rhos(j))
      enddo
      WRITE(10)
c-----fci for the first ion component
      do j=1,nsteps_freq
        fci= charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))
      enddo
      WRITE(10)
c-----lh
      do j=1,nsteps_freq 
        fpi=charge(2)*fpes(j)/dsqrt(dmas(2))
        fci= charge(2)*fces(j)/dmas(2)
        flh=dsqrt(fci**2+fpi**2*fces(j)**2/(fces(j)**2+fpes(j)**2))
        write(10) real(rs(j)),real(flh),real(rhos(j)) 
      enddo
      WRITE(10)

c-----2*fci for the first ion component
      do j=1,nsteps_freq
        fci=2*charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))  
      enddo

      WRITE(10)
c-----3*fci for the first ion component
      do j=1,nsteps_freq
        fci=3*charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))  
      enddo
      WRITE(10)

c-----4*fci for the first ion component
      do j=1,nsteps_freq
        fci=4*charge(2)*fces(j)/dmas(2)
        write(10) real(rs(j)),real(fci),real(rhos(j))  
      enddo
      WRITE(10)

      CLOSE(10)

      return
      end
c----------------subroutine mk_gronetwo1-------------*
c Prepares for file drawonet.in for xdraw            *
c INPUT: from common block  'onetwo.i'              *
c****************************************************
      SUBROUTINE mk_gronetwo_1
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'onetwo.i'

      if(myrank.ne.0) return

      OPEN(11,file='onetwo1.bin',form='unformatted')
cBH110331      OPEN(21,file='onetwo1.dat')
20    format(11e12.3)

c      write(*,*)'mk_gronetwo_1 NR=',NR

cBH110313      WRITE(21,*)'i rho powden spower powden_e powden_i powden_cl
cBH110313     *curden_par curden_onetwo curden_tor cur_den_pol currden '

      hrho=1.d0/(NR-1)
 
      DO 10 I=1,NR-1
         rho=hrho*(I-1+0.5d0)
         WRITE(11) REAL(rho),REAL(powden(I)),
     &   REAL(spower(I)),
     &   REAL(powden_e(I)),REAL(powden_i(I)),REAL(powden_cl(I)),   
     &   REAL(s_cur_den_parallel(I)),REAL(s_cur_den_onetwo(I)),
     &   REAL(s_cur_den_toroidal(I)),REAL(s_cur_den_poloidal(I)),
     &   REAL(currden(I))
 
cBH110313         WRITE(21,20)rho,powden(I),
cBH110313     &   spower(I),
cBH110313     &   powden_e(I),powden_i(I),powden_cl(I),   
cBH110313     &   s_cur_den_parallel(I),s_cur_den_onetwo(I),
cBH110313     &   s_cur_den_toroidal(I),s_cur_den_poloidal(I),
cBH110313     &   currden(I)
        
c         WRITE(*,*)i,rho,powden(I),spower(I),
c     &   powden_e(I),powden_i(I),powden_cl(I),   
c     &   s_cur_den_parallel(I),s_cur_den_onetwo(I),
c     &   s_cur_den_toroidal(I),s_cur_den_poloidal(I),
c     &   currden(I)

10    CONTINUE
      WRITE(11)
     
      CLOSE(11)
cBH110313      CLOSE(21)

      write(*,*)'mk_graph.f: mk_gronetwo1 created onewto1.bin'
    
      RETURN
      END



      subroutine MK_graph_chamber_wall
c---------------------------------------------------
c     writes file wall.bin
c     with chamber wall coordinates
c----------------------------------------------------
      implicit none
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one_nml.i'
      include 'gr.i'
      include 'write.i'
      include 'fourb.i'
       
      integer i,j,m

      if(myrank.ne.0) return

      open(1,file='wall.bin',form='unformatted')
        
      DO J=1,3 
       DO I=1,NP+1
         WRITE(1) REAL(AR(J,I)),REAL(AZ(J,I))
       enddo
      enddo
      WRITE(1) 

      DO J=4,NL
       DO I=1,NP+1
         WRITE(1) REAL(AR(J,I)),REAL(AZ(J,I))
       enddo
      enddo
      WRITE(1) 

c      write(*,*)'mk_graph.f n_wall_add',n_wall_add

      if (n_wall.gt.1) then
         do i=1,n_wall  
           WRITE(1) REAL(r_wall(i)*100.0),REAL(z_wall(i)*100.0)
         enddo
         WRITE(1)
         do m=0,max_limiters          
            do i=1,n_wall_add(m)  
              WRITE(1) REAL(r_wall_add(i,m)*100.0),
     &                 REAL(z_wall_add(i,m)*100.0)
            enddo !i
            WRITE(1)
         enddo !m
      endif

      CLOSE(1)
      
      return
      end

      subroutine plot_cold_n_perp_omega_npar(z,r,phi,cnpar,
     &ratio_freq_min,ratio_freq_max,n_plot_freq)
c     
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'
c-----input
      real*8 
     &z,r,phi,         ! space coordinates
     & ratio_freq_min, ! frequency_min/frqncy
     & ratio_freq_max,  ! frequency_max/frqncy
     & cnpar           ! N_parallel

      integer n_plot_freq !number of points of frequency mesh
c-----locals
      real*8 freq_loc,freq_min,freq_max,step,df, cnpar_loc,
     &eps,g,eta,det,yi,xi,a2,a0,n_z,n_x

      complex*16 im_one,det_x,det_y,det_z,e_x,e_y,e_z
      real*8     det_xy,det_yz,det_xz,det_max,emod,emod2

      integer j

      real*8, dimension(1:nbulk) :: v_loc,w_loc
      real*8, dimension(1:n_plot_freq) :: freq_ar,n_perp2_p,n_perp2_m
      real*8, dimension(1:n_plot_freq) :: cnpar_ar
      complex*16, dimension(1:n_plot_freq) :: Ex_p,Ey_p,Ez_p,E_long_p
      complex*16, dimension(1:n_plot_freq) :: Ex_m,Ey_m,Ez_m,E_long_m
      real*8, dimension(1:n_plot_freq) :: temp

      integer n,i,n_rung
c-----external
      real*8 b,x,y

      if(myrank.ne.0) return
      
      im_one=dcmplx(0.d9,1.d0)

      freq_min=ratio_freq_min*frqncy
      freq_max=ratio_freq_max*frqncy

      step=(freq_max-freq_min)/(n_plot_freq-1)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  v_loc(i)=v(i)
	  w_loc(i)=w(i)
      enddo

      do n=1,n_plot_freq

         freq_loc=  freq_min+step*(n-1)  
         freq_ar(n)= freq_loc
         df=frqncy/freq_loc 

         do i=1,nbulk
       	    v(i)=v_loc(i)*df* df
	    w(i)=w_loc(i)*df
         enddo !i
       
c         cnpar_loc=cnpar*df
         cnpar_loc=cnpar
         cnpar_ar(n)=cnpar_loc

         eps=1.d0
         g=0.d0
         eta=1.d0

         do i=1,nbulk
            xi=x(z,r,phi,i)
            yi=y(z,r,phi,i)

       	    eps=eps-xi/(1.d0-yi**2)
            eta=eta-xi

            if (i.eq.1) then
               g=g+xi*yi/(1.d0-yi**2)
            else
	       g=g-xi*yi/(1.d0-yi**2)
            endif  
         enddo !i

c--------dispersion relation
c        eps*N_perp**4-a2*N_perp**2+a0=0

         a2=(eps+eta)*(eps-cnpar_loc**2)-g*g
         a0=eta*((eps-cnpar_loc**2)**2-g*g)

         det=a2*a2-4.d0*eps*a0
  
         if (det.ge.0) then
           n_perp2_p(n)=(a2+dsqrt(det))/(2.d0*eps)
           n_perp2_m(n)=(a2-dsqrt(det))/(2.d0*eps)
         else
           n_perp2_p(n)=0.d0
           n_perp2_m(n)=0.d0
         endif

c        write(*,*)'n,eps*n_perp2_p(n)**2-a2*n_perp2_p(n)+a0',
c     &              n,eps*n_perp2_p(n)**2-a2*n_perp2_p(n)+a0
c         write(*,*)'n,n_perp2_p(n)',n,n_perp2_p(n)

cc         write(*,*)'n,eps*n_perp2_m(n)**2-a2*n_perp2_m(n)+a0',
cc     &              n,eps*n_perp2_m(n)**2-a2*n_perp2_m(n)+a0

c------------------------------------------------------------
c        polarizatiion
c
c        (eps-N_z^2)E_x + igE_y                + N_xN_zE_z      =0
c
c         - igE_x       + (eps-N_x^2-N_z^2)E_y                  =0
c
c         N_zN_xE_x     +                        (eta-N_x^2)E_z =0
c
c-------------------------------------------------------------------
         
          if (n_perp2_p(n).ge.0) then
             n_z=cnpar_loc
             n_x=dsqrt(n_perp2_p(n))

             det_xy=(eps-n_z**2)*(eps-n_x**2-n_z**2)-g**2
             det_yz=(eps-n_x**2-n_z**2)*(eta-n_x**2)    
             det_xz=(eps-n_z**2)*(eta-n_x**2)-n_x**2*n_z**2    
            
             det_max=dabs(det_yz)
             j=1

c             write(*,*)'n_z,n_x,det_xy,det_yz,det_xz,det_max',
c     &                  n_z,n_x,det_xy,det_yz,det_xz,det_max

             if(dabs(det_yz).gt.det_max) then
                det_max=dabs(det_yz)
                j=2
             endif
             if(dabs(det_xz).gt.det_max) then
                det_max=dabs(det_xz)
                j=3
             endif

             if (det_max.gt.0) then
                n_rung=2
             else
                n_rung=1
             endif  

c             write(*,*)'j,n_rung',j,n_rung

             if(j.eq.1) then
               e_z=1.d0
               det_x=-n_z*n_x*(eps-n_x**2-n_z**2)
               det_y= im_one*g*n_z*n_x
               e_x=det_x/det_xy
               e_y=det_y/det_xy
             endif

             if(j.eq.2) then
               e_x=1.d0
               det_y=im_one*g*(eta-n_x**2)
               det_z=-(eps-n_x**2-n_z**2)*n_x*n_z 
               e_y=det_y/det_yz
               e_z=det_z/det_yz
             endif

             if(j.eq.3) then
               e_y=1.d0
               det_x=im_one*g*(eta-n_x**2)
               det_z=im_one*g*n_x*n_z 
               e_x=det_x/det_xz
               e_z=det_z/det_xz
             endif
             emod2=dreal(e_x*dconjg(e_x)+e_y*dconjg(e_y)+
     &                   e_z*dconjg(e_z))
             emod=dsqrt(emod2)
             Ex_p(n)=e_x/emod
             Ey_p(n)=e_y/emod
             Ez_p(n)=e_z/emod
cSAP110117
c             E_long_p(n)=(Ez_p(n)*cnpar_loc+Ex_p(n)*n_perp2_p(n))/
c     &                   dsqrt(cnpar_loc**2+n_perp2_p(n)**2)
             E_long_p(n)=(Ez_p(n)*cnpar_loc+Ex_p(n)*n_x)/
     &                   dsqrt(cnpar_loc**2+n_perp2_p(n))
          else
             Ex_p(n)=(0.0,0.0)
             Ey_p(n)=(0.0,0.0)
             Ez_p(n)=(0.0,0.0)
             E_long_p(n)=(0.0,0.0)
          endif
c           write(*,*)'Ex N_p', Ex_p(n)
c          write(*,*)'Ey N_p', Ey_p(n)
c          write(*,*)'Ez N_p', Ez_p(n)
c          write(*,*)'E_long N_p', E_long_p(n)
c-----------------------------------------      
c          write(*,*)'n_perp2_m(n)',n_perp2_m(n)

          if (n_perp2_m(n).ge.0) then
             n_z=cnpar_loc
             n_x=dsqrt(n_perp2_m(n))

             det_xy=(eps-n_z**2)*(eps-n_x**2-n_z**2)-g**2
             det_yz=(eps-n_x**2-n_z**2)*(eta-n_x**2)    
             det_xz=(eps-n_z**2)*(eta-n_x**2)-n_x**2*n_z**2    
            
             det_max=dabs(det_yz)
           
             j=1
             if(dabs(det_yz).gt.det_max) then
                det_max=dabs(det_yz)
                j=2
             endif
             if(dabs(det_xz).gt.det_max) then
                det_max=dabs(det_xz)
                j=3
             endif

c             write(*,*)'det_max,j',det_max,j

             if (det_max.gt.0) then
                n_rung=2
             else
                n_rung=1
             endif  

             if(j.eq.1) then
               e_z=1.d0
               det_x=-n_z*n_x*(eps-n_x**2-n_z**2)
               det_y= im_one*g*n_z*n_x
               e_x=det_x/det_xy
               e_y=det_y/det_xy
c               write(*,*)'j,det_xy,e_z,e_x,e_y',j,det_xy,e_z,e_x,e_y
             endif

             if(j.eq.2) then
               e_x=1.d0
               det_y=im_one*g*(eta-n_x**2)
               det_z=-(eps-n_x**2-n_z**2)*n_x*n_z 
               e_y=det_y/det_yz
               e_z=det_z/det_yz
c               write(*,*)'j,det_yz,e_z,e_x,e_y',j,det_yz,e_z,e_x,e_y
             endif

             if(j.eq.3) then
               e_y=1.d0
               det_x=im_one*g*(eta-n_x**2)
               det_z=im_one*g*n_x*n_z 
               e_x=det_x/det_xz
               e_z=det_z/det_xz
c               write(*,*)'j,det_xz,e_z,e_x,e_y',j,det_xz,e_z,e_x,e_y
             endif
             emod2=dreal(e_x*dconjg(e_x)+e_y*dconjg(e_y)+
     &                   e_z*dconjg(e_z))
             emod=dsqrt(emod2)
             Ex_m(n)=e_x/emod
             Ey_m(n)=e_y/emod
             Ez_m(n)=e_z/emod
cSAP110118
c             E_long_m(n)=(Ez_m(n)*cnpar_loc+Ex_m(n)*n_perp2_m(n))/
c     &                   dsqrt(cnpar_loc**2+n_perp2_m(n)**2)
             E_long_m(n)=(Ez_m(n)*cnpar_loc+Ex_m(n)*n_x)/
     &                   dsqrt(cnpar_loc**2+n_perp2_m(n))

          else
             Ex_m(n)=(0.0,0.0)
             Ey_m(n)=(0.0,0.0)
             Ez_m(n)=(0.0,0.0)
             E_long_m(n)=(0.0,0.0)
          endif
c          write(*,*)'Ex N_m', Ex_m(n)
c          write(*,*)'Ey N_m', Ey_m(n)
c          write(*,*)'Ez N_m', Ez_m(n)
c          write(*,*)'E_long N_m', E_long_m(n)
c---------------------------------------------------     
      enddo !n

      do i=1,nbulk
       	 v(i)=v_loc(i)
	 w(i)=w_loc(i)
      enddo !i   

      call plot1dt(freq_ar,cnpar_ar,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','N_par')

       call plot1dt(freq_ar,n_perp2_p,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','N_perp**2_p')
 
      do n=1,n_plot_freq
         temp(n)=dreal(Ex_p(n))
         write(*,*)'n,Ex_p(n),temp(n)',n,Ex_p(n),temp(n)
      enddo

      write(*,*)'n_perp2_p irayl'

      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ReE_x at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dimag(Ex_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ImE_x at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dreal(Ey_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ReE_y at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dimag(Ey_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ImE_y at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dreal(Ez_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ReE_z at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dimag(Ez_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ImE_z at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dreal( E_long_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','Re E_long_p at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dimag( E_long_p(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','Im E_long_p at N_perp_p')

c-----------------------------------------------------
      call plot1dt(freq_ar,n_perp2_m,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','N_perp**2_m')

      do n=1,n_plot_freq
         temp(n)=dreal(Ex_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ReE_x at N_perp_m')

      do n=1,n_plot_freq
         temp(n)=dimag(Ex_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ImE_x at N_perp_m')

      do n=1,n_plot_freq
         temp(n)=dreal(Ey_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ReE_y at N_perp_m')

      do n=1,n_plot_freq
         temp(n)=dimag(Ey_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ImE_y at N_perp_m')

      do n=1,n_plot_freq
         temp(n)=dreal(Ez_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ReE_z at N_perp_m')

      do n=1,n_plot_freq
         temp(n)=dimag(Ez_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','ImE_z at N_perp_m')

      do n=1,n_plot_freq
         temp(n)=dreal( E_long_m(n))
      enddo
      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','Re E_long_m at N_perp_p')

      do n=1,n_plot_freq
         temp(n)=dimag( E_long_m(n))
      enddo

      call plot1dt(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'freq GHZ','Im E_long_m at N_perp_p')

      return
      end

      

      
      subroutine plot_lsc_powers
c-----plot old and new power along ray at one LSC power iteration
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c      include 'onetwo.i'   

c-----locals
      real*4
     &delpwr_min,delpwr_max, 
     &s_pol_min,s_pol_max  !YuP: arguments of PGPLOT subs are real*4
      integer is,iray

      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal s_poloidal along rays
c------------------------------------------------------------------
      delpwr_min=0.0 
      delpwr_max=0.0

      s_pol_min=0.0
      s_pol_max=0.0

c      write(*,*)'in plot_lsc_powers nrayl',nrayl
      do iray=1,nrayl
        do is=1,nrayelt_nc(iray)
         if(delpwr_nc(is,iray).gt.delpwr_max)  
     &     delpwr_max=delpwr_nc(is,iray)

         if(ws_nc(is,iray).gt.s_pol_max) 
     &     s_pol_max=ws_nc(is,iray)
        enddo       
      enddo

c      write(*,*)'in plot_lsc_powers delpwr_max,s_pol_max',
c     & delpwr_max,s_pol_max


C***** Plot  dimensionless scale
c      write(*,*)' in plot_lsc_powers 1'  

      CALL PGSLW(4)
      CALL PGSCH(R41P5)
      CALL PGENV(s_pol_min,s_pol_max,delpwr_min,delpwr_max,0,1)

c      write(*,*)' in plot_lsc_powers 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='pol. dist. (cm)'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in plot_lsc_powers 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'delpwr'
     &  ,titlestr(1:titlelen))
        
       
c      write(*,*)' in plot_lsc_powers 4'
  
 
C***  note on screen, BLACK=0 means background color
C***  WHITE=1 means default foreground color
c      CALL PGSLS(1)
c      CALL PGSCI(1)
c      CALL PGSLW(6) 
c      CALL PGLINE(nsteps_freq,rs,fces)

C***  plot delpwr 
     
      DO iray=1,nrayl

c         write(*,*)' in plot_lsc_powers 5 iray,nrayelt_nc(iray)',
c     &               iray,nrayelt_nc(iray)
  
 
         allocate(newy(1:nrayelt_nc(iray)),STAT=istat)
         !write(*,*)' in plot_lsc_powers 5_1 istat',istat
         !write(*,*)' nrayelt_nc(iray)=',nrayelt_nc(iray) !YuP: should be at least 2
     
         call r4bcast(newy,zero,SIZE(newy))
         !write(*,*)' in plot_lsc_powers 6'
  
 
         allocate(newx(1:nrayelt_nc(iray)),STAT=istat)
         call r4bcast(newx,zero,SIZE(newx))
!         write(*,*)'in  plot_lsc_powers aft newx istat',
!     &              istat 

c         write(*,*)' in plot_lsc_powers 7'
  
         DO is=1,nrayelt_nc(iray)
            newy(is)=delpwr_nc(is,iray)  
            newx(is)=ws_nc(is,iray)  
         ENDDO

         !write(*,*)' in plot_lsc_powers 8'
  
         CALL PGSLS(2)
         CALL PGSCI(BLUE)
         CALL PGLINE(nrayelt_nc(iray),newx,newy)
c         print *,'BLUE: new delpwr'

         !write(*,*)' in plot_lsc_powers 9'
  
         DO is=1,nrayelt_nc(iray)
            newy(is)=delpwr_nc_old_lsc(is,iray)    
         ENDDO

   
         !write(*,*)' in plot_lsc_powers 10'
  
         CALL PGSLS(3)
         CALL PGSCI(RED)
         CALL PGLINE(nrayelt_nc(iray),newx,newy)

         !write(*,*)' in plot_lsc_powers 11'
  
c         print *,'RED: old delpwr'

         CALL PGSLS(1)  
         CALL PGSCI(GREEN)
         newx(1)=0.
         newx(2)=0.
         newy(1)=delpwr_min
         newy(2)=delpwr_max
         CALL PGLINE(2,newx,newy)

         deallocate (newx,STAT=istat)
         deallocate (newy,STAT=istat)

      ENDDO

c        CALL PGSLS(4)
c        CALL PGSCI(GREEN)
c        newx(1)=raxis
c        newx(2)=raxis
c        newy(1)=ymin
c        newy(2)=ymax
c        CALL PGLINE(2,newx,newy)

C***** magnetic axis
C** Now add labels   :BEED TO ADJUST FOR .not.ARIES (BH)
c        CALL PGSCI(WHITE)
c        CALL PGTEXT(4.0,31.0,'f\dce\u')
c        CALL PGTEXT(3.8,70.0,'2 f\dce\u')
c        CALL PGTEXT(3.8,110.0,'3 f\dce\u')
c        CALL PGSCI(RED)
c        CALL PGTEXT(1.4,45.0,'f\dpe\u')
c        CALL PGSCI(BLUE)
c        CALL PGTEXT(2.0,165.0,'f\duh\u')
c        CALL PGSCI(GREEN)
c        CALL PGTEXT(4.75,10.0,'R\dm\u')
c        call PGEND



      return
      end

      subroutine plot_lsc_d_ql(irho)
c-----plot old and new QL coefficient 
c      at given small radiual point irho
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c      include 'onetwo.i'   
c-----input
      integer irho !number of radial point 1=<irho<NR-1

c-----locals
      real*4
     &d_ql_lsc_lh_ar_min,d_ql_lsc_lh_ar_max,
     &v_par_dt_min,v_par_dt_max

      integer nv
      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal s_poloidal along rays
c------------------------------------------------------------------
      d_ql_lsc_lh_ar_min=0.d0 
      d_ql_lsc_lh_ar_max=0.d0

      v_par_dt_min=0.d0
      v_par_dt_max=0.d0
            
      
      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         if(d_ql_lsc_lh_ar(irho,nv).gt.d_ql_lsc_lh_ar_max)  
     &     d_ql_lsc_lh_ar_max=d_ql_lsc_lh_ar(irho,nv)

       if(d_ql_lsc_lh_ar(irho,nv).lt.d_ql_lsc_lh_ar_min)  
     &     d_ql_lsc_lh_ar_min=d_ql_lsc_lh_ar(irho,nv)

       if(v_par_mesh_lsc(nv).gt.v_par_dt_max) 
     &     v_par_dt_max=v_par_mesh_lsc(nv)

       if(v_par_mesh_lsc(nv).lt.v_par_dt_min) 
     &     v_par_dt_min=v_par_mesh_lsc(nv)
      enddo

c      write(*,*)'plot_lsc_d_ql d_ql_lsc_lh_ar_min,d_ql_lsc_lh_ar_max'
c     &          ,d_ql_lsc_lh_ar_min,d_ql_lsc_lh_ar_max

      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         if(d_ql_lsc_lh_ar_old(irho,nv).gt.d_ql_lsc_lh_ar_max)  
     &     d_ql_lsc_lh_ar_max=d_ql_lsc_lh_ar_old(irho,nv)

       if(d_ql_lsc_lh_ar_old(irho,nv).lt.d_ql_lsc_lh_ar_min)  
     &     d_ql_lsc_lh_ar_min=d_ql_lsc_lh_ar_old(irho,nv)

      
      enddo

C***** Plot  dimensionless scale
      write(*,*)'in plot_lsc_d_ql 1  v_par_dt_min,v_par_dt_max',
     &           v_par_dt_min,v_par_dt_max
      write(*,*)'plot_lsc_d_ql d_ql_lsc_lh_ar_min,d_ql_lsc_lh_ar_max'
     &          ,d_ql_lsc_lh_ar_min,d_ql_lsc_lh_ar_max

      CALL PGSLW(4)
      CALL PGSCH(R41P5)
      CALL PGENV(v_par_dt_min,v_par_dt_max,
     &           d_ql_lsc_lh_ar_min,d_ql_lsc_lh_ar_max,0,1)

c      write(*,*)' in plot_lsc__d_ql 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='v par dev vt '
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in plot_lsc_d_ql 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'D_QL'
     &  ,titlestr(1:titlelen))
        

       
c      write(*,*)' in plot_lsc__d_ql 4'
  
 
C***  note on screen, BLACK=0 means background color
C***  WHITE=1 means default foreground color
c      CALL PGSLS(1)
c      CALL PGSCI(1)
c       CALL PGSLW(6) 
c      CALL PGLINE(nsteps_freq,rs,fces)

C***  plot d_ql 
         
c      write(*,*)' in plot_lsc__d_ql 5  bef allocate newy'
  
      allocate(newy(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
c      write(*,*)'in plot_lsc_d_ql 5_1 istat',istat
     
      call r4bcast(newy,zero,SIZE(newy))
        
c      write(*,*)' in plot_lsc__d_ql 6'
  
      allocate(newx(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))

c      write(*,*)'in  plot_lsc__d_ql aft newx istat',istat 

c      write(*,*)' in plot_lsc__d_ql 7'
  
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=d_ql_lsc_lh_ar(irho,nv)  
         newx(nv)=v_par_mesh_lsc(nv) 
c         write(*,*)'id_ql_lsc_lh_ar irho,nv,newy(nv)',irho,nv,newy(nv)       
      ENDDO

c      write(*,*)' in plot_lsc__d_ql 8'
  
      CALL PGSLS(2)
      CALL PGSCI(BLUE)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)

c      print *,'BLUE: new d_ql'

c      write(*,*)' in plot_lsc_d_ql 9'
  
   
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=d_ql_lsc_lh_ar_old(irho,nv)
c         write(*,*)'id_ql_lsc_lh_old irho,nv,newy(nv)',irho,nv,newy(nv)       
      ENDDO

c      write(*,*)' in plot_lsc__d_ql 10'
  
      CALL PGSLS(3)
      CALL PGSCI(RED)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)

c      write(*,*)' in plot_lsc__d_ql 11'
  
c      print *,'RED: old d_ql'

      
c        CALL PGSLS(4)
c        CALL PGSCI(GREEN)
c        newx(1)=raxis
c        newx(2)=raxis
c        newy(1)=ymin
c        newy(2)=ymax
c        CALL PGLINE(2,newx,newy)

c      CALL PGSLS(4)
c      CALL PGSCI(BLACK)
c      newx(1)=0.
c      newx(2)=0.
c      newy(1)=d_ql_lsc_lh_ar_min
c      newy(2)=d_ql_lsc_lh_ar_max
c      CALL PGLINE(2,newx,newy)

       CALL PGSLS(1)  
       CALL PGSCI(GREEN)
       newx(1)=0.
       newx(2)=0.
       newy(1)=d_ql_lsc_lh_ar_min
       newy(2)=d_ql_lsc_lh_ar_max
       CALL PGLINE(2,newx,newy)


      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)
     
      return
      end
      

      subroutine plot_log_fe(irho)
c-----plot old and new log_e [fe(v_par**2*sign(v_par))] 
c      at given small radiual point irho
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c      include 'onetwo.i'   
c-----input
      integer irho !number of radial point 1=<irho<NR-1

c-----locals
      real*4
     &Lnfe_min,Lnfe_max,
     &v_par_dt_min,v_par_dt_max,  
     &sv_par_dt_min,sv_par_dt_max,sign
    

      integer nv
      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      write(*,*)'in subroutine plot_log_fe irho=',irho

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal s_poloidal along rays
c------------------------------------------------------------------
      Lnfe_min=0.d0 
      Lnfe_max=0.d0

      v_par_dt_min=0.
      v_par_dt_max=0.
                   
      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh

c         write(*,*)'nv,fe_lsc(irho,nv)',nv,fe_lsc(irho,nv)
c         write(*,*)'v_par_mesh_lsc(nv)',v_par_mesh_lsc(nv)
         if(dlog(abs(fe_lsc(irho,nv))).gt. Lnfe_max)  
     &      Lnfe_max=dlog(abs(fe_lsc(irho,nv)))

         if(dlog(abs(fe_lsc(irho,nv))).lt. Lnfe_min)  
     &      Lnfe_min=dlog(abs(fe_lsc(irho,nv)))

c        write(*,*)'Lnfe_max,Lnfe_min',Lnfe_max,Lnfe_min

       if(v_par_mesh_lsc(nv).gt.v_par_dt_max) 
     &     v_par_dt_max=v_par_mesh_lsc(nv)

       if(v_par_mesh_lsc(nv).lt.v_par_dt_min) 
     &     v_par_dt_min=v_par_mesh_lsc(nv)
 
      enddo

      if(v_par_dt_min.ge.0.d0) then 
        sign=1.d0
      else
        sign=-1.d0
      endif
      sv_par_dt_min=v_par_dt_min**2*sign

      if(v_par_dt_max.ge.0.d0) then 
        sign=1.d0
      else
        sign=-1.d0
      endif
      sv_par_dt_max=v_par_dt_max**2*sign

c      write(*,*)' plot_log_fe Lnfe_min,Lnfe_max',
c     &                         Lnfe_min,Lnfe_max

      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh

         if(dlog(abs(fe_lsc_old(irho,nv))).gt. Lnfe_max)  
     &      Lnfe_max=dlog(abs(fe_lsc_old(irho,nv)))


         if(dlog(abs(fe_lsc_old(irho,nv))).lt. Lnfe_min)  
     &      Lnfe_min=dlog(abs(fe_lsc_old(irho,nv)))
            
      enddo

C***** Plot  dimensionless scale
c      write(*,*)'in plot_log_fe 1  v_par_dt_min,v_par_dt_max',
c     &           v_par_dt_min,v_par_dt_max
c      write(*,*)'in plot_log_fe 1  sv_par_dt_min,sv_par_dt_max',
c     &           sv_par_dt_min,sv_par_dt_max
c      write(*,*)'plot_log_fe  Lnfe_min Lnfe_max',Lnfe_min,Lnfe_max

    
      CALL PGSLW(4)
      CALL PGSCH(R41P5)
c      CALL PGENV(v_par_dt_min,v_par_dt_max,
c    &           Lnfe_min,Lnfe_max,0,1)
      CALL PGENV(sv_par_dt_min,sv_par_dt_max,
     &           Lnfe_min,Lnfe_max,0,1)

c      write(*,*)' in plot_log_fe 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='(v_par/vt)**2*sign(v_par)'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in plot_log_fe 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'log_e[fe]'
     &  ,titlestr(1:titlelen))
        

       
c      write(*,*)' in plot_log_fe 4'
  
 
C***  note on screen, BLACK=0 means background color
C***  WHITE=1 means default foreground color
      CALL PGSLS(1)
      CALL PGSCI(1)
      CALL PGSLW(6) 
c      CALL PGLINE(nsteps_freq,rs,fces)

C***  plot d_ql 
         
c      write(*,*)' in plot_log_fe 5  bef allocate newy'
  
      allocate(newy(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
 
c      write(*,*)'in plot_log_fel 5_1 istat',istat
     
      call r4bcast(newy,zero,SIZE(newy))
        
c      write(*,*)' in plot_log_fe 6'
  
      allocate(newx(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))

c      write(*,*)'in plot_log_fe aft newx istat',istat 

c      write(*,*)' in plot_log_fe 7'
  
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=log(abs(fe_lsc(irho,nv)))
c         newx(nv)=v_par_mesh_lsc(nv)
      
         if(v_par_mesh_lsc(nv).ge.0.d0) then 
           sign=1.d0
         else
           sign=-1.d0
         endif
         newx(nv)=v_par_mesh_lsc(nv)**2*sign

      ENDDO

c      write(*,*)' in plot_log_fe 8'
  
      CALL PGSLS(2)
      CALL PGSCI(BLUE)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)
c      print *,'BLUE: new plot_log_fe'

c      write(*,*)' in plot_log_fe 9'
  
   
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=log(abs(fe_lsc_old(irho,nv)))
c         write(*,*)'log_fe_old irho,nv,newy(nv)',irho,nv,newy(nv)       
      ENDDO

c      write(*,*)' in plot_log_fe 10'
  
      CALL PGSLS(3)
      CALL PGSCI(RED)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)

c      write(*,*)' in plot_log_fe 11'
  
c      print *,'RED: old log_fe'

    
      CALL PGSLS(1)
      CALL PGSCI(GREEN)
c      newx(1)=v_par_dt_min
c      newx(2)=v_par_dt_min
      newx(1)=sv_par_dt_min
      newx(2)=sv_par_dt_min
      newy(1)=Lnfe_min
      newy(2)=Lnfe_max
      CALL PGLINE(2,newx,newy)


      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)
      

c      write(*,*)'end plot_log_fe'

      return
      end
      
     
     
      subroutine plot_integral_lsc(irho)
c-----plot old and new log_e [fe(v_par**2*sign(v_par))] 
c      at given small radiual point irho
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c      include 'onetwo.i'   
c-----input
      integer irho !number of radial point 1=<irho<NR-1

c-----locals
      real*4
     &Lnfe_min,Lnfe_max,
     &v_par_dt_min,v_par_dt_max,  
     &sv_par_dt_min,sv_par_dt_max,sign,
     &integral_min,integral_max

      integer nv
      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal s_poloidal along rays
c------------------------------------------------------------------
      Lnfe_min=0.d0 
      Lnfe_max=0.d0

      v_par_dt_min=0.
      v_par_dt_max=0.
      integral_min=0.
      integral_max=0.            
      
      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         if(dlog(abs(fe_lsc(irho,nv))).gt. Lnfe_max)  
     &      Lnfe_max=dlog(abs(fe_lsc(irho,nv)))

         if(dlog(abs(fe_lsc(irho,nv))).lt. Lnfe_min)  
     &      Lnfe_min=dlog(abs(fe_lsc(irho,nv)))
       
       if(v_par_mesh_lsc(nv).gt.v_par_dt_max) 
     &     v_par_dt_max=v_par_mesh_lsc(nv)

       if(v_par_mesh_lsc(nv).lt.v_par_dt_min) 
     &     v_par_dt_min=v_par_mesh_lsc(nv)
 
c       write(*,*)'irho,nv,integral_fe_lsc(irho,nv)',
c     &            irho,nv,integral_fe_lsc(irho,nv)

       if(integral_fe_lsc(irho,nv).gt.integral_max) 
     &     integral_max=integral_fe_lsc(irho,nv)

       if(integral_fe_lsc(irho,nv).lt.integral_min) 
     &     integral_min=integral_fe_lsc(irho,nv)
      enddo

      if(v_par_dt_min.ge.0.d0) then 
        sign=1.d0
      else
        sign=-1.d0
      endif
      sv_par_dt_min=v_par_dt_min**2*sign

      if(v_par_dt_max.ge.0.d0) then 
        sign=1.d0
      else
        sign=-1.d0
      endif
      sv_par_dt_max=v_par_dt_max**2*sign

c      write(*,*)'v_par_dt_min,v_par_dt_max',
c     &           v_par_dt_min,v_par_dt_max
 
c     write(*,*)'sv_par_dt_min,sv_par_dt_max',
c     &          sv_par_dt_min,sv_par_dt_max

c      write(*,*)' plot_log_fe Lnfe_min,Lnfe_max',
c     &                        Lnfe_min,Lnfe_max

c      write(*,*)' plot_log_fe integral_min,integral_max',
c     &                        integral_min,integral_max


      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh

         if(dlog(abs(fe_lsc_old(irho,nv))).gt. Lnfe_max)  
     &      Lnfe_max=dlog(abs(fe_lsc_old(irho,nv)))


         if(dlog(abs(fe_lsc_old(irho,nv))).lt. Lnfe_min)  
     &      Lnfe_min=dlog(abs(fe_lsc_old(irho,nv)))
            
      enddo

C***** Plot  dimensionless scale
c      write(*,*)'in plot_integra1  v_par_dt_min,v_par_dt_max',
c     &           v_par_dt_min,v_par_dt_max
c      write(*,*)'in plot_integral 1  sv_par_dt_min,sv_par_dt_max',
c     &           sv_par_dt_min,sv_par_dt_max
c      write(*,*)'plot_log_fe  integral_min,integral_max',
c     &          integral_min,integral_max

    
      CALL PGSLW(4)
      CALL PGSCH(R41P5)
c      CALL PGENV(v_par_dt_min,v_par_dt_max,
c    &           Lnfe_min,Lnfe_max,0,1)
c      CALL PGENV(sv_par_dt_min,sv_par_dt_max,
c     &           Lnfe_min,Lnfe_max,0,1)
      CALL PGENV(sv_par_dt_min,sv_par_dt_max,
     &           integral_min, integral_max,0,1)

c      write(*,*)' in plot_integral 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='(v_par/vt)**2*sign(v_par)'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in plot_integral 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'integral for log_e[fe]'
     &  ,titlestr(1:titlelen))
        

       
c      write(*,*)' in plot_integral 4'
  
 
C***  note on screen, BLACK=0 means background color
C***  WHITE=1 means default foreground color
      CALL PGSLS(1)
      CALL PGSCI(1)
      CALL PGSLW(6) 
c      CALL PGLINE(nsteps_freq,rs,fces)

C***  plot d_ql 
         
c      write(*,*)' in plot_integral 5  bef allocate newy'
  
      allocate(newy(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
 
c      write(*,*)'in plot_integral 5_1 istat',istat
     
      call r4bcast(newy,zero,SIZE(newy))
        
c      write(*,*)' in plot_integral 6'
  
      allocate(newx(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))

c      write(*,*)'in plot_integral aft newx istat',istat 

c      write(*,*)' in plot_integral 7'
  
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=integral_fe_lsc(irho,nv)
c         newx(nv)=v_par_mesh_lsc(nv)
      
         if(v_par_mesh_lsc(nv).ge.0.d0) then 
           sign=1.d0
         else
           sign=-1.d0
         endif
         newx(nv)=v_par_mesh_lsc(nv)**2*sign

      ENDDO

c      write(*,*)' in plot_integral nv_lsc_ql_lh,irho',nv_lsc_ql_lh,irho
c      write(*,*)'integral newx',newx
c      write(*,*)'integral newy',newy

      CALL PGSLS(2)
      CALL PGSCI(BLUE)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)
c      print *,'BLUE: new plot_integral '

c      write(*,*)' in plot_integral 9'
  
   
c      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
c         newy(nv)=log(abs(fe_lsc_old(irho,nv)))
c         write(*,*)'log_fe_old irho,nv,newy(nv)',irho,nv,newy(nv)       
c      ENDDO

c      write(*,*)' in plot_integral 10'
  
c      CALL PGSLS(3)
c      CALL PGSCI(RED)
c      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)

c      write(*,*)' in plot_integral 11'
  
c      print *,'RED: old log_integral'

    
      CALL PGSLS(1)
      CALL PGSCI(GREEN)
c      newx(1)=v_par_dt_min
c      newx(2)=v_par_dt_min
      newx(1)=sv_par_dt_min
      newx(2)=sv_par_dt_min
c      newy(1)=Lnfe_min
c      newy(2)=Lnfe_max
      newy(1)=integral_min
      newy(2)=integral_max
      CALL PGLINE(2,newx,newy)


      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)
      

      return
      end
      
     
    
      subroutine plot_lsc_powdens_rho
c-----plot power density versus small radius rho
c     power and CD density powden(i),currden(i):rho_bin_center(i)
c     i=1,NR-1
c
c     power and CD density calculted using QL diffusion and nonmaxwellian fe
c     power_dens_watt_m3_TSC_1D(i),CD_dens_no_Edc_a_m2_TSC_1D(i):
c     rho_bin_center_lsc(i) i=1,n_psi_TSC
c
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
      include 'onetwo.i'   

c-----locals
      real*8 hrho

      real
     &powdens_min,powdens_max, 
     &CDdens_min,CDdens_max,
     &powdens_f_min,powdens_f_max, 
     &CDdens_f_min,CDdens_f_max,
     &rho_min,rho_max,y_min,y_max
      integer i

      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal values
c------------------------------------------------------------------
      powdens_min=1.e30
      powdens_max=-1.e30
      CDdens_min=1.e30
      CDdens_max=-1.e30
      
      powdens_f_min=1.e30
      powdens_f_max=-1.e30
      CDdens_f_min=1.e30
      CDdens_f_max=-1.e30

      hrho=1.d0/(NR-1)

      do i=1,NR-1        
         rho_bin_center(i)=hrho*(i-0.5d0)

c         write(*,*)'in plot_lsc_powdens_rho i,rho_bin_center(i)'//
c     &   ',powden(i)',i,rho_bin_center(i),powden(i)

         if(powden(i).gt.powdens_max)  
     &      powdens_max=powden(i)
         if(powden(i).lt.powdens_min)  
     &      powdens_min=powden(i)

c         if(currden(i).gt.CDdens_max)  
c     &       CDdens_max=currden(i)
c         if(currden(i).lt. CDdens_min)  
c     &       CDdens_min=currden(i)

         if(s_cur_den_parallel(i).gt.CDdens_max)  
     &       CDdens_max=s_cur_den_parallel(i)
         if(s_cur_den_parallel(i).lt. CDdens_min)  
     &       CDdens_min=s_cur_den_parallel(i)
        
      enddo

c      write(*,*)'in plot_lsc_powdens_rho powdens_min,powdens_max',
c     &                                   powdens_min,powdens_max

      do i=1,n_psi_TSC 

c        write(*,*)'in plot_lsc_powdens_rho i,rho_bin_center_lsc(i)'//
c     &  ',power_dens_watt_m3_TSC_1D(i)*10',
c     &  i,rho_bin_center_lsc(i),power_dens_watt_m3_TSC_1D(i)*10

         if(power_dens_watt_m3_TSC_1D(i)*10.gt.powdens_f_max)  
     &      powdens_f_max=power_dens_watt_m3_TSC_1D(i)*10
         if(power_dens_watt_m3_TSC_1D(i)*10.lt.powdens_f_min)  
     &      powdens_f_min=power_dens_watt_m3_TSC_1D(i)*10

c         if( CD_dens_no_Edc_a_m2_TSC_1D(i)*1.d-4.gt.CDdens_f_max)  
c     &       CDdens_f_max=CD_dens_no_Edc_a_m2_TSC_1D(i)*1.d-4
c         if( CD_dens_no_Edc_a_m2_TSC_1D(i)*1.d04.lt. CDdens_f_min)  
c     &       CDdens_f_min=CD_dens_no_Edc_a_m2_TSC_1D(i)*1.d-4          

         if( CD_dens_small_Edc_a_m2_TSC_1D(i)*1.d-4.gt.CDdens_f_max)  
     &       CDdens_f_max=CD_dens_small_Edc_a_m2_TSC_1D(i)*1.d-4
         if( CD_dens_small_Edc_a_m2_TSC_1D(i)*1.d04.lt. CDdens_f_min)  
     &       CDdens_f_min=CD_dens_small_Edc_a_m2_TSC_1D(i)*1.d-4     
     
      enddo
c      write(*,*)'in plot_lsc_powdens_rho powdens_f_min,powdens_f_max',
c     &                                   powdens_f_min,powdens_f_max

c      write(*,*)'in plot_lsc_powdens_rho'

      rho_min=0.
      rho_max=1.
      y_max=max(powdens_max,powdens_f_max)
      y_min=min(powdens_min,powdens_f_min)

c      write(*,*)'in plot_lsc_powdens_rho y_min,y_max',y_min,y_max
C***** Plot  dimensionless scale
   
      CALL PGSLW(4)
      CALL PGSCH(R41P5)
      CALL PGENV(rho_min,rho_max, y_min, y_max,0,1)

c      write(*,*)' in  plot_lsc_powdens_rho 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='small radius, blue ImNperp, red fe'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in  plot_lsc_powdens_rho 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'powdens egr/s/cm**3'
     &  ,titlestr(1:titlelen))
        
       
c      write(*,*)' in plot_lsc_powdens_rho 4'
  
 
C***  note on screen, BLACK=0 means background color
C***  WHITE=1 means default foreground color
c      CALL PGSLS(1)
c      CALL PGSCI(1)
c      CALL PGSLW(6) 
c      CALL PGLINE(nsteps_freq,rs,fces)

C***  plot delpwr 
     
      
c      write(*,*)'in  plot_lsc_powdens_rho 5 '
  
 
      allocate(newy(NR-1),STAT=istat)
c      write(*,*)' in plot_lsc_powers 5_1 istat',istat
     
      call r4bcast(newy,zero,SIZE(newy))
        
c      write(*,*)' in  plot_lsc_powdens_rho 6'
  
      allocate(newx(NR-1),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))
c      write(*,*)'in plot_lsc_powdens_rho aft newx istat',istat 

  
      DO i=1,NR-1
         newy(i)=powden(i)
         newx(i)=rho_bin_center(i)
      ENDDO

c      write(*,*)' in plot_lsc_powdens_rho 8'
  
      CALL PGSLS(2)
      CALL PGSCI(BLUE)
      CALL PGLINE(NR-1,newx,newy)
c      print *,'BLUE: powdens '
      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)


c-------------------------------------------------------------

      allocate(newy(n_psi_TSC),STAT=istat)
      write(*,*)' in plot_lsc_powers 5_1 istat',istat
     
      call r4bcast(newy,zero,SIZE(newy))
        
c      write(*,*)' in  plot_lsc_powdens_rho 6'
  
      allocate(newx(n_psi_TSC),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))
c      write(*,*)'in plot_lsc_powdens_rho aft newx istat',istat 

c      write(*,*)' in plot_lsc_powdens_rho 9'
  
      DO i=1,n_psi_TSC
         newx(i)=rho_bin_center_lsc(i)
         newy(i)=power_dens_watt_m3_TSC_1D(i)*10
      ENDDO
c      write(*,*)' in plot_lsc_powdens_rho 10'
  
      CALL PGSLS(3)
      CALL PGSCI(RED)
      CALL PGLINE(n_psi_TSC,newx,newy)

c      write(*,*)' in plot_lsc_powdens 11'
  
c      print *,'RED: powdens from fe'

      CALL PGSLS(1)  
      CALL PGSCI(GREEN)
c      CALL PGSCI(BLACK)
      newx(1)=0.
      newx(2)=0.
      newy(1)=y_min
      newy(2)=y_max
      CALL PGLINE(2,newx,newy)

c
c      write(*,*)' in  plot_lsc_powdens_rho 3'
      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)
c--------------------------------------------------------------  
c     Current density plot
c------------------------------------------------------------- 
      y_max=max(CDdens_max,CDdens_f_max)
      y_min=min(CDdens_min,CDdens_f_min)
C***** Plot  dimensionless scale
   
      CALL PGSLW(4)
      CALL PGSCH(R41P5)
      CALL PGENV(rho_min,rho_max, y_min, y_max,0,1)
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='small radius, blue efficiency, red fe'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)
  
      CALL PGLAB(xtitle(1:xtitlelen),'curdens A/cm**2'
     &  ,titlestr(1:titlelen))

      allocate(newy(NR-1),STAT=istat)        
      call r4bcast(newy,zero,SIZE(newy))
      allocate(newx(NR-1),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))

      DO i=1,NR-1
c         newy(i)=currden(i) 
         newy(i)=s_cur_den_parallel(i)
         newx(i)=rho_bin_center(i)
      ENDDO
       
      CALL PGSLS(2)
      CALL PGSCI(BLUE)
      CALL PGLINE(NR-1,newx,newy)
c      print *,'BLUE: currden '

      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)


c-------------------------------------------------------------

      allocate(newy(n_psi_TSC),STAT=istat)     
      call r4bcast(newy,zero,SIZE(newy))
         
      allocate(newx(n_psi_TSC),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))
      DO i=1,n_psi_TSC
         newx(i)=rho_bin_center_lsc(i)
         newy(i)=CD_dens_no_Edc_a_m2_TSC_1D(i)*1.d-4
      ENDDO
c      write(*,*)' in plot_lsc_powdens_rho 10'
  
c      CALL PGSLS(3)
c      CALL PGSCI(RED)
c      CALL PGLINE(n_psi_TSC,newx,newy)
  
c      print *,'RED: currdens from fe'

      CALL PGSLS(1)  
      CALL PGSCI(GREEN)
c      CALL PGSCI(BLACK)
      newx(1)=0.
      newx(2)=0.
      newy(1)=y_min
      newy(2)=y_max
      CALL PGLINE(2,newx,newy)
 
      DO i=1,n_psi_TSC
         newx(i)=rho_bin_center_lsc(i)
         newy(i)=CD_dens_small_Edc_a_m2_TSC_1D(i)*1.d-4
      ENDDO

      CALL PGSLS(3)
c      CALL PGSCI(CYAN)
      CALL PGSCI(RED)
      CALL PGLINE(n_psi_TSC,newx,newy)
c      print *,'CYAN: currdens small E_DC'
c     print *,'RED: currdens small E_DC'

      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)
c------------------------------------------------------------
    
      


      return
      end

  

      subroutine MK_GRAPH_LSC(nray)
c-----creates lsc.bin file
      IMPLICIT none 
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'writencdf.i' 
c-----input
      integer nray !total number of rays

c-----local
      integer iray,is

      if(myrank.ne.0) return

      open(1,file='lsc.bin',form='unformatted')

      do iray=1,nray
        do is=1,nrayelt_nc(iray)
           WRITE(1) 
     &     REAL(ws_nc(is,iray)),REAL(delpwr_nc(is,iray))
        enddo 
        WRITE(1)
      enddo

      close(1)

      return
      end


      SUBROUTINE MK_GRAPH_delpwr
c-----prepare delp[wr.bin file for xdraw plot
      IMPLICIT double precision (a-h,o-z)

      include 'pgconst.i'
      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'writencdf.i' 

      if(myrank.ne.0) return

      OPEN(1,file='delpwr.bin',form='unformatted')
      do iray=1,nrayl
         do is=1,nrayelt_nc(iray)
            WRITE(1) REAL(delpwr_nc(is,iray)),REAL(ws_nc(is,iray))
         enddo
         WRITE(1)
      enddo
      CLOSE(1)

      return
      end
       
      subroutine plot_lsc_powers_1
c-----plot old and new power along ray at one LSC power iteration
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c-----locals
      real
     &delpwr_min,delpwr_max, 
     &s_pol_min,s_pol_max  
      integer is,iray

      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:  
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal s_poloidal along rays
c------------------------------------------------------------------
      delpwr_min=0.d0 
      delpwr_max=0.d0

      s_pol_min=0.d0
      s_pol_max=0.d0

      write(*,*)'in plot_lsc_powers_1 nrayl',nrayl
      do iray=1,nrayl
        do is=1,nrayelt_nc(iray)
         if(delpwr_nc(is,iray).gt.delpwr_max)  
     &     delpwr_max=delpwr_nc(is,iray)

         if(ws_nc(is,iray).gt.s_pol_max) 
     &     s_pol_max=ws_nc(is,iray)
        enddo       
      enddo

      write(*,*)'in plot_lsc_powers_1 delpwr_max,s_pol_max',
     & delpwr_max,s_pol_max

C***** Plot  dimensionless scale
c      write(*,*)' in plot_lsc_powers 1'  

      CALL PGSLW(4)
      CALL PGSCH(R41P5)
      CALL PGENV(s_pol_min,s_pol_max,delpwr_min,delpwr_max,0,1)

c      write(*,*)' in plot_lsc_powers 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='pol. dist. (cm)'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in plot_lsc_powers 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'delpwr'
     &  ,titlestr(1:titlelen))
        
       
c      write(*,*)' in plot_lsc_powers 4'
  
c ***  plot delpwr 
     
      DO iray=1,nrayl

         write(*,*)' in plot_lsc_powers_1 5 iray,nrayelt_nc(iray)',
     &               iray,nrayelt_nc(iray)
  
 
         allocate(newy(1:nrayelt_nc(iray)),STAT=istat)
c         write(*,*)' in plot_lsc_powers 5_1 istat',istat
     
         call r4bcast(newy,zero,SIZE(newy))
c         write(*,*)' in plot_lsc_powers 6'
  
 
         allocate(newx(1:nrayelt_nc(iray)),STAT=istat)
         call r4bcast(newx,zero,SIZE(newx))
c         write(*,*)'in  plot_lsc_powers aft newx istat',
c     &              istat 

c         write(*,*)' in plot_lsc_powers 7'
  
         DO is=1,nrayelt_nc(iray)
            newy(is)=delpwr_nc(is,iray)  
            newx(is)=ws_nc(is,iray)  
            write(*,*)'iray,is,newx(is),newy(is)',
     &                 iray,is,newx(is),newy(is)
         ENDDO
c         write(*,*)' in plot_lsc_powers_1 8'
  
         CALL PGSLS(2)
         CALL PGSCI(BLUE)
         CALL PGLINE(nrayelt_nc,newx,newy)
c         print *,'BLUE: new delpwr'

c         write(*,*)' in plot_lsc_powers 9'
         DO is=1,nrayelt_nc(iray)
            newy(is)=delpwr_nc_old_lsc(is,iray)   
            write(*,*)'iray,is,newx(is),newy_old(is)',
     &                 iray,is,newx(is),newy(is) 
         ENDDO
 
         CALL PGSLS(3)
         CALL PGSCI(RED)
         CALL PGLINE(nrayelt_nc,newx,newy)

        CALL PGSLS(1)  
         CALL PGSCI(GREEN)
         newx(1)=0.
         newx(2)=0.
         newy(1)=delpwr_min
         newy(2)=delpwr_max
         CALL PGLINE(2,newx,newy)

         deallocate (newx,STAT=istat)
         deallocate (newy,STAT=istat)


      ENDDO
c         write(*,*)' in plot_lsc_powers 11'
  
c         print *,'RED: old delpwr'

c         CALL PGSLS(1)  
c         CALL PGSCI(GREEN)
c         newx(1)=0.
c         newx(2)=0.
c         newy(1)=delpwr_min
c         newy(2)=delpwr_max
c         CALL PGLINE(2,newx,newy)

      DO iray=1,nrayl
         write(*,*)'iray,nrayelt_nc(iray)',iray,nrayelt_nc(iray)

         allocate(newy(1:nrayelt_nc(iray)),STAT=istat)
         call r4bcast(newy,zero,SIZE(newy))
         allocate(newx(1:nrayelt_nc(iray)),STAT=istat)
         call r4bcast(newx,zero,SIZE(newx))
         
         DO is=1,nrayelt_nc(iray)
            newy(is)=delpwr_nc(is,iray)  
            newx(is)=ws_nc(is,iray)
            write(*,*)'5,is,newx(is),newy(is)',is,newx(is),newy(is)
         enddo

         CALL PGSLW(4)
         CALL PGSCH(R41P5)
         CALL PGENV(s_pol_min,s_pol_max,delpwr_min,delpwr_max,0,1)

         CALL PGLAB(xtitle(1:xtitlelen),'delpwr'
     &  ,titlestr(1:titlelen))

         CALL PGSLS(2)
         CALL PGSCI(BLUE)
         CALL PGLINE(nrayelt_nc(iray),newx,newy)

         DO is=1,nrayelt_nc(iray)
            newy(is)=delpwr_nc_old_lsc(is,iray)
            write(*,*)'5old ,is,newx(is),newy(is)',is,newx(is),newy(is)
         ENDDO
 
         CALL PGSLS(3)
         CALL PGSCI(RED)
         CALL PGLINE(nrayelt_nc(iray),newx,newy)

         CALL PGSLS(1)  
         CALL PGSCI(GREEN)
         newx(1)=0.
         newx(2)=0.
         newy(1)=delpwr_min
         newy(2)=delpwr_max
         CALL PGLINE(2,newx,newy)

         deallocate (newx,STAT=istat)
         deallocate (newy,STAT=istat)

      ENDDO

      return
      end


      subroutine plot_cold_n_perp_omega_npar_iray(z,r,phi,cnpar,
     &ratio_freq_min,ratio_freq_max,n_plot_freq)
c,iray)  
c-----It plots two roots N perpendicular^2 versus frequency
c     of the cold plasma dispersion relation
c     in the form:  eps*N_perp**4-a2*N_perp**2+a0=0
c     at given N parallel=cnpar
c     in the space point (z,r,phi),
c     for ray wth the nmuber irayl givern in writnetcdf.i
c
c     The roots hve the form
c     N_perp2_p=(a2+dsqrt(det))/(2.d0*eps)
c     N_perp2_m=(a2-dsqrt(det))/(2.d0*eps)
c     Here det=a2*a2-4.d0*eps*a0
c
c     For each root it plots 
c     complex electric fied polarisation versus frequency
c     (Ex,Ey,Ez) in Stix frame and
c     the electric field along wave vector  E_long.
c
c     frequency is the interval:
c     frequency_min < frequency/frqncy <frequency_max
c     frqncy is the wave frequency given in 'one_nml.i'
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'one.i'

      include 'writencdf.i'
c-----input
      real*8 
     &z,r,phi,          ! space coordinates
     & ratio_freq_min,  ! frequency_min/frqncy
     & ratio_freq_max,  ! frequency_max/frqncy
     & cnpar            ! N_parallel
      integer
c     & irayl,            ! number of ray
     & n_plot_freq      ! number of points of frequency mesh

c-----locals
      real*8 freq_loc,freq_min,freq_max,step,df, cnpar_loc,
     &eps,g,eta,det,yi,xi,a2,a0,n_z,n_x

      complex*16 im_one,det_x,det_y,det_z,e_x,e_y,e_z
      real*8     det_xy,det_yz,det_xz,det_max,emod,emod2

      integer j

      real*8, dimension(1:nbulk) :: v_loc,w_loc
      real*8, dimension(1:n_plot_freq) :: freq_ar,n_perp2_p,n_perp2_m
      real*8, dimension(1:n_plot_freq) :: cnpar_ar
      complex*16, dimension(1:n_plot_freq) :: Ex_p,Ey_p,Ez_p,E_long_p
      complex*16, dimension(1:n_plot_freq) :: Ex_m,Ey_m,Ez_m,E_long_m
      real*8, dimension(1:n_plot_freq) :: temp

      integer n,i,n_rung

      integer  n_param                    !number of input parameters
      parameter  (n_param=2) 
      character*(10) name_param(n_param)   !names  of the input parameters
      REAL*4 param(n_param)                 !values of the input parameters

c-----external
      real*8 b,x,y

      if(myrank.ne.0) return
      
      im_one=dcmplx(0.d9,1.d0)

      freq_min=ratio_freq_min*frqncy
      freq_max=ratio_freq_max*frqncy

      step=(freq_max-freq_min)/(n_plot_freq-1)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  v_loc(i)=v(i)
	  w_loc(i)=w(i)
      enddo

      do n=1,n_plot_freq

         freq_loc=  freq_min+step*(n-1)  
         freq_ar(n)= freq_loc
         df=frqncy/freq_loc 

         do i=1,nbulk
       	    v(i)=v_loc(i)*df* df
	    w(i)=w_loc(i)*df
         enddo !i
       
         cnpar_loc=cnpar
         cnpar_ar(n)=cnpar_loc

         eps=1.d0
         g=0.d0
         eta=1.d0

         do i=1,nbulk
            xi=x(z,r,phi,i)
            yi=y(z,r,phi,i)

       	    eps=eps-xi/(1.d0-yi**2)
            eta=eta-xi

            if (i.eq.1) then
               g=g+xi*yi/(1.d0-yi**2)
            else
	       g=g-xi*yi/(1.d0-yi**2)
            endif  
         enddo !i


c--------dispersion relation
c        eps*N_perp**4-a2*N_perp**2+a0=0

         a2=(eps+eta)*(eps-cnpar_loc**2)-g*g
         a0=eta*((eps-cnpar_loc**2)**2-g*g)

         det=a2*a2-4.d0*eps*a0
  
         if (det.ge.0) then
           n_perp2_p(n)=(a2+dsqrt(det))/(2.d0*eps)
           n_perp2_m(n)=(a2-dsqrt(det))/(2.d0*eps)
         else
           n_perp2_p(n)=0.d0
           n_perp2_m(n)=0.d0
         endif

c         write(*,*)'n,eps*n_perp2_p(n)**2-a2*n_perp2_p(n)+a0',
c     &              n,eps*n_perp2_p(n)**2-a2*n_perp2_p(n)+a0
c         write(*,*)'iray.n,n_perp2_p(n)',n,n_perp2_p(n)
c------------------------------------------------------------
c        polarizatiion
c
c        (eps-N_z^2)E_x + igE_y                + N_xN_zE_z      =0
c
c         - igE_x       + (eps-N_x^2-N_z^2)E_y                  =0
c
c         N_zN_xE_x     +                        (eta-N_x^2)E_z =0
c
c-------------------------------------------------------------------
         
          if (n_perp2_p(n).ge.0) then
             n_z=cnpar_loc
             n_x=dsqrt(n_perp2_p(n))

             det_xy=(eps-n_z**2)*(eps-n_x**2-n_z**2)-g**2
             det_yz=(eps-n_x**2-n_z**2)*(eta-n_x**2)    
             det_xz=(eps-n_z**2)*(eta-n_x**2)-n_x**2*n_z**2    
            
             det_max=dabs(det_yz)
             j=1

c             write(*,*)'n_z,n_x,det_xy,det_yz,det_xz,det_max',
c     &                  n_z,n_x,det_xy,det_yz,det_xz,det_max


             if(dabs(det_yz).gt.det_max) then
                det_max=dabs(det_yz)
                j=2
             endif
             if(dabs(det_xz).gt.det_max) then
                det_max=dabs(det_xz)
                j=3
             endif

             if (det_max.gt.0) then
                n_rung=2
             else
                n_rung=1
             endif  

c             write(*,*)'j,n_rung',j,n_rung
             if(j.eq.1) then
               e_z=1.d0
               det_x=-n_z*n_x*(eps-n_x**2-n_z**2)
               det_y= im_one*g*n_z*n_x
               e_x=det_x/det_xy
               e_y=det_y/det_xy
             endif

             if(j.eq.2) then
               e_x=1.d0
               det_y=im_one*g*(eta-n_x**2)
               det_z=-(eps-n_x**2-n_z**2)*n_x*n_z 
               e_y=det_y/det_yz
               e_z=det_z/det_yz
             endif


             if(j.eq.3) then
               e_y=1.d0
               det_x=im_one*g*(eta-n_x**2)
               det_z=im_one*g*n_x*n_z 
               e_x=det_x/det_xz
               e_z=det_z/det_xz
             endif
             emod2=dreal(e_x*dconjg(e_x)+e_y*dconjg(e_y)+
     &                   e_z*dconjg(e_z))
             emod=dsqrt(emod2)
             Ex_p(n)=e_x/emod
             Ey_p(n)=e_y/emod
             Ez_p(n)=e_z/emod
             E_long_p(n)=(Ez_p(n)*cnpar_loc+
     &                   Ex_p(n)*n_x)/
     &                   dsqrt(cnpar_loc**2+n_perp2_p(n))
          else
             Ex_p(n)=(0.0,0.0)
             Ey_p(n)=(0.0,0.0)
             Ez_p(n)=(0.0,0.0)
             E_long_p(n)=(0.0,0.0)
          endif
c           write(*,*)'Ex N_p', Ex_p(n)
c          write(*,*)'Ey N_p', Ey_p(n)
c          write(*,*)'Ez N_p', Ez_p(n)
c          write(*,*)'E_long N_p', E_long_p(n)
c-----------------------------------------      
c          write(*,*)'n_perp2_m(n)',n_perp2_m(n)

          if (n_perp2_m(n).ge.0) then
             n_z=cnpar_loc
             n_x=dsqrt(n_perp2_m(n))

             det_xy=(eps-n_z**2)*(eps-n_x**2-n_z**2)-g**2
             det_yz=(eps-n_x**2-n_z**2)*(eta-n_x**2)    
             det_xz=(eps-n_z**2)*(eta-n_x**2)-n_x**2*n_z**2    
            
             det_max=dabs(det_yz)
           
             j=1
             if(dabs(det_yz).gt.det_max) then
                det_max=dabs(det_yz)
                j=2
             endif
             if(dabs(det_xz).gt.det_max) then
                det_max=dabs(det_xz)
                j=3
             endif

             if (det_max.gt.0) then
                n_rung=2
             else
                n_rung=1
             endif  

             if(j.eq.1) then
               e_z=1.d0
               det_x=-n_z*n_x*(eps-n_x**2-n_z**2)
               det_y= im_one*g*n_z*n_x
               e_x=det_x/det_xy
               e_y=det_y/det_xy
c               write(*,*)'j,det_xy,e_z,e_x,e_y',j,det_xy,e_z,e_x,e_y
             endif

             if(j.eq.2) then
               e_x=1.d0
               det_y=im_one*g*(eta-n_x**2)
               det_z=-(eps-n_x**2-n_z**2)*n_x*n_z 
               e_y=det_y/det_yz
               e_z=det_z/det_yz
c               write(*,*)'j,det_yz,e_z,e_x,e_y',j,det_yz,e_z,e_x,e_y
             endif


             if(j.eq.3) then
               e_y=1.d0
               det_x=im_one*g*(eta-n_x**2)
               det_z=im_one*g*n_x*n_z 
               e_x=det_x/det_xz
               e_z=det_z/det_xz
c               write(*,*)'j,det_xz,e_z,e_x,e_y',j,det_xz,e_z,e_x,e_y
             endif
             emod2=dreal(e_x*dconjg(e_x)+e_y*dconjg(e_y)+
     &                   e_z*dconjg(e_z))
             emod=dsqrt(emod2)
             Ex_m(n)=e_x/emod
             Ey_m(n)=e_y/emod
             Ez_m(n)=e_z/emod
             E_long_m(n)=(Ez_m(n)*cnpar_loc+Ex_m(n)*n_x)/
     &                   dsqrt(cnpar_loc**2+n_perp2_m(n))

          else
             Ex_m(n)=(0.0,0.0)
             Ey_m(n)=(0.0,0.0)
             Ez_m(n)=(0.0,0.0)
             E_long_m(n)=(0.0,0.0)
          endif
c          write(*,*)'Ex N_m', Ex_m(n)
c          write(*,*)'Ey N_m', Ey_m(n)
c          write(*,*)'Ez N_m', Ez_m(n)
c          write(*,*)'E_long N_m', E_long_m(n)
c---------------------------------------------------     
      enddo !n

      name_param(1)='iray'
      name_param(2)='N parallel'
      param(1)=irayl+1
      param(2)=cnpar

      do i=1,nbulk
       	 v(i)=v_loc(i)
	 w(i)=w_loc(i)
      enddo !i   

c----plots for n_perp2_p ---------------------------------------------
      call plot1dt_param(freq_ar,n_perp2_p,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0',
     &'freq GHZ','N_perp**2_p',n_param,name_param,
     &param)

c      write(*,*)'n_perp2_p irayl'

      do n=1,n_plot_freq
         temp(n)=dreal(Ex_p(n))
c          write(*,*)'n,Ex_p(n),temp(n)',n,Ex_p(n),temp(n)
      enddo  
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ReEx_p',n_param,name_param,
     &param)
 

      do n=1,n_plot_freq
         temp(n)=dimag(Ex_p(n))
      enddo
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ImEx_p',n_param,name_param,
     &param)
      
      do n=1,n_plot_freq
         temp(n)=dreal(Ey_p(n))
      enddo
c
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ReEy_p',n_param,name_param,
     &param)

      do n=1,n_plot_freq
         temp(n)=dimag(Ey_p(n))
      enddo 
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c    & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ImEy_p',n_param,name_param,
     &param)
    
      do n=1,n_plot_freq
         temp(n)=dreal(Ez_p(n))
      enddo
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ReEz_p',n_param,name_param,
     &param)
    
      do n=1,n_plot_freq
         temp(n)=dimag(Ez_p(n))
      enddo 
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ImEz_p',n_param,name_param,
     &param)
      
      do n=1,n_plot_freq
         temp(n)=dreal( E_long_p(n))
      enddo 
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','Re E_long_p at N_perp_p',n_param,name_param,
     &param)
      

      do n=1,n_plot_freq
         temp(n)=dimag( E_long_p(n))
      enddo
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_p=a2+dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','Im E_long_p at N_perp_p',n_param,name_param,
     &param)
   
 10   continue !end plots for N_perp2_p

c----plots for n_perp2_m ---------------------------------------------

      call plot1dt_param(freq_ar,n_perp2_m,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','N_perp**2_m',n_param,name_param,
     &param)

      do n=1,n_plot_freq
         temp(n)=dreal(Ex_m(n))
      enddo  
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ReEx_m',n_param,name_param,
     &param)
 

      do n=1,n_plot_freq
         temp(n)=dimag(Ex_m(n))
      enddo
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2+-sqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ImEx_m',n_param,name_param,
     &param)
      
      do n=1,n_plot_freq
         temp(n)=dreal(Ey_m(n))
      enddo
c
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ReEy_m',n_param,name_param,
     &param)

      do n=1,n_plot_freq
         temp(n)=dimag(Ey_m(n))
      enddo 
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ImEy_m',n_param,name_param,
     &param)
    
      do n=1,n_plot_freq
         temp(n)=dreal(Ez_m(n))
      enddo
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ReEz_m',n_param,name_param,
     &param)
    
      do n=1,n_plot_freq
         temp(n)=dimag(Ez_m(n))
      enddo 
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','ImEz_m',n_param,name_param,
     &param)
      
      do n=1,n_plot_freq
         temp(n)=dreal( E_long_m(n))
      enddo 
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','Re E_long_m at N_perp_m',n_param,name_param,
     &param)
      

      do n=1,n_plot_freq
         temp(n)=dimag( E_long_m(n))
      enddo
      call plot1dt_param(freq_ar,temp,0,0,n_plot_freq,1,
     &'linlin',0.d0,0.d0,
     &'N_perp^2_m=a2-dsqrt(det))/(2.d0*eps)',
c     & Disp. relation: eps*N_perp^4-a2*N_perp^2+a0=0'.
     &'freq GHZ','Im E_long_m at N_perp_m',n_param,name_param,
     &param)
   
 20   continue !end plots for N_perp2_m
c-----------------------------------------------------

      return
      end

   
   

 

c----------------subroutine MK_GRAPH----------------
c new version 110211
c Firstly plots tokamak cross-sections after that rays 
c
c Prepares for file drawgenr.in for xdraw              
c INPUT: from common blocks 'gr' and 'write'        
c        iray -number of the ray at antenna        
c        n_wall is a number of points in arrays
c                which set the wall form
c*****************************************************
cSAP090313
c      SUBROUTINE MK_GRAPH(iray,nray,ifreq,nfreq,nbulk)
      SUBROUTINE MK_GRAPH(iray,nray,ifreq)
      IMPLICIT double precision (a-h,o-z)
      include 'pgconst.i'
      INTEGER I,J

c-----input
      integer iray,nray,ifreq

      integer isave

cSm061209
      integer isave_freq

      save isave

cSm061209
      save isave_freq

      include 'param.i'

      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'gr.i'
      include 'write.i'
cSAP090313
      include 'one.i'
      include 'fourb.i'
      include 'three.i'

      data isave /0/
      data isave_freq/0/

      if(myrank.ne.0) return
      
      if (isave.lt.1) then
        isave=iray
      endif

cSm061229
      if(isave_freq.lt.1) then
        isave_freq=ifreq
      endif
      
c      write(*,*)'!!!!!!in mk_graph iray,nrayelt,isave,isave_freq',
c     + iray,nrayelt,isave,isave_freq

      if((iray.eq.1).or.(isave.eq.iray)) then
c         CALL ASSIGN("assign -F f77 -N ieee u:82",ier)
cBH111031         write(*,*)'mk_grph before open 82 and 81'
         write(*,*)'mk_grph before open 82'
         OPEN(82,file='genray.bin',form='unformatted')
cBH111031         open(81,file='genray.doc')
         open(85,file='absorp.bin',form='unformatted')
         open(86,file='absorp.doc')
         open(75,file='efield.bin',form='unformatted')
         open(76,file='eps_r.bin',form='unformatted')
         open(77,file='eps_i.bin',form='unformatted')
         open(78,file='em_res.bin',form='unformatted')
         open(79,file='cd.bin',form='unformatted')
      end if
cSm030508
c      if (nrayelt.eq.0)then
c      write(*,*)'nrayelt,iray,isave,ifreq',nrayelt,iray,isave,ifreq
      if (nrayelt.lt.2)then
        if((iray.gt.isave).or.(ifreq.gt.1)) then
c          write(*,*)'before goto 70'
          goto 70  
        else
          !YuP[07-2017] added:
          nrayelt=1 !Set nrayelt to 1, to avoid out-of-bounds problem,
          ! in case of nrayelt=0
        endif
        
      endif
     
c      write(*,*)'iray,isave,ifreq,isave_freq',
c     &            iray,isave,ifreq,isave_freq

c      if (((iray.gt.isave).or.(ifreq.gt.1)).and.
      if (((iray.gt.isave).or.(ifreq.gt.1)).or.
     &    (ifreq.gt.isave_freq))then
c        data for noncentral antenna rays      
c         write(*,*)'before goto 60' 
         goto 60
      endif

c-------------------------------------------------------------------
c     data for the central antenna ray
c-------------------------------------------------------------------
 1    format(15(1pe10.3))
 7    format(7(1pe10.3))
 2    format(/)
      write(*,*)'data for central ray nrayelt',nrayelt

cSAP110211
c-----tokamak cross-section  
      DO 30 J=1,3 
       DO 20 I=1,NP+1
         WRITE(82) REAL(AR(J,I)),REAL(AZ(J,I)),REAL(XT(J,I)),
     1   REAL(YT(J,I)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

cBH111031        WRITE(81,1)AR(J,I),AZ(J,I),XT(J,I),
cBH111031     1   YT(J,I),
cBH111031     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
cBH111031     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
cBH111031     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt),
cBH111031     5   wal_emis(nrayelt),wj_emis(nrayelt)
20     CONTINUE
       WRITE(82)
cBH111031       WRITE(81,2)
30    CONTINUE

      DO 50 J=4,NL
       DO 40 I=1,NP+1
         WRITE(82) REAL(AR(J,I)),REAL(AZ(J,I)),
     1   REAL(XT(3,NP+1)),REAL(YT(3,NP+1)),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
    
cBH111031         WRITE(81,1)AR(J,I),AZ(J,I),
cBH111031     1   XT(3,NP+1),YT(3,NP+1),
cBH111031     2   ws(nrayelt),delpwr(nrayelt),rez(nrayelt),
cBH111031     3   spsi(nrayelt),wnpar(nrayelt),wnper(nrayelt),
cBH111031     4   salphal(nrayelt),wye(nrayelt),wyi(nrayelt),
cBH111031     5   wal_emis(nrayelt),wj_emis(nrayelt)
40     CONTINUE
       WRITE(82)
cBH111031       WRITE(81,2)
50    CONTINUE
     
cSAP090313 write wall coordinates

      if (n_wall.gt.0) then
         do i=1,n_wall  
         WRITE(82) 
     1    REAL(r_wall(i)*100.0),REAL(z_wall(i)*100.0),
     &    REAL(xarr(nrayelt)),REAL(yarr(nrayelt)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3    REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5    real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         enddo

         WRITE(82)

      endif

      if (n_wall.gt.1) then        
         do m=0,max_limiters          
            do i=1,n_wall_add(m)           
               WRITE(82) 
     &         REAL(r_wall_add(i,m)*100.0),
     &         REAL(z_wall_add(i,m)*100.0),
     &         REAL(xarr(nrayelt)),REAL(yarr(nrayelt)),
     2    REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4    REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5         real(wal_emis(nrayelt)),real(wj_emis(nrayelt))  
           enddo ! i=1,n_wall_add(m)       
 
           WRITE(82)
         enddo !m
      endif
         
cSAP090413     
      p=pi/180.d0
      do m=1,max_limiters
         r_min=1.d0
         do i=1,n_limiter(m)
            if(r_min.gt.r_limiter(i,m)) r_min=r_limiter(i,m)
         enddo
         r_max=rr(nxeqd)*100.d0
         r_min=r_min*100.d0
         write(*,*)'degree phi_limiter(1,m),phi_limiter(2,m)',
     &              phi_limiter(1,m),phi_limiter(2,m)
         write(*,*)'radian phi_limiter(1,m)*p,phi_limiter(2,m)*p',
     &              phi_limiter(1,m)*p,phi_limiter(2,m)*p

         x_lim_min=r_min*dcos(phi_limiter(1,m)*p)
         x_lim_max=r_max*dcos(phi_limiter(1,m)*p)
         y_lim_min=r_min*dsin(phi_limiter(1,m)*p)
         y_lim_max=r_max*dsin(phi_limiter(1,m)*p)
 
         write(*,*)'mk_graph r_min,r_max',r_min,r_max
         write(*,*)'x_lim_min,y_lim_min',x_lim_min,y_lim_min
         write(*,*)'x_lim_max,y_lim_max',x_lim_max,y_lim_max

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_min),REAL(y_lim_min),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_max),REAL(y_lim_max),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         WRITE(82)

         x_lim_min=r_min*dcos(phi_limiter(2,m)*p)
         x_lim_max=r_max*dcos(phi_limiter(2,m)*p)
         y_lim_min=r_min*dsin(phi_limiter(2,m)*p)
         y_lim_max=r_max*dsin(phi_limiter(2,m)*p)

         write(*,*)'x_lim_min,y_lim_min',x_lim_min,y_lim_min
         write(*,*)'x_lim_max,y_lim_max',x_lim_max,y_lim_max

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_min),REAL(y_lim_min),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))

         WRITE(82) REAL(AR(1,1)),REAL(AZ(1,1)),
     1   REAL(x_lim_max),REAL(y_lim_max),
     2   REAL(ws(nrayelt)),REAL(delpwr(nrayelt)),REAL(rez(nrayelt)),
     3   REAL(spsi(nrayelt)),REAL(wnpar(nrayelt)),REAL(wnper(nrayelt)),
     4   REAL(salphal(nrayelt)),real(wye(nrayelt)),real(wyi(nrayelt)),
     5   real(wal_emis(nrayelt)),real(wj_emis(nrayelt))
         WRITE(82)
      enddo !m
c-------------------------------------------------------------------          
c     end data for tokamak cross - section
c---------------------------------------------------------------------
      if(nrayelt.eq.0) nrayelt=1  ! YuP[07-2017] added: 
         ! could happen if the root Nper at the launch was not found.
         ! Set nrayelt to 1, to avoid out-of-bounds problem.

      DO 10 I=1,nrayelt
cSm070720
        if (nrayelt.eq.1) goto 11

c        write(*,*)'2before 82 I',I
c        write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),wye(I),wyi(I)'
c        write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
c     3  REAL(wyi(i))

c        write(*,*)'I,wal_emis(I)',I,wal_emis(I)
        WRITE(82) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),REAL(wye(I)),
     3  REAL(wyi(I)),real(wal_emis(I)),real(wj_emis(I))

cBH111031        WRITE(81,1)wr(I),wz(I),xarr(I),yarr(I),
cBH111031     1  ws(I),delpwr(I),rez(I),spsi(I),
cBH111031     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I),
cBH111031     3  wal_emis(I),wj_emis(I)
         
c       for ioin ion resonance for nbulk >2

c        write(*,*)'nbulk,I,wyi(I),wxi2(I),wyi2(I),wxi(I)',
c     &             nbulk,I,wyi(I),wxi2(I),wyi2(I),wxi(I)

c        if (nbulk.gt.2) then
        if (nbulk.gt.3) then
           p12=wyi(I)*wxi2(I)/(wyi2(I)*wxi(I))

           write(*,*)'p12',p12

           del=1.d0/(1.d0+p12)
           yii=wyi(I)*wyi2(I)*((1.d0-del)*wyi(I)+del*(wyi2(I)))/
     .          ((1.d0-del)*wyi2(I)+del*(wyi(I)))   
        else
           p12=0.d0
           del=1.d0
           yii=0.d0
        endif

c        write(*,*)'mk_graph bef 85'
                       
        WRITE(85)REAL(ws(I)),REAL(delpwr(I)),REAL(spsi(I)),REAL(wye(I)),
     3  REAL(wyi(I)),REAL(wyi2(I)),REAL(dsqrt(yii))

c        write(*,*)'mk_graph bef 86'

        WRITE(86,7)ws(I),delpwr(I),spsi(I),wye(I),
     3  wyi(I),wyi2(I),dsqrt(dabs(yii))

c        write(*,*)'mk_graph bef 75'

        write(75)real(ws(I)),
     +   real(dreal(cwexde(I))),real(dimag(cwexde(I))),
     +   real(dreal(cweyde(I))),real(dimag(cweyde(I))),
     +   real(dreal(cwezde(I))),real(dimag(cwezde(I))),
     +   real(fluxn(I)),REAL(spsi(I))

c        write(*,*)'mk_graph bef 76'

        write(76)real(ws(I)),
     &  real(dreal(w_ceps(1,1,I))),real(dreal(w_ceps(1,2,I))),
     &  real(dreal(w_ceps(1,3,I))),
     &  real(dreal(w_ceps(2,1,I))),real(dreal(w_ceps(2,2,I))),
     &  real(dreal(w_ceps(2,3,I))),
     &  real(dreal(w_ceps(3,1,I))),real(dreal(w_ceps(3,2,I))),
     &  real(dreal(w_ceps(3,3,I)))


c        write(*,*)'mk_graph bef 77'

        write(77)real(ws(I)),
     &  real(dimag(w_ceps(1,1,I))),real(dimag(w_ceps(1,2,I))),
     &  real(dimag(w_ceps(1,3,I))),
     &  real(dimag(w_ceps(2,1,I))),real(dimag(w_ceps(2,2,I))),
     &  real(dimag(w_ceps(2,3,I))),
     &  real(dimag(w_ceps(3,1,I))),real(dimag(w_ceps(3,2,I))),
     &  real(dimag(w_ceps(3,3,I)))

c        write(*,*)'mk_graph bef 78'

        WRITE(78) REAL(wr(I)),REAL(wye(I)),
     &  real(wp_perpmax_dmvt(I,2)),
     &  real(wp_parmin_dmvt(I,2)),real(wp_parmax_dmvt(I,2)),
     &  real(wp_0_dmvt(I,2)),real(wdel_p_dmvt(I,2)),
     &  REAL(wnpar(I))

c      write(*,*)'mk_graph i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)',
c     &i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)

c       write(*,*)'mk_graph bef 79'

       WRITE(79) REAL(ws(I)),REAL(wr(I)),REAL(spsi(I)),REAL(wye(I)),
     &  real(eff(I)),REAL(wnpar(I)),real(delpow_e_ar(I)),
     &  real(delcur_par_ar(I)),real(wtheta_pol(I))
     

10    CONTINUE
cyup      write(*,*)'after 10'
      WRITE(82)
      write(85)
      write(75)
      write(78)
      write(79)
cBH111031      WRITE(81,2)

cSm070720
 11   continue
 
c-------------------------------------------------------------------          
      goto 70
c  end data for the central ray
c---------------------------------------------------------------------

60    continue
c      write(*,*)'!!!!!!in mk_graph 60 iray,nrayelt',iray,nrayelt
c--------------------------------------------------------------
c        data for noncentral antenna rays
c--------------------------------------------------------------
cSm070729
      if(nrayelt.eq.0) nrayelt=1  ! YuP[07-2017] added: 
         ! could happen if the root Nper at the launch was not found.
         ! Set nrayelt to 1, to avoid out-of-bounds problem.

      if (nrayelt.eq.1) goto 62

      DO 61 I=1,nrayelt

c       write(*,*)'befor 82 i',i
c       write(*,*)'REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I))'
c       write(*,*)REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
c     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
c     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I))


       WRITE(82) REAL(wr(I)),REAL(wz(I)),REAL(xarr(I)),REAL(yarr(I)),
     1  REAL(ws(I)),REAL(delpwr(I)),REAL(rez(I)),REAL(spsi(I)),
     2  REAL(wnpar(I)),REAL(wnper(I)),REAL(salphal(I)),real(wye(I)),
     3  REAL(wyi(I)),real(wal_emis(I)),real(wj_emis(I))

cBH111031       WRITE(81,1)wr(I),wz(I),xarr(I),yarr(I),
cBH111031     1  ws(I),delpwr(I),rez(I),spsi(I),
cBH111031     2  wnpar(I),wnper(I),salphal(I),wye(I),wyi(I),
cBH111031     3  wal_emis(I),wj_emis(I)

        write(75)real(ws(I)),
     +   real(dreal(cwexde(I))),real(dimag(cwexde(I))),
     +   real(dreal(cweyde(I))),real(dimag(cweyde(I))),
     +   real(dreal(cwezde(I))),real(dimag(cwezde(I))),
     +   real(fluxn(I)),REAL(spsi(I))

        WRITE(78) REAL(wr(I)),REAL(wye(I)),
     &  real(wp_perpmax_dmvt(I,2)),
     &  real(wp_parmin_dmvt(I,2)),real(wp_parmax_dmvt(I,2)),
     &  real(wp_0_dmvt(I,2)),real(wdel_p_dmvt(I,2)),
     &  REAL(wnpar(I))

c      write(*,*)'mk_graph i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)',
c     &i,wp_0_dmvt(I,2),wdel_p_dmvt(I,2)

       WRITE(79) REAL(ws(I)),REAL(wr(I)),REAL(spsi(I)),REAL(wye(I)),
     &  real(eff(I)),REAL(wnpar(I)),real(delpow_e_ar(I)),
     &  real(delcur_par_ar(I)),real(wtheta_pol(I))
61    CONTINUE
c      write(*,*)'after 61'
      WRITE(82)
cBH111031      WRITE(81,2)
      write(75)
      WRITE(78)
      write(79) 

cSm070720 
 62   CONTINUE

c  end data for the noncentral rays
c---------------------------------------------------------------------

70    continue
c      if (iray.eq.nray) then
       
c      write(*,*)'mk_graph iray,nraym,ifreq,nfreq',
c     & iray,nray,ifreq,nfreq
      if ((iray.eq.nray).and.(ifreq.eq.nfreq)) then
c       write(*,*)'mk_graph before close 82'
       CLOSE(82)
cBH111031       CLOSE(81)
       close(85)
       close(86)
       close(75)
       close(76)
       close(77)
       close(78)
       close(79)
c       write(*,*)'mk_graph after close 82 and 81'
      endif
      RETURN
      END

      subroutine plot_log_fe_param(irho)
c-----plot  log_e [fe(v_par**2*sign(v_par))] 
c     at given small radiual point irho
      implicit none

      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c      include 'onetwo.i'   
c-----input
      integer irho !number of radial point 1=<irho<NR-1

c-----locals
      real*4
     &Lnfe_min,Lnfe_max,
     &v_par_dt_min,v_par_dt_max,  
     &sv_par_dt_min,sv_par_dt_max,sign
    

      integer nv
      integer maxsteps
      real*4, pointer :: newx(:),newy(:)
      real*4  zero
      integer istat

      integer xtitlelen,titlelen,PGOPEN
      CHARACTER filename*40, line*1000, xtitle*100, titlestr*200
C Color index:
      INTEGER BLACK, WHITE, RED, GREEN, BLUE, CYAN, MAGENT, YELLOW
      PARAMETER (BLACK=0)
      PARAMETER (WHITE=1)
      PARAMETER (RED=2)
      PARAMETER (GREEN=3)
      PARAMETER (BLUE=4)
      PARAMETER (CYAN=5)
      PARAMETER (MAGENT=6)
      PARAMETER (YELLOW=7)
C Line style:
      INTEGER FULL, DASHED, DOTDSH, DOTTED, FANCY
      PARAMETER (FULL=1)
      PARAMETER (DASHED=2)
      PARAMETER (DOTDSH=3)
      PARAMETER (DOTTED=4)
      PARAMETER (FANCY=5)
C Character font:
      INTEGER NORMAL, ROMAN, ITALIC, SCRIPT
      PARAMETER (NORMAL=1)
      PARAMETER (ROMAN=2)
      PARAMETER (ITALIC=3)
      PARAMETER (SCRIPT=4)
C Fill-area style:
      INTEGER SOLID, HOLLOW
      PARAMETER (SOLID=1)
      PARAMETER (HOLLOW=2)

      character*4 name_param 
      REAL*4 RILIN
      character*72 text_param      
c      write(*,*)'in subroutine plot_log_fe irho=',irho

      if(myrank.ne.0) return
      
      zero=0.
c-----plot to plot.ps file

c-----------------------------------------------------------------
c     calculate maximal power and maximal s_poloidal along rays
c------------------------------------------------------------------
      Lnfe_min=0.d0 
      Lnfe_max=0.d0

      v_par_dt_min=0.
      v_par_dt_max=0.
                   
      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh

c         write(*,*)'nv,fe_lsc(irho,nv)',nv,fe_lsc(irho,nv)
c         write(*,*)'v_par_mesh_lsc(nv)',v_par_mesh_lsc(nv)
         if(dlog(abs(fe_lsc(irho,nv))).gt. Lnfe_max)  
     &      Lnfe_max=dlog(abs(fe_lsc(irho,nv)))

         if(dlog(abs(fe_lsc(irho,nv))).lt. Lnfe_min)  
     &      Lnfe_min=dlog(abs(fe_lsc(irho,nv)))

c        write(*,*)'Lnfe_max,Lnfe_min',Lnfe_max,Lnfe_min

       if(v_par_mesh_lsc(nv).gt.v_par_dt_max) 
     &     v_par_dt_max=v_par_mesh_lsc(nv)

       if(v_par_mesh_lsc(nv).lt.v_par_dt_min) 
     &     v_par_dt_min=v_par_mesh_lsc(nv)
 
      enddo

      if(v_par_dt_min.ge.0.d0) then 
        sign=1.d0
      else
        sign=-1.d0
      endif
      sv_par_dt_min=v_par_dt_min**2*sign

      if(v_par_dt_max.ge.0.d0) then 
        sign=1.d0
      else
        sign=-1.d0
      endif
      sv_par_dt_max=v_par_dt_max**2*sign

c      write(*,*)' plot_log_fe Lnfe_min,Lnfe_max',
c     &                         Lnfe_min,Lnfe_max

     
C***** Plot  dimensionless scale
c      write(*,*)'in plot_log_fe 1  v_par_dt_min,v_par_dt_max',
c     &           v_par_dt_min,v_par_dt_max
c      write(*,*)'in plot_log_fe 1  sv_par_dt_min,sv_par_dt_max',
c     &           sv_par_dt_min,sv_par_dt_max
c      write(*,*)'plot_log_fe  Lnfe_min Lnfe_max',Lnfe_min,Lnfe_max

    
      CALL PGSLW(4)
      CALL PGSCH(R41P5)
c      CALL PGENV(v_par_dt_min,v_par_dt_max,
c    &           Lnfe_min,Lnfe_max,0,1)
      CALL PGENV(sv_par_dt_min,sv_par_dt_max,
     &           Lnfe_min,Lnfe_max,0,1)

c      write(*,*)' in plot_log_fe 2'
  
c     print *,'Past PGENV '
C***  Sets limit for plot.
      xtitle='(v_par/vt)**2*sign(v_par)'
      xtitlelen=LEN_TRIM(xtitle)
      titlestr=''
      titlelen=LEN_TRIM(titlestr)

c      write(*,*)' in plot_log_fe 3'
  
      CALL PGLAB(xtitle(1:xtitlelen),'log_e[fe]'
     &  ,titlestr(1:titlelen))
        

       
c      write(*,*)' in plot_log_fe 4'
  
 
C***  note on screen, BLACK=0 means background color
C***  WHITE=1 means default foreground color
      CALL PGSLS(1)
      CALL PGSCI(1)
      CALL PGSLW(6) 
c      CALL PGLINE(nsteps_freq,rs,fces)

C***  plot d_ql 
         
c      write(*,*)' in plot_log_fe 5  bef allocate newy'
  
      allocate(newy(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
 
c      write(*,*)'in plot_log_fel 5_1 istat',istat
     
      call r4bcast(newy,zero,SIZE(newy))
        
c      write(*,*)' in plot_log_fe 6'
  
      allocate(newx(-nv_lsc_ql_lh:nv_lsc_ql_lh),STAT=istat)
      call r4bcast(newx,zero,SIZE(newx))

c      write(*,*)'in plot_log_fe aft newx istat',istat 

c      write(*,*)' in plot_log_fe 7'
  
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=log(abs(fe_lsc(irho,nv)))
c         newx(nv)=v_par_mesh_lsc(nv)
      
         if(v_par_mesh_lsc(nv).ge.0.d0) then 
           sign=1.d0
         else
           sign=-1.d0
         endif
         newx(nv)=v_par_mesh_lsc(nv)**2*sign

      ENDDO

c      write(*,*)' in plot_log_fe 8'
  
      CALL PGSLS(2)
      CALL PGSCI(BLUE)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)
c      print *,'BLUE: new plot_log_fe'

c      write(*,*)' in plot_log_fe 9'
  
   
      DO nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv)=log(abs(fe_lsc_old(irho,nv)))
c         write(*,*)'log_fe_old irho,nv,newy(nv)',irho,nv,newy(nv)       
      ENDDO

c      write(*,*)' in plot_log_fe 10'
  
      CALL PGSLS(3)
      CALL PGSCI(RED)
      CALL PGLINE(2*nv_lsc_ql_lh+1,newx,newy)

c      write(*,*)' in plot_log_fe 11'
  
c      print *,'RED: old log_fe'

    
      CALL PGSLS(1)
      CALL PGSCI(GREEN)
c      newx(1)=v_par_dt_min
c      newx(2)=v_par_dt_min
      newx(1)=sv_par_dt_min
      newx(2)=sv_par_dt_min
      newy(1)=Lnfe_min
      newy(2)=Lnfe_max
      CALL PGLINE(2,newx,newy)


      deallocate (newx,STAT=istat)
      deallocate (newy,STAT=istat)

 570  format(1X,A,1pe10.3)
      rilin=8.
      rilin=1.
      name_param='irho'
      write(text_param,570)name_param,real(irho)

      CALL PGMTXT('B',RILIN,-R4P1,R40,text_param)  

c      write(*,*)'end plot_log_fe'

      return
      end

      
      subroutine plot_log_fe_vpar2_param(irho)
      implicit none 
      include 'pgconst.i'
      include 'param.i'
      
      integer myrank !In serial run: myrank=0; In MPI run: myrank=rank
      common/mpimyrank/myrank !In serial run: myrank=0; In MPI run: myrank=rank
      
      include 'grill.i'
      include 'one.i'                 
      include 'lsc_approach.i'
      include 'write.i'
      include 'writencdf.i'
c-----input
      integer irho !number of radial point 1=<irho<NR-1
c-----locals
      integer n_param                    !number of parameters
      parameter  (n_param=2)
      
      character*20 name_param(n_param)  !names  of parameters
      REAL*4 param(n_param)                      !values of parameters  

      real*8, pointer :: newx(:),newy(:) 
      real*8  zero,sign
      integer istat
      integer nv
 
      if(myrank.ne.0) return
      
      zero=0.
      allocate(newy(1:2*nv_lsc_ql_lh+1),STAT=istat)
c      call r4bcast(newy,zero,SIZE(newy))
      call bcast(newy,zero,SIZE(newy))
      allocate(newx(1:2*nv_lsc_ql_lh+1),STAT=istat)
c      call r4bcast(newx,zero,SIZE(newx))
      call bcast(newy,zero,SIZE(newy))
      do nv=-nv_lsc_ql_lh,nv_lsc_ql_lh
         newy(nv+nv_lsc_ql_lh+1)=log(abs(fe_lsc(irho,nv)))

         if(v_par_mesh_lsc(nv).ge.0.d0) then 
           sign=1.d0
         else
           sign=-1.d0
         endif

         newx(nv+nv_lsc_ql_lh+1)=v_par_mesh_lsc(nv)**2*sign
      enddo

      name_param(1)='irho'
      name_param(2)='rho_bin_center_lsc'
      param(1)=irho
      param(2)=rho_bin_center_lsc(irho)
      call plot1dt_param(newx,newy,0,0,2*nv_lsc_ql_lh+1,1,
     &'linlin',0.d0,0.d0,
     &'ln fe',
     &'v_par**2*sign(vpar)','ln fe',n_param,name_param,
     &param)

      return
      end
