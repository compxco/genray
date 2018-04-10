

c   routine name   - mmdei
c
c-----------------------------------------------------------------------
c
c   computer            - vax/double
c
c   latest revision     - november 1, 1979
c
c   purpose             - exponential integrals
c
c   usage               - function mmdei (iopt,arg,ier)
c
c   arguments    mmdei  - output value of the integral. mmdei must be
c                           typed appropriately in the calling program.
c                           (see precision/hardware section.)
c                iopt   - input option.
c                         for iopt = 1, the integral (from -infinity to
c                           arg) of exp(t)/t dt will be evaluated if arg
c                           is greater than 0. if arg is less than 0.0,
c                           (-1)*the integral (from -arg to infinity)
c                           of exp(-t)/t dt will be evaluated.
c                         for iopt = 2, the integral (from arg to
c                           infinity) of exp(-t)/t dt will be evaluated.
c                           arg must be greater than 0.
c                         for iopt = 3, exp(-arg)*the integral (from
c                           -infinity to arg) of exp(t)/t dt will be
c                           evaluated if arg is greater than 0. if arg
c                           is less than 0.0, exp(-arg)*(-1)*the
c                           integral (from -arg to infinity) of
c                           exp(-t)/t dt will be evaluated.
c                arg    - input parameter. see iopt description.
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129 indicates that iopt was less than
c                             1 or greater than 3. mmdei is set to
c                             machine infinity.
c                           ier = 130 indicates that arg was equal to
c                             0.0. mmdei is set to machine infinity if
c                             iopt = 2 and negative machine infinity if
c                             iopt = 1 or 3.
c                           ier = 131 indicates that an overflow would
c                             have occurred. if iopt = 1, mmdei is set
c                             to machine infinity.
c                           ier = 132 indicates that an underflow would
c                             have occurred. mmdei is set to 0.0 if
c                             iopt = 1 or 2.
c                         warning with fix
c                           ier = 69 indicates that arg was negative
c                             for iopt = 2. calculation continues using
c                             abs(arg).
c
c   precision/hardware  - double/h32,h36
c                       - single/h48,h60
c
c   reqd. imsl routines - uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c-----------------------------------------------------------------------
c
      double precision function mmdei (iopt,arg,ier)
c                                  specifications for arguments
      implicit none

      integer            iopt,ier
      double precision   arg
c                                  specifications for local variables
      integer i
      integer iend,jend,kend,iendm1,iendp1,kendm1,j,jendp1
	character*6 name
      double precision   a(6),b(6),c(8),d(8),e(8),f(8),p1(9),q1(9),
     *                   p2(9),q2(8),p3(10),q3(9),p4(10),q4(9),p0(6),
     *                   q0(6),px(9),qx(9),frac,sump,sumq,t,w,x,x0,
     *                   xx0,x01,x02,xmx0,y,dexp40,xinf,xmax,xmin,
     *                   six,twelve,two,three,zero,one,half,twent4,
     *                   four,forty
      double precision   sumb,sumc,sumd,sume,sumf
      data               dexp40/.2353852668370200d18/
      data               xinf/1.700000000d+38/ 
      data               xmax/88.02969193111305d0/
      data               xmin/-.1750475447540571d03/
      data               x0/.3725074107813666d00/
      data               x01/.3725074107805994d00/
      data               x02/.7671772501993940d-12/
      data               a(1)/  -.5772156649015328d 00/
      data               a(2)/   .7541643136630163d 00/
      data               a(3)/   .1298492329273731d 00/
      data               a(4)/   .2406813556839774d-01/
      data               a(5)/   .1320843092096093d-02/
      data               a(6)/   .6577393997532639d-04/
      data               b(1)/   .1000000000000000d 01/
      data               b(2)/   .4258991938115897d 00/
      data               b(3)/   .7977947184102281d-01/
      data               b(4)/   .8302084760987714d-02/
      data               b(5)/   .4864271383930161d-03/
      data               b(6)/   .1306551958228487d-04/
      data               c(1)/   .8677459548384432d-07/
      data               c(2)/   .9999955193013902d 00/
      data               c(3)/   .1184831055549458d 02/
      data               c(4)/   .4559306442533897d 02/
      data               c(5)/   .6992794512910029d 02/
      data               c(6)/   .4252020347688406d 02/
      data               c(7)/   .8836718088038437d 01/
      data               c(8)/   .4013776649406646d 00/
      data               d(1)/   .1000000000000000d 01/
      data               d(2)/   .1284819353791566d 02/
      data               d(3)/   .5644335695618032d 02/
      data               d(4)/   .1066451837699138d 03/
      data               d(5)/   .8973110971252896d 02/
      data               d(6)/   .3149718491704406d 02/
      data               d(7)/   .3795590037621223d 01/
      data               d(8)/   .9088045691888691d-01/
      data               e(1)/  -.9999999999999731d 00/
      data               e(2)/  -.3440619950066848d 02/
      data               e(3)/  -.4275326712019885d 03/
      data               e(4)/  -.2396019432474904d 04/
      data               e(5)/  -.6168852100554763d 04/
      data               e(6)/  -.6576096987480218d 04/
      data               e(7)/  -.2106077371426331d 04/
      data               e(8)/  -.1489908499729481d 02/
      data               f(1)/   .1000000000000000d 01/
      data               f(2)/   .3640619950064596d 02/
      data               f(3)/   .4943450702099040d 03/
      data               f(4)/   .3190272374895431d 04/
      data               f(5)/   .1033707530858408d 05/
      data               f(6)/   .1632414535577834d 05/
      data               f(7)/   .1114977528710966d 05/
      data               f(8)/   .2378138991021601d 04/
      data               p0(1)/   .1525388359511120d 03/
      data               p0(2)/   .3402682862739600d 03/
      data               p0(3)/   .2597386446160079d 03/
      data               p0(4)/   .7787096586760712d 02/
      data               p0(5)/   .7460510544921461d 01/
      data               p0(6)/   .5574718225325585d-01/
      data               q0(1)/   .1525388359511120d 03/
      data               q0(2)/   .4165377042495159d 03/
      data               q0(3)/   .4171612180903951d 03/
      data               q0(4)/   .1857403824840772d 03/
      data               q0(5)/   .3490362129565328d 02/
      data               q0(6)/   .2000000000000000d 01/
      data               p1(1)/   .5531977362081977d 01/
      data               p1(2)/   .2063133336244559d 03/
      data               p1(3)/   .1427234409068234d 05/
      data               p1(4)/   .3685738950923286d 05/
      data               p1(5)/   .4493416458218790d 07/
      data               p1(6)/  -.1684497821007958d 07/
      data               p1(7)/   .3529608047950282d 09/
      data               p1(8)/  -.1251949974431755d 09/
      data               p1(9)/   .2997849734461850d 10/
      data               q1(1)/   .2562890625000000d 02/
      data               q1(2)/  -.1512618411191135d 04/
      data               q1(3)/   .4268855000903744d 05/
      data               q1(4)/  -.7478772860127960d 06/
      data               q1(5)/   .8857915400539992d 07/
      data               q1(6)/  -.7249035719651191d 08/
      data               q1(7)/   .4014138914734781d 09/
      data               q1(8)/  -.1398366755614922d 10/
      data               q1(9)/   .1279632488038080d 10/
      data               p2(1)/  -.2469409834483613d 01/
      data               p2(2)/  -.3677831134783113d 02/
      data               p2(3)/   .2327302338390390d 02/
      data               p2(4)/   .7894722092944569d 01/
      data               p2(5)/  -.1941329675144305d 02/
      data               p2(6)/   .5886582407532809d 01/
      data               p2(7)/   .4181024225628565d 01/
      data               p2(8)/   .5731167057445080d 01/
      data               p2(9)/   .9989576665165515d 00/
      data               q2(1)/   .2639830073180245d 01/
      data               q2(2)/   .9654052174292799d 03/
      data               q2(3)/  -.8387670841896405d 01/
      data               q2(4)/   .3172794892543692d 03/
      data               q2(5)/   .5231655687345586d 02/
      data               q2(6)/   .3413652125243753d 03/
      data               q2(7)/  -.1991496002312352d 03/
      data               q2(8)/   .1146252532490162d 01/
      data               p3(1)/  -.1647721172463462d 01/
      data               p3(2)/  -.1860092121726437d 02/
      data               p3(3)/  -.1000641913989284d 02/
      data               p3(4)/  -.2105740799548040d 02/
      data               p3(5)/  -.9134835699998741d 00/
      data               p3(6)/  -.3323612579343960d 02/
      data               p3(7)/   .2495487730402057d 02/
      data               p3(8)/   .2652575818452800d 02/
      data               p3(9)/  -.1845086232391277d 01/
      data               p3(10)/   .9999933106160568d 00/
      data               q3(1)/   .9792403599217288d 02/
      data               q3(2)/   .6403800405352414d 02/
      data               q3(3)/   .5994932325667409d 02/
      data               q3(4)/   .2538819315630708d 03/
      data               q3(5)/   .4429413178337928d 02/
      data               q3(6)/   .1192832423968600d 04/
      data               q3(7)/   .1991004470817741d 03/
      data               q3(8)/  -.1093556195391090d 02/
      data               q3(9)/   .1001533852045341d 01/
      data               p4(1)/   .1753388012654660d 03/
      data               p4(2)/  -.2231276707776324d 03/
      data               p4(3)/  -.1819496649298688d 02/
      data               p4(4)/  -.2797985286243052d 02/
      data               p4(5)/  -.7631477016202536d 01/
      data               p4(6)/  -.1528566236369296d 02/
      data               p4(7)/  -.7068109778950293d 01/
      data               p4(8)/  -.5000066404131309d 01/
      data               p4(9)/  -.3000000003209813d 01/
      data               p4(10)/   .1000000000000010d 01/
      data               q4(1)/   .3978459771674147d 05/
      data               q4(2)/   .3972771091004144d 01/
      data               q4(3)/   .1377903902357480d 03/
      data               q4(4)/   .1171792205020864d 03/
      data               q4(5)/   .7048318471804246d 02/
      data               q4(6)/  -.1201877635471546d 02/
      data               q4(7)/  -.7992435957763395d 01/
      data               q4(8)/  -.2999998940403249d 01/
      data               q4(9)/   .1999999999990480d 01/
      data               six/6.d0/,twelve/12.d0/,three/3.d0/,two/2.d0/,
     *                   one/1.d0/,half/.5d0/,twent4/24.d0/,four/4.d0/,
     *                   forty/40.d0/,zero/0.d0/
c                                  first executable statement

      iend = 8
      jend = 9
      kend = 6
      iendm1 = iend-1
      iendp1 = iend+1
      jendp1 = jend+1
      kendm1 = kend-1
      x = arg
      ier = 0
      if (iopt.lt.1.or.iopt.gt.3) go to 125
      go to (5,100,5), iopt
c                                  iopt = 1 or 3
    5 if (x) 60,115,10
   10 if (x.ge.twelve) go to 35
      if (x.ge.six) go to 25
c                                  x greater than or equal to 0 and
c                                  less than 6 - rational approximation
c                                  used is expressed in terms of
c                                  chebyshev polynomials to improve
c                                  conditioning
      t = x+x
      t = t/three-two
      px(1) = zero
      qx(1) = zero
      px(2) = p1(1)
      qx(2) = q1(1)
      do 15 i=2,iend
         px(i+1) = t*px(i)-px(i-1)+p1(i)
         qx(i+1) = t*qx(i)-qx(i-1)+q1(i)
   15 continue
      sump = half*t*px(iendp1)-px(iend)+p1(iendp1)
      sumq = half*t*qx(iendp1)-qx(iend)+q1(iendp1)
      frac = sump/sumq
      xmx0 = (x-x01)-x02
      if (dabs(xmx0).lt.0.037d0) go to 20
      xx0 = x/x0
      mmdei = dlog(xx0)+xmx0*frac
      if (iopt.eq.3) mmdei = dexp(-x)*mmdei
      go to 9005
c                                  evaluate approximation for ln(x/x0)
c                                  for x close to x0
   20 y = xmx0/x0
      sump = ((((p0(6)*y+p0(5))*y+p0(4))*y+p0(3))*y+p0(2))*y+p0(1)
      sumq = ((((q0(6)*y+q0(5))*y+q0(4))*y+q0(3))*y+q0(2))*y+q0(1)
      mmdei = (sump/(sumq*x0)+frac)*xmx0
      if (iopt.eq.3) mmdei = dexp(-x)*mmdei
      go to 9005
c                                  x greater than or equal to 6 and
c                                  less than 12
   25 frac = zero
      do 30 i=1,iend
         frac = q2(i)/(p2(i)+x+frac)
   30 continue
      mmdei = (p2(iendp1)+frac)/x
      if (iopt.ne.3) mmdei = mmdei*dexp(x)
      go to 9005
c                                  x greater than or equal to 12 and
c                                  less than 24
   35 if (x.ge.twent4) go to 45
      frac = zero
      do 40 i=1,jend
         frac = q3(i)/(p3(i)+x+frac)
   40 continue
      mmdei = (p3(jendp1)+frac)/x
      if (iopt.ne.3) mmdei = mmdei*dexp(x)
      go to 9005
c                                  x greater than or equal to 24
   45 if ((x.ge.xmax).and.(iopt.lt.3)) go to 110
      y = one/x
      frac = zero
      do 50 i=1,jend
         frac = q4(i)/(p4(i)+x+frac)
   50 continue
      frac = p4(jendp1)+frac
      mmdei = y+y*y*frac
      if (iopt.eq.3) go to 9005
      if (x.gt.170.0d0) go to 55
      mmdei = mmdei*dexp(x)
      go to 9005
c                                  calculation reformulated to avoid
c                                  premature overflow
   55 mmdei = (mmdei*dexp(x-forty))*dexp40
      go to 9005
c                                  original x was negative.
   60 y = -x
   65 w = one/y
      if (y.gt.four) go to 85
      if (y.gt.one) go to 75
c                                  -x greater than 0 and less than or
c                                  equal to 1
      sumb = b(kend)
      do 70 i=1,kendm1
         j = kend-i
         sumb = (sumb*y)+b(j)
   70 continue
      mmdei = dlog(y)-(((((a(6)*y+a(5))*y+a(4))*y+a(3))*y+a(2))*y+a(1))
     1/sumb
      if (iopt.eq.3) mmdei = mmdei*dexp(y)
      go to 95
c                                  -x greater than -1 and less than or
c                                  equal to 4
   75 sumc = c(iend)
      sumd = d(iend)
      do 80 i=1,iendm1
         j = iend-i
         sumc = (sumc*w)+c(j)
         sumd = (sumd*w)+d(j)
   80 continue
      mmdei = -sumc/sumd
      if (iopt.eq.3) go to 9005
      mmdei = mmdei*dexp(-y)
      go to 95
c                                  -x greater than 4
   85 if ((-dabs(x).lt.xmin).and.(iopt.lt.3)) go to 105
      sume = e(iend)
      sumf = f(iend)
      do 90 i=1,iendm1
         j = iend-i
         sume = (sume*w)+e(j)
         sumf = (sumf*w)+f(j)
   90 continue
      mmdei = -w*(1.0d0+w*sume/sumf)
      if (iopt.eq.3) go to 9005
      mmdei = mmdei*dexp(-y)
   95 if (iopt.eq.2) mmdei = -mmdei
      if (ier.eq.69) go to 9000
      go to 9005
  100 y = x
      if (y) 120,115,65
c                                  terminal error - arg is less than
c                                  xmin causing underflow
  105 mmdei = zero
      ier = 132
      go to 9000
c                                  terminal error - x is greater than
c                                  xmax causing overflow
  110 mmdei = xinf
      ier = 131
      go to 9000
c                                  terminal error - arg = 0
  115 mmdei = -xinf
      if (iopt.eq.2) mmdei = -mmdei
      ier = 130
      go to 9000
c                                  warning with fix - arg is less than
c                                  0.0 for iopt = 2
  120 ier = 69
      go to 60
c                                  terminal error - iopt is out of range
  125 mmdei = xinf
      ier = 129
 9000 continue

cSm990901
  	name='mmdei ' 
	call uertst(ier,name)
c     call uertst (ier,6hmmdei )


 9005 return
      end
