c     deck main
$	debug
      program scat2
c***********************************************************************
c                         *                   *                        *
c                         *  programme scat2  *                        *
c                         *                   *                        *
c                         *********************                        *
c                                                                      *
c     o.bersillon          juillet 1977                                *
c                  rev.1   juillet 1979                                *
c                                                                      *
c***********************************************************************
      parameter(noe=200,nol=40)
      real mi,mt
      dimension en(noe)
      dimension ami(6),azi(6),asi(6)
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/poten/a(6),pot(5,5),r(6),beta,ewmax
      common/tce/tc(noe,nol)
      common/xs/se(noe),sr(noe),st(noe)
      common/elb/elab(noe)
c
      data ami/1.008665,1.007825,2.014102,3.016050,3.016030,4.002603/
      data asi/0.5,0.5,1.0,0.5,0.5,0.0/
      data azi/0.,1.,1.,1.,2.,2./
c
      ie=5
      is=7
      is1=11
      is2=12
      is3=8
      open(unit=ie,file='scat.in',status='old')
      open(unit=is,file='scat.out',status='unknown')
      open(unit=14,file='sc.tr',status='unknown')
      open(unit=15,file='sc.cs',status='unknown')
      open(unit=16,file='sc.eps',status='unknown')
      open(unit=17,file='sc.dt',status='unknown')
      open(unit=18,file='sc.ass',status='unknown')
      open(unit=21,file='p.pc',status='unknown')
      open(unit=22,file='q.pc',status='unknown')
      open(unit=23,file='p.pv',status='unknown')
      open(unit=24,file='q.pv',status='unknown')
      open(unit=25,file='ele.out',status='unknown')
      open(unit=26,file='delta.out',status='unknown')
      open(unit=27,file='defasagem.out',status='unknown')
      open(unit=28,file='Re_delta.out',status='unknown')
      open(unit=29,file='Im_delta.out',status='unknown')
      open(unit=is1,status='scratch')
      open(unit=is2,status='scratch')
      open(unit=is3,file='coftr.dat',status='unknown')
      key=0
      pi=3.141592654
      do 7 n=1,noe
      do 7 j=1,nol
      tc(n,j)=0.0
    7 continue
c
      read(ie,1) ipr,ida,iba,ipu
c     ****
c     si ipr = 1 ecriture des t(l,j) for gnash or stapre code
c
c            = 0 pas d'ecriture des t(l,j)
c     si ida = 1 calcul de la distribution angulaire de la diffusion
c                elastique potentielle  ( angles equirepartis )
c            = 2 idem  ( cosinus des angles equirepartis )
c            = 0 pas de calcul
c     si iba = 1 ecriture sur tape07 des parametres, sections efficaces
c
c            = 0 pas d'ecriture
c     si ipu = 1 punch des resultats  (desativado 24/10/84)
c            = 0 pas de punch         (desativado 24/10/84)
c
      if(ida.ne.0) call fact
      if(ida.ne.0) call preang(ida,na)
  10  ne=0
      nel=1
 1000 continue
      read(ie,1) me
      ne=ne+me
c     ****
      read(ie,5) (en(j),j=nel,ne)
c     ****
c     me = nombre d'energies a calculer
c     en(j) = energies (mev)
c             si en(1) > 0. energies centre de masse
c                      < 0. energies laboratoire
      do 11 j=nel,ne
      elab(j)=abs(en(j))
   11 ein(j)=abs(en(j))
c
   20 continue
      read(ie,1) izt,imt
c     ****
c     izt = numero atomique de la cible
c     imt = nombre de masse de la cible
      zt=float(izt)
      mt=float(imt)
   30 read(ie,1) ip,ipot,key
c     ****
c     ip = 1 neutron         ipot = 1 --> wilmore - hodgson
c                                   2 --> beccheti - greenless
c                                   3 --> ferrer - carlson - rapaport
c                                   4 --> bersillon - cindro
c                                   5 --> madland (actinides)
c          2 proton          ipot = 1 --> perey
c                                   2 --> beccheti - greenless
c          3 deuteron        ipot = 1 --> lohr - haeberli
c                                   2 --> perey
c          4 triton          ipot = 1 --> beccheti - greenless
c          5 helium-3        ipot = 1 --> beccheti - greenless
c          6 alpha           ipot = 1 --> potentiel moyen
c
c     si ipot = 0 , les parametres sont lus
c
      mi=ami(iabs(ip))
      if(ip.lt.0) mi=float(ifix(mi))
      ip=iabs(ip)
      zi=azi(ip)
      si=asi(ip)
      ipl =ifix(2.0*si+1.001)
      amu=mt/(mi+mt)
c
c     transformation des energies laboratoire dans le centre de masse
c                                                               if any
      if(en(1).gt.0.0) go to 18
      do 17 j=nel,ne
   17 ein(j)=abs(en(j))*amu
   18 continue
c
      do 19 i=1,5
      a(i)=0.0
      r(i)=0.0
      do 19 j=1,5
      pot(i,j)=0.0
   19 continue
      a(6)=0.0
      r(6)=0.0
      beta=0.0
      ewmax=0.0
c
      if(ipot.eq.0) go to 21
      call syspot(ip,iabs(ipot))
      go to 22

c
   21 continue
      read(ie,3) r(1),a(1),(pot(1,i),i=1,5),beta
c     ****
c-----potentiel reel : woods - saxon
c     r(1) = rayon (fm)
c     a(1) = diffusivite (fm)
c     pot(1,i) = parametres de la profondeur du puits
c       v = pot(1,1) + pot(1,2) * e + pot(1,3) *e*e + pot(1,4) * ln(e) +
c           pot(1,5) * sqrt(e)
c     beta = nonlocality range
c            si beta.ne.0 le potentiel imaginaire est purement
c                                                     de surface (d-w-s)
c
c      read(ie,3) r(2),a(2),(pot(2,i),i=1,5),a(5)
       read(ie,300) r(2),a(2),(pot(2,i),i=1,5),a(5)
300    format(2f10.5,5e10.5,f6.3)
c     ****
c-----potentiel imaginaire de surface
c     r(2) = rayon (fm)    si > 0.  derivee de woods - saxon
c                          si < 0.  gaussien
c     a(2) = diffusivite (fm)
c            peut dependre de l'energie suivant la relation
c                  a(2) = a(2) + a(5) * e
c     pot(2,i) = parametres de la profondeur du puits
c      wd = pot(2,1) + pot(2,2) * e + pot(2,3) *e*e + pot(2,4) * ln(e) +
c           + pot(2,5) * sqrt(e)
c								    
      read(ie,3) r(3),a(3),(pot(3,i),i=1,5)
c     ****
c-----potentiel imaginaire volume
c     r(3) = rayon (fm)
c     a(3) = diffusivite (fm)
c     pot(3,i) = parametres de la profondeur du puits
c      wv = pot(3,1) + pot(3,2) * e + pot(3,3) *e*e + pot(3,4) * ln(e) +
c           + pot(3,5) * sqrt(e)
c
      read(ie,3) r(4),a(4),(pot(4,i),i=1,5)
c     ****
c-----potentiel spin - orbite
c     r(4) = rayon (fm)
c     a(4) = diffusivite (fm)
c     pot(4,i) = parametres de la profondeur du puits
c     vso = pot(4,1) + pot(4,2) * e + pot(4,3) *e*e + pot(4,4) * ln(e) +
c           + pot(4,5) * sqrt(e)
c
      read(ie,3) r(5),ewmax
c     ****
c-----r(5) = rayon coulombien
c     ewmax = energie a partir de laquelle les profondeurs imaginaires
c             sont constantes
c
   22 if(ipot.gt.0) go to 23
      read(ie,3) r(6),a(6),(pot(5,i),i=1,5)
c
   23 read(ie,1) isuit
c     ****
c     isuit = 0  sortie
c             1  nouveau cas complet
c             2  on conserve la grille en energie
c             3  on ne change que les potentiels
      isuit=isuit+1
c
      call pripot(iabs(ipot))
c
      do 100 n=nel,ne
      e1=ein(n)
      if(ipr.eq.1) write(is,4)
      if(ipr.eq.2) write(is,4)
c
      call scat(lmax,ipr)
           
c
      if(ipl.eq.1) call spin0 (n,lmax,ipr)
      if(ipl.eq.2) call spin05(n,lmax,ipr)
      if(ipl.eq.3) call spin1 (n,lmax,ipr)
c
c     if(ida.ne.0.and.ip.eq.1) call shapel(lmax,ida,na,ipl)
      if(ida.ne.0) call shapec(lmax,ida,na,ipl,ipu,elab(n))
c
      if(ipu.ne.0) write(is3,6) ein(n),se(n),sr(n),st(n)
c
c     if(ipr.eq.1) call tpun(ip,me,ipl,0,lmax)
  100 continue
      if(isuit.eq.5) then
                      nel=nel+me
                     go to 1000
                     endif
c
      call pritc(ne)
      if(ipr.ne.0) call pritd(ne,ip,iba)
c
c     if(ipr.eq.1) call tpun(ip,ne,ipl,1,lmax)
      do 24 n=1,noe
      do 24 j=1,nol
  24  tc(n,j)=0.0
c
      go to (9000,10,20,30),isuit
c
c     formats
c
    1 format(10i3)
    3 format(7f10.5,f6.3)
    4 format(1h1)
    5 format(6f10.5)
    6 format(1p6f10.5)
c
 9000 stop
      end
c     deck cgamma
      complex function cgamma(z)
c***********************************************************************
c                                                                      *
c     calculates complex gamma function                                *
c     z must be declared complex in the calling program                *
c                                                                      *
c     reference    y.l.luke                                            *
c                  the special functions and their approximations      *
c                  vol.2, academic press, new york and london          *
c                  (1969) p.304-305                                    *
c                                                                      *
c***********************************************************************
      complex   h,s,u,v,z,piu
      dimension g(16)
c
      double precision pi,g,con1
      data pi/3.141592653589793d+00/
      data g/
     1 41.624436916439068d+00,-51.224241022374774d+00,
     2 11.338755813488977d+00, -0.747732687772388d+00,
     3 +0.8782877493061d-02, -0.1899030264d-05,
     4 +0.1946335d-08, -0.199345d-09, +0.8433d-11,
     5 +0.1486d-11, -0.806d-12, +0.293d-12,
     6 -0.102d-12, +0.37d-13, -0.14d-13,
     7 +0.6d-14/
      data con1/2.506628274631001d+00/
      u=z
      x=real(u)
      if(x .ge. 1.0) go to 3
      if(x .ge. 0.0) go to 2
      v=1.0-u
      l=1
      go to 11
    2 v=u+1.0
      l=2
      go to 11
    3 v=u
      l=3
c
   11 h=1.0
      s=g(1)
      do 1 k=2,16
      fk=k-2
      fk1=fk+1.0
      h=((v-fk1)/(v+fk))*h
    1 s=s+g(k)*h
      h=v+4.5
      cgamma=con1*cexp((v-0.5)*clog(h)-h)*s
c
      go to (21,22,23),l
c
   21 piu=pi*u
      cgamma=pi/(csin(piu)*cgamma)
      return
c
   22 cgamma=cgamma/u
   23 return
c
      end
c     deck cleb
      function cleb(aj1,aj2,aj3,am1,am2,am3)
c***********************************************************************
c     calcul des coefficients de clebsch-gordan                        *
c                                                                      *
c     attention :  cg(j1,j2,j3,;m1,m2,m3) = (-1)**(j1+j2-m3)*          *
c                                           3-j(j1,j2,j3;m1,m2,-m3)    *
c                                           ---                        *
c                                                                      *
c     d'apres john.g.wills     ornl-tm-1949 (august 1967)              *
c                           et comp.phys.comm. 2(1971)381              *
c                                                                      *
c     o.bersillon     aout 1977                                        *
c                                                                      *
c***********************************************************************
      dimension i(11)
      common/facto/g(101)
      equivalence (i(1),i1),(i(2),i2),(i(3),i3),(i(4),i4),(i(5),i5),
     *(i(6),i6),(i(7),i7),(i(8),i8),(i(9),i9),(i(10),i10),(i(11),i11)
c
      is=6
      im=101
      cleb=0.0
c     convert the arguments to integer
      j1=int(2.0*aj1+0.001)
      j2=int(2.0*aj2+0.001)
      j3=int(2.0*aj3+0.001)
      m1=int(2.0*am1+sign(0.001,am1))
      m2=int(2.0*am2+sign(0.001,am2))
      m3=int(2.0*am3+sign(0.001,am3))
c     test m1 + m2 = m3
      if(m1+m2-m3)300,40,300
c     test table size
   40 i(10)=(j1+j2+j3)/2+2
      n=i(10)
      i(11)=j3+2
      if(i(10)-im)70,70,50
   50 write(is,60) i(10),im,aj1,aj2,aj3,am1,am2,am3
      go to 300
   60 format(1h ,11htable size ,2i5,6f5.1)
   70 i(1)=j1+j2-j3
      i(2)=j2+j3-j1
      i(3)=j3+j1-j2
      i(4)=j1-m1
      i(5)=j1+m1
      i(6)=j2-m2
      i(7)=j2+m2
      i(8)=j3-m3
      i(9)=j3+m3
c     check i(j) = even, triangular inequality, m less than j,
c     find no of terms
      do 110 j=1,9
      k=i(j)/2
      if(i(j)-2*k)300,80,300
   80 if(k)300,90,90
   90 if(k-n)100,110,110
  100 n=k
  110 i(j)=k+1
      if(m3)115,400,115
  115 il=0
      la=i1-i5
      lb=i1-i6
      if(il-la)120,130,130
  120 il=la
  130 if(il-lb)140,145,145
  140 il=lb
c     form coefficient of sum
  145 c=(g(i11)-g(i11-1)+g(i1)+g(i2)+g(i3) -g(i10)+g(i4)+g(i5)+g(i6)+
     *g(i7)+g(i8)+g(i9))/2.0
      j1=i1-il
      j2=i4-il
      j3=i7-il
      m1=il+1
      m2=il-la+1
      m3=il-lb+1
      c=c-g(j1)-g(j2)-g(j3)-g(m1)-g(m2)-g(m3)
      c=exp(c)
      if((il-2*(il/2)).ne.0) c=-c
      if(n)300,150,160
  150 cleb=c
      go to 300
c     form sum
  160 a=j1-1
      b=j2-1
      h=j3-1
      d=m1
      e=m2
      f=m3
      s=1.0
      q=n-1
      do 170 j=1,n
      t=(a-q)/(d+q)*(b-q)/(e+q)*(h-q)/(f+q)
      s=1.0-s*t
      q=q-1.0
  170 continue
      cleb=c*s
  300 return
c     special formula for m3 = 0 and m1 = 0 or 1/2
  400 k=i10/2
      if(i10-2*k)410,420,410
  410 k=1
      go to 430
  420 k=0
  430 if(m1)115,440,460
  440 l=0
      if(k)300,480,300
  460 if(m1-1)115,470,115
  470 l=1
  480 x=l
      m=i3+(i1+k+1)/2-l
      m1=i10/2+k
      m2=i4+i5
      m3=i6+i7
      j1=(i1+1-k)/2
      j2=(i2+1+k-l)/2
      j3=(i3+1+k-l)/2
      cleb=          exp((g(i11)-g(i11-1)+g(i1)+g(i2)+g(i3)-g(i10))/2.+
     *g(m1)-g(j1)-g(j2)-g(j3)+x*(g(3)-(g(m2)-g(m2-1)+g(m3)-g(m3-1))/2.))
      if((m-2*(m/2)).ne.0) cleb=-cleb
      go to 300
      end
c     deck fact
      subroutine fact
c***********************************************************************
c     calcul des logarithmes des factorielles                          *
c***********************************************************************
      common/facto/g(101)
c
      im=101
      g(1)=0.0
      g(2)=0.0
      do 1 j=3,im
      x=float(j-1)
    1 g(j)=g(j-1)+alog(x)
      return
      end
c     deck integ
      subroutine integ(h,fl,sl,spo,jpl)
c***********************************************************************
c     integration de l'equation de schroedinger par la methode         *
c     de numerov matricielle                                           *
c                                                                      *
c     voir      a.c.allison                                            *
c     journal of computational physics  6(1970)378-391                 *
c***********************************************************************
      parameter(npt=200,npt5=npt+5)
      real mi,mt
      dimension y(2,npt5),z(2,npt5),u(npt5),w(npt5)
      dimension f1(2,2),f2(2,2),f3(2,2),s1(2),s2(2)
      double precision f1,f2,f3,s1,s2,det,y,z
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/inout/ie,is,is1,is2
      common/potn/poti(npt5),potr(npt5),vso(npt5),popi(npt5)
      common/psi/psir,psirp,psii,psiip,dwr,dwi
      double precision psir,psirp,psii,psiip,dwr,dwi
      double precision dw1,dw2,dw,ddw,fl4,d23,h212,r,r2
      data d23/0.666666666666667d+00/
c
c
c     formules d'integration
c
c               2         -1           2               2
c     y   = (i+h *f   /12)  *((2*i-10*h *f /12)*y -(i+h *f   /12)*y   ))
c      n+1         n+1                    n      n        n-1      n-1
c
c     f1 = f(n-1) * h**2/12
c     f2 = f( n ) * h**2/12
c     f3 = f(n+1) * h**2/12
c     s1 = i + f1
c     s2 = 2*i - 10*f2
c     det = determinant ( i + f3 )
c
      h212=h*h/(12.0d+00)
      npt3=npt+3
c
      do 1 n=1,npt3
      u(n)=potr(n)+spo*vso(n)
      w(n)=poti(n)
    1 continue
c
c     y(1,n) = partie   reelle   de la fonction d'onde
c     y(2,n) = partie imaginaire de la fonction d'onde
c
c     conditions initiales
      y(1,1)=0.0d+00
      y(2,1)=0.0d+00
      y(1,2)=1.0d-20
      y(2,2)=1.0d-20
c
      if(jpl.ne.1) then
            z(1,1)=0.0d+00
            z(2,1)=0.0d+00
            z(1,2)=1.0d-20
            z(2,2)=1.0d-20
            endif
c
      fl4=4.0d+00*fl
      dwr=0.0d+00
      dwi=0.0d+00
      dw=d23
      ddw=-0.5d+00*dw
c
      f1(1,1)=0.0d+00
      f1(1,2)=0.0d+00
      f1(2,1)=0.0d+00
      f1(2,2)=0.0d+00
c
      f2(1,1)=-1.0d+00
      f2(1,2)= 0.0d+00
      f2(2,1)= 0.0d+00
      f2(2,2)=-1.0d+00
c
      do 10 n=3,npt5
      m=n-2
      r=float(m)*h
      r2=r*r
c
      f3(1,1)=(ak2-sl/r2-u(m))*h212
      f3(1,2)=( +w(m)        )*h212
      f3(2,1)=( -w(m)        )*h212
      f3(2,2)=(ak2-sl/r2-u(m))*h212
c
      s1(1)=((1.d+00)+f1(1,1))*y(1,n-2)+f1(1,2)*y(2,n-2)
      s1(2)=((1.d+00)+f1(2,2))*y(2,n-2)+f1(2,1)*y(1,n-2)
c
      s2(1)=((2.d+00)-(10.d+00)*f2(1,1))*y(1,n-1)
     1      -(10.d+00)*f2(1,2)*y(2,n-1)
      s2(2)=((2.d+00)-(10.d+00)*f2(2,2))*y(2,n-1)
     1     -(10.d+00)*f2(2,1)*y(1,n-1)
c
      det=(f3(1,1)+(1.d+00))*(f3(2,2)+(1.d+00))-f3(1,2)*f3(2,1)
c
      y(1,n)=(f3(2,2)+(1.d+00))*(s2(1)-s1(1))-f3(1,2)*(s2(2)-s1(2))
      y(2,n)=(f3(1,1)+(1.d+00))*(s2(2)-s1(2))-f3(2,1)*(s2(1)-s1(1))
      y(1,n)=y(1,n)/det
      y(2,n)=y(2,n)/det
c
      do 3 i=1,2
      do 3 j=1,2
      f1(i,j)=f2(i,j)
      f2(i,j)=f3(i,j)
    3 continue
c
      if(jpl.eq.1) go to 5
      z(1,n)=y(1,n)
      z(2,n)=y(2,n)
      go to 10
c
    5 if(n.lt.4) go to 10
      nm=n-1
      dy1=y(1,n)-y(1,m)
      dy2=y(2,n)-y(2,m)
      dz1=z(1,n)-z(1,m)
      dz2=z(2,n)-z(2,m)
      dw1=dy1*z(1,nm)-dy2*z(2,nm)-y(1,nm)*dz1+y(2,nm)*dz2
      dw1=dw1+fl4*(y(1,nm)*z(1,nm)-y(2,nm)*z(2,nm))/(m-1)
      dw2=dy1*z(2,nm)+dy2*z(1,nm)-y(1,nm)*dz2-y(2,nm)*dz1
      dw2=dw2+fl4*(y(1,nm)*z(2,nm)+y(2,nm)*z(1,nm))/(m-1)
      dwr=dwr+dw*popi(nm)*dw1
      dwi=dwi+dw*popi(nm)*dw2
      dw=dw+ddw
      ddw=-ddw
   10 continue
c
c     calcul des derivees
c
      n=npt+2
      psir=y(1,n)
      psii=y(2,n)
      psirp=(y(1,n+3)-y(1,n-3)+(9.d+00)*(y(1,n-2)-y(1,n+2))
     *       +(45.d+00)*(y(1,n+1)-y(1,n-1)))/ ((60.d+00)*h)
      psiip=(y(2,n+3)-y(2,n-3)+(9.d+00)*(y(2,n-2)-y(2,n+2))
     *       +(45.d+00)*(y(2,n+1)-y(2,n-1)))/ ((60.d+00)*h)
c
      return
      end
c     deck preang
      subroutine preang(ida,na)
c***********************************************************************
c     calcul des valeurs des polynomes de legendre                     *
c     si ida = 1 pour des angles equirepartis entre 0 et 180 par       *
c                pas de 2.5                                            *
c            = 2 pour des angles dont les cosinus sont equirepartis    *
c                entre -1 et 1 par pas de 0.02                         *
c***********************************************************************
      parameter(nol=40,noa=101)
      common/angles/a(noa),c(noa),pl(nol,noa),pl1(nol,noa),pl2(nol,noa)
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
c
      rad=180./pi
c
      if(ida.eq.2) go to 20
c
c     calculate cosines of equally spaced angles
      na=min0(73,noa)
      d=180./float(na-1)
      do 1 i=1,na
      a(i)=float(i-1)*d
    1 c(i)=cos(a(i)/rad)
      go to 30
c
c     calculate equally spaced cosines
   20 na=noa
      d=2.0/float(na-1)
      do 25 i=1,na
   25 c(i)=1.0-float(i-1)*d
   30 continue
c
      do 50 i=1,na
      x=c(i)
      pl(1,i)=1.0
      pl(2,i)=x
      pl(3,i)=1.5*x*x-0.5
      pl1(1,i)=0.0
      pl1(2,i)=sqrt(1.0-x*x)
      pl1(3,i)=3.0*x*pl1(2,i)
      pl2(1,i)=0.
      pl2(2,i)=0.
      pl2(3,i)=3.*(1.-x*x)
      do 40 l=4,nol
      xl=float(l-1)
      pl(l,i)=((2.0*xl-1.0)*x*pl(l-1,i)-(xl-1.0)*pl(l-2,i))/xl
      pl1(l,i)=((2.0*xl-1.0)*x*pl1(l-1,i)-xl*pl1(l-2,i))/(xl-1.0)
      pl2(l,i)=((2.0*xl-1.0)*x*pl2(l-1,i)-(xl+1.0)*pl2(l-2,i))/(xl-2.0)
   40 continue
   50 continue
c
      return
      end
c     deck pripot
      subroutine pripot(ipot)
c***********************************************************************
c     ecriture des parametres du potentiel optique                     *
c***********************************************************************
      real mi,mt
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/inout/ie,is,is1,is2
      common/poten/a(6),pot(5,5),r(6),beta,ewmax
c
      write(is,1)
      write(is,10)
      write(is,2) zi,mi,ipot,zt,mt
      write(is,3)
      write(is,4) (pot(1,i),i=1,5),r(1),a(1)
      if(beta.gt.0.0) write(is,11) beta
      if(r(2).ge.0.0) write(is,5)
      if(r(2).lt.0.0) write(is,6)
      r2=abs(r(2))
      if(key.ne.0) then
      write(is,16)pot(2,1),pot(2,2),pot(3,4),pot(3,4),pot(2,3),r2,a(2)
      else
      write(is,4)(pot(2,i),i=1,5),r2,a(2)
      endif
      if(a(5).ne.0.0) write(is,13) a(5)
      write(is,7)
      if(key.ne.0) then
      write(is,17)pot(3,1),pot(3,2),pot(3,4),pot(3,4),pot(3,3),r(3),a(3)
      else
      write(is,4)(pot(3,i),i=1,5),r(3),a(3)
      endif
      if(ewmax.ne.0.0) write(is,14) ewmax
      write(is,8)
      write(is,4) (pot(4,i),i=1,5),r(4),a(4)
      if(a(6).gt.0.001) then
            write(is,15)
            write(is,18) (pot(5,i),i=1,5),r(6),a(6)
            endif
      if(r(5).eq.0.0.and.zi.gt.0.0) write(is,12)
      write(is,9) r(5)
c
c     formats
c
    1 format(1h1,14x,80htransmission coefficients calculated from the fo
     *llowing optical model parameters)
    2 format(1h ,34x,6hcharge,29x,4hmass,//,15x,10hprojectile,10x,f6.1,
     *20x,1pe13.6,20x,22hparametrization ipot =,i2,
     */,15x,6htarget,14x,0pf6.1,20x,1pe13.6,//)
    3 format(1h ,15hsaxon real well,/)
    4 format(1h ,10x,4hv = ,f9.4,3h + ,f9.4,6h *e + ,f9.4,9h * e*e + ,
     +f9.4,10h * ln(e) +,/,10x,f9.4,10h * sqrt(e),10x,4hr = ,f8.4,10x,
     *4ha = ,f8.4,//)
    5 format(1h ,31hsaxon derivative imaginary well,/)
    6 format(1h ,23hgaussian imaginary well,/)
    7 format(1h ,29hsaxon imaginary well (volume),/)
    8 format(1h ,12hspin - orbit,/)
    9 format(1h ,17hcoulomb radius = ,f9.4,//)
   10 format(1h ,14x,80(1h*),///)
   11 format(1h ,105x,7hbeta = ,f8.4,//)
   12 format(1h ,20x,38herror : r(5) must be different from 0.,5x,
     *60(1h*))
   13 format(1h ,116x,2h+ ,f8.4,4h * e,//)
   14 format(1h ,46himaginary depth and radius are constant above ,
     *f8.4,4h mev,/)
c   15 format(1h ,21hparity-violating term,/)
   15 format(1h ,33hparity-violating term *2*10**7/hc,/)
   16 format(11x,'v =( ',f9.4,' + ',f9.4,' * ( e -',f9.4,' )**2 )/',/,
     +14x,'(( e -',f9.4,')**2 + ',f9.4,' **2 )**2 ',/,
     +38x,' r =',f9.4,10x,'a ='
     +,f9.4)
   17 format(11x,'v =( ',f9.4,' + ',f9.4,' * ( e -',f9.4,' )**2 )/',/,
     +14x,'(( e -',f9.4,')**2 + ',f9.4,' **2 )',/,
     +38x,' r =',f9.4,10x,'a = ',f9.4)
   18 format(1h ,10x,4hv = ,g9.3,3h + ,g9.3,6h *e + ,g9.3,9h * e*e + ,
     +g9.3,10h * ln(e) +,/,10x,g9.3,10h * sqrt(e),10x,4hr = ,f8.4,10x,
     *4ha = ,f8.4,//)
      return
      end
c     deck pritc
      subroutine pritc(ne)
c***********************************************************************
c     ecriture du tableau recapitulatif des coefficients de            *
c     transmission (moyennes sur les j)                                *
c     et de la section efficace de formation du noyau compose          *
c***********************************************************************
      parameter(noe=200,nol=40)
      real mi,mt
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/tce/tc(noe,nol)
      common/xs/se(noe),sr(noe),st(noe)
c
      write(is,14) mi,zi,mt,zt
c
c      determination du lmax pour l'impression des t(l)
      nec=1
      if(ein(ne).gt.ein(1)) nec=ne
      do 1 l=1,nol
      if(tc(nec,l).lt.0.5e-06) go to 2
    1 continue
    2 lmax=l-1
      lmax=min0(lmax,22)
      kf=1
      if(lmax.gt.11) kf=2
c
      do 10 ke=1,kf
      lmin=0
      lpmax=10
      if(ke.eq.2) go to 11
      write(is,15) (l,l=lmin,lpmax)
      go to 12
   11 lmin=11
      lpmax=lmax-1
      write(is,16) (l,l=lmin,lpmax)
      write(is,17)
   12 continue
      lmin=lmin+1
      lpmax=lpmax+1
      do 20 m=1,ne
      if(ke.eq.1) write(is,21) ein(m),sr(m),(tc(m,l),l=lmin,lpmax)
      if(ke.eq.2) write(is,22) ein(m),(tc(m,l),l=lmin,lpmax)
   20 continue
      write(is,30)
   10 continue
c
c     formats
c
   14 format(1h1,29x,14hprojectile  a=,f4.0,3x,2hz=,f4.0,10x,9hcible  a=
     *,f5.0,3x,2hz=,f5.0,//)
   15 format(1h ,7henergie,2x,8hsigma r.,1x,11(5x,i2,3x),/)
   16 format(1h ,7henergie,11x,11(5x,i2,3x),/)
   17 format(56x)
   21 format(1h ,f6.3,2x,f9.2,1x,11(2x,f8.6))
   22 format(1h ,f6.3,12x,11(2x,f8.6))
   30 format(//)
c
      return
      end
c     deck pritd
      subroutine pritd(ne,ip,iba)
c***********************************************************************
c     ecriture du tableau recapitulatif des coefficients de            *
c     transmission (moyennes sur les j)                                *
c     this subroutine prepares the table of transmission coefficients  *
c     for stapre code                                                  *
c***********************************************************************
      parameter(noe=200,nol=40)
      character*10 aa(6)
      dimension sl(nol)
      real mi,mt
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/iw,is,is1,is2
      common/tce/tc(noe,nol)
      common/xs/se(noe),sr(noe),st(noe)
      common/elb/ elab(noe)
      data aa/' neutron',' proton',' deuteron',' triton',
     *' he-3',' alpha'/
c
c     for interpolation, for each l order, stapre requires transmission
c     coefficents at a minimum of 4 energies. if less than 4 energies
c     output is not possible.
c
      if(iba.ne.0) go to 80
      if(ne.lt.4) go to 60
c
c     determine number of l values to output. output all l orders
c     for which 4 or more energies have transmission coefficients of
c     0.5e-06 or more. assume that energies are in ascending order,
c     transmission coefficents increase with increasing energy and
c     decrease with increasing l value. find highest l value that
c     has 4 or more energies with coefficents of 0.5e-06 or more
c     (do not output coefficents for more than nol l values).
c
      do 5 ie=1,ne
      if(tc(ie,1).ge.0.5e-6) go to 7
 5    continue
      ie=ne
 7    neb=ie
      net=ne-neb+1
      do 20 l=1,nol
      do 10 ie=neb,ne
      if(tc(ie,l).ge.0.5e-06) go to 15
   10 continue
      go to 25
   15 ke=ne-ie+1
      if(ke.lt.4) go to 25
   20 continue
      l=41
   25 lout=l-1
c
c     return if no output.
c
      if(lout.le.0) go to 70
c
c     identify projectile and target.
c
c     write(6,6000) mi,za,mt,zt
c
c     print number of l values for which coefficients will be output
c
      write(7,6030) aa(ip),net,lout
      write(7,6020) (ein(i),i=neb,ne)
c
c     set up loop over l orders.
c
      do 50 l=1,lout
c
c     find lowest energy with transmission coefficents of 0.5e-06
c     or more (insure there are at least 4 energies per l order).
c
      do 30 ie=neb,ne
      if(tc(ie,l).ge.0.5e-06) go to 40
   30 continue
      ie=ne
   40 ke=ne-ie+1
      if(ke.lt.4) ke=4
      ie=ne-ke+1
c
c     output index of the first energy and then the coefficients.
c
      ieb=ie-neb+1
      write(7,6010) ieb
   50 write(7,6020)         (tc(i,l),i=ie,ne)
c     write(7,6030)
      return
c
c     stapre output not possible. print warning message.
c
   60 ke=ne
   70 write(7,6040) ke
      return
 80   sfac=4000.*atan(1.)*(mi+mt)/(4.784468*mi*mt)
      cml=(mi+mt)/mt
      do 110 i=1,ne
      sfe=sfac/ein(i)
      el=cml*ein(i)
      do 90 l=1,nol
      lmax=l
      sl(l)=sfe*(2.*l-1.)*tc(i,l)
      if(sl(l).lt.1.0e-6) go to 100
 90   continue
 100  lmax=max0(1,lmax-1)
      write(7,6050) lmax,el,lmax,aa(ip)
      write(7,6020) (sl(l),l=1,lmax)
 110  continue
      return
c
c     formats
c
 6000 format(1h1,31htransm. coef. for stapre prgr. ,
     *14hprojectile  a=,f4.0,3x,2hz=,f4.0,10x,9htarget a=
     *,f5.0,3x,2hz=,f5.0,//)
 6010 format(i4)
 6020 format(6f12.7)
 6030 format(42x,a10,12x,2i4)
 6040 format(//48h output for stapre is not possible. transmission,
     1 35h coefficents are only available at ,i3,10h energies./
     2 41h stapre requires a minimum of 4 energies.)
 6050 format(i10,f10.3,i10,30x,a10)
      end
c     deck racah
      function racah(a,b,c,d,e,f)
c***********************************************************************
c     calcul des coefficients de racah  w(a,b,c,d;e,f)                 *
c                                                                      *
c     attention :   w(a,b,c,d;e,f) = (-1)**(a+b+c+d)*6-j(a,b,e;d,c,f)  *
c                                                    ---               *
c                                                                      *
c     d'apres john.g.wills     ornl-tm-1949 (august 1967)              *
c                           et comp.phys.comm. 2(1971)381              *
c                                                                      *
c     o.bersillon     aout 1977                                        *
c                                                                      *
c***********************************************************************
      dimension i(16)
      common/facto/g(101)
      equivalence (i(1),i1),(i(2),i2),(i(3),i3),(i(4),i4),(i(5),i5),
     *(i(6),i6),(i(7),i7),(i(8),i8),(i(9),i9),(i(10),i10),(i(11),i11),
     *(i(12),i12),(i(13),i13),(i(14),i14),(i(15),i15),(i(16),i16)
c
      racah=0.0
c     convert arguments to integer and make usefull combinations
      ja=int(2.0*a+0.001)
      jb=int(2.0*b+0.001)
      jc=int(2.0*c+0.001)
      jd=int(2.0*d+0.001)
      je=int(2.0*e+0.001)
      jf=int(2.0*f+0.001)
      i1 =ja+jb-je
      i2 =jb+je-ja
      i3 =je+ja-jb
      i4 =jc+jd-je
      i5 =jd+je-jc
      i6 =je+jc-jd
      i7 =ja+jc-jf
      i8 =jc+jf-ja
      i9 =jf+ja-jc
      i10=jb+jd-jf
      i11=jd+jf-jb
      i12=jf+jb-jd
      i13=ja+jb+je
      i14=jc+jd+je
      i15=ja+jc+jf
      i16=jb+jd+jf
c     check triangular inequalities, find no. of terms in sum,
c     div. i by 2
      n=i16
      do 80 j=1,12
      k=i(j)/2
      if(i(j)-2*k)300,40,300
   40 if(k)300,50,50
   50 if(k-n)60,70,70
   60 n=k
   70 i(j)=k+1
   80 continue
c     find minimum value of summation index
      il=0
      do 100 j=13,16
      i(j)=i(j)/2
      if(il-i(j))90,100,100
   90 il=i(j)
  100 continue
      j1=il-i13+1
      j2=il-i14+1
      j3=il-i15+1
      j4=il-i16+1
      j5=i13+i4-il
      j6=i15+i5-il
      j7=i16+i6-il
      h= -         exp((g(i1)+g(i2)+g(i3)-g(i13+2)+g(i4)+g(i5)+g(i6)-
     *g(i14+2)+g(i7)+g(i8)+g(i9)-g(i15+2)+g(i10)+g(i11)+g(i12)-g(i16+2))
     */2.0+g(il+2)-g(j1)-g(j2)-g(j3)-g(j4)-g(j5)-g(j6)-g(j7))
      if((j5-2*(j5/2)).ne.0) h=-h
      if(n)300,110,120
  110 racah=h
      go to 300
  120 s=1.0
      q=n-1
      p=il+2
      r=j1
      o=j2
      v=j3
      w=j4
      x=j5-1
      y=j6-1
      z=j7-1
      do 130 j=1,n
      t=(p+q)/(r+q)*(x-q)/(o+q)*(y-q)/(v+q)*(z-q)/(w+q)
      s=1.0-s*t
      q=q-1.0
  130 continue
      racah=h*s
  300 return
      end
c     deck rcwfn
      subroutine rcwfn(rhs,ets,minl,maxl,fc,fcp,gc,gcp,accur,step)
c***********************************************************************
c     coulomb wave functions calculated at r = rho                     *
c     by the continued-fraction method of j.w.steed                    *
c     minl, maxl are actual l-values                                   *
c                                                                      *
c     see a.r.barnett, d.h.feng, j.w.steed and l.j.b.goldfarb          *
c     computer physics communications  8 (1974) 377-395                *
c***********************************************************************
      parameter(nol=40,nol1=nol+1)
      implicit double precision(a-h,o-z)
      double precision k,k1,k2,k3,k4,m1,m2,m3,m4
      real   step,accur,ets,rhs
      dimension fc(nol1),fcp(nol1),gc(nol1),gcp(nol1)
      data gpmax/1.0d+60/,one/1.0d+00/,hundr/0.1d+02/,c1/1.0d-12/
      data c2/1.0d-06/,half/0.5d+00/,two/0.2d+01/,zer/0.0d+00/
      data six/0.6d+01/,twous/0.2d+05/,fsous/0.46d+05/
      data otwo/0.2d+00/,thrn/0.999d+03/,oone/0.1d-02/
      pace=dble(step)
      acc=dble(accur)
      rho=dble(rhs)
      eta=dble(ets)
      if(pace.lt.hundr) pace=hundr
      if(acc.lt.c1.or.acc.gt.c2) acc=c2
      r=rho
      ktr=1
      lmax=maxl
      lmin1=minl+1
      xll1=  dble(minl*lmin1)
      eta2=eta*eta
      turn=eta+dsqrt(eta2+xll1)
      if(r.lt.turn.and.dabs(eta).ge.c2) ktr=-1
      ktrp=ktr
      go to 2
    1 r=turn
      tf=f
      tfp=fp
      lmax=minl
      ktrp=1
    2 etar=eta*r
      rho2=r*r
      pl=  dble(lmax+1)
      pmx=pl+half
c     continued fraction for fp(maxl)/f(maxl) , xl is f , xlprime is fp
      fp=eta/pl+pl/r
      dk=etar*two
      del=zer
      d=zer
      f=one
      k=(pl*pl-pl+etar)*(two*pl-one)
      if((pl*pl+pl+etar).ne.zer) go to 3
      r=r+c2
      go to 2
    3 h=(pl*pl+eta2)*(one-pl*pl)*rho2
      k=k+dk+pl*pl*six
      d=one/(d*h+k)
      del=del*(d*k-one)
      if(pl.lt.pmx) del=-r*(pl*pl+eta2)*(pl+one)*d/pl
      pl=pl+one
      fp=fp+del
      if(d.lt.0.0) f=-f
      if(pl.gt.twous) go to 11
      if(dabs(del/fp).ge.acc) go to 3
      fp=f*fp
      if(lmax.eq.minl) go to 5
      fc(lmax+1)=f
      fcp(lmax+1)=fp
c     downward recursion to minl for f and fp, arrays gc,gcp are storage
      l=lmax
      do 4 lp=lmin1,lmax
      pl=  dble(l)
      gc(l+1)=eta/pl+pl/r
      gcp(l+1)=dsqrt(eta2+pl*pl)/pl
      fc(l)=(gc(l+1)*fc(l+1)+fcp(l+1))/gcp(l+1)
      fcp(l)=gc(l+1)*fc(l)-gcp(l+1)*fc(l+1)
    4 l=l-1
      f=fc(lmin1)
      fp=fcp(lmin1)
    5 if(ktrp.eq.-1) go to 1
c     repeat for r = turn if rho lt turn
c     now obtain p + i.q for minl from continued fraction (32)
c     real arithmetic to facilitate conversion to ibm using real*8
      p=zer
      q=r-eta
      pl=zer
      ar=-(eta2+xll1)
      ai=eta
      br=two*q
      bi=two
      wi=two*eta
      dr=br/(br*br+bi*bi)
      di=-bi/(br*br+bi*bi)
      dp=-(ar*di+ai*dr)
      dq=ar*dr-ai*di
    6 p=p+dp
      q=q+dq
      pl=pl+two
      ar=ar+pl
      ai=ai+wi
      bi=bi+two
      d=ar*dr-ai*di+br
      di=ai*dr+ar*di+bi
      t=one/(d*d+di*di)
      dr=d*t
      di=-t*di
      h=br*dr-bi*di-one
      k=bi*dr+br*di
      t=dp*h-dq*k
      dq=dp*k+dq*h
      dp=t
      if(pl.gt.fsous) go to 11
      cnt=(dabs(dp)+dabs(dq))-((dabs(p)+dabs(q))*acc)
      if(cnt) 66,66,6
   66 p=p/r
      q=q/r
c     solve for fp,g,gp and normalise f at l = minl
      g=(fp-p*f)/q
      gp=p*g-q*f
      if(dabs(g).lt.1.0d100) then
        w=one/dsqrt(fp*g-f*gp)
       else
        w=one/(dsqrt(g)*dsqrt(fp-(f/g)*gp))
       endif
      g=w*g
      gp=w*gp
      if(ktr.eq.1) go to 8
      f=tf
      fp=tfp
      lmax=maxl
c     runge-kutta integration of g(minl) and gp(minl) inwards from turn
c
      if(rho.lt.(otwo*turn)) pace=thrn
      r3=0.33333333333333333d+00
      h=(rho-turn)/(pace+one)
      h2=half*h
      i2=idint(pace+oone)
      etah=eta*h
      h2ll=h2*xll1
      s=(etah+h2ll/r)/r-h2
    7 rh2=r+h2
      t=(etah+h2ll/rh2)/rh2-h2
      k1=h2*gp
      m1=s*g
      k2=h2*(gp+m1)
      m2=t*(g+k1)
      k3=h*(gp+m2)
      m3=t*(g+k2)
      m3=m3+m3
      k4=h2*(gp+m3)
      rh=r+h
      s=(etah+h2ll/rh)/rh-h2
      m4=s*(g+k3)
      g=g+(k1+k2+k2+k3+k4)*r3
      gp=gp+(m1+m2+m2+m3+m4)*r3
      r=rh
      i2=i2-1
      gpg=gp
      if(dabs(gpg).gt.gpmax) go to 11
      if(i2.ge.0) go to 7
      w=one/(fp*g-f*gp)
c     upward recursion from gc(minl) and gcp(minl),stored values are r,s
c     renormalise fc,fcp for each l-value
    8 gc(lmin1)=g
      gcp(lmin1)=gp
      if(lmax.eq.minl) go to 10
      do 9 l=lmin1,lmax
      t=gc(l+1)
      gc(l+1)=(gc(l)*gc(l+1)-gcp(l))/gcp(l+1)
      gcp(l+1)=gc(l)*gcp(l+1)-gc(l+1)*t
      fc(l+1)=w*fc(l+1)
    9 fcp(l+1)=w*fcp(l+1)
      fc(lmin1)=fc(lmin1)*w
      fcp(lmin1)=fcp(lmin1)*w
      go to 12
   10 fc(lmin1)=w*f
      fcp(lmin1)=w*fp
      go to 12
   11 w=zer
      g=zer
      gp=zer
      go to 8
   12 return
      end
c     deck scat
      subroutine scat(lmax,ipr)
c***********************************************************************
c     calcul des coefficients de transmission t(l,j)                   *
c     et des amplitudes de diffusion eta(l,j)                          *
c***********************************************************************
      parameter(noe=200,nol=40,nol1=nol+1)
      parameter(npt=200,npt5=npt+5)
      real mi,mt,mt3,mu
      double precision a1,a2,a3,a4,a5,ar,ai,t5,t6,t7,t8
      double precision fc,fcp,gc,gcp
      dimension a(6)
      dimension u1(7),y1(7)
      dimension rv(6),pote(5)
      dimension fc(nol1),fcp(nol1),gc(nol1),gcp(nol1)
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/poten/aa(6),pot(5,5),rr(6),beta,ewmax
      common/potn/poti(npt5),potr(npt5),vso(npt5),popi(npt5)
      common/psi/psir,psirp,psii,psiip,dwr,dwi
      double precision psir,psirp,psii,psiip,dwr,dwi
      common/tlj/br(3,nol),bi(3,nol),t(3,nol)
      data hc7o2/9.866e-6/
c
      epstl=1.0e-06
      fws=0.
      fwv=0.
c
      do 1 j=1,nol
      do 1 i=1,3
      br(i,j)=0.0
      bi(i,j)=0.0
      t(i,j)=0.0
    1 continue
c
c     constantes
c
      mu=mi*mt/(mi+mt)
      el=e1*(mi+mt)/mt
      w2=4.784468*mu/100.
      ak2=w2*e1
      ak=sqrt(ak2)
      zz=zi*zt
      eta=0.15748603*zz*sqrt(mu/e1)
             write(7,*)zi,zt,mu,e1,eta
      mt3=mt**0.333333
c
c     mu = masse reduite
c     ak2 = k**2
c     ak = k
c     w2 = 2*amu/(hb**2)
c
c     rayons reel
c
      do 2 i=1,6
    2 rv(i)=abs(rr(i))*mt3
c
c     profondeurs reelles
c
c*****
      do 5 i=1,4
      a(i)=aa(i)
      if(a(i).lt.0.001) a(i)=1.0
    5 continue
      a(6)=aa(6)
      if(a(6).lt.0.001) a(6)=1.0
c*****
c
c      key=0        normal potentials
c      otherwise    parametrized dispersion terms used
c
      if(key.eq.0) go to 12
      do 3 i=1,4,3
      pote(i)=pot(i,1)+pot(i,2)*e1+pot(i,3)*e1*e1+
     +  pot(i,4)*alog(e1)+pot(i,5)*sqrt(e1)
      if(pote(i).lt.0) pote(i)=0.0
    3 continue
      pote(5)=pot(5,1)+pot(5,2)*e1+pot(5,3)*e1*e1+
     +  pot(5,4)*alog(e1)+pot(5,5)*sqrt(e1)
      pote(5)=hc7o2*pote(5)
      def=pot(3,4)-e1
      def2=def*def
      eo2=pot(2,3)**2
      dw=1./(eo2+def2)
      fws=0.5*(pot(2,1)+pot(2,2)*eo2)*dw
      dw=dw*dw
      fws=def*(fws/eo2+(pot(2,1)-pot(2,2)*eo2)*dw)/pot(2,3)
      pote(2)=(pot(2,1)+pot(2,2)*def2)*dw
      eo2=pot(3,3)**2
      dw=1./(eo2+def2)
      fwv=(pot(3,1)-pot(3,2)*eo2)*def*dw/pot(3,3)
      pote(3)=(pot(3,1)+pot(3,2)*def2)*dw
      go to 8
c
c*****
c
   12 do 11 i=1,4
      pote(i)=pot(i,1)+pot(i,2)*el+pot(i,3)*el*el+pot(i,4)*alog(el)+
     + pot(i,5)*sqrt(el)
      if(pote(i).lt.0.) pote(i)=0.0
   11 continue
      pote(5)=pot(5,1)+pot(5,2)*el+pot(5,3)*el*el+pot(5,4)*alog(el)+
     + pot(5,5)*sqrt(el)
      pote(5)=hc7o2*pote(5)
c
c     diffusivites
c
      a(2)=a(2)+aa(5)*el
c
c     potentiels imaginaires constants
c
      if(ewmax*(el-ewmax))8,8,6
    6 continue
      do 7 i=2,3
      pote(i)=pot(i,1)+pot(i,2)*ewmax+pot(i,3)*ewmax*ewmax
     *                +pot(i,4)*alog(ewmax)+pot(i,5)*sqrt(ewmax)
      if(pote(i).lt.0.0) pote(i)=0.0
    7 continue
c
      a(2)=aa(2)
      if(a(2).lt.0.001) a(2)=1.0
      a(2)= a(2)+aa(5)*ewmax
c
c
c     rayon de raccordement
c
    8 r1=rv(1)+7.0*a(1)
      r2=rv(2)+7.0*a(2)
      r3=rv(3)+7.0*a(3)
      rm=1.5*amax1(r1,r2,r3)
      rho=ak*rm
c
      if(rv(6)+7.0*a(6).gt.rm.and.aa(6).gt.0.001) write(is,9)
    9 format(//,' WARNING -- parity-violating DWBA integrals extend ',
     1   'beyond integration limit.')
c
c     fonctions de coulomb au point de raccordement
c
      accur=1.0e-14
      step=999.0
      call rcwfn(rho,eta,0,nol,fc,fcp,gc,gcp,accur,step)
c
c     calcul des potentiels
c
      npt3=npt+3
      h=rm/float(npt)
c     h = pas d'integration
c
      vcl=ak*eta/w2
      if(beta.le.0.0) go to 15
      pote(3)=0.0
      c2=beta*beta/16.
      c1=4.*w2*c2
      res=4.*c2*ak2
      d1=exp(res)
   15 continue
c
      t1=1.0/exp(rv(1)/a(1))
      t2=1.0/exp(rv(2)/a(2))
      t3=1.0/exp(rv(3)/a(3))
      t4=1.0/exp(rv(4)/a(4))
      t5=1.0/exp(rv(6)/a(6))
c
      dt1=exp(h/a(1))
      dt2=exp(h/a(2))
      dt3=exp(h/a(3))
      dt4=exp(h/a(4))
      dt5=exp(h/a(6))
c
      do 100 i=1,npt3
      r=float(i)*h
c
c     potentiel reel
      t1=t1*dt1
      potr(i)=+pote(1)/(1.0+t1)
c
c     potentiel imaginaire de volume
      t3=t3*dt3
      poti(i)=+pote(3)/(1.0+t3)
      potr(i)=potr(i)+fwv/(1.0+t3)
      if(rr(2))19,19,25
   25 continue
c
c     + potentiel imaginaire de surface ( derivee de w - s )
      t2=t2*dt2
      poti(i)=poti(i)+4.0*pote(2)*t2/((1.0+t2)**2)
      potr(i)=potr(i)+4.0*fws*t2/((1.0+t2)**2)
      go to 20
c
c     + potentiel imaginaire de surface ( gaussien)
   19 yy=-(((r-rv(2))/a(2))**2)
      if(yy.lt.-600.) yy=-600.
      poti(i)=poti(i)+pote(2)*exp(yy)
   20 continue
      if(beta)23,23,24
c
c     potentiel "equivalent local" d'un potentiel non-local
   24 continue
      p1=c2/(potr(i)*potr(i)+poti(i)*poti(i))
      p4=t1/(1.+t1)
      p5=t2/(1.+t2)
      p2=-potr(i)*p4/a(1)
      p3=poti(i)*(1.-2.*p5)/a(2)
      p4= p2    *(1.-2.*p4)/a(1)
      p5=poti(i)*(1.-6.*p5*(1.-p5))/(a(2)*a(2))
      u2=p1*((potr(i)*p2+poti(i)*p3)*2./r+potr(i)*p4+poti(i)*p5)
      y2=p1*((potr(i)*p3-poti(i)*p2)*2./r+potr(i)*p5-poti(i)*p4)
      y1(1)=poti(i)/(d1*d1+2.*d1*c1*potr(i))
      u1(1)=(potr(i)+poti(i)*c1*y1(1))/(d1+c1*potr(i))
      do 21 k=1,6
      p1=c1*y1(k)-y2
      p2=sin(p1)
      p1=cos(p1)
      p3=1./(d1*exp(c1*u1(k)-u2))
      u1(k+1)=(potr(i)*p1+poti(i)*p2)*p3
      y1(k+1)=(poti(i)*p1-potr(i)*p2)*p3
   21 continue
      p1=u1(7)-2.*u1(6)+u1(5)
      if(p1)31,32,31
   32 potr(i)=u1(7)
      go to 33
   31 potr(i)=u1(7)-((u1(7)-u1(6))**2)/p1
   33 p2=y1(7)-2.*y1(6)+y1(5)
      if(p2)35,34,35
   34 poti(i)=y1(7)
      go to 23
   35 poti(i)=y1(7)-((y1(7)-y1(6))**2)/p2
   23 continue
c
c     potentiel reel + potentiel coulombien
      if(zz)14,18,14
   14 if(r-rv(5))17,16,16
   16 potr(i)=potr(i)-2.0*vcl/r
      go to 18
   17 potr(i)=potr(i)-vcl*(3.0-(r/rv(5))**2)/rv(5)
   18 continue
c
c     potentiel spin - orbite
      t4=t4*dt4
      vso(i)=+2.043655*pote(4)*t4/(a(4)*r*((1.0+t4)**2))
c
c     parity - violating term
      t5=t5*dt5
      popi(i)=pote(5)/(1.+t5)
c
      potr(i)=-potr(i)*w2
      poti(i)=-poti(i)*w2
      vso(i) =-vso(i) *w2
      popi(i)= popi(i)*w2
  100 continue
c
      ipl=ifix(2.0*si+1.001)
c     ipl = 2*s+1
      do 200 l=1,nol
      lmax=l
      fl=float(l-1)
      sl=fl*(fl+1.)
      fj=fl-si-1.
      do 190 j=1,ipl
      fj=fj+1.
      if(fj.lt.abs(fl-si)) go to 190
      spo=fj*(fj+1.)-fl*(fl+1.)-si*(si+1.)
c
      call integ(h,fl,sl,spo,j)
c
c     psir  =               partie   reelle   de la fonction d'onde int.
c     psirp = derivee de la partie   reelle   de la fonction d'onde int.
c     psii  =               partie imaginaire de la fonction d'onde int.
c     psiip = derivee de la partie imaginaire de la fonction d'onde int.
      t5=fc (l)
      t6=fcp(l)*ak
      t7=gc (l)
      t8=gcp(l)*ak
      a1=-(t5*psirp-t6*psir+t7*psiip-t8*psii)
      a2=-(t5*psiip-t6*psii-t7*psirp+t8*psir)
      a3=+(t5*psirp-t6*psir-t7*psiip+t8*psii)
      a4=+(t5*psiip-t6*psii+t7*psirp-t8*psir)
      a5=a1*a1+a2*a2
c     a1 = partie   reelle   du denominateur
c     a2 = partie imaginaire du denominateur
c     a3 = partie   reelle   du numerateur
c     a4 = partie imaginaire du numerateur
c     a5 = module du denominateur
      etr=sngl(1.0d+00-(a3*a1+a2*a4)/a5)
      eti=sngl((a4*a1-a2*a3)/a5)
c     etr = 1.0 - partie   reelle   de l'amplitude de diffusion
c     eti =       partie imaginaire de l'amplitude de diffusion
                  write(6,*)l,etr,eti
      br(j,l)=etr
      bi(j,l)=eti
      t (j,l)=sngl(1.0d+00-(a3*a3+a4*a4)/a5)
c
      if(ipl.ne.2) go to 190
      a1=t7*etr+t5*eti
      a2=t5*(2.-etr)+t7*eti
      a5=2.*(psir*psir+psii*psii)
      if(j.eq.ipl) then
        ar=(a2*psir+a1*psii)/a5
        ai=(a1*psir-a2*psii)/a5
       else
        a3=(a2*psir+a1*psii)/a5
        a4=(a1*psir-a2*psii)/a5
        a1=a3*ar-a4*ai
        a2=a3*ai+a4*ar
        br(3,l-1)=-2.*sngl(a2*dwr+a1*dwi)/ak
        bi(3,l-1)=2.*sngl(a1*dwr-a2*dwi)/ak
        etr=br(3,l-1)**2+bi(3,l-1)**2
        t(1,l)=t(1,l)-etr
        t(2,l-1)=t(2,l-1)-etr
        t(3,l-1)=sqrt(etr)
       endif
c
  190 continue
      if(br(ipl,l))205,205,195
  195 if(abs(t(ipl,l)/t(ipl,1))-epstl)210,200,200
  200 continue
      go to 210
c
  205 do 206 j=1,3
      br(j,lmax)=0.0
      bi(j,lmax)=0.0
      t (j,lmax)=0.0
  206 continue
      lmax=lmax-1
  210 lm1=lmax-1
      if(ipr.gt.0) write(is,10) ak,eta,rm,h,npt,lm1
   10 format(1h ,3x,3hk =,1pe12.5,5x,5heta =,e12.5,5x,4hrm =,e12.5,5x,
     *4hdr =,e12.5,5x,i4,7h points,10x,6hlmax =,i3,///)
c
      return
      end
c     deck shapec
      subroutine shapec(lmax,ida,na,ipl,ipu,el)
c***********************************************************************
c     shape elastic differential cross-section                         *
c     for charged particles                                            *
c***********************************************************************
      parameter(nol=40,noa=101)
      complex cgamma
      complex ai,unc,z,zeroc
      complex csigl(nol),etac(3,nol),fc(noa)
      complex c1(nol),c2(nol),c3(nol),c4(nol),c5(nol)
      complex z1,z2,z3,z4,z5
      dimension sigl(nol),dsig(noa),dsigr(noa)
      dimension rap(noa),polar(noa),ppolar(noa),qpolar(noa)
      common/angles/ang(noa),cang(noa),pl(nol,noa),pl1(nol,noa),
     1                                                 pl2(nol,noa)
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/inout/ie,is,is1,is2
      common/tlj/etar(3,nol),etai(3,nol),t(3,nol)
      data ai/(0.0,1.0)/
      data unc/(1.0,0.0)/
      data zeroc/(0.0,0.0)/
      data sq18/0.35355339/
c
      is3=7
      ak=sqrt(ak2)
c
      do 10 i=1,na
      rap(i)=0.0
      dsig(i)=0.0
      dsigr(i)=1.0
      polar(i)=0.0
      fc(i)=zeroc
   10 continue
      rap(1)=1.0
c
      do 20 l=1,lmax
      csigl(l)=unc
      do 15 i=1,ipl
      etac(i,l)=cmplx(1.0-etar(i,l),etai(i,l))
	    15 continue
      if(ipl.eq.2) etac(3,l)=cmplx(etar(3,l),etai(3,l))
   20 continue
c
      if(eta)25,90,25
c
c     calcil des dephasages coulombiens
   25 continue
      z=cmplx(1.0,eta)
      z=cgamma(z)
      sig0=aimag(clog(z))
	               write(6,*)sig0
c
      som=0.0
      do 40 l=1,lmax
      al=float(l)
      sigl(l)=sig0+som
      csigl(l)=cexp(2.0*ai*sigl(l))
      som=som+atan2(eta,al)                                                                                                                                                     ,al)
   40 continue
c
c     calcul de l'amplitude de diffusion coulombienne
      do 50 i=1,na
      arg=0.5*(1.0-cang(i))
      if(arg)49,50,49
   49 continue
      z=-ai*eta*alog(arg)+2.0*ai*sig0
      fc(i)=-eta*cexp(z)/(2.0*arg)
      dsigr(i)=cabs(fc(i))**2/ak2
      dsigr(i)=10.0*dsigr(i)
   50 continue
c
   90 go to (100,200,300),ipl
c
c     spin 0 particles
c     ****************
  100 continue
      do 110 l=1,lmax
      al=float(l-1)
      c1(l)=csigl(l)*(2.0*al+1.0)*(1.0-etac(1,l))
  110 continue
c
      do 130 i=1,na
      z1=zeroc
      do 120 l=1,lmax
      z1=z1+c1(l)*pl(l,i)
  120 continue
      z1=(fc(i)+ai*z1/2.0)/ak
      dsig(i)=cabs(z1)**2
      dsig(i)=10.0*dsig(i)
      if(i-1)130,130,129
  129 rap(i)=dsig(i)/dsigr(i)
  130 continue
      go to 1000
c
c     spin 1/2 particles
c     ******************
  200 continue
c
      ell=1000.*el
c
      x3=0.0
      do 210 l=1,lmax
	   write(6,*)l,clog(etac(2,l))
      al=float(l-1)
*************************  15/10/2007	 ************************
	write(25,*)al
	write(26,*)clog(etac(2,l))
c	write(28,*)dreal(clog(etac(2,l)))
c	write(29,*)dimag(clog(etac(2,l)))
	write(27,*)al, clog(etac(2,l))
***************************************************************
      c1(l)=csigl(l)*((al+1.0)*(1.0-etac(2,l))+al*(1.0-etac(1,l)))
      c2(l)=csigl(l)*(etac(1,l)-etac(2,l))
      c3(l)=cexp(ai*(sigl(l)+sigl(l+1)))*(al+1.0)*etac(3,l)
c
c      if(l.eq.1) then
c        c3(l)=5.1e-9*(0.,1.)/(ell/30.0-(1.,-0.5))/sqrt(1.0e3*ell)
c       else
c        c3(l)=0.0d0
c       endif
c
  210 continue
c
c        c3(1)=c3(1)+5.1e-9*(0.,1.)/(ell/30.0-(1.,-0.5))/sqrt(1.0e3*ell)
c
      z2=2.*c2(2)/c1(1)
      z3=2.*c3(1)/c1(1)
      write(21,'(e15.8,5x,e15.8)') ell,aimag(z2)
      write(22,'(e15.8,5x,e15.8)') ell,-real(z2)
      write(23,'(e15.8,5x,e15.8)') ell,real(z3)
      write(24,'(e15.8,5x,e15.8)') ell,aimag(z3)
      write(is,409) -ai*z2,z3
c      
      do 230 i=1,na
      z1=zeroc
      z2=zeroc
      z3=zeroc
      do 220 l=1,lmax
      z1=z1+c1(l)*pl(l,i)
      z2=z2+c2(l)*pl1(l,i)
      z3=z3+c3(l)*(pl(l,i)+pl(l+1,i))
  220 continue
      z1=(fc(i)+ai*z1/2.0)/ak
      z2=-z2/(2.0*ak)
      z3=-ai*z3/(2.0*ak)
      x3=amax1(x3,cabs(z3))
      dsig(i)=cabs(z1)**2+cabs(z2)**2
      dsig(i)=10.0*dsig(i)
c  sign of polar(i) changed to agree with Rose (??) or Wolfenstein (??)
      polar(i)=-20.0*real(z2*conjg(z1))/dsig(i)
      ppolar(i)=-20.0*real(z3*conjg(z1))/dsig(i)
      qpolar(i)=-20.0*aimag(z3*conjg(z1))/dsig(i)
      if(i-1)230,230,229
  229 rap(i)=dsig(i)/dsigr(i)
  230 continue
      go to 1000
c
c     spin 1 particles
c     ****************
  300 continue
      do 310 l=1,lmax
      al=float(l-1)
      al2=2.*al+1.
      fl=al*(al+1.)
      c1(l)=csigl(l)*(al2-(al+1.)*etac(3,l)-al*etac(1,l))
      c2(l)=csigl(l)*(al2*(2.-etac(2,l))-(al+2.)*etac(3,l)
     1                                       -(al-1.)*etac(1,l))
      c3(l)=csigl(l)*(etac(1,l)-etac(3,l))
      c4(l)=csigl(l)*((fl+al)*etac(3,l)-al2*etac(2,l)
     1                            -(fl-al-1.)*etac(1,l))/fl
      c5(l)=csigl(l)*(al2*etac(2,l)-al*etac(3,l)-(al+1.)*etac(1,l))/fl
  310 continue
      do 330 i=1,na
      z1=zeroc
      z2=zeroc
      z3=zeroc
      z4=zeroc
      z5=zeroc
      do 320 l=1,lmax
      z1=z1+c1(l)*pl(l,i)
      z2=z2+c2(l)*pl(l,i)
      z3=z3+c3(l)*pl1(l,i)
      z4=z4+c4(l)*pl1(l,i)
      z5=z5+c5(l)*pl2(l,i)
  320 continue
      z1=(fc(i)-0.5*ai*z1)/ak
      z2=(fc(i)-0.25*ai*z2)/ak
      z3=-ai*sq18*z3/ak
      z4=-ai*sq18*z4/ak
      z5=-0.25*ai*z5/ak
      dsig(i)=cabs(z2)**2+cabs(z3)**2+cabs(z4)**2+cabs(z5)**2
      dsig(i)=10.*(cabs(z1)**2+2.*dsig(i))/3.
      polar(i)=aimag(z1*conjg(z3)-z2*conjg(z4)-z4*conjg(z5))
      polar(i)=80.*sq18*polar(i)/(3.*dsig(i))
      if(i-1) 330,330,329
  329 rap(i)=dsig(i)/dsigr(i)
  330 continue
      go to 1000
c
 1000 continue
c
c
c     ecriture de la distribution angulaire
      if(x3.lt.1.0e-20) then
        if(ida.eq.1) write(is,410)
        if(ida.eq.2) write(is,411)
       else
        if(ida.eq.1) write(is,412)
        if(ida.eq.2) write(is,413)
       endif
      naa=na/2+1
      do 420 i=1,naa
      imin=2*(i-1)+1
      imax=2*i
      if(imax.gt.na) imax=na
      if(x3.lt.1.0d-20) then
        if(ida.eq.1)write(is,422)( ang(j),dsig(j),dsigr(j),rap(j),
     *    polar(j),j=imin,imax)
        if(ida.eq.2)write(is,422)(cang(j),dsig(j),dsigr(j),rap(j),
     *    polar(j),j=imin,imax)
       else
        if(ida.eq.1)write(is,423)( ang(j),dsig(j),dsigr(j),rap(j),
     *    polar(j),ppolar(j),qpolar(j),j=imin,imax)
        if(ida.eq.2)write(is,423)(cang(j),dsig(j),dsigr(j),rap(j),
     *    polar(j),ppolar(j),qpolar(j),j=imin,imax)
      endif
  420 continue
c
      if(ipu.eq.0) go to 999
      do 450 i=1,na
      write(is3,451) ang(i),dsig(i),dsigr(i),rap(i),polar(i)
  450 continue
c     formats
c
  409 format(1h1,'  P0(pc)=',e12.5,'  Q0(pc)=',e12.5,'  P0(pv)=', 
     * e12.5,'  P0(pv)=',e12.5,//)
  410 format(1h1,46x,40hshape elastic differential cross-section,/,1h ,
     * 46x,40(1h*),///,1h ,2(4x,9h     teta,3x,10hdsig(teta),2x,
     * 11hdsigr(teta),3x,10hdsig/dsigr,1x,12hpolarisation),/)
  411 format(1h1,46x,40hshape elastic differential cross-section,/,1h ,
     * 46x,40(1h*),///,1h ,2(4x,9hcos(teta),3x,10hdsig(teta),2x,
     * 11hdsigr(teta),3x,10hdsig/dsigr,1x,12hpolarisation),/) 
  412 format(1h1,46x,40hshape elastic differential cross-section,/,1h ,
     * 46x,40(1h*),///,1h ,4x,9h     teta,3x,10hdsig(teta),2x,
     * 11hdsigr(teta),3x,10hdsig/dsigr,1x,12hpolarisation,1x,
     * 12hpar-viol pol,3x, 10hpar-viol q,/)
  413 format(1h1,46x,40hshape elastic differential cross-section,/,1h ,
     * 46x,40(1h*),///,1h ,4x,9hcos(teta),3x,10hdsig(teta),2x,
     * 11hdsigr(teta),3x,10hdsig/dsigr,1x,12hpolarisation,1x,
     * 12hpar-viol pol,3x, 10hpar-viol q,/)
  422 format(1h ,1p,5e13.5,2x,5e13.5)
  423 format(1h ,1p,7e13.5)
  451 format(1p6e12.5)
c
  999 return
      end
c     deck shapel
      subroutine shapel(lmax,ida,na,ipl)
c***********************************************************************
c     shape elastic differential cross-section                         *
c     for neutrons                                                     *
c***********************************************************************
      parameter(nol=40,noa=101)
      dimension bl(nol),da(noa)
      common/angles/a(noa),c(noa),pl(nol,noa),pl1(nol,noa),pl2(nol,noa)
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/inout/ie,is,is1,is2
      common/tlj/br(3,nol),bi(3,nol),t(3,nol)
c
      rad=180./pi
      lmax2=2*lmax-1
c
      do 100 il=1,lmax2
      al=float(il-1)
      som=0.0
c
      do 90 il1=1,lmax
      al1=float(il1-1)
c
      do 80 il2=1,lmax
c     test   l + l1 + l2  pair
      k1=il+il1+il2
      k2=k1/2
      if(k1-2*k2)50,80,50
   50 continue
      al2=float(il2-1)
      aj1=al1-si-1.0
c
      do 70 ij1=1,ipl
      aj1=aj1+1.0
      if(aj1.lt.abs(al1-si)) go to 70
      aj2=al2-si-1.0
c
      do 60 ij2=1,ipl
      aj2=aj2+1.0
      if(aj2.lt.abs(al2-si)) go to 60
      z=sqrt((2.0*al1+1.0)*(2.0*al2+1.0)*(2.0*aj1+1.0)*(2.0*aj2+1.0))*
     *cleb(al1,al2,al,0.0,0.0,0.0)*racah(al1,aj1,al2,aj2,si,al)
c     z = coefficient de blatt et biedenharn, voir valeurs tabulees
c         l.c.biedenharn  ornl-1501 (1953)
      z2=z*z
      s=br(ij1,il1)*br(ij2,il2) + bi(ij1,il1)*bi(ij2,il2)
      som=som+z2*s
   60 continue
   70 continue
   80 continue
   90 continue
      bl(il)=10.0*som/(8.0*ak2)
  100 continue
c
c     distribution angulaire
      do 210 i=1,na
      som=0.0
      do 200 l=1,lmax2
  200 som=som+bl(l)*pl(l,i)
      da(i)=som
  210 continue
c
c     calcul de l'integrale de la d.a.
      som=0.0
      if(ida.eq.2) go to 290
      do 270 i=2,na
  270 som=som+0.5*(da(i)*sin(a(i)/rad)+da(i-1)*sin(a(i-1)/rad))
      som=som*2.5/rad
      go to 310
  290 continue
      do 300 i=1,na
  300 som=som+da(i)
      som=som-0.5*(da(1)+da(na))
      som=som*0.02
  310 som=2.0*pi*som
c
c     ecriture de la distribution angulaire
      if(ida.eq.1) write(is,225)
      if(ida.eq.2) write(is,221)
      naa=na/4+1
      do 220 i=1,naa
      imin=4*(i-1)+1
      imax=4*i
      if(imax.gt.na) imax=na
      if(ida.eq.1) write(is,222) (a(j),da(j),j=imin,imax)
      if(ida.eq.2) write(is,222) (c(j),da(j),j=imin,imax)
  220 continue
c
      b1=4.0*pi*bl(1)
      write(is,240) som,b1
c
c     ecriture des coefficients du developpement en polynome de
c     legendre  ( definition de endf )
      do 400 l=2,lmax2
      il=l-1
      al=float(il)
      bl(l)=bl(l)/((2.0*al+1.0)*bl(1))
  400 continue
      bl(1)=1.0
      write(is,401)
      nl=lmax2/5
      if((lmax2-5*nl).gt.0) nl=nl+1
      do 410 l=1,nl
      lmi=5*(l-1)
      lma=5*l-1
      if((lma+1).gt.lmax2) lma=lmax2-1
      write(is,402) (lm,bl(lm+1),lm=lmi,lma)
  410 continue
c
c     formats
c
  221 format(1h1,46x,40hshape elastic differential cross-section,/,1h ,
     *46x,40(1h*),///,1h ,5x,4(9hcos(teta),2x,15hd.sigma/d.omega,6x),/)
  222 format(1h ,5x,4(1pe12.5,2x,e12.5,6x))
  225 format(1h1,46x,40hshape elastic differential cross-section,/,1h ,
     *46x,40(1h*),///,1h ,5x,4(9h    teta ,2x,15hd.sigma/d.omega,6x),/)
  240 format(//,1h ,10hintegral =,1pe12.5,3h mb,10x,7hbl(0) =,e12.5,
     *3h mb,//)
  401 format(///,1h ,7x,5(2x,1hl,7x,5hbl(l),10x),/)
  402 format(1h ,7x,5(i3,3x,1pe14.7,5x))
      return
      end
c     deck spin0
      subroutine spin0(n,lmax,ipr)
c***********************************************************************
c     ecriture des coefficients de transmission t(l,j)                 *
c     calcul des sections efficaces - formation du noyau compose       *
c                                     shape elastique                  *
c                                     totale                           *
c     spin de la particule incidente = 0                               *
c***********************************************************************
      parameter(noe=200,nol=40)
      real mi,mt
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/tce/tc(noe,nol)
      common/tlj/br(3,nol),bi(3,nol),t(3,nol)
      common/xs/se(noe),sr(noe),st(noe)
c
      se(n)=0.0
      sr(n)=0.0
      st(n)=0.0
c
      if(ipr.ne.0) write(is,5)
c
      do 10 l=1,lmax
      k=l-1
      al2=2.0*float(l-1)+1.0
      tm=t(1,l)
      tc(n,l)=tm
      sr(n)=sr(n)+al2*(2.0*br(1,l)-br(1,l)*br(1,l)-bi(1,l)*bi(1,l))
      se(n)=se(n)+al2*(br(1,l)*br(1,l)+bi(1,l)*bi(1,l))
      st(n)=st(n)+al2*2.0*br(1,l)
      if(ipr.eq.0) go to 10
      write(is,20) k,tm,br(1,l),bi(1,l)
   10 continue
c     section efficaces en mb
      se(n)=10.0*se(n)*pi/ak2
      sr(n)=10.0*sr(n)*pi/ak2
      st(n)=10.0*st(n)*pi/ak2
c
c     if(ipr.eq.0) go to 200
      el=e1*(mi+mt)/mt
      write(is,30) e1,sr(n),el,se(n),st(n)
c     formats
c
c
    5 format(1h ,3x,1hl,7x,5htc(l),3x,9h1 - eta r,7x,5heta i,/)
   20 format(1h ,i4,1p3e12.4)
   30 format(//,1h ,f8.4,7h mev cm,5x,30hcompound nucleus cross section,
     *1pe14.7,/,1h ,0pf8.4,8h mev lab,7x,27hshape elastic cross section,
     *1pe14.7,/,1h ,31x,19htotal cross section,e14.7,3h mb,//)
c
  200 return
      end
c     deck spin05
      subroutine spin05(n,lmax,ipr)
c***********************************************************************
c     ecriture des coefficients de transmission t(l,j)                 *
c     calcul des sections efficaces - formation du noyau compose       *
c                                     shape elastique                  *
c                                     totale                           *
c     spin de la particule incidente = 1/2                             *
c***********************************************************************
      parameter(noe=200,nol=40)
      real mi,mt
      dimension pe(noe),pr(noe),pt(noe)
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/poten/a(6),pot(5,5),r(6),beta,ewmax
      common/tce/tc(noe,nol)
      common/tlj/br(3,nol),bi(3,nol),t(3,nol)
      common/xs/se(noe),sr(noe),st(noe)
c
      se(n)=0.0
      sr(n)=0.0
      st(n)=0.0
c
      pe(n)=0.0
      pr(n)=0.0
      pt(n)=0.0
c
c   Used to calculate the j=1/2 cross sections
c      br(2,2)=0.
c      bi(2,2)=0.
c      br(3,2)=0.
c      bi(3,2)=0.
c
c   To calculate just the j=1/2 l=1 cross section
c      br(2,1)=0.
c      bi(2,1)=0.
c
      if(ipr.eq.0) go to 8
      if(a(6).lt.0.001) then
            write(is,5)
            else
            write(is,7)
            endif
c
c    Used to calculate the j=1/2 cross sections
c    8 do 10 l=1,2
c
    8 do 10 l=1,lmax
      k=l-1
      al=float(l-1)
      al1=al+1.0
      al2=2.0*al+1.0
      tm=(al1*t(2,l)+al*t(1,l))/al2
      tc(n,l)=tm
      sr(n)=sr(n)+al1*(2.0*br(2,l)-br(2,l)*br(2,l)-bi(2,l)*bi(2,l))
     *           +al *(2.0*br(1,l)-br(1,l)*br(1,l)-bi(1,l)*bi(1,l))
     *           -(al2+1.0)*t(3,l)*t(3,l)
      se(n)=se(n)+al1*(br(2,l)*br(2,l)+bi(2,l)*bi(2,l))
     *           +al *(br(1,l)*br(1,l)+bi(1,l)*bi(1,l))
     *           +(al2+1.0)*t(3,l)*t(3,l)
      st(n)=st(n)+al1*(2.0*br(2,l))
     *           +al *(2.0*br(1,l))
      pr(n)=pr(n)+(al2+1.0)*real(br(3,l)*(2.0-br(2,l)-br(1,l+1))
     *                         +bi(3,l)*(bi(2,l)+bi(1,l+1)))
      pe(n)=pe(n)+(al2+1.0)*real(br(3,l)*(br(2,l)+br(1,l+1))
     *                          -bi(3,l)*(bi(2,l)+bi(1,l+1)))
      pt(n)=pt(n)+(al2+1.0)*2.0*br(3,l)
      if(ipr.eq.0) go to 10
      b1=(al1*br(2,l) + al*br(1,l))/al2
      b2=(al1*bi(2,l) + al*bi(1,l))/al2
      if(a(6).lt.0.001) then
            write(is,20) k,tm,b1,b2,(t(i,l),br(i,l),bi(i,l),i=1,2)
            else
            write(is,20) k,(t(i,l),br(i,l),bi(i,l),i=1,3)
            endif
   10 continue
   11 continue
c     section efficaces en mb
      pe(n)=pe(n)/se(n)
      pr(n)=pr(n)/sr(n)
      pt(n)=pt(n)/st(n)
c Used to calculate j=1/2 cross sections
c      pe(n)=10.0*pe(n)*pi/ak2
c      pr(n)=10.0*pr(n)*pi/ak2
c      pt(n)=10.0*pt(n)*pi/ak2
      se(n)=10.0*se(n)*pi/ak2
      sr(n)=10.0*sr(n)*pi/ak2
      st(n)=10.0*st(n)*pi/ak2
c
      if(e1.lt.0.1) then
            erg=1.0e+06*e1
            s0=tc(n,1)/(2.0*pi*sqrt(erg))
            r2=r(1)*r(1)*(mt**0.666666667)
            p1=(ak2*r2)/(1.0+ak2*r2)
            s1=tc(n,2)/(2.0*pi*p1*sqrt(erg))
            rp=sqrt(se(n)/(40.0*pi))
            pp=2.*br(3,1)/br(1,2)
            write(18,800) mt,s0,s1,pp
            else
            pp=0.
            endif
c
c     if(ipr.eq.0) go to 200
      el=e1*(mi+mt)/mt
      write(14,800) el,t(2,1),t(1,2),t(2,2)
      write(15,800) el,se(n),sr(n),st(n)
c      write(16,800) el,pe(n),pr(n),pt(n)
      write(16,800) el,st(n)*pt(n)
      write(17,800) el,t(3,1),t(3,2),pp
      if(a(6).lt.0.001) then
            write(is,30) e1,sr(n),el,se(n),st(n)
            if(e1.lt.0.1) write(is,31) s0,s1,rp
            else
            write(is,40) e1,sr(n),pr(n),el,se(n),pe(n),st(n),pt(n)
            if(e1.lt.0.1) write(is,32) s0,s1,rp,pp
            endif
 800  format(1pe12.4,2x,3e12.4)
c
c     formats
c
    5 format(1h ,3x,1hl,7x,5htc(l),3x,9h1 - eta r,7x,5heta i,9x,
     *10ht(l,l-1/2),3x,9h1 - eta r,7x,5heta i,9x,10ht(l,l+1/2),3x,
     *9h1 - eta r,7x,5heta i,/)
    7 format(1h ,3x,1hl,2x,10ht(l,l-1/2),3x,9h1 - eta r,7x,5heta i,
     *9x,10ht(l,l+1/2),3x,9h1 - eta r,7x,5heta i,10x,9hdt(l,l+1),7x,
     *5heta r,7x,5heta i,/)
   20 format(1h ,i4,1p3e12.4,2(7x,3e12.4))
   30 format(//,1h ,f10.6,7h mev cm,5x,30hcompound nucleus cross section
     *,1pe14.7,/,1h ,0pf10.6,8h mev lab,7x,27hshape elastic cross sectio
     *n,1pe14.7,/,1h ,33x,19htotal cross section,e14.7,3h mb,//)
   31 format(1h ,29x,23hstrength functions   s0,1pe14.7,/,1h ,50x,2hs1,
     *e14.7,/,1h ,35x,17hscattering radius,e14.7,//)
   32 format(1h ,29x,23hstrength functions   s0,1pe14.7,/,1h ,50x,2hs1,
     *e14.7,/,1h ,35x,17hscattering radius,e14.7,/,1h ,38x,
     *14hp1/2 asymmetry,e14.7,//)
   40 format(//,1h ,f10.6,7h mev cm,5x,30hcompound nucleus cross section
     *,1pe14.7,3h mb,5x,'absorption asymmetry ',e14.7,/,1h ,0pf10.6,
     *8h mev lab,7x,27hshape elastic cross section,1pe14.7,3h mb,8x,
     *'elastic asymmetry ',e14.7,/,1h ,33x,19htotal cross section,
     *e14.7,3h mb,10x,'total asymmetry ',e14.7,//)
c
  200 return
      end
c     deck spin1
      subroutine spin1(n,lmax,ipr)
c***********************************************************************
c     ecriture des coefficients de transmission t(l,j)                 *
c     calcul des sections efficaces - formation du noyau compose       *
c                                     shape elastique                  *
c                                     totale                           *
c     spin de la particule incidente = 1                               *
c***********************************************************************
      parameter(noe=200,nol=40)
      real mi,mt
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/tce/tc(noe,nol)
      common/tlj/br(3,nol),bi(3,nol),t(3,nol)
      common/xs/se(noe),sr(noe),st(noe)
c
      se(n)=0.0
      sr(n)=0.0
      st(n)=0.0
c
      if(ipr.ne.0) write(is,5)
c
      do 10 l=1,lmax
      k=l-1
      al=float(l-1)
      al1=2.0*al-1.0
      al2=2.0*al+1.0
      al3=2.0*al+3.0
      tm=(al3*t(3,l)+al2*t(2,l)+al1*t(1,l))/(3.0*al2)
      tc(n,l)=tm
      sr(n)=sr(n)+al3*(2.0*br(3,l)-br(3,l)*br(3,l)-bi(3,l)*bi(3,l))
     *           +al2*(2.0*br(2,l)-br(2,l)*br(2,l)-bi(2,l)*bi(2,l))
     *           +al1*(2.0*br(1,l)-br(1,l)*br(1,l)-bi(1,l)*bi(1,l))
      se(n)=se(n)+al3*(br(3,l)*br(3,l)+bi(3,l)*bi(3,l))
     *           +al2*(br(2,l)*br(2,l)+bi(2,l)*bi(2,l))
     *           +al1*(br(1,l)*br(1,l)+bi(1,l)*bi(1,l))
      st(n)=st(n)+al3*2.0*br(3,l)
     *           +al2*2.0*br(2,l)
     *           +al1*2.0*br(1,l)
      if(ipr.eq.0) go to 10
      b1=(al3*br(3,l)+al2*br(2,l)+al1*br(1,l))/(3.0*al2)
      b2=(al3*bi(3,l)+al2*bi(2,l)+al1*bi(1,l))/(3.0*al2)
      write(is,20) k,tm,b1,b2,(t(i,l),br(i,l),bi(i,l),i=1,3)
   10 continue
      se(n)=se(n)/3.0
      sr(n)=sr(n)/3.0
      st(n)=st(n)/3.0
c     section efficaces en mb
      se(n)=10.0*se(n)*pi/ak2
      sr(n)=10.0*sr(n)*pi/ak2
      st(n)=10.0*st(n)*pi/ak2
c
c     if(ipr.eq.0) go to 200
      el=e1*(mi+mt)/mt
      write(is,30) e1,sr(n),el,se(n),st(n)
c
c     formats
c
    5 format(1h ,1x,1hl,7x,5htc(l),3x,9h1 - eta r,7x,5heta i,3x,
     *8ht(l,l-1),1x,9h1 - eta r,5x,5heta i,5x,6ht(l,l),1x,9h1 - eta r,
     *5x,5heta i,3x,8ht(l,l+1),1x,9h1 - eta r,5x,5heta i,/)
   20 format(1h ,i2,1p3e12.4,3(1x,3e10.3))
   30 format(//,1h ,f8.4,7h mev cm,5x,30hcompound nucleus cross section,
     *1pe14.7,/,1h ,0pf8.4,8h mev lab,7x,27hshape elastic cross section,
     *1pe14.7,/,1h ,31x,19htotal cross section,e14.7,3h mb,//)
c
  200 return
      end
c     deck syspot
      subroutine syspot(ip,ipot)
c***********************************************************************
c     potentiels extraits de la compilation de c.m.perey et f.g.perey  *
c     atomic data and nuclear data tables   17 (1976) 1-101            *
c***********************************************************************
      real mt,mt2,mt3,nmzsa
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/poten/a(6),pot(5,5),r(6),beta,ewmax
c
      mt2=mt*mt
      mt3=mt2*mt
      nmzsa=(mt-2.0*zt)/mt
c
      go to (10,20,30,40,50,60),ip
c
   10 go to (11,12,13,14,15),ipot
c     neutrons
c     ********
   11 continue
c     parametres de wilmore - hodgson                 40 < a
c                                                          e < # 10
      r(1)=1.322-7.6e-04*mt+4.0e-06*mt2-8.0e-09*mt3
      a(1)=0.66
      pot(1,1)=47.01
      pot(1,2)=-0.267
      pot(1,3)=-0.00118
      r(2)=1.266-3.7e-04*mt+2.0e-06*mt2-4.0e-09*mt3
      a(2)=0.48
      pot(2,1)=9.52
      pot(2,2)=-0.053
      r(4)=r(1)
      a(4)=a(1)
      pot(4,1)=7.0
      go to 100
   12 continue
c     parametres de beccheti - greenless              40 < a
c                                                   # 10 < e < 50
      r(1)=1.17
      a(1)=0.75
      pot(1,1)=56.3-24.0*nmzsa
      pot(1,2)=-0.32
      r(2)=1.26
      a(2)=0.58
      pot(2,1)=13.0-12.0*nmzsa
      pot(2,2)=-0.25
      r(3)=1.26
      a(3)=0.58
      pot(3,1)=-1.56
      pot(3,2)=0.22
      r(4)=1.01
      a(4)=0.75
      pot(4,1)=6.2
      go to 100
   13 continue
c     parametres de ferrer et al.                     24 < a < 209
c     nucl.phys. a275(1977)325-341                         e = 11
c
      r(1)=1.27
      a(1)=0.71
      pot(1,1)=47.14-22.50*nmzsa
      r(2)=1.27
      a(2)=0.434
      pot(2,1)=12.16-2.03*nmzsa
      r(4)=1.08
      a(4)=0.71
      pot(4,1)=4.55
      go to 100
   14 continue
c     parametres bersillon cindro
c     contribution to the 5th international symposium on
c     interactions of fast neutrons with nuclei,
c     gaussig, ddr, 17-21 nov 1975
c
      r(1)=1.182+1.93e-04*mt
      a(1)=0.65
      pot(1,1)=71.
      r(2)=1.21
      a(2)=0.47
      pot(2,1)=7.
      pot(2,2)=0.4
      r(4)=r(1)
      a(4)=a(1)
      pot(4,1)=7.
      beta=0.85
      go to 100
   15 continue
c     parametres de madland                      actinides
c     harwell conference                         e < 10 mev
c     september 25-29, 1978
c     -----  temporary values  -----
c
      r(1)=1.264
      a(1)=0.612
      pot(1,1)=50.378-27.073*nmzsa
      pot(1,2)=-0.354
      r(2)=1.256
      a(2)=0.553
      a(5)=0.0144
      pot(2,1)=9.265-12.666*nmzsa
      pot(2,2)=-0.232
      pot(2,3)=+0.03318
      r(4)=1.01
      a(4)=0.75
      pot(4,1)=6.2
      ewmax=10.
      go to 100
c
   20 go to (21,22),ipot
c     protons
c     *******
   21 continue
c     parametres de perey                        30 < a < 100
c                                                     e < 20
      r(1)=1.25
      a(1)=0.65
      pot(1,1)=53.3+27.0*nmzsa+0.4*zt/(mt**0.333333)
      pot(1,2)=-0.55
      r(2)=1.25
      a(2)=0.47
      pot(2,1)=13.5
      r(4)=1.25
      a(4)=0.47
      pot(4,1)=7.5
      r(5)=1.25
      go to 100
   22 continue
c     parametres de beccheti - greenless          a < 40
c                                            20 < e < 50
      r(1)=1.17
      a(1)=0.75
      pot(1,1)=54.0+24.0*nmzsa+0.4*zt/(mt**0.333333)
      pot(1,2)=-0.32
      r(2)=1.32
      a(2)=0.51+0.7*nmzsa
      pot(2,1)=11.8+12.0*nmzsa
      pot(2,2)=-0.25
      r(3)=1.32
      a(3)=0.51+0.7*nmzsa
      pot(3,1)=-2.7
      pot(3,2)=0.22
      r(4)=1.01
      a(4)=0.75
      pot(4,1)=6.2
      r(5)=1.25
      go to 100
c
   30 go to (31,32),ipot
c     deuterons
c     *********
   31 continue
c     parametres de lohr - haeberli               40 < a
c                                                  8 < e < 13
      r(1)=1.05
      a(1)=0.86
      pot(1,1)=91.13+2.2*zt/(mt**0.333333)
      r(2)=1.43
      a(2)=0.5+0.013*(mt**0.666667)
      pot(2,1)=218.0/(mt**0.666667)
      r(4)=0.75
      a(4)=0.50
      pot(4,1)=7.0
      r(5)=1.3
      go to 100
   32 continue
c     parametres de perey                             12 < e < 25
      r(1)=1.15
      a(1)=0.81
      pot(1,1)=81.0+2.0*zt/(mt**0.333333)
      pot(1,2)=-0.22
      r(2)=1.34
      a(2)=0.68
      pot(2,1)=14.4
      pot(2,2)=0.24
      r(5)=1.15
      go to 100
c
   40 continue
c     tritons
c     *******
c     parametres de beccheti - greenless              40 < a
c                                                          e < 40
      r(1)=1.20
      a(1)=0.72
      pot(1,1)=165.0-6.4*nmzsa
      pot(1,2)=-0.17
      r(3)=1.40
      a(3)=0.84
      pot(3,1)=46.0-110.0*nmzsa
      pot(3,2)=-0.33
      r(4)=1.20
      a(4)=0.72
      pot(4,1)=2.5
      r(5)=1.30
      go to 100
c
   50 continue
c     helium-3
c     ********
c     parametres de beccheti - greenless              40 < a
c                                                          e < 40
      r(1)=1.20
      a(1)=0.72
      pot(1,1)=151.9+50.0*nmzsa
      pot(1,2)=-0.17
      r(3)=1.40
      a(3)=0.88
      pot(3,1)=41.7-44.0*nmzsa
      pot(3,2)=-0.33
      r(4)=1.20
      a(4)=0.72
      pot(4,1)=2.5
      r(5)=1.30
      go to 100
c
   60 go to (61,62),ipot
c     alphas
c     ******
   61 continue
c     parametres moyens
c     mac fadden et satchler  nucl.phys. 84(1966)177
c
      r(1)=1.40
      a(1)=0.52
      pot(1,1)=185.
      r(3)=1.40
      a(3)=0.52
      pot(3,1)=25.
      r(5)=1.40
      go to 100
   62 continue
c     mac fadden and satchler    nea parameters  13/10/83
c
      r(1)=1.442
      a(1)=0.520
      pot(1,1)=164.7
      r(3)=1.442
      a(3)=0.520
      pot(3,1)=22.4
      r(5)=1.25
      go to 100
c
  100 return
      end
c     deck tpun
      subroutine tpun(ip,ne,ipl,it,lmax)
c***********************************************************************
c     write the transm. coef. for gnash code on coftr                  *
c***********************************************************************
      parameter(noe=200,nol=40)
      real mi,mt
      common/const/mi,si,zi,mt,zt,pi,ak2,eta,key
      common/ener/e1,ein(noe)
      common/inout/ie,is,is1,is2
      common/poten/a(6),pot(5,5),r(6),beta,ewmax
      common/tlj/br(3,nol),bi(3,nol),t(3,nol)
      common/xs/se(noe),sr(noe),st(noe)
      character*10 aa(6)
      data aa/' neutron',' proton',' deuteron',' triton',
     *' he-3',' alpha'/
c
      nl=min0(lmax,nol)
      k=0
      ip1=6/ipl
c
      if(it.ne.0) go to 100
c
      do 10 j=1,nol,ip1
      ju=j+ip1-1
      if(ju.gt.nol) ju=nol
      write(is2,15) ((t(k,l),k=1,ipl),l=j,ju)
   10 continue
      go to 200
c
 100  continue
      rewind is2
      net=ne+1
      nlp=nl*ipl
      write(7,16) aa(ip),net,nlp,k
c
      do 130 i=1,ne,6
      k=k+1
      iu=i+5
      if(iu.gt.ne) iu=ne
      write(7,17) (ein(j),j=i,iu),k
  130 continue
      nup=ip1*((nl-1)/ip1+1)
c
      do 140 i=1,ne
      do 140 j=1,nol,ip1
      ju=j+ip1-1
      if(ju.gt.nol) ju=nol
      read(is2,15) ((t(k,l),k=1,ipl),l=j,ju)
      if(ju.gt.nup) go to 140
      write(7,17) ((t(k,l),k=1,ipl),l=j,ju),i
  140 continue
c
      rewind is2
c
c     formats
c
   15 format(1p6e12.5)
   16 format(42x,a10,12x,2i4,i8)
   17 format(1p6e12.5,0p,i8)
c
  200 return
      end
