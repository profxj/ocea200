      program srnchn
c
c --- this is a package of routines for a finite-element solution of
c --- the quasi-geostrophic equations,as described in "a baroclinic
c --- quasi-geostrophic open ocean model" (1981) by r.n. miller,
c --- a.r. robinson, and d.b. haidvogel of harvard.
c
c --- bottom topography (for slopes less than o(r0)) and top vertical
c --- velocity have been implemented (10/15/82).
c
c --- a driving package is required to provide initial and boundary
c --- values of the streamfunction, vorticity, and top and bottom
c --- density variations; as well as the topography. this package
c --- is called infld.
c
c --- a variety of options are available:
c
c        ifdif = 1/0     finite difference/collocation in depth
c        ifpert = 1/0    calculate perturbation fields (yes/no)
c        ifrst = 1/0     restart at time rdt (yes/no)
c        iftop = 1/0     use top density information (yes/no)
c        ifbot = 1/0     use bottom density information (yes/no)
c        iftvv = 1/0     use top vertical velocity (yes/no)
c        ifrel = 1/0     use bottom relief (yes/no)
c        ifsbl = 1/0     use surface boundary layer
c        ifeva = 1/0     execute eva preprocessor
c
	include 'prob2.pert'
c
      integer ldv, ldz, maxn, maxnev, maxncv
c
c
c.....The following parameter statement is for qg normal mode on day 5.
c	parameter(kpnt=18)
c.....The following paremeter statment is for the fg mode in energy on day5
 	parameter(kpnt=9)
      parameter(maxn=(m0-1)*n0*k0, maxnev=100,
     &   maxncv=40,ldv=maxn,ldz=maxn)
c      parameter(maxn=(m0-1)*n0*k0, maxnev=998,
c     &   maxncv=1000,ldv=maxn,ldz=maxn)
c
      real vss(ldv,maxncv), zz(ldz, maxnev),
     &     workl(3*maxncv*maxncv+6*maxncv),
     &     workd(3*maxn), dsrn(maxnev), resid(maxn)
     &     ,d(maxncv,3),workev(3*maxncv)
        real tmpwrk(maxn),zzstm(ldz),tmpwrk1(maxn)
      logical select(maxncv)
      integer iparam(11), ipntr(14)
c
      real amm1(mn0,k0),amm2(mn0,k0)
c
      character bmat*1, which*2
      integer   ido, npnt, nev, ncv, lworkl, info, j
      logical   rvec,first
      real      tol, sigma, snrm2
c
c
      character*4 restid,runid
      character*70 idstrg
c>>>>>>>>>>>
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/grate/gre(mn0),grip(mn0),tezm(mn0),aipzm(mn0)
c
c@@@@@@@@@THE FOLLOWING COMMON BLOCK IS ONLY IN CLINIC!!!!
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/egl2/pl,akl,alls,allv
c
	common/mmde/jmode
c
	common/nrest/irst
c
	common/adrst/ityrst
c
	common/cnorm/jnorm
c
	common/jctrl/ktim,icstart,iczero
c
	common/gcord/ir1,ir2,jr1,jr2
c
        common/ccnt/callcnt
c
        common/nmmt/nmit,nmiter,idone
c
        common/comm1/vmat(k0,k0),vinv(k0,k0)
c
        common/amiv/amx(m0-1,m0-1),amy(n0,n0)
c
        real work(m0-1),det(2),z(m0-1),worky(n0),zy(n0)
        integer ipvt(m0-1),ipvty(n0)
c
c>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/count/icnt,itcnt,ibcnt
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),topo(mn0),
     2     topd(mn0),botd(mn0)
c
	common/atop/radt(mn0),s2t(mn0),st(mn0),fyt(mx0,4)
c
	common/abot/radb(mn0),rad1b(mn0),rad2b(mn0),s2b(mn0),sb(mn0)
     1       ,fyb(mx0,4)
c
c
c
	common/adm/tiam(k0,k0),taml(k0,k0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
c
	common/grad/gomeg(mn0,k0),gphi(mn0,k0),ggam(mn0),
     1     gdel(mn0),dm(mn0,k0),dp(mn0,k0),dg(mn0),dd(mn0)
     2    ,cgf
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/steps/sumt,sumb,cost
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/ctest/ctol
c
	common/swoi/iter
c
	common/toff/tadj0
c
	common/istar/icall
c
        common/iwfq/ibgw,iclw
c
	common/bdata/scrp(n0,k0,2),vcrp(n0,k0,2)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/lan1/maxt,ifts
c
c
        common/newbg/amms(mn0,k0),ammv(mn0,k0)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      character*80 ident,fniin,fniorr,fname
c
      character*80 ritzfil,vecfil
c
      character charcen*9,charpre*5,charvec*11
c
      dimension pdt(4),nplev(4),iplev(4,k02),iplot(k02),lplot(k0)
      dimension therm(m20,k0),dxxpsi(m20),dyypsi(m20)
	real crap(mn0)
	real xx(m0)
	real bigs(mnbig,k0),bigv(mnbig,k0)
     &       ,bigt(mnbig),bigb(mnbig)
c
	integer ints(kpnt),jnts(kpnt),ipnts(kpnt)
c
        real tbs(78,49),tbv(78,49)
c
      equivalence (nprnt, nplot)
c
c --- rdt = time at which to restart
c --- iin = unit from which to read job specifications
c --- iout = unit on which to write job results
c --- idff = unit on which to write full field diagnostics
c --- idif = unit on which to write interior field diagnostics
c --- iorr = unit from which to read restart data
c --- iowr = unit on which to write restart data
c --- iplt = unit on which to write plot data
c --- iced = no of times diagnostic routine has been called
c --- iwrite = write title for plot data (0=yes/1=no)
c
      data rdt/1000./,runid/'    ','    ','    ','    ','    ','    '/
      data restid/'    ','    ','    ','    ','    ','    ','    ',
     +            '    ','    '/
      data iin/3/,iout/4/,idff/7/,idif/8/,iorr/9/,iowr/10/,iplt/11/
      data iced/0/, iwrite/0/, idiag/0/, irest/0/, ifndd/0/
c     data icnt,itcnt,ibcnt/0,0,0/
c
c
c......The following data are for the qg normal mode on day 5.....
c	data ints/24,22,18,18,25,35,41,39,45,55,48,45,55,60,67,
c    &            74,38,41/
c	data jnts/40,36,29,23,26,29,23,15,8,12,18,30,33,24,18,22,
c    &            38,42/
c
c......The following data are for the fg mode in energy on day 5....
 	data ints/37,37,38,39,42,47,43,46,36/
 	data jnts/27,22,18,13,10,7,28,31,31/
c
c.....initialise filter counters...
c
      icnt=0
      itcnt=0
      ibcnt=0
c
      charpre='ritz.'
      charvec='srnvec_big.'
c
c---   open units for input/output and read * file names
c
c      ident    used to generate the file names for the
c	        listing file (*.out), full field diagnostics
c	        (*.dff), interior field diagnostics (*.dif),
c	        restart output file (*.wrs), and plot output
c		file (*.plt)
c
c	fniin	the parameter input file
c
c	fniorr  the input restart file
c
c
c       write(idstrg,
c    +  '(''@(#)version 1.12 of the harvard open ocean model'')')
c       write(*,'('' enter identifier for file names    '',$)')
c	read(*,'(a)')ident
c       write(*,'('' enter name of input parameter file '',$)')
c	read(*,'(a)')fniin
c       write(*,'('' enter name of restart file         '',$)')
c	read(*,'(a)')fniorr
c
c.....open unit 99 and use to write out the qg model fields for the
c     adjoint calculation......
c
c>>>>>>>>>>>>
c
c        call lib$get_symbol('gsfil',charcen)
c
c        ritzfil=charpre//charcen
c        vecfil=charvec//charcen
c
c        print *,'rtz=',ritzfil(1:lenstr(ritzfil,80))
c        print *,'vec=',vecfil(1:lenstr(vecfil,80))
c
c
        open(unit=8,file='prob2.input',status='old')
c
c       open(unit=99,file='qgfield',form='unformatted')
c       open(unit=95,file='qginf',form='unformatted')
        open(unit=3,file='prob2.dat3')
        open(unit=5,file='chan.dat5')
        open(unit=4,file='chan.out4')
c        open(unit=6,file='prob2.out6')
c       open(unit=10,file='restart',form='unformatted')
        open(unit=23,file='ftnmplt.dat',form='unformatted')
        open(unit=43,file='ftnm.dat',form='unformatted')
c
c
c>>>>>>>>>>>>>
c
c.......read in the adjoint model parameters.......
c
	read(5,*)sigs
	write(6,*)' enter weighting for vorticity observations'
	read(5,*)sigv
	write(6,*)' tolerance test on cost func.? 1/0 (yes/no)'
	read(5,*)itol
	if(itol.eq.0)
     1  write(6,*)' maximum number of iterations?'
	if(itol.eq.0)read(5,*)numi
	write(6,*)' enter tolerance for adjoint
     1  model assimilation convergence'
	read(5,*)toler
	write(6,*)' enter the initial step size estimate'
	read(5,*)stepi
	write(6,*)' enter icg: 0 for fletcher-reeves scheme'
	write(6,*)'            1 for polak-ribiere scheme  '
	read(5,*)icg
	write(6,*)'which scheme do you want?'
	write(6,*)'1=conjugate gradient, 2=scaled c-g,'
	read(5,*)isem
	write(6,*)' use steepest descent every nsk iterations?
     1    1/0 (yes/no)'
	read(5,*)idis
	write(6,*)' enter iteration interval between conjugate
     1  gradient and steepest descent methods'
	read(5,*)nsk
	write(6,*)' employ conjugacy test also? 1/0 (yes/no)'
	read(5,*)itest
	write(6,*)' enter tolerance for conjugacy test'
	read(5,*)ctol
	write(6,*)'starting time relative control run (days)'
	read(5,*)tadj0
	write(6,*)'adjoint filter parameters: iadn,iadf,iadc'
	read(5,*)iadn,iadf,iadc
c
	write(6,*)'sigs= ',sigs,'  sigv= ',sigv,'  toler= ',toler
	write(6,*)'initial step size= ',stepi
	if(isem.eq.1)write(6,*)'conjugate gradient used'
	if(isem.eq.2)write(6,*)'scaled conjugate gradient used'
	if(isem.eq.1.or.isem.eq.2.and.icg.eq.0)
     1   write(6,*)'fletcher-reeves scheme selected'
	if(isem.eq.1.or.isem.eq.2.and.icg.eq.1)
     1   write(6,*)'polak-ribiere scheme used'
	if(idis.eq.1)
     1   write(6,*)'steepest descent used every ',nsk,' timesteps'
	if(idis.eq.1.and.itest.eq.1)
     1    write(6,*)'conjugacy test employed also for switching to steepest
     2    descent, with tolerance= ',ctol
	write(6,*)'starting time relative to control run = ',tadj0
c
C>>>>>>>>>>>
c
	if(jmode.ne.2.and.ityrst.eq.1)then
	write(6,*)'ITYRST MUST BE ZERO FOR OPTIMAL PERTURBATION RUNS!'
	stop
	endif
c
C>>>>>>>>>>>
c
c
c
c     iflen=lnblk(ident,80)
c     open(unit=iin,file=fniin,form='formatted')
c     fname=ident(:iflen)//'.out'
c     open(unit=iout,file=fname,form='formatted')
c     fname=ident(:iflen)//'.dff'
c     open(unit=idff,file=fname,form='formatted')
c     fname=ident(:iflen)//'.dif'
c     open(unit=idif,file=fname,form='formatted')
c     fname=ident(:iflen)//'.plt'
c     open(unit=iplt,file=fname,form='unformatted')
c
c......initialise the adjoint pointer variables......
c
	iter=1
	ipass=2
	iflag=0
c
c......set icall = 1. this is a parameter for infld....
c
        icall=1
c
c
c --- initialize timing and get specs for this run
c --- an end code read in estart terminates processing
c
   10 continue
      icalc=0
      ntitl=0
c      call timstp(runid)
c
c --- read and write info on data cards
c
      write(4, *)idstrg
      write(4, 600)runid
c     write(4,*)' file identifier:    ',a,/,
c    +             ' parameter file :    ',a,/,
c    +             ' restart file   :    ',a,/,ident,fniin,fniorr
c
c
    1 continue
      read(iin, *, end=120)icard
      if(icard.le.0 .or. icard.gt.19) go to 132
      go to(101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 
     +  111, 112, 113, 114, 115, 116, 117, 118), icard
c
  101 continue
      ntitl=ntitl+1
      if(ntitl.gt.2)ntitl=2
      read(iin, 501, end=130)(titl(ititl,ntitl),ititl=1,20)
      write(4, 601)(titl(ititl,ntitl),ititl=1,20)
      go to 1
c
  102 continue
      read(iin, *, end=130)alpha, beta, rf
c
        read(8,*)beta
c
      write(4, 602)alpha, beta, rf
      go to 1
c
  103 continue
      read(iin, *, end=130)m, n, kz
      write(4, 603)m, n, kz
      kzp1=kz+1
      go to 1
c
  104 continue
      read(iin, *, end=130)(hz(i), i=1, kz)
      write(4, 604)kz
      write(4,1234)
     1  (hz(i), i=1, kz)
 1234 format(10x,4(f9.3,2x))
      go to 1
c
c --- note that sigz(1) is no longer scsig and the
c --- addition of sigz(1) and sigz(kz+1)
c
  105 continue
      read(iin, *, end=130)scsig, (sigz(i), i=1, kzp1)
      write(4,605)scsig
      write(4,1235)
     1   (sigz(i),i=1,kzp1)
 1235 format(10x,4(f9.3,2x))
      go to 1
c
  106 continue
      read(iin, *,end=130)iltord,iltfrq,iltcnt
      write(4,606)iltord,iltfrq,iltcnt
      go to 1
c
c --- real array for driver
c
  107 continue
      read(iin,*,end=130)(pram(i),i=1,16)
      write(4,607)(pram(i),i=1,16)
      go to 1
  108 continue
      read(iin,*,end=130)dt,tstart,tmax,xbasin
      write(4,608)dt,tstart,tmax,xbasin
      go to 1
c
c --- switches for finite diff./colloc., f-/beta-plane, diagnostics,
c --- top and bottom density, top vertical velocity, and bottom
c --- topography.
c
  109 continue
      read(iin,*,end=130)ifdiff,ifpert,ifrst,iftop,ifbot,
     +                   iftvv,ifrel,ifsbl,ifeva
      if(iftvv.ne.0)iftop=1
      if(ifrel.ne.0)ifbot=1
      write(4,609)ifdiff,ifpert,ifrst,iftop,ifbot,
     +                   iftvv,ifrel,ifsbl,ifeva
      if(ifrst.eq.0)go to 1
c     open(unit=iorr,file=fniorr,form='unformatted',status='old')
c      read(iorr,end=190)restid,titl
      write(4,618)restid,titl
      go to 1
c
c --- if collocation used in depth, read in eigenstuff:
c
  110 continue
      read(iin,*,end=130)(eigval(i),i=1,kz)
      do 1110 i=1,kz
        read(iin,*,end=130)(amat(j1,i),j1=1,kz)
 1110 continue
      go to 1
c
c --- parameters for diagnostic options
c
  111 continue
      read(iin,*,end=130)nddt,(ddt(i),i=1,iabs(nddt))
      if(nddt.lt.0)ifndd=1
      nddt=iabs(nddt)
      if(nddt.gt.0)go to 1111
      write(4,1611)
      goto 1
 1111 continue
      write(4,610)
      do 2111 i=1,nddt
        write(4,611)ddt(i)
 2111 continue
      go to 1
c
c --- parameters for array printing options
c
  112 continue
      read(iin,*,end=130)nadt,(adt(i),i=1,nadt)
      if(nadt.gt.0)go to 1112
      write(4,1612)
      go to 1
 1112 continue
      write(4,612)
      read(iin,*,end=130)(nalev(i),(ialev(i,j),j=1,nalev(i)),i=1,nadt)
      do 2112 i=1,nadt
        write(4,613)adt(i),(ialev(i,j),j=1,nalev(i))
 2112 continue
      go to 1
c
c --- parameters for array plotting options
c
  113 continue
      read(iin,*,end=130)npdt,(pdt(i),i=1,npdt)
      if(npdt.gt.0)go to 1113
      write(4,1614)
      go to 1
 1113 continue
      write(4,614)
      read(iin,*,end=130)(nplev(i),(iplev(i,j),j=1,nplev(i)),i=1,npdt)
      do 2113 i=1,npdt
        write(4,613)pdt(i),(iplev(i,j),j=1,nplev(i))
 2113 continue
      go to 1
c
c --- rdt = time interval at which to print restart information
c
  114 continue
      read(iin,*,end=130)rdt
      write(4,616)rdt
      go to 1
c
  115 continue
      read(iin,*,end=130)theta
      write(4,617)theta
      theta=theta*(3.1415926/180.)
      goto 1
c
  116 continue
      read(iin,*)rlat0,rlng0,t0,v0,dhor,ht,time0,buoysur
      go to 1
c
  117 continue
      read(iin,*)tsurf, (xc(i),i=1,kzp1)
      write(4,623)tsurf
      write(4,1236)
     1     (xc(i),i=1,kzp1)
 1236 format(5(e13.6,',  '),/)
      go to 1
  118 continue
      read(iin,*)salsur, (eta(i),i=1,kzp1)
      write(4,624)salsur
      write(4,1237)
     1    (xc(i),i=1,kzp1)
 1237 format(5(e13.6,',  '),/)
      go to 1
c
c
  120 continue
      write(4,121)
  121 format(' eof encountered reading icard')
      stop
c
  130 continue
      write(4,131)
  131 format(' eof encountered reading input variable')
      stop
c
c --- come here after reading input parameters
c
  132 continue
      if(rdt.le.(tmax-tstart+dt/2))then
c       fname=ident(:iflen)//'.wrs'
c       open(unit=iowr,file=fname,form='unformatted',status='new')
      endif
c
c
      if(t0.eq.0.0)go to 136
      write(4,625)rlat0,rlng0,t0,v0,dhor,ht,time0,buoysur
c
      pi=acos(-1.0)
      f0star=4*pi*sin(rlat0*pi/180.)/86164.
      betstr=4*pi*cos(rlat0*pi/180.)/86164./6356912.
      r0=1.0/(t0*f0star)
      astar=t0*v0/dhor
      bstar=t0*betstr*dhor
      write(4,626)f0star,betstr,r0,astar,bstar
c
	icstart=1./dt
	icstart=icstart*iday
c	print *,'icstart=',icstart,' icalc=',icalc
c
      dt=dt*86400./t0
      tstart=(tstart-time0)*86400./t0
      tmax=(tmax-time0)*86400./t0
c
      do 135 i=1,4
        ddt(i)=ddt(i)*86400./t0
        adt(i)=adt(i)*86400./t0
        pdt(i)=pdt(i)*86400./t0
  135 continue
c
      rdt=rdt*86400./t0
      go to 140
c
c --- translation for printed output if not dimensional input
c
  136 continue
      time0=0.0
      t0=86400.0
      ht=1.0
      dhor=1.0
      v0=1.0
c
  140 continue
c
c --- set up the stratification and associated matrices and constants
c
      call set3d(kz,ifdiff)
c
c --- constants for common/parm/
c
      mm1=m-1
      mm2=m-2
      mp1=m+1
      msq=m*n
c
      nm1=n-1
      nm2=n-2
      np1=n+1
c
      hy=xbasin/float(m-1)
      rea=hy*float(nm1*mm1)
      kzm1=kz-1
      kzp2=kz+2
c
c......read in the observational data......
c
c	call getob
c
c
c.......set up mass matrices amx and amy and find inverses......
c
        do 2030 j=1,m-1
        do 2030 i=1,m-1
 2030   amx(i,j)=0.0
        do 2031 j=1,n
        do 2031 i=1,n
 2031   amy(i,j)=0.0
c
        do 2020 i=1,m-1
 2020   amx(i,i)=4.0
        do 2021 i=1,m-2
 2021   amx(i,i+1)=1.0
        do 2022 i=2,m-1
 2022   amx(i,i-1)=1.0
c
        amx(1,m-1)=1.0
        amx(m-1,1)=1.0
c
        do 2023 i=1,n
 2023   amy(i,i)=4.0
        do 2024 i=1,n-1
 2024   amy(i,i+1)=1.0
        do 2025 i=2,n
 2025   amy(i,i-1)=1.0
        amy(1,1)=2.0
        amy(n,n)=2.0
c
c       do 2026 i=1,m-1
c       write(6,2027)(amx(i,j),j=1,m-1)
c2027   format(8f10.7)
c2026   continue
c
c       do 2028 i=1,n
c2028   write(6,2027)(amy(i,j),j=1,n)
c
        call sgeco(amx,m0-1,m0-1,ipvt,rcond,z)
c
        print *,'rcond for amx =',rcond
c
        call sgedi(amx,m0-1,m0-1,ipvt,det,work,01)
c
        call sgeco(amy,n0,n0,ipvty,rcond,zy)
c
        print *,'rcond for amy =',rcond
c
        call sgedi(amy,n0,n0,ipvty,det,worky,01)
c
c --- bottom topography and surface density variation is only
c --- implemented for finite difference method
c
      if(iftop.eq.0 .and. ifbot.eq.0)goto 146
      if(ifdiff.eq.0)goto 142
      do 141 i=1,kz
        topfct(i)=0.
        botfct(i)=0.
  141 continue
      topfct(1)=1./hz(1)
      botfct(kz)=1./hz(kz)
      goto 146
  142 continue
      write(4,143)
  143 format('  non-uniform surface or bottom density implemented',
     +  /,'  only for finite differences.')
      stop
  146 continue
c
c --- level printing and plotting information. allows for printing
c --- and plotting of density fields.
c
      do 151 na=1,nadt
        do 148 i=1,kzp2
          iprnt(i)=0
  148   continue
        do 149 i=1,nalev(na)
          lp=ialev(na,i)
          if(lp.eq.kzp1)lp=kzp2
          if(lp.eq.0)lp=kzp1
          iprnt(lp)=1
  149   continue
        do 150 i=1,kzp2
          ialev(na,i)=iprnt(i)
  150   continue
  151 continue
c
      do 155 np=1,npdt
        do 152 i=1,kzp2
          iplot(i)=0
  152   continue
        do 153 i=1,nplev(np)
          lp=iplev(np,i)
          if(lp.eq.kzp1)lp=kzp2
          if(lp.eq.0)lp=kzp1
          iplot(lp)=1
  153   continue
        do 154 i=1,kzp2
          iplev(np,i)=iplot(i)
  154   continue
  155 continue
c
c --- array 'mast'
c
      mast(1)=1
      mast(2)=m
      mast(3)=-1
      mast(4)=-m
c
c --- array 'nast'
c
      nast(1)=m
      nast(2)=-1
      nast(3)=-m
      nast(4)=1
c
c --- array 'icorn'
c
      icorn(1)=1
      icorn(2)=m
      icorn(3)=msq
      icorn(4)=msq-m+1
c
c --- array 'ncorn'
c
      ncorn(1)=m+2
      ncorn(2)=2*m-1
      ncorn(3)=msq-m-1
      ncorn(4)=msq-2*m+2
c
c --- coordinate values for 'latitude'
c
C>>>>>>>>>>>>>
	ir1=1
	ir2=m
	jr1=1
	jr2=n
	ybasin=hy*float(n-1)
	xxpc=0.5*hy*float(mbig-1)
	yypc=0.5*hy*float(nbig-1)
	xcl=xxpc-hy*float(ir1-1)
	ycl=yypc-hy*float(jr1-1)
      do 159 j=1,n
        yy(j)=hy*float(j-1)-ycl
  159 continue
	do 1590 i=1,m
	xx(i)=hy*float(i-1)-xcl
 1590   continue
c
c --- define the meridional coordinate under model rotation
c
      ip=0
      do 1592 j=1,n
        do 1591 i=1,m
          ip=ip+1
          pvort(ip)=(yy(j)*cos(theta)+xx(i)*sin(theta))*beta
 1591   continue
 1592 continue
C>>>>>>>>>>>>>>>
c
c
c --- constants for common/solver/
c
      vbm(1)=4.0
      vbn(1)=2.0
      do 160 i=3,mm1
        vam(i-1)=1./vbm(i-2)
        vbm(i-1)=4.-vam(i-1)
  160 continue
      vam(mm1)=1./vbm(mm1-1)
      vbm(mm1)=vbm(1)-vam(mm1)
      do 165 i=3,n
        van(i-1)=1./vbn(i-2)
        vbn(i-1)=4.-van(i-1)
  165 continue
      van(n)=1./vbn(n-1)
      vbn(n)=vbn(1)-van(n)
c
c.....determine the maximum timestep...
c
      maxt=(tmax-tstart)/dt
	print *,'maxt=',maxt
      maxt1=maxt
c
c
c
      icalc=0
c
        call cjet(u0,bstm,bvrt)
c
       call thermvt(bstm,therm,kz,msq)
c
        do 7544 k=1,kz
        do 7544 ip=1,msq
 7544   bvrt(ip,k)=bvrt(ip,k)+therm(ip,k)
c
	print *,'bstm=',(bstm(i,1),i=1,msq-m,m)
c	print *,'bvrt=',bvrt
c
c
      npnt=(m-1)*n*kz
        print *,'npnt=',npnt
      nev = 6
c      nev = maxnev
      ncv=maxncv
      mnaupd=3
c
        if(nev.gt.maxnev)then
        print *,'nev exceeds maxnev!!!'
        stop
        endif
        if(ncv.gt.maxncv)then
        print *,'ncv exceeds maxncv!!!'
        stop
        endif
c
      bmat = 'I'
      which = 'LM'
      lworkl=3*ncv**2+6*ncv
c      tol = 0.0
      tol = 1.e-6
      info = 0
      ishfts = 1
      nb     = 1
      mode= 1
      maxitr= 250
      iparam(1) = ishfts
      iparam(2) = 1
      iparam(3) = maxitr
      iparam(4) = nb
c      iparam(6) = 1
      iparam(7)= mode
      ido = 0
c
        nmit=0
c
c
 1010 continue
c
      call snaupd ( ido, bmat, npnt, which, nev, tol,
     &              resid, ncv, vss, ldv,
     &              iparam, ipntr, workd, workl, lworkl, info )
c
       print *,'snaupd done: nev=',nev,' nconv=',iparam(5)
c
      if (ido .eq. -1 .or. ido .eq. 1) then
c
        call matveca (npnt, workd(ipntr(1)), workd(ipntr(2)))
c
c       call solvem (npnt,tmpwrk,workd(ipntr(2)))
c
c       else if (ido.eq.2)then
c
c       call matvecm(npnt,workd(ipntr(1)),workd(ipntr(2)))
c
        goto 1010
c
       endif
c
        if(info.lt.0)then
c
         print *,'Error with snaupd, info=',info
         print *,'check documentation of snaupd'
c
        else
c
        rvec=.true.
c
        call sneupd(rvec,'A',select,d,d(1,2),vss,ldv,
     &   sigmar,sigmai,workev,bmat,npnt,which,nev,tol,
     &   resid,ncv,vss,ldv,iparam,ipntr,workd,workl,
     &   lworkl,ierr)
c
          if(ierr.ne.0)then
c
            print *,'Error in sneupd, info=',ierr
            print *,'Check documentation for sneupd'
c
          else
c
             first=.true.
             nconv=iparam(5)
c
c             write(35,*)'Ritz values and residuals'
c
             write(23)nconv,u0,beta
             write(23)bstm
             write(23)bvrt
c
             do j=1,nconv
c
               if(d(j,2).eq.0.0)then
c
c..........pure real case.....
c
                 nmit=1
                 write(23)d(j,1),d(j,2),d(j,3)
                 call matveca(npnt,vss(1,j),tmpwrk)
                 nmit=0
                 call saxpy(npnt,-d(j,1),vss(1,j),1,tmpwrk,1)
                 d(j,3)=snrm2(npnt,tmpwrk,1)
	         write(43)d(j,1),d(j,2),d(j,3)
	         write(43)(vss(i,j),i=1,npnt)
c
                else if (first) then
c
c..........complex case......
c
                 nmit=1
                 write(23)d(j,1),d(j,2),d(j,3)
                 call matveca(npnt,vss(1,j),tmpwrk)
                 call saxpy(npnt,-d(j,1),vss(1,j),1,tmpwrk,1)
                 call saxpy(npnt,d(j,2),vss(1,j+1),1,tmpwrk,1)
                 d(j,3)=snrm2(npnt,tmpwrk,1)
                 nmit=1
                 write(23)d(j+1,1),d(j+1,2),d(j+1,3)
                 call matveca(npnt,vss(1,j+1),tmpwrk)
                 call saxpy(npnt,-d(j,2),vss(1,j),1,tmpwrk,1)
                 call saxpy(npnt,-d(j,1),vss(1,j+1),1,tmpwrk,1)
                 d(j,3)=slapy2(d(j,3),snrm2(npnt,tmpwrk,1))
                 d(j+1,3)=d(j,3)
	         write(43)d(j,1),d(j,2),d(j,3)
	         write(43)(vss(i,j),i=1,npnt)
	         write(43)d(j+1,1),d(j+1,2),d(j+1,3)
	         write(43)(vss(i,j+1),i=1,npnt)
                 first=.false.
                else
                 first=.true.
                endif
c
c                write(35,*)d(j,1),d(j,2),d(j,3)
c
            enddo
c
          call smout(6,nconv,3,d,maxncv,-6,
     &      'Ritz values (Real,Imag) and direct residuals')
c
c
        endif
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit
     &                  Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
c
         print *, ' '
         print *, '_NDRV1 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', npnt
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
c
      end if
c
	print *,'npnt=',npnt
c
c
c --- input format statements for data cards
c
  500 format(40i2)
  501 format(20a4)
  502 format(8f10.0,/,8f10.0)
  503 format(i5,4f10.0)
  504 format(11f8.0)
  505 format(16i5)
  506 format(f10.2)
  507 format(12f10.0)
c
c --- output format statements for data cards
c
  600 format(' ',6a4)
  601 format(/,20a4)
  602 format(/'ratio of characteristic to advective time scales=',f8.4,
     + //,'ratio of characteristic to rossby wave time scales=',f8.4,//,
     +    'bottom friction=',f8.4)
  603 format(/'number of horizontal grid points, m=',i3,/,
     +        ' number of merdional grid points, n=',i3,/,
     +        '      number of vertical levels, kz=',i3)
  604 format(/'no of layers= ',i3,'  layer depths : ')

  605 format(/'stratification scale = ',1pe11.4,/,' sigma values:')

  606 format(/,'shapiro filter:  order=',i3,5x,'frequency=',i3,
     +       5x,'every ',i3,' time steps.')
  607 format(/,'array pram:  ',4(1pg13.6,1x),/,13x,4(g13.6,1x),/,
     +       13x,4(g13.6,1x),/,13x,4(g13.6,1x))
  608 format(/'dt= ',f8.5,'   tstart= ',f12.5,'   tmax= ',f12.5,
     +         //'xbasin= ',f8.5)
  609 format(/,'switches:                 diff/colloc: ',i1,
     +       /,'                               ifpert: ',i1,
     +       /,'                                ifrst: ',i1,
     +       /,'              non-uniform top density: ',i1,
     +       /,'           non-uniform bottom density: ',i1,
     +       /,'                top vertical velocity: ',i1,
     +       /,'                         bottom slope: ',i1,
     +       /,'                         writing data: ',i1,
     +       /,'               surface boundary layer: ',i1,
     +       /,'                     eva preprocessor: ',i1)
  610 format(/'diag printing:'//,4x,'interval')
  611 format(3x,f10.6)
 1611 format(/'no diagnostic printing')
  612 format(/'array printing:'//,4x,'interval',6x,'levels')
 1612 format(/'no array printing')
  613 format(3x,f10.6,4x,9i3)
  614 format(/'array plotting:'//,4x,'interval',6x,'levels')
 1614 format(/'no array plotting')
  615 format(/'       plotting interval    no of levels   ',40i3)
  616 format(/'interval at which to write restart information:',1pe10.3)
  617 format(/' domain is rotated',f12.3,' degrees cck')
  618 format(/,9a4,/,' ',20a4)
  619 format(/' dt specified by input and specified by restart data',
     +     /,'  file differ by more than .01%.  odt,dt:',1pe15.8,', ',
     +       e15.8)
  621 format(/' no complete sets of data.')
  622 format(/' unable to find set of data, whose time is within dt',
     + /,'  of the time specified.  tstart, tx:',1pe15.8,', ',e15.8)
  623 format(/' vertical temperature profile ',//,
     +       ' surface temperature: ',1pe13.6,' deg c.',
     +      /' vertical derivative of temperature: ( deg c. / meter)')
  624 format(/' vertical salinity profile ',//,
     +       ' surface salinity: ',1pe13.6,' .001',
     +      /' vertical derivative of salinity: ( .001 / meter)')
  625 format(/' **** run includes dimensional external parameters ***',
     +     //,' the units of time specifed in the input data file are',
     +     /,'  in days.  ',
     +     /,'  warning : the input alpha and  beta should be checked',
     +     /,'  with those computed from the dimensional input.',//,
     +     //,' center grid latitude, longitude: ',2(f12.4),
     +             '(deg north and west)',
     +     //,' t0 ',e14.6,' s, v0',e14.6,' m/s,',/,
     +        '  l ',e14.6,' m,  h',e14.6,' m',
     +     //,' time offset ',0pf10.2,' (day), surface buoyancy '
     +         ,1pe14.6,
     +              ' g/cm/cm/cm')
  626 format(/' parameters determined from the dimensional input',/,
     +       /'     f0* ',1pe14.6,'(/s), beta* ',e14.6,'(/s/m)',
     +       /'     r0* ',e14.6,
     +       /'     a*  ',e14.6,'         b* ',e14.6)
c
c --- format statements for array printing
c
  700 format(/'level ',i2)
  701 format(/,30x,'streamfunction',5x,5htime=,f10.4,
     +' levels',9(1x,i2))
  702 format(/,30x,'vorticity',5x,5htime=,f10.4,' levels',9(1x,i2))
  703 format(/,30x,'perturbation streamfunction',5x,5htime=,f10.4,
     +' levels',9(1x,i2))
  704 format(/,30x,'perturbation vorticity',5x,5htime=,f10.4,
     +' levels',9(1x,i2))
  705 format(/,30x,'top density',8x,5htime=,f10.4)
  706 format(/,30x,'bottom density',5x,5htime=,f10.4)
  707 format('scale = ',1pe6.0)
c
c --- formats for diagnostic printing.
c
  800 format(/'level ',
     +             40hicalc   t    i(z)      i(q)    i(nrg)   ,
     +             36hmax(v)   di(nrg)  maxd(s)  rmsd(s)  ,
     +             43hrmsd(z)  std(s)   ratio   abs corr ten corr)
  801 format(6a4,9a4,/,' ',20a4,/,' ',20a4,//,
     !                                ' diagnostics for the full field')
  802 format(6a4,9a4,/,' ',20a4,/,' ',20a4,//,
     !                            ' diagnostics for the interior field')
c
      stop
      end
c
      subroutine arhs(s,prv)
c
c --- this subroutine calculates the right-hand-side (rhs) of
c --- the finite-element vorticity equation
c
          include 'prob2.pert'
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      dimension s(1),prv(1)
c
      do ip=1,msq
      fq=-0.5*s(ip)
      s(ip)=1.5*s(ip)+prv(ip)
      prv(ip)=fq
      enddo
c
      fac=-36*dt/(hy**2)
      do 10 ip=1,msq
        s(ip)=fac*s(ip)
   10 continue
c
      return
      end
c
c
      subroutine arhs2(s,prv)
c
c --- this subroutine calculates the right-hand-side (rhs) of
c --- the finite-element vorticity equation
c
          include 'prob2.pert'
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      dimension s(1),prv(1)
c
      do ip=1,msq
      s(ip)=s(ip)+prv(ip)
      enddo
c
      fac=-36*dt/(hy**2)
      do 10 ip=1,msq
        s(ip)=fac*s(ip)
        prv(ip)=0.0
   10 continue
c
      return
      end
c
c
c
      subroutine adjacv(psi,q1,r,m,n)
      dimension psi(m,n),q(m,n),r(m,n),q1(m,n)
c
c --- subroutine to compute the arakawa jacobian j(psi,q)=r
c
      do j=1,n
      do i=1,m
      q(i,j)=0.0
      enddo
      enddo
c
c --- interior
c
      m1=m-1
      n1=n-1
      do 15 i=1,m
        do 10 j=2,n1
          im1=i-1
          ip1=i+1
          if(i.eq.1)im1=m-1
          if(i.eq.m)ip1=2
c
          q(im1,j-1)=q(im1,j-1)+
     &          r(i,j)*(psi(im1,j)-psi(i,j-1))
c
          q(i,j-1)=q(i,j-1)+
     &     r(i,j)*(psi(im1,j-1)+psi(im1,j)-psi(ip1,j-1)-psi(ip1,j))
c
          q(ip1,j-1)=q(ip1,j-1)+
     &      r(i,j)*(psi(i,j-1)-psi(ip1,j))
c
          q(im1,j)=q(im1,j)+
     &      r(i,j)*(psi(im1,j+1)+psi(i,j+1)-psi(im1,j-1)-psi(i,j-1))
c
          q(ip1,j)=q(ip1,j)+
     &      r(i,j)*(psi(i,j-1)+psi(ip1,j-1)-psi(i,j+1)-psi(ip1,j+1))
c
          q(im1,j+1)=q(im1,j+1)+
     &      r(i,j)*(psi(i,j+1)-psi(im1,j))
c
          q(i,j+1)=q(i,j+1)+
     &      r(i,j)*(psi(ip1,j)+psi(ip1,j+1)-psi(im1,j)-psi(im1,j+1))
c
          q(ip1,j+1)=q(ip1,j+1)+
     &      r(i,j)*(psi(ip1,j)-psi(i,j+1))
c
   10   continue
   15 continue
c
c --- left edge
c
      do 20 i=1,m
        im1=i-1
        ip1=i+1
        if(i.eq.1)im1=m-1
        if(i.eq.m)ip1=2
c
        q(im1,1)=q(im1,1)+
     &   r(i,1)*(psi(im1,2)+psi(i,2)-2.*psi(i,1))
c
        q(i,1)=q(i,1)+
     &     r(i,1)*(2.*(psi(im1,1)-psi(ip1,1)))
c
        q(ip1,1)=q(ip1,1)+
     &     r(i,1)*(2.*psi(i,1)-psi(i,2)-psi(ip1,2))
c
        q(im1,2)=q(im1,2)+
     &     r(i,1)*(psi(i,2)-psi(im1,1))
c
        q(i,2)=q(i,2)+
     &     r(i,1)*(psi(ip1,1)+psi(ip1,2)-psi(im1,1)-psi(im1,2))
c
        q(ip1,2)=q(ip1,2)+
     &     r(i,1)*(psi(ip1,1)-psi(i,2))
c
   20 continue
c
c --- bottom
c
c      do 30 i=2,n1
c        im1=i-1
c        ip1=i+1
c
c        q(m,im1)=q(m,im1)+
c     &    r(m,i)*(psi(m-1,im1)+psi(m-1,i)-2.*psi(m,i))
c
c        q(m,i)=q(m,i)+
c     &     r(m,i)*(2.*(psi(m,im1)-psi(m,ip1)))
c
c        q(m,ip1)=q(m,ip1)+
c     &     r(m,i)*(2.*psi(m,i)-psi(m-1,i)-psi(m-1,ip1))
c
c        q(m-1,im1)=q(m-1,im1)+
c     &     r(m,i)*(psi(m-1,i)-psi(m,im1))
c
c        q(m-1,i)=q(m-1,i)+
c     &     r(m,i)*(psi(m,ip1)+psi(m-1,ip1)-psi(m,im1)-psi(m-1,im1))
c
c        q(m-1,ip1)=q(m-1,ip1)+
c     &     r(m,i)*(psi(m,ip1)-psi(m-1,i))
c
c   30 continue
c
c --- right edge
c
      do 40 i=1,m
        im1=i-1
        ip1=i+1
        if(i.eq.1)im1=m-1
        if(i.eq.m)ip1=2
c
        q(ip1,n)=q(ip1,n)+
     &    r(i,n)*(psi(ip1,n1)+psi(i,n1)-2.*psi(i,n))
c
        q(i,n)=q(i,n)+
     &     r(i,n)*(2.*(psi(ip1,n)-psi(im1,n)))
c
        q(im1,n)=q(im1,n)+
     &     r(i,n)*(2.*psi(i,n)-psi(i,n1)-psi(im1,n1))
c
        q(ip1,n1)=q(ip1,n1)+
     &     r(i,n)*(psi(i,n1)-psi(ip1,n))
c
        q(i,n1)=q(i,n1)+
     &     r(i,n)*(psi(im1,n)+psi(im1,n1)-psi(ip1,n)-psi(ip1,n1))
c
        q(im1,n1)=q(im1,n1)+
     &     r(i,n)*(psi(im1,n)-psi(i,n1))
c
   40 continue
c
c
c --- top
c
c      do 50 i=2,n1
c        im1=i-1
c        ip1=i+1
c
c        q(1,ip1)=q(1,ip1)+
c     &    r(1,i)*(psi(2,ip1)+psi(2,i)-2.*psi(1,i))
c
c        q(1,i)=q(1,i)+
c     &     r(1,i)*(2.*(psi(1,ip1)-psi(1,im1)))
c
c        q(1,im1)=q(1,im1)+
c     &     r(1,i)*(2.*psi(1,i)-psi(2,i)-psi(2,im1))
c
c        q(2,ip1)=q(2,ip1)+
c     &     r(1,i)*(psi(2,i)-psi(1,ip1))
c
c        q(2,i)=q(2,i)+
c     &     r(1,i)*(psi(1,im1)+psi(2,im1)-psi(1,ip1)-psi(2,ip1))
c
c        q(2,im1)=q(2,im1)+
c     &     r(1,i)*(psi(1,im1)-psi(2,i))
c
c   50 continue
c
c --- corners
c
c
c      q(1,1)=q(1,1)+
c     &   r(1,1)*(2.*(psi(1,2)-psi(2,1)))
c
c      q(2,1)=q(2,1)+
c     &   r(1,1)*(2.*psi(1,1)-psi(2,2)-psi(1,2))
c
c      q(2,2)=q(2,2)+
c     &   r(1,1)*(psi(2,1)-psi(1,2))
c
c      q(1,2)=q(1,2)+
c     &   r(1,1)*(psi(2,1)+psi(2,2)-2.*psi(1,1))
c
c
c
c
c      q(m,1)=q(m,1)+
c     &   r(m,1)*(2.*(psi(m-1,1)-psi(m,2)))
c
c      q(m,2)=q(m,2)+
c     &   r(m,1)*(2.*psi(m,1)-psi(m-1,2)-psi(m-1,1))
c
c      q(m-1,2)=q(m-1,2)+
c     &   r(m,1)*(psi(m,2)-psi(m-1,1))
c
c      q(m-1,1)=q(m-1,1)+
c     &   r(m,1)*(psi(m,2)+psi(m-1,2)-2.*psi(m,1))
c
c
c
c      q(m,n)=q(m,n)+
c     &   r(m,n)*(2.*(psi(m,n1)-psi(m-1,n)))
c
c      q(m-1,n)=q(m-1,n)+
c     &   r(m,n)*(2.*psi(m,n)-psi(m-1,n1)-psi(m,n1))
c
c      q(m-1,n1)=q(m-1,n1)+
c     &   r(m,n)*(psi(m-1,n)-psi(m,n1))
c
c      q(m,n1)=q(m,n1)+
c     &   r(m,n)*(psi(m-1,n)+psi(m-1,n1)-2.*psi(m,n))
c
c
c
c      q(1,n)=q(1,n)+
c     &  r(1,n)*(2.*(psi(2,n)-psi(1,n1)))
c
c      q(1,n1)=q(1,n1)+
c     &   r(1,n)*(2.*psi(1,n)-psi(2,n1)-psi(2,n))
c
c      q(2,n1)=q(2,n1)+
c     &   r(1,n)*(psi(1,n1)-psi(2,n))
c
c      q(2,n)=q(2,n)+
c     &   r(1,n)*(psi(1,n1)+psi(2,n1)-2.*psi(1,n))
c
c
      do j=1,n
      do i=1,m
      q1(i,j)=q1(i,j)+q(i,j)/12.
      enddo
      enddo
c
      return
      end
c
c
      subroutine adjacs(psi1,q,r,m,n)
      dimension psi(m,n),q(m,n),r(m,n),psi1(m,n)
c
c --- subroutine to compute the arakawa jacobian j(psi,q)=r
c
      do j=1,n
      do i=1,m
      psi(i,j)=0.0
      enddo
      enddo
c
c --- interior
c
      m1=m-1
      n1=n-1
      do 15 i=1,m
        do 10 j=2,n1
          im1=i-1
          ip1=i+1
          if(i.eq.1)im1=m-1
          if(i.eq.m)ip1=2
c
          psi(im1,j)=psi(im1,j)+q(im1,j-1)*r(i,j)
          psi(i,j-1)=psi(i,j-1)-q(im1,j-1)*r(i,j)
c
          psi(im1,j-1)=psi(im1,j-1)+q(i,j-1)*r(i,j)
          psi(im1,j)=psi(im1,j)+q(i,j-1)*r(i,j)
          psi(ip1,j-1)=psi(ip1,j-1)-q(i,j-1)*r(i,j)
          psi(ip1,j)=psi(ip1,j)-q(i,j-1)*r(i,j)
c
          psi(i,j-1)=psi(i,j-1)+q(ip1,j-1)*r(i,j)
          psi(ip1,j)=psi(ip1,j)-q(ip1,j-1)*r(i,j)
c
          psi(im1,j+1)=psi(im1,j+1)+q(im1,j)*r(i,j)
          psi(i,j+1)=psi(i,j+1)+q(im1,j)*r(i,j)
          psi(im1,j-1)=psi(im1,j-1)-q(im1,j)*r(i,j)
          psi(i,j-1)=psi(i,j-1)-q(im1,j)*r(i,j)
c
          psi(i,j-1)=psi(i,j-1)+q(ip1,j)*r(i,j)
          psi(ip1,j-1)=psi(ip1,j-1)+q(ip1,j)*r(i,j)
          psi(i,j+1)=psi(i,j+1)-q(ip1,j)*r(i,j)
          psi(ip1,j+1)=psi(ip1,j+1)-q(ip1,j)*r(i,j)
c
          psi(i,j+1)=psi(i,j+1)+q(im1,j+1)*r(i,j)
          psi(im1,j)=psi(im1,j)-q(im1,j+1)*r(i,j)
c
          psi(ip1,j)=psi(ip1,j)+q(i,j+1)*r(i,j)
          psi(ip1,j+1)=psi(ip1,j+1)+q(i,j+1)*r(i,j)
          psi(im1,j)=psi(im1,j)-q(i,j+1)*r(i,j)
          psi(im1,j+1)=psi(im1,j+1)-q(i,j+1)*r(i,j)
c
          psi(ip1,j)=psi(ip1,j)+q(ip1,j+1)*r(i,j)
          psi(i,j+1)=psi(i,j+1)-q(ip1,j+1)*r(i,j)
c
   10   continue
   15 continue
c
c --- left edge
c
      do 20 i=1,m
        im1=i-1
        ip1=i+1
        if(i.eq.1)im1=m-1
        if(i.eq.m)ip1=2
c
        psi(im1,2)=psi(im1,2)+q(im1,1)*r(i,1)
        psi(i,2)=psi(i,2)+q(im1,1)*r(i,1)
        psi(i,1)=psi(i,1)-2.*q(im1,1)*r(i,1)
c
        psi(im1,1)=psi(im1,1)+2.*q(i,1)*r(i,1)
        psi(ip1,1)=psi(ip1,1)-2.*q(i,1)*r(i,1)
c
        psi(i,1)=psi(i,1)+2.*q(ip1,1)*r(i,1)
        psi(i,2)=psi(i,2)-q(ip1,1)*r(i,1)
        psi(ip1,2)=psi(ip1,2)-q(ip1,1)*r(i,1)
c
        psi(i,2)=psi(i,2)+q(im1,2)*r(i,1)
        psi(im1,1)=psi(im1,1)-q(im1,2)*r(i,1)
c
        psi(ip1,1)=psi(ip1,1)+q(i,2)*r(i,1)
        psi(ip1,2)=psi(ip1,2)+q(i,2)*r(i,1)
        psi(im1,1)=psi(im1,1)-q(i,2)*r(i,1)
        psi(im1,2)=psi(im1,2)-q(i,2)*r(i,1)
c
        psi(ip1,1)=psi(ip1,1)+q(ip1,2)*r(i,1)
        psi(i,2)=psi(i,2)-q(ip1,2)*r(i,1)
c
   20 continue
c
c --- bottom
c
c      do 30 i=2,n1
c        im1=i-1
c	ip1=i+1
c
c        psi(m-1,im1)=psi(m-1,im1)+q(m,im1)*r(m,i)
c        psi(m-1,i)=psi(m-1,i)+q(m,im1)*r(m,i)
c        psi(m,i)=psi(m,i)-2.*q(m,im1)*r(m,i)
c
c        psi(m,im1)=psi(m,im1)+2.*q(m,i)*r(m,i)
c        psi(m,ip1)=psi(m,ip1)-2.*q(m,i)*r(m,i)
c
c        psi(m,i)=psi(m,i)+2.*q(m,ip1)*r(m,i)
c        psi(m-1,i)=psi(m-1,i)-q(m,ip1)*r(m,i)
c        psi(m-1,ip1)=psi(m-1,ip1)-q(m,ip1)*r(m,i)
c
c        psi(m-1,i)=psi(m-1,i)+q(m-1,im1)*r(m,i)
c        psi(m,im1)=psi(m,im1)-q(m-1,im1)*r(m,i)
c
c        psi(m,ip1)=psi(m,ip1)+q(m-1,i)*r(m,i)
c        psi(m-1,ip1)=psi(m-1,ip1)+q(m-1,i)*r(m,i)
c        psi(m,im1)=psi(m,im1)-q(m-1,i)*r(m,i)
c        psi(m-1,im1)=psi(m-1,im1)-q(m-1,i)*r(m,i)
c
c        psi(m,ip1)=psi(m,ip1)+q(m-1,ip1)*r(m,i)
c        psi(m-1,i)=psi(m-1,i)-q(m-1,ip1)*r(m,i)
c
c
c   30 continue
c
c --- right edge
c
      do 40 i=1,m
       im1=i-1
       ip1=i+1
       if(i.eq.1)im1=m-1
       if(i.eq.m)ip1=2
c
        psi(ip1,n1)=psi(ip1,n1)+q(ip1,n)*r(i,n)
        psi(i,n1)=psi(i,n1)+q(ip1,n)*r(i,n)
        psi(i,n)=psi(i,n)-2.*q(ip1,n)*r(i,n)
c
        psi(ip1,n)=psi(ip1,n)+2.*q(i,n)*r(i,n)
        psi(im1,n)=psi(im1,n)-2.*q(i,n)*r(i,n)
c
        psi(i,n)=psi(i,n)+2.*q(im1,n)*r(i,n)
        psi(i,n1)=psi(i,n1)-q(im1,n)*r(i,n)
        psi(im1,n1)=psi(im1,n1)-q(im1,n)*r(i,n)
c
        psi(i,n1)=psi(i,n1)+q(ip1,n1)*r(i,n)
        psi(ip1,n)=psi(ip1,n)-q(ip1,n1)*r(i,n)
c
        psi(im1,n)=psi(im1,n)+q(i,n1)*r(i,n)
        psi(im1,n1)=psi(im1,n1)+q(i,n1)*r(i,n)
        psi(ip1,n)=psi(ip1,n)-q(i,n1)*r(i,n)
        psi(ip1,n1)=psi(ip1,n1)-q(i,n1)*r(i,n)
c
        psi(im1,n)=psi(im1,n)+q(im1,n1)*r(i,n)
        psi(i,n1)=psi(i,n1)-q(im1,n1)*r(i,n)
c
   40 continue
c
c --- top
c
c      do 50 i=2,n1
c       im1=i-1
c       ip1=i+1
c
c        psi(2,ip1)=psi(2,ip1)+q(1,ip1)*r(1,i)
c        psi(2,i)=psi(2,i)+q(1,ip1)*r(1,i)
c        psi(1,i)=psi(1,i)-2.*q(1,ip1)*r(1,i)
c
c        psi(1,ip1)=psi(1,ip1)+2.*q(1,i)*r(1,i)
c        psi(1,im1)=psi(1,im1)-2.*q(1,i)*r(1,i)
c
c        psi(1,i)=psi(1,i)+2.*q(1,im1)*r(1,i)
c        psi(2,i)=psi(2,i)-q(1,im1)*r(1,i)
c        psi(2,im1)=psi(2,im1)-q(1,im1)*r(1,i)
c
c        psi(2,i)=psi(2,i)+q(2,ip1)*r(1,i)
c        psi(1,ip1)=psi(1,ip1)-q(2,ip1)*r(1,i)
c
c        psi(1,im1)=psi(1,im1)+q(2,i)*r(1,i)
c        psi(2,im1)=psi(2,im1)+q(2,i)*r(1,i)
c        psi(1,ip1)=psi(1,ip1)-q(2,i)*r(1,i)
c        psi(2,ip1)=psi(2,ip1)-q(2,i)*r(1,i)
c
c        psi(1,im1)=psi(1,im1)+q(2,im1)*r(1,i)
c        psi(2,i)=psi(2,i)-q(2,im1)*r(1,i)
c
c   50 continue
c
c
c --- corners
c
c
c      psi(1,2)=psi(1,2)+2.*q(1,1)*r(1,1)
c      psi(2,1)=psi(2,1)-2.*q(1,1)*r(1,1)
c
c      psi(1,1)=psi(1,1)+2.*q(2,1)*r(1,1)
c      psi(2,2)=psi(2,2)-q(2,1)*r(1,1)
c      psi(1,2)=psi(1,2)-q(2,1)*r(1,1)
c
c      psi(2,1)=psi(2,1)+q(2,2)*r(1,1)
c      psi(1,2)=psi(1,2)-q(2,2)*r(1,1)
c
c      psi(2,1)=psi(2,1)+q(1,2)*r(1,1)
c      psi(2,2)=psi(2,2)+q(1,2)*r(1,1)
c      psi(1,1)=psi(1,1)-2.*q(1,2)*r(1,1)
c
c
c      psi(m-1,1)=psi(m-1,1)+2.*q(m,1)*r(m,1)
c      psi(m,2)=psi(m,2)-2.*q(m,1)*r(m,1)
c
c      psi(m,1)=psi(m,1)+2.*q(m,2)*r(m,1)
c      psi(m-1,2)=psi(m-1,2)-q(m,2)*r(m,1)
c      psi(m-1,1)=psi(m-1,1)-q(m,2)*r(m,1)
c
c      psi(m,2)=psi(m,2)+q(m-1,2)*r(m,1)
c      psi(m-1,1)=psi(m-1,1)-q(m-1,2)*r(m,1)
c
c      psi(m,2)=psi(m,2)+q(m-1,1)*r(m,1)
c      psi(m-1,2)=psi(m-1,2)+q(m-1,1)*r(m,1)
c      psi(m,1)=psi(m,1)-2.*q(m-1,1)*r(m,1)
c
c
c
c      psi(m,n1)=psi(m,n1)+2.*q(m,n)*r(m,n)
c      psi(m-1,n)=psi(m-1,n)-2.*q(m,n)*r(m,n)
c
c      psi(m,n)=psi(m,n)+2.*q(m-1,n)*r(m,n)
c      psi(m-1,n1)=psi(m-1,n1)-q(m-1,n)*r(m,n)
c      psi(m,n1)=psi(m,n1)-q(m-1,n)*r(m,n)
c
c      psi(m-1,n)=psi(m-1,n)+q(m-1,n1)*r(m,n)
c      psi(m,n1)=psi(m,n1)-q(m-1,n1)*r(m,n)
c
c      psi(m-1,n)=psi(m-1,n)+q(m,n1)*r(m,n)
c      psi(m-1,n1)=psi(m-1,n1)+q(m,n1)*r(m,n)
c      psi(m,n)=psi(m,n)-2.*q(m,n1)*r(m,n)
c
c
c
c      psi(2,n)=psi(2,n)+2.*q(1,n)*r(1,n)
c      psi(1,n1)=psi(1,n1)-2.*q(1,n)*r(1,n)
c
c      psi(1,n)=psi(1,n)+2.*q(1,n1)*r(1,n)
c      psi(2,n1)=psi(2,n1)-q(1,n1)*r(1,n)
c      psi(2,n)=psi(2,n)-q(1,n1)*r(1,n)
c
c      psi(1,n1)=psi(1,n1)+q(2,n1)*r(1,n)
c      psi(2,n)=psi(2,n)-q(2,n1)*r(1,n)
c
c      psi(1,n1)=psi(1,n1)+q(2,n)*r(1,n)
c      psi(2,n1)=psi(2,n1)+q(2,n)*r(1,n)
c      psi(1,n)=psi(1,n)-2.*q(2,n)*r(1,n)
c
c
        do j=1,n
        do i=1,m
        psi1(i,j)=psi1(i,j)+psi(i,j)/12.
        enddo
        enddo
c
      return
      end
c
c
      subroutine shpiro(g,m,nord)
c
c --- g = portion of a larger array to be filtered
c --- m = length of this array
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter it using subroutine filter
c
      dimension g(1)
      mm1=m-1
      mm2=m-2
      do 100 kord=1,nord
        do 10 i=1,mm1
          g(i)=g(i+1)-g(i)
   10   continue
        g(m)=-2.*g(mm1)
        do 20 ii=1,mm2
          i=m-ii
          g(i)=g(i)-g(i-1)
   20   continue
        g(1)=2.*g(1)
  100 continue
      return
      end
c
c
      subroutine ashpro(g,m,nord)
c
c --- g = portion of a larger array to be filtered
c --- m = length of this array
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter it using subroutine filter
c
      dimension g(1)
      mm1=m-1
      mm2=m-2
      do 100 kord=1,nord
c
        g(1)=2.*g(1)
c
        do 20 ii=mm2,1,-1
          i=m-ii
          g(i-1)=g(i-1)-g(i)
   20   continue
c
        g(mm1)=g(mm1)-2.*g(m)
        g(m)=0.0
c
        do 10 i=mm1,1,-1
          g(i+1)=g(i+1)+g(i)
          g(i)=-g(i)
   10   continue
c
  100 continue
      return
      end
c
	subroutine scalar
c
        include 'prob2.pert'
c
	common/egy/peint,akeint
	common/l2n/al2s,al2v
	common/egl2/pl,akl,alls,allv
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
	common/swoi/iter
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	call energy
	call l2norm
c
        max=(tmax-tstart)/dt
	if(icalc.ne.0.and.icalc.ne.max)goto 10
	pl=0.
	akl=0.
	alls=0.
	allv=0.
 10     continue
c
	write(89,1)icalc,peint,akeint,al2s,al2v
 1      format(i5,1p,4e10.3)
	if(icalc.eq.0.or.icalc.eq.max)return
	denom=(dt*(peint+pl))
	if(denom.gt.0.)grp=2.*(peint-pl)/denom
	denom=(dt*(akeint+akl))
	if(denom.gt.0.)grk=2.*(akeint-akl)/denom
	denom=(dt*(al2s+alls))
	if(denom.gt.0.)gr2s=2.*(al2s-alls)/denom
	denom=(dt*(al2v+allv))
	if(denom.gt.0.)gr2v=2.*(al2v-allv)/denom
	write(88,1)icalc,grp,grk,gr2s,gr2v
	pl=peint
	akl=akeint
	alls=al2s
	allv=al2v
c
	return
	end
c
c
c
      subroutine msolve
c
c --- this subroutine solves the adams-bashforth finite-element
c --- equations for the interior of the field
c
c --- array r (entering) = right-hand side of eq's
c --- array r (leaving) = solution
c
          include 'prob2.pert'
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      dimension aux(mx0)
      dimension va(mx0),vb(mx0)
c
      do 5 j=1,m
	va(j)=vam(j)
	vb(j)=vbm(j)
    5 continue
      do 40 j=2,nm1
        i1=(j-1)*m+2
        aux(1)=r(i1)
        do 10 i=2,mm2
          aux(i)=r(i1+i-1)-va(i)*aux(i-1)
   10   continue
        aux(mm2)=aux(mm2)/vb(mm2)
        do 20 i=2,mm2
          ii=mm1-i
          aux(ii)=(aux(ii)-aux(ii+1))/vb(ii)
   20   continue
        do 30 i=1,mm2
          r1(i1+i-1)=aux(i)
   30   continue
   40 continue
c
      do 45 j=1,n
	va(j)=van(j)
	vb(j)=vbn(j)
   45 continue
      do 80 i=2,mm1
        i1=m+i
        aux(1)=r1(i1)
        do 50 j=2,nm2
          aux(j)=r1(i1+(j-1)*m)-va(j)*aux(j-1)
   50   continue
        aux(nm2)=aux(nm2)/vb(nm2)
        do 60 j=2,nm2
          jj=nm1-j
          aux(jj)=(aux(jj)-aux(jj+1))/vb(jj)
   60   continue
        do 70 j=1,nm2
          ip=i1+(j-1)*m
          r(ip)=aux(j)
   70   continue
   80 continue
      return
      end
c
	subroutine adjmod
c
c --- this is a package of routines for a finite-element solution of
c --- the adjoint quasi-geostrophic equations.
c
	include 'prob2.pert'
c
c>>>>>>>>>>>>
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topd(mn0),botd(mn0),rlf(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c
	common/egl2/pl,akl,alls,allv
c
	common/cnorm/jnorm
c
	common/mmde/jmode
c
	common/nrest/irst
c
	common/adrst/ityrst
c
	common/jctrl/ktim,icstart,iczero
c
	common/gcord/ir1,ir2,jr1,jr2
c
	common/bdata/scrp(n0,k0,2),vcrp(n0,k0,2)
c
c
c>>>>>>>>>>>>>>>
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),relief(mn0),
     2     topden(mn0),botden(mn0)
c
	common/atop/radt(mn0),s2t(mn0),st(mn0),fyt(mx0,4)
c
	common/abot/radb(mn0),rad1b(mn0),rad2b(mn0),s2b(mn0),sb(mn0)
     1       ,fyb(mx0,4)
c
c
        common/newbg/amms(mn0,k0),ammv(mn0,k0)
c
c
	common/adm/tiam(k0,k0),taml(k0,k0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
      common/zerot/iday0,imth0,iyear0,ihr0,mint0
c
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/steps/sumt,sumb,cost
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/swoi/iter
c
	common/icm/intc
c
	common/diag/asm(mn0,k0)
c
c
        common/iwfq/ibgw,iclw
c
	common/lan1/maxt,ifts
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/count/icnt,itcnt,ibcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      character*4 restid,runid
      character*70 idstrg
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
       real savstm1(mn0,k0),savvrt1(mn0,k0),savbot1(mn0),savtop1(mn0)
       real savstm2(mn0,k0),savvrt2(mn0,k0),savbot2(mn0),savtop2(mn0)
c
	write(6,*)'entering adjoint model'
c
c
c.......increment icalc for tmax+dt
c
	tim=tim+dt
c
        tm1=tim+dt
c
c.......clear prvjac.......
c
        do k=1,kz
        do i=1,msq
        prvjac(i,k)=0.0
        enddo
        enddo
c
c        maxrec=((maxt-1)/ibgw+1)*(2*kz+2)
c
c>>>>>>>>>>>>>>>
c
      print *,'maxt in adjoint=',maxt
c
      do 9999 itime=maxt,1,-1
c
      icalc=itime
c
      tim=tim-dt
c
c        print *,'icalc=',icalc,' bstm=',(bstm(i,1),i=781,785)
c
c
      call alop3d
c
        if(mod(itime-1,8).eq.0)then
c        write(53)stream
c        write(53)vort
        endif
c
c
c --- return to top of main loop
c
 9999 continue
c
	write(6,*)'exit adjoint model: all done'
c
	return
	end
c
c
      subroutine helm4(klevel,xlam)
c
c --- subroutine to solve the poisson equation.
c
c --- this is a fourth-order poisson solver. to maintain
c --- accuracy for the vorticity equation, the solution of the
c --- poisson equation for the streamfunction must also be of
c --- fourth order. this is done using the method of deferred
c --- corrections. see haidvogel,schulman, and robinson p 12.
c
          include 'prob2.pert'
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common /cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      hx2=hy*hy
      cd=8.-xlam*hx2
      ce=xlam*xlam*hx2/12.
      one12=1./12
      cc2=1./(6.*hx2)
c
c --- call ncar poisson solver package
c
      dbasin=hy*nm1
      call pwscrt(0,0.,xbasin,mm1,1,bda,bdb,0.,dbasin,nm1,1,bdc,bdd
     +,xlam,stream(1,klevel),m,1.,ierror,r1)
c
      return
c
      end
c
      function ncheck (n)
      np = n
  101 k = np/2
      if (2*k.ne.np) go to 102
      np = k
      go to 101
  102 if (np.ne.1) go to 103
      ncheck = 0
      return
  103 k = np/3
      if (3*k.ne.np) go to 104
      np = k
      go to 103
  104 k = np/5
      if (5*k.ne.np) go to 105
      np = k
      go to 104
  105 if (np.ne.1) go to 106
      ncheck = 1
      return
  106 ncheck = 2
      return
      end
      subroutine trdi (n,a,b,c,y,x,ks,work)
      dimension       a(n)       ,b(n)       ,c(n)       ,y(n)       ,
     1                x(n)       ,work(1)
  103 format(45h0array dimension n not greater than 2 in trdi)
  104 format('0hks not equal to 0 or 1 in call to trdi')
      if (n .gt. 2) go to 101
      print 103
      return
  101 continue
      if (ks.eq.0 .or. ks.eq.1) go to 102
      print 104
      return
  102 continue
      call trdslv (n,a,b,c,y,x,ks,work(1),work(n))
      return
      end
c
c
        subroutine adtherm(psi,tvort,kz,msq)
          include 'prob2.pert'
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
        dimension tvort(mn0,k0),psi(mn0,k0),vtm(k0,k0)
c
       do 5 k=1,kz
       do 6 kk=1,kz
       vtm(k,kk)=0.0
 6    continue
 5    continue
c
c
      do 20 i=1,kz
        ip=0
        do 19 j=1,kz
          tmp=0.0
          do 18 k=1,kz
            tmp=tmp+amat(i,k)*eigval(k)*ainv(k,j)
   18     continue
          vtm(i,j)=tmp
   19   continue
   20 continue
c
c --- now calculate the thermal vorticity
c
       do 500 k=1,kz
        do 170 kk=1,kz
          tmp=vtm(k,kk)
          do 165 ip=1,msq
            psi(ip,kk)=psi(ip,kk)+tmp*tvort(ip,k)
  165     continue
  170   continue
 500    continue
c
c
        return
        end
c
c
c
      subroutine afiltr(zz,m,n,nord)
c
c --- zz = field to be filtered
c ---  m = no of points along a side of the array
c --- nord = order of the filter
c
          include 'prob2.pert'
      dimension zz(1),g(mx0+10)
      if(m.gt.mx0+10.or.n.gt.mx0+10)stop'filerr'
      mm1=m-1
      iodev=(nord+1)/2 -nord/2
      fac=-1.+2.*float(iodev)
      fac=fac/2.**(2*nord)
c
      do 200 i=1,m
c
        do 150 j=1,n
          irow=(j-1)*m+i
          g(j)=fac*zz(irow)
  150   continue
c
        call ashpro(g,n,nord)
c
        do 110 j=1,n
          zz((j-1)*m+i)=zz((j-1)*m+i)+g(j)
  110   continue
c
  200 continue
c
      do 100 j=1,n
        irow=(j-1)*m
c
        do 50 i=1,m
          g(i)=fac*zz(irow+i)
   50   continue
c
        call ashprox(g,m,nord)
c
        do 10 i=1,m
          zz(irow+i)=zz(irow+i)+g(i)
   10   continue
c
  100 continue
      return
      end
c
      subroutine eigvec (nm,n,wr,wi,z,select,k,zr,zi,ierr)
      integer  nm,n,k,ierr
      real wr(n),wi(n),z(nm,1),zr(n),zi(n)
      logical select(n)
      ierr = 0
      if (select(k)) go to 10
      ierr = 33
      print 1002
 1002 format(' error in eigenvector selection')
      do 5 i = 1,n
      zr(i) = 0.
      zi(i) = 0.
    5 continue
      return
   10 continue
      sign = 1.
      kvec = 0
      do 30 i = 1,k
      if (.not. select(i)) go to 30
      if (wi(i) .eq. 0.) go to 20
      kvec = kvec + 2
      sign = 1.
      if (i .eq. 1) go to 30
      if (wi(i) .ne. (-wi(i-1)) .or. .not. select(i-1)) go to 30
      kvec = kvec - 2
      sign = -1.
      go to 30
   20 kvec = kvec + 1
   30 continue
      if (wi(k) .eq. 0.) go to 50
      do 40 i = 1,n
      zr(i) = z(i,kvec-1) 
      zi(i) = z(i,kvec)*sign
   40 continue
      return
   50 continue
      do 60 i = 1,n
      zr(i) = z(i,kvec)
      zi(i) = 0.
   60 continue
      return
      end
c
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,
     1        igh,its,low,mp2,enm2,ierr
      real h(nm,n),wr(n),wi(n),z(nm,n)
      real p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,machep
      logical notlas
      complex z3
      machep = 2.**(-47)
      ierr = 0
      norm = 0.0
      k = 1
      do 50 i = 1, n
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0
   50 continue
      en = igh
      t = 0.0
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0) s = norm
         if (abs(h(l,l-1)) .le. machep * s) go to 100
   80 continue
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (its .eq. 30) go to 1000 
      if (its .ne. 10 .and. its .ne. 20) go to 130
      t = t + x
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75 * s
      y = x
      w = -0.4375 * s * s
  130 its = its + 1
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. machep * abs(p)
     1    * (abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))) go to 150
  140 continue
  150 mp2 = m + 2
      do 160 i = mp2, en
         h(i,i-2) = 0.0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0
  160 continue
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         do 210 j = k, n
           p = h(k,j) + q * h(k+1,j)
           if (.not. notlas) go to 200
           p = p + r * h(k+2,j)
           h(k+2,j) = h(k+2,j) - p * zz
  200      h(k+1,j) = h(k+1,j) - p * y
           h(k,j) = h(k,j) - p * x
  210    continue
         j = min0(en,k+3)
         do 230 i = 1, j
           p = x * h(i,k) + y * h(i,k+1)
           if (.not. notlas) go to 220
           p = p + zz * h(i,k+2)
           h(i,k+2) = h(i,k+2) - p * r
  220      h(i,k+1) = h(i,k+1) - p * q
           h(i,k) = h(i,k) - p
  230    continue
         do 250 i = low, igh
           p = x * z(i,k) + y * z(i,k+1)
           if (.not. notlas) go to 240
           p = p + zz * z(i,k+2)
           z(i,k+2) = z(i,k+2) - p * r
  240      z(i,k+1) = z(i,k+1) - p * q
           z(i,k) = z(i,k) - p
  250    continue
  260 continue
      go to 70
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0
      en = na
      go to 60
  280 p = (y - x) / 2.0
      q = p * p + w
      zz = sqrt(abs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0) go to 320
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0) wr(en) = x - w / zz
      wi(na) = 0.0
      wi(en) = 0.0
      x = h(en,na)
      r = sqrt(x*x+zz*zz)
      p = x / r
      q = zz / r
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
      go to 330
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
  340 if (norm .eq. 0.0) go to 1001
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
  600    m = en
         h(en,en) = 1.0
         if (na .eq. 0) go to 800
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = h(i,en)
            if (m .gt. na) go to 620
            do 610 j = m, na
  610       r = r + h(i,j) * h(j,en)
  620       if (wi(i) .ge. 0.0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0) go to 640
            t = w
            if (w .eq. 0.0) t = machep * norm
            h(i,en) = -r / t
            go to 700
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 700
  650       h(i+1,en) = (-s - y * t) / zz
  700    continue
         go to 800
  710    m = na
         if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    z3 = cmplx(0.0,-h(na,en)) / cmplx(h(na,na)-p,q)
         h(na,na) = real(z3)
         h(na,en) = aimag(z3)
  730    h(en,na) = 0.0
         h(en,en) = 1.0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
         do 790 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0
            sa = h(i,en)
            do 760 j = m, na
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
            if (wi(i) .ge. 0.0) go to 770
            zz = w
            r = ra
            s = sa
            go to 790
  770       m = i
            if (wi(i) .ne. 0.0) go to 780
            z3 = cmplx(-ra,-sa) / cmplx(w,q)
            h(i,na) = real(z3)
            h(i,en) = aimag(z3)
            go to 790
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0 * q
            if (vr .eq. 0.0 .and. vi .eq. 0.0) vr = machep * norm
     1       * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
            z3 = cmplx(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra) / cmplx(vr,vi)
            h(i,na) = real(z3)
            h(i,en) = aimag(z3)
            if (abs(x) .le. abs(zz) + abs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       z3 = cmplx(-r-y*h(i,na),-s-y*h(i,en)) / cmplx(zz,q)
            h(i+1,na) = real(z3)
            h(i+1,en) = aimag(z3)
  790    continue
  800 continue
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
         do 820 j = i, n
  820    z(i,j) = h(i,j)
  840 continue
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
         do 880 i = low, igh
            zz = 0.0
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
            z(i,j) = zz
  880 continue
      go to 1001
 1000 ierr = en
 1001 return
      end
c
c
      subroutine ngenrl(s,zp,klevel)
c
          include 'prob2.pert'
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/bsolv/ab(mx0),bb(mx0),cb(mx0),workb(mnwork)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
      dimension s(1),zp(1)
      dimension aac(5),fx(mx0),fy(mx0),inout(mx0)
      dimension inflow(ibp),ipos(ibp)
c
c.....clear inflow and ipos.....
c
	do 9000 ip=1,ibp
	ipos(ip)=1
 9000   inflow(ip)=0
c
      do 5 iside=1,4
        ib=icorn(iside)
        ip=ncorn(iside)
        aac(iside)=r(ib)-r(ip)
    5 continue
      aac(5)=aac(1)
c
	ict=0
c
      do 100 iside=1,4
	ict=ict+1
        i1=icorn(iside)
        i2=i1+mast(iside)
        i3=i2+mast(iside)
c
c --- determine inflow points
c
	ist=mast(iside)
	if(ist.eq.1.or.ist.eq.-1)mn2=mm2
	if(ist.eq.m.or.ist.eq.-m)mn2=nm2
        do 10 ii=1,mn2
          inout(ii)=0
          if(s(i3).gt.s(i1)) inout(ii)=1
          i1=i2
          i2=i3
          i3=i3+mast(iside)
   10   continue
        ib=icorn(iside)+mast(iside)
	ipos(ict)=ib
        if(inout(1).eq.1) go to 15
        i5=ib+nast(iside)
        i6=i5+mast(iside)
        fy(1)=r(ib)-4.*r(i5)-r(i6)-0.5*aac(iside)
        bb(1)=7.
        cb(1)=2.
        go to 20
   15   continue
        fy(1)=r1(ib)-zp(ib)
        bb(1)=1.
        cb(1)=0.
	inflow(ict)=1
   20   continue
	if(ist.eq.1.or.ist.eq.-1)mn3=m-3
	if(ist.eq.m.or.ist.eq.-m)mn3=n-3
        do 40 iz=2,mn3
	ict=ict+1
          ib=ib+mast(iside)
	  ipos(ict)=ib
          if(inout(iz).eq.1) go to 35
          i5=ib+nast(iside)
          i4=i5-mast(iside)
          i6=i5+mast(iside)
          fy(iz)=r(ib)-4.*r(i5)-r(i4)-r(i6)
          bb(iz)=8.
          ab(iz)=2.
          cb(iz)=2.
          go to 40
   35     continue
          fy(iz)=r1(ib)-zp(ib)
          bb(iz)=1.
          ab(iz)=0.
          cb(iz)=0.
	  inflow(ict)=1
   40   continue
	ict=ict+1
        ib=ib+mast(iside)
	ipos(ict)=ib
        if(inout(mn2).eq.1) go to 45
        i5=ib+nast(iside)
        i4=i5-mast(iside)
        fy(mn2)=r(ib)-4.*r(i5)-r(i4)-0.5*aac(iside+1)
        bb(mn2)=7.
        ab(mn2)=2.
        go to 50
   45   continue
        fy(mn2)=r1(ib)-zp(ib)
        bb(mn2)=1.
        ab(mn2)=0.
	inflow(ict)=1
   50   continue
c
c --- call the ncar tri-diagonal matrix inverter
c
        call trdi(mn2,ab,bb,cb,fy,fx,0,workb)
c
        ib=icorn(iside)
        do 60 iz=1,mn2
          ib=ib+mast(iside)
          r(ib)=fx(iz)
   60   continue
  100 continue
      do 120 iside=1,4
      ict=ict+1
        ib=icorn(iside)
	ipos(ict)=ib
        i1=ib+nast(iside)
        i2=ib+mast(iside)
        if(s(i2).gt.s(i1)) go to 110
        r(ib)=0.25*aac(iside)-0.5*(r(i1)+r(i2))
        go to 120
  110   continue
        r(ib)=r1(ib)-zp(ib)
	inflow(ict)=1
  120 continue
c
	if(ict.gt.ibp)write(6,*)'error in ngenrl'
	if(ict.gt.ibp)stop
c
c
c.....write out inflow.....
c
c	if(ipass.eq.2)write(95)icalc,klevel,inflow,ipos
c
      return
      end
c
      subroutine trid (kstart,kstop,ising,m,mm1,a,b,c,y,twocos,d,w)
      dimension        a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     1                 twocos(1)  ,d(1)       ,w(1)
c
c uro is a machine dependent constant.  it is the smallest floating
c point number such that 1 + uro is greater than 1 in the
c arithmetic of this computer.
c
      uro = 1./2.**20
  101     do 106 k=kstart,kstop
          x = -twocos(k)
          d(1) = c(1)/(b(1)+x)
          y(1) = y(1)/(b(1)+x)
              do 103 i=2,m
              z = b(i)+x-a(i)*d(i-1)
              d(i) = c(i)/z
  102         y(i) = (y(i)-a(i)*y(i-1))/z
  103         continue
              do 105 ip=1,mm1
              i = m-ip
  104         y(i) = y(i)-d(i)*y(i+1)
  105         continue
  106     continue
  107 if (ising.eq.1) return
      d(1) = c(1)/(b(1)-2.)
      y(1) = y(1)/(b(1)-2.)
          do 109 i=2,mm1
          z = b(i)-2.-a(i)*d(i-1)
          d(i) = c(i)/z
  108     y(i) = (y(i)-a(i)*y(i-1))/z
  109     continue
      z = b(m)-2.-a(m)*d(m-1)
      x = y(m)-a(m)*y(m-1)
      if (abs(z).gt.(30.*m*abs(a(m)))*uro) go to 111
  110 y(m) = 0.
      go to 112
  111 y(m) = x/z
  112     do 114 ip=1,mm1
          i = m-ip
  113     y(i) = y(i)-d(i)*y(i+1)
  114     continue
      return
      end
c
c
c
c
      subroutine balanc(nm,n,a,low,igh,scale)
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      real a(nm,n),scale(n)
      real c,f,g,r,s,b2,radix
      logical noconv
      radix = 2.
      b2 = radix * radix
      k = 1
      l = n
      go to 100
   20 scale(m) = j
      if (j .eq. m) go to 50
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
   50 go to (80,130), iexc
   80 if (l .eq. 1) go to 280
      l = l - 1
  100 do 120 jj = 1, l
         j = l + 1 - jj
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.0) go to 120
  110    continue
         m = l
         iexc = 1
         go to 20
  120 continue
      go to 140
  130 k = k + 1
  140 do 170 j = k, l
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0) go to 170
  150    continue 
         m = k
         iexc = 2
         go to 20
  170 continue
      do 180 i = k, l
  180 scale(i) = 1.0
  190 noconv = .false.
      do 270 i = k, l
         c = 0.0
         r = 0.0
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(a(j,i)) 
            r = r + abs(a(i,j))
  200    continue
         g = r / radix
         f = 1.0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
  240    if ((c + r) / f .ge. 0.95 * s) go to 270
         g = 1.0 / f
         scale(i) = scale(i) * f 
         noconv = .true.
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
  270 continue
      if (noconv) go to 190
  280 low = k
      igh = l
      return
      end
c
      subroutine elmhes(nm,n,low,igh,a,int)
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real a(nm,n)
      real x,y
      integer int(igh)
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0
         i = m
         do 100 j = m, igh
            if (abs(a(j,mm1)) .le. abs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
         int(m) = i
         if (i .eq. m) go to 130
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
  130    if (x .eq. 0.0) go to 180
         mp1 = m + 1
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0) go to 160
            y = y / x
            a(i,mm1) = y
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
  160    continue
  180 continue
  200 return
      end
c
      function ierinv (n,na,nv)
  103 format(23h0* n .lt. 1 in invmtx *)
  104 format(24h0* na .lt. n in invmtx *)
  105 format(24h0* nv .lt. n in invmtx *)
      ierinv = 0
      if (n .ge. 1) go to 101
      ierinv = 34
      print 103
      return
  101 if (na .ge. n) go to 102
      ierinv = 35
      print 104
      return
  102 if (nv .ge. n) return
      ierinv = 36
      print 105
      return
      end
c
c
c
      subroutine nlop3d
c
c --- subroutine nlop3d is the main computation loop for this
c --- code.  it advances the vorticity and solves for the new
c --- streamfunction fields.
c
          include 'prob2.pert'
c
c>>>>>>>>>>>>
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c
c>>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/dynam/alpha,beta,rf
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +  topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +  prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/pert/stf(nosp),vtf(nosp)
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/steps/sumt,sumb,cost
c
c
c
c ---advance the vorticity
c
c.....skip nvortq if icalc=0.....
c
      if(icalc.eq.0)goto 8080
      do 100 k=1,kz
        call nvortq(k)
  100 continue
 8080 continue
c
c       print *,'done nvortq'
c
c --- solve for three-dimensional streamfunction
c
c
c --- solve helmholtz equation level by level.
c
      do 300 k=1,kz
c
c --- first use r to hold the transformed interior
c       vorticity field for this level.
c
        do j=1,n
        do i=1,m
          ip=(j-1)*m+i
          r(ip)=0.0
          do k1=1,kz
            r(ip)=r(ip)+ainv(k,k1)*vort(ip,k1)
          enddo
        enddo
        enddo
c
c
        do 290 ip=1,msq
 290    stream(ip,k)=r(ip)
c
c
c.......impose boundary conditions on r......
c
        do 2503 ip=1,m
 2503   stream(ip,k)=0.0
c
        do 2504 ip=msq-m+1,msq
 2504   stream(ip,k)=0.0
c
c<<<<<<<<<>>>>>>>>>>>
c
        mb=0
c
      dbasin=hy*nm1
      call pwscrt(0,0.,xbasin,mm1,mb,bda,bdb,0.,dbasin,nm1,1,bdc,bdd
     +,eigval(k),stream(1,k),m,1.,ierror,r1)
c<<<<<>>>>>>>>
c
  300 continue
c
c --- new streamfunction--point by point
c
      do 350 ip=1,msq
        do 340 k=1,kz
           r(k)=stream(ip,k)
           stream(ip,k)=0.
  340   continue
c
        do 345 k=1,kz
        do 345 j=1,kz
           stream(ip,k)=stream(ip,k)+amat(k,j)*r(j)
  345   continue
  350 continue
c
c
      return
      end
c
      subroutine balbak(nm,n,low,igh,scale,m,z)
      integer i,j,k,m,n,ii,nm,igh,low
      real scale(n),z(nm,m)
      real s
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
      do 110 i = low, igh
         s = scale(i)
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
  110 continue
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
  140 continue
  200 return
      end
c
      subroutine eltran(nm,n,low,igh,a,int,z)
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real a(nm,igh),z(nm,n)
      integer int(igh)
      do 80 i = 1, n
         do 60 j = 1, n
   60    z(i,j) = 0.0
         z(i,i) = 1.0
   80 continue
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
         i = int(mp)
         if (i .eq. mp) go to 140
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0
  130    continue
         z(i,mp) = 1.0
  140 continue
  200 return
      end
c
c
c
c
c
      subroutine vwrite(vector, n, r1, r2, r3, r4, r5, iunit)
      dimension vector(n)
      write(iunit,1)n,r1,r2,r3,r4,r5
      write(iunit,2)vector
 1    format(i5,1p,5e10.3)
 2    format(1p,8e10.3)
      return
      end
c
c
c
	subroutine bener
c
        include 'prob2.pert'
c
	common/egy/peint,akeint
	common/l2n/al2s,al2v
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
c
	real pe(mn0,k01),ke(mn0,k0),pex(n0,k01),kex(n0,k0),pez(k01),
     1       kez(k0)
c
	total=0.0
	do 1000 k=1,kz
	do 1000 j=1,n
	do 1000 i=1,m-1
	ip=(j-1)*m+i
c1000   total=total+stream(ip,k)*vort(ip,k)*hz(k)
 1000   total=total+stream(ip,k)*vort(ip,k)*hz(k)
c
	akeint=-0.5*total
        peint=0.0
c
	return
	end
c
      subroutine infld(r,k,t,ityp,ifbnd)
c
c --- this subroutine provides the streamfunction, vorticity, density,
c --- topography, and determines when to write data for a baroclinic
c --- simulation. ityp determines the values returned.
c
c
c --- the argument list for this subroutine is as follows:
c
c ---      r : an (mxm) array in which the values are returned to
c              the main routine
c 
c ---      k : the level for which the computation is being done
c
c ---      t : the time at which the computation is being done
c
c ---   ityp : the indicator for which computation is to be done
c
c              = 1 : streamfunction
c              = 2 : vorticity
c              = 3 : verification streamfunction
c              = 4 : verification vorticity
c              = 5 : top vertical velocity (w/r0 = w1)
c              = 6 : bottom slope (zb = -h + r0*eta(x,y))
c              = 7 : sigma*dpsi/dz values for the top
c              = 8 : sigma*dpsi/dz values for the bottom
c              = 9 : write unformatted fields
c
c ---  ifbnd : calculate boundary values only (1/0=yes/no)
c
c
c --- the real array pram(16) in common /user/ is used to provide
c --- parameters for the calculations. current definitions are 
c --- as follows:
c                    pram(2) = yes/no read the vort from input file
c                    pram(3) = 1/0 calc therm+sst vort/don't calculate
c                    pram(5) = iskip, no of days to skip in data set
c                    pram(6) = day1, time to skip to first data set
c                    pram(9) = r0, the rossby number
c                    pram(10)= d,horizontal length scale of motion
c                    pram(11)= lz, average depth of the basin
c                    pram(12)= ht, vertical length of motion
c                              (typically, height of thermocline)
c                    pram(16)= ifpers (1/0 = persis/bench)
c
      include 'prob2.pert'
c
      common/user/pram(16)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +  ifsbl,ifeva
      common/field/psi(mn0,k0),zeta(mn0,k0),tlast,day1,tkeep,
     +  curpsi(mn0,k0),curzet(mn0,k0),psiprv(mn0,k0),zetprv(mn0,k0),
     +  topprv(mn0),topcur(mn0),topdens(mn0)
      common/verify/verpsi(mn0,k0),verzeta(mn0,k0)
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n, nm1,nm2,np1
      common/istar/icall
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c
	common/swoi/iter
c
	common/perst/psf(mn0,k0),pvr(mn0,k0)
c
	common/bdata/scrp(n0,k0,2),vcrp(n0,k0,2)
c
      dimension r(1),topog(mn0)
c      real lz
c      character*80 fname33,fname13
c@@@@@@@@@@@@@@@@@
c
c.......read stream and vort boundary conditions from file
c
	goto 8888
	rewind 21
c	rewind 22
	read(21,*)icrap
c	read(22,*)icrap
	do 3000 klev=1,kz
	read(21,3001)(scrp(j,klev,1),j=1,n)
	read(21,3001)(vcrp(j,klev,1),j=1,n)
c	read(22,3001)(scrp(j,klev,2),j=1,n)
c	read(22,3001)(vcrp(j,klev,2),j=1,n)
 3000   continue
 3001   format(1p,8e10.3)
c
c	if(icalc.eq.32)then
c	print *,'scrp=',scrp
c	print *,'vcrp=',vcrp
c	endif
 8888   continue
c
c
c@@@@@@@@@@@@@@@@@@@@
c
c*********************
c......skip to end if iter > 1.....
c
c	if(iter.gt.1)goto 9123
c
c**********************
c
c
c
c
c --- get the proper field according to the value of ityp
c
c --- streamfunction
c
  100 continue
      if(ityp.ne.1)goto 200
      if(ifbnd.ne.0)goto 153
c
c......add random noise to the streamfunction initial condition...
c
	print *,'m=',m,' n=',n
	s=1
	is=s
 	do 154 j=1,n
 	do 154 i=1,m
 	ip=(j-1)*m+i
	r(ip)=0.0
	r(ip)=1.e-5*gasdev(is)
c 	r(ip)=1.e-5*ranf()
c	print *,'r=',r(ip)
  154 continue
c       r(2230)=1.e-10
c
 153    continue
	do 150 j=1,n,(n-1)
	do 150 i=1,m
	ip=(j-1)*m+i
 150    r(ip)=0.0
c
 	do 151 i=1,m,(m-1)
 	do 151 j=2,(n-1)
 	ip=(j-1)*m+i
 151    r(ip)=0.0
c@@@@@@@@@@@@@@@@@@
c
c......set western and eastern boundary condition.......
c
	do 155 j=1,n
	ip=(j-1)*m+1
	ipe=j*m
c	r(ipe)=scrp(j,k,2)
c155    r(ip)=1.e-6*bstm(ip,k)
 155    r(ip)=0.0
c
c@@@@@@@@@@@@@@@@@@@
      return
c
c --- vorticity
c
  200 continue
      if(ityp.ne.2)goto 300
      if(ifbnd.ne.0)goto 253
	do 254 ip=1,msq
        r(ip)=0.0
  254 continue
	return
c
 253  continue
	do 250 j=1,n,(n-1)
	do 250 i=1,m
	ip=(j-1)*m+i
 250    r(ip)=0.0
c
 	do 251 i=1,m,(m-1)
 	do 251 j=2,(n-1)
 	ip=(j-1)*m+i
 251    r(ip)=0.0
c@@@@@@@@@@@@@@@@@@
c
c......set western and eastern boundary condition.......
c
	do 255 j=1,n
	ip=(j-1)*m+1
	ipe=j*m
c	r(ipe)=vcrp(j,k,2)
c255    r(ip)=1.e-6*bvrt(ip,k)
 255    r(ip)=0.0
c
c@@@@@@@@@@@@@@@@@@@
      return
c
c --- verification streamfunction
c
  300 continue
c
c........set top density......
c
      if(ityp.ne.7)goto 3002
      if(ifbnd.ne.0)goto 2503
	do 2504 ip=1,msq
        r(ip)=0.0
 2504  continue
	return
c
 2503  continue
	do 2500 j=1,n,(n-1)
	do 2500 i=1,m
	ip=(j-1)*m+i
 2500    r(ip)=0.0
c
 	do 2501 i=1,m,(m-1)
 	do 2501 j=2,(n-1)
 	ip=(j-1)*m+i
 2501    r(ip)=0.0
c@@@@@@@@@@@@@@@@@@
c
c......set western and eastern boundary condition.......
c
	do 2505 j=1,n
	ip=(j-1)*m+1
	ipe=j*m
c2505    r(ip)=1.e-6*btopd(ip)
 2505    r(ip)=0.0
c
c@@@@@@@@@@@@@@@@@@@
      return
c
 3002   continue
c
c........set bottom density......
c
      if(ityp.ne.8)goto 3003
      if(ifbnd.ne.0)goto 2703
	do 2704 ip=1,msq
        r(ip)=0.0
 2704  continue
	return
c
 2703  continue
	do 2700 j=1,n,(n-1)
	do 2700 i=1,m
	ip=(j-1)*m+i
 2700    r(ip)=0.0
c
 	do 2701 i=1,m,(m-1)
 	do 2701 j=2,(n-1)
 	ip=(j-1)*m+i
 2701    r(ip)=0.0
c@@@@@@@@@@@@@@@@@@
c
c......set western and eastern boundary condition.......
c
	do 2705 j=1,n
	ip=(j-1)*m+1
	ipe=j*m
c2705    r(ip)=1.e-6*bbotd(ip)
 2705    r(ip)=0.0
c
c@@@@@@@@@@@@@@@@@@@
      return
c
 3003   continue
c
c
	return
	end
c
      subroutine poinit (nperod,n,iex2,iex3,iex5,n2pw,n2p3p,n2p3p5,ks3,
     1                   ks5,twocos)
      dimension        twocos(1)
c
      pi = acos(-1.0)
  101 np = n
      if (nperod.eq.1) np = np+1
      if (nperod.eq.3) np = np-1
      iex2 = 0
  102 k = np/2
      if (2*k.ne.np) go to 103
      iex2 = iex2+1
      np = k
      go to 102
  103 iex3 = 0
  104 k = np/3
      if (3*k.ne.np) go to 105
      iex3 = iex3+1
      np = k
      go to 104
  105 iex5 = 0
  106 k = np/5
      if (5*k.ne.np) go to 107
      iex5 = iex5+1
      np = k
      go to 106
  107 continue
      n2pw = 2**iex2
      n2p3p = n2pw*(3**iex3)
      n2p3p5 = n2p3p*(5**iex5)
      np = 1
      twocos(1) = 0.
      k = 1
      if (iex2.eq.0) go to 113
      l = iex2
      if (iex3+iex5.ne.0) go to 108
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 142
  108     do 112 kount=1,l
          np = 2*np
              do 110 i=1,np
  109         twocos(k+i) = 2.*cos((i-.5)*pi/np)
  110         continue
  111     k = k+np
  112     continue
  113 if (iex3.eq.0) go to 118
      l = iex3
      if (iex5.ne.0) go to 114
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 122
  114     do 117 kount=1,l
          np = 3*np
              do 116 i=1,np
  115         twocos(k+i) = 2.*cos((i-.5)*pi/np)
  116         continue
          k = k+np
  117     continue
  118 l = iex5
      if (nperod.eq.1) l = l-1
      if (l.le.0) go to 122
          do 121 kount=1,l
          np = 5*np
              do 120 i=1,np
  119         twocos(k+i) = 2.*cos((i-.5)*pi/np)
  120         continue
          k = k+np
  121     continue
  122 if (nperod.eq.1.or.nperod.eq.2) go to 128
      if (nperod.eq.0) go to 125
      npt2 = 2*np
          do 124 i=1,npt2
  123     twocos(k+i) = 2.*cos(i*pi/np)
  124     continue
      k = k+npt2
      go to 128
  125     do 127 i=1,np
  126     twocos(k+i) = 2.*cos(2.*i*pi/np)
  127     continue
      k = k+np
  128 np = n2p3p5
      if (iex5.eq.0) go to 137
      ks5 = k
          do 136 kount=1,iex5
          np = np/5
          npt2 = 2*np
              do 130 i=1,npt2
  129         twocos(k+i) = 2.*cos((i-.5)*pi/npt2)
  130         continue
          k = k+npt2
              do 132 i=1,np
              twocos(k+4*i-3) = 2.*cos((i-.8)*pi/np)
              twocos(k+4*i-2) = 2.*cos((i-.6)*pi/np)
              twocos(k+4*i-1) = 2.*cos((i-.4)*pi/np)
  131         twocos(k+4*i) = 2.*cos((i-.2)*pi/np)
  132         continue
          k = k+4*np
              do 134 i=1,np
              twocos(k+2*i-1) = 2.*cos((3*i-2)*pi/(3*np))
  133         twocos(k+2*i) = 2.*cos((3*i-1)*pi/(3*np))
  134         continue
  135     k = k+2*np
  136     continue
  137 if (iex3.eq.0) go to 142
      ks3 = k
          do 141 kount=1,iex3
          np = np/3
              do 139 i=1,np
              twocos(k+2*i-1) = 2.*cos((3*i-2)*pi/(3*np))
  138         twocos(k+2*i) = 2.*cos((3*i-1)*pi/(3*np))
  139         continue
  140     k = k+2*np
  141     continue
  142 return
      end
c
c
c
c
c
      subroutine ahelm4(klevel,xlam)
c
c --- subroutine to solve the poisson equation.
c
c --- this is a fourth-order poisson solver. to maintain
c --- accuracy for the vorticity equation, the solution of the
c --- poisson equation for the streamfunction must also be of
c --- fourth order. this is done using the method of deferred
c --- corrections. see haidvogel,schulman, and robinson p 12.
c
          include 'prob2.pert'
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common /cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      hx2=hy*hy
      cd=8.-xlam*hx2
      ce=xlam*xlam*hx2/12.
      one12=1./12
      cc2=1./(6.*hx2)
c
c --- call ncar poisson solver package
c
      dbasin=hy*nm1
      call apwscrt(0,0.,xbasin,mm1,1,bda,bdb,0.,dbasin,nm1,1,bdc,bdd
     +,xlam,stream(1,klevel),m,1.,ierror,r1,r2)
c
      return
      end
c
	subroutine bl2n
c
        include 'prob2.pert'
c
	common/egy/peint,akeint
	common/l2n/al2s,al2v
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
c
	real sx(n0,k0),vx(n0,k0),sz(k0),vz(k0)
c
c.......integrate norms over whole domain.....
c
	sums=0.0
	sumv=0.0
	do k=1,kz
	do j=1,n
	do i=1,m
	ip=(j-1)*m+i
	sums=sums+stream(ip,k)**2*hz(k)
        sumv=sumv+vort(ip,k)**2*hz(k)
        enddo
        enddo
        enddo
c
	al2s=sums
	al2v=sumv
c
	return
	end
c
      subroutine initial(t)
c
c --- it call getfld for the first and second day of data
c --- preparing the arrays for interp
c
          include 'prob2.pert'
      common/user/pram(16)
      common/field/psi(mn0,k0),zeta(mn0,k0),tlast,day1,tkeep,
     +  curpsi(mn0,k0),curzet(mn0,k0),psiprv(mn0,k0),zetprv(mn0,k0),
     +  topprv(mn0),topcur(mn0),topdens(mn0)
c
c
c --- skip in the data set
c
        iskip=pram(5)
        if(iskip.eq.0) go to 100
        do 99 i=1,iskip
        call getfld(psiprv,zetprv,topprv,t)
   99   continue
100     continue
        call getfld(psiprv,zetprv,topprv,t)
        call getfld(curpsi,curzet,topcur,t)
        tlast=pram(6)
      return
      end
c
      subroutine pois (iflag,nperod,n,mperod,m,a,b,c,idimy,y,w)
      external         trid       ,tridp
      dimension        w(1)       ,b(1)       ,a(1)       ,c(1)       ,
     1                 y(idimy,1)
      iwdim1 = 5.25*n+1
      iwdim2 = iwdim1+m
      iwdim3 = iwdim2+m
      iwdim4 = iwdim3+m
      iwdim5 = iwdim4+m
          do 102 i=1,m
          a(i) = -a(i)
          c(i) = -c(i)
  101     w(iwdim5+i-1) = -b(i)+2.
  102     continue
      if (mperod.eq.0) go to 103
      call poisgn (nperod,n,iflag,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),trid)
      go to 104
  103 call poisgn (nperod,n,iflag,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),tridp)
  104     do 106 i=1,m
          a(i) = -a(i)
  105     c(i) = -c(i)
  106     continue
  107 return
      end
c
c
c
      subroutine sort(x,amat,icd)
          include 'prob2.pert'
      dimension x(icd),amat(k0,k0)
c
c --- routine to sort the floating point array x from greatest to
c --- smallest and arrange the columns of the matrix amat in the same
c --- order as the array x
c
      if(icd.le.1)return
      k=1
   10 continue
      if(k.eq.0) goto 35
      k=0
      do 30 i=2,icd
        if(x(i-1).ge.x(i)) goto 30
        k=1
        temp=x(i-1)
        x(i-1)=x(i)
        x(i)=temp
        do 25 j=1,icd
          temp=amat(j,i-1)
          amat(j,i-1)=amat(j,i)
          amat(j,i)=temp
   25   continue
   30 continue
      goto 10
   35 continue
      return
      end
c
c
      subroutine poisgn (nperod,n,iflag,m,ba,bb,bc,idimq,q,twocos,b,t,
     1                   d,w,tri)
	external tri
      dimension        ba(1)      ,bb(1)      ,bc(1)      ,q(idimq,1) ,
     1                 twocos(1)  ,b(1)       ,t(1)       ,d(1)       ,
     2                 w(1)
      if (iflag.ne.0) go to 101
      call poinit (nperod,n,iex2,iex3,iex5,n2pw,n2p3p,n2p3p5,ks3,ks5,
     1             twocos)
  101 mm1 = m-1
          do 106 j=1,n
              do 103 i=1,m
  102         b(i) = -q(i,j)
  103         continue
          call tri (1,1,1,m,mm1,ba,bb,bc,b,twocos,d,w)
              do 105 i=1,m
  104         q(i,j) = b(i)
  105         continue
  106     continue
      np = 1
  107 if (iex2.eq.0) go to 118
      l = iex2
      if (iex3+iex5.ne.0) go to 108
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 156
  108     do 117 kount=1,l
          k = np
          np = 2*np
          k2 = np-1
          k3 = np
          k4 = 2*np-1
          jstart = np
          if (nperod.eq.3) jstart = 1
              do 116 j=jstart,n,np
              jm1 = j-k
              jp1 = j+k
              if (j.eq.1) jm1 = jp1
              if (j.ne.n) go to 109
              jp1 = jm1
              if (nperod.eq.0) jp1 = k
  109             do 111 i=1,m
  110             b(i) = q(i,jm1)+q(i,jp1)
  111             continue
              call tri (k,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 113 i=1,m
                  t(i) = q(i,j)+b(i)
  112             b(i) = t(i)
  113             continue
              call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 115 i=1,m
  114             q(i,j) = t(i)+2.*b(i)
  115             continue
  116         continue
  117     continue
  118 if (iex3.eq.0) go to 134
      l = iex3
      if (iex5.ne.0) go to 119
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 156
  119 k2 = np-1
          do 133 kount=1,l
          k = np
          np = 3*np
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
          k4 = k2+np
          jstart = np
          if (nperod.eq.3) jstart = 1
              do 132 j=jstart,n,np
              if (j.ne.1) go to 120
              jm1 = j+k
              jm2 = jm1+k
              go to 122
  120         jm1 = j-k
              jm2 = jm1-k
              if (j.ne.n) go to 122
              if (nperod.eq.0) go to 121
              jp1 = jm1
              jp2 = jm2
              go to 123
  121         jp1 = k
              jp2 = jp1+k
              go to 123
  122         jp1 = j+k
              jp2 = jp1+k
  123             do 125 i=1,m
  124             b(i) = 2.*q(i,j)+q(i,jm2)+q(i,jp2)
  125             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 127 i=1,m
                  t(i) = b(i)+q(i,jm1)+q(i,jp1)
  126             b(i) = t(i)
  127             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 129 i=1,m
                  q(i,j) = q(i,j)+b(i)
  128             b(i) = t(i)
  129             continue
              call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 131 i=1,m
  130             q(i,j) = q(i,j)+3.*b(i)
  131             continue
  132         continue
  133     continue
  134 l = iex5
  135 if (nperod.eq.1) l = l-1
      if (l.le.0) go to 156
      k2 = (n2pw+n2p3p)/2-1
          do 155 kount=1,l
          k = np
          np = 5*np
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
          k4 = k2+np
          jstart = np
          if (nperod.eq.3) jstart = 1
              do 154 j=jstart,n,np
              if (j.ne.1) go to 136
              jm1 = j+k
              jm2 = jm1+k
              jm3 = jm2+k
              jm4 = jm3+k
              go to 138
  136         jm1 = j-k
              jm2 = jm1-k
              jm3 = jm2-k
              jm4 = jm3-k
              if (j.ne.n) go to 138
              if (nperod.eq.0) go to 137
              jp1 = jm1
              jp2 = jm2
              jp3 = jm3
              jp4 = jm4
              go to 139
  137         jp1 = k
              jp2 = jp1+k
              jp3 = jp2+k
              jp4 = jp3+k
              go to 139
  138         jp1 = j+k
              jp2 = jp1+k
              jp3 = jp2+k
              jp4 = jp3+k
  139             do 141 i=1,m
  140             b(i) = 6.*q(i,j)+4.*(q(i,jm2)+q(i,jp2))+q(i,jm4)+
     1                   q(i,jp4)
  141             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 143 i=1,m
  142             b(i) = b(i)+3.*(q(i,jm1)+q(i,jp1))+q(i,jm3)+q(i,jp3)
  143             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 145 i=1,m
                  t(i) = b(i)
  144             b(i) = 2.*q(i,j)+q(i,jm2)+q(i,jp2)+b(i)
  145             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 147 i=1,m
  146             b(i) = b(i)+q(i,jm1)+q(i,jp1)
  147             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 149 i=1,m
                  temp = b(i)+q(i,j)
                  b(i) = 4.*q(i,j)+3.*(q(i,jm2)+q(i,jp2))+q(i,jm4)+
     1                   q(i,jp4)-t(i)
  148             q(i,j) = temp
  149             continue
              call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 151 i=1,m
  150             b(i) = b(i)+2.*(q(i,jm1)+q(i,jp1))+q(i,jm3)+q(i,jp3)
  151             continue
              call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 153 i=1,m
  152             q(i,j) = q(i,j)+5.*b(i)
  153             continue
  154         continue
  155     continue
  156 if (nperod.eq.1.or.nperod.eq.2) go to 171
      if (nperod.eq.0) go to 161
          do 158 i=1,m
  157     b(i) = 2.*q(i,1)
  158     continue
      call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
          do 160 i=1,m
          q(i,n) = q(i,n)+b(i)
  159     b(i) = 4.*q(i,n)
  160     continue
      call tri (k4+1,k4+2*np-1,2,m,mm1,ba,bb,bc,b,twocos,d,w)
      go to 164
  161     do 163 i=1,m
  162     b(i) = 2.*q(i,n)
  163     continue
      call tri (k4+1,k4+np-1,2,m,mm1,ba,bb,bc,b,twocos,d,w)
  164     do 166 i=1,m
  165     q(i,n) = q(i,n)+b(i)
  166     continue
      if (nperod.eq.0) go to 171
          do 168 i=1,m
  167     b(i) = 2.*q(i,n)
  168     continue
      call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
          do 170 i=1,m
  169     q(i,1) = q(i,1)+b(i)
  170     continue
  171 continue
      if (iex5.eq.0) go to 209
      np = n2p3p5
      k8 = ks5
      if (nperod.eq.1) k3 = k3+np/5
          do 208 kount=1,iex5
          k = np
          np = np/5
          k4 = k3-1
          k3 = k3-np
          k1 = k8+1
          k2 = k8+2*np
          k5 = k2+1
          k6 = k2+4*np
          k7 = k6+1
          k8 = k6+2*np
          jstart = np
          if (nperod.eq.3) jstart = 1+np
              do 207 j=jstart,n,k
              jm1 = j-np
              jp1 = j+np
              jp2 = jp1+np
              jp3 = jp2+np
              jp4 = jp3+np
              if (jm1.ne.0) go to 177
              if (nperod.eq.0) go to 174
                  do 173 i=1,m
  172             b(i) = q(i,jp1)
  173             continue
              go to 180
  174             do 176 i=1,m
  175             b(i) = q(i,jp1)+q(i,n)
  176             continue
              go to 180
  177             do 179 i=1,m
  178             b(i) = q(i,jp1)+q(i,jm1)
  179             continue
  180         call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 182 i=1,m
                  q(i,j) = q(i,j)+b(i)
  181             b(i) = q(i,jp1)+q(i,jp3)
  182             continue
              call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 184 i=1,m
  183             q(i,jp2) = q(i,jp2)+b(i)
  184             continue
              if (jp4.gt.n) go to 187
                  do 186 i=1,m
  185             b(i) = 2.*q(i,jp2)+q(i,j)+q(i,jp4)
  186             continue
              go to 190
  187             do 189 i=1,m
  188             b(i) = 2.*q(i,jp2)+q(i,j)
  189             continue
  190         call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 192 i=1,m
                  q(i,jp2) = q(i,jp2)+b(i)
  191             b(i) = q(i,jp2)+q(i,j)
  192             continue
              call tri (k5,k6,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 194 i=1,m
                  q(i,jp2) = q(i,jp2)+b(i)
  193             b(i) = q(i,j)+q(i,jp2)
  194             continue
              call tri (k7,k8,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 196 i=1,m
                  q(i,j) = q(i,j)+b(i)
  195             b(i) = q(i,j)+q(i,jp2)
  196             continue
              call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 198 i=1,m
  197             q(i,jp1) = q(i,jp1)+b(i)
  198             continue
              if (jp4.gt.n) go to 201
                  do 200 i=1,m
  199             b(i) = q(i,jp2)+q(i,jp4)
  200             continue
              go to 204
  201             do 203 i=1,m
  202             b(i) = q(i,jp2)
  203             continue
  204         call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 206 i=1,m
  205             q(i,jp3) = q(i,jp3)+b(i)
  206             continue
  207         continue
  208     continue
  209 if (iex3.eq.0) go to 239
      np = n2p3p
      k2 = ks3
      if (nperod.eq.1.and.iex5.eq.0) k3 = k3+np/3
          do 238 kount=1,iex3
          k = np
          np = np/3
          k4 = k3-1
          k3 = k3-np
          k1 = k2+1
          k2 = k2+2*np
          jstart = np
          if (nperod.eq.3) jstart = np+1
              do 237 j=jstart,n,k
              jm1 = j-np
              jp1 = j+np
              jp2 = jp1+np
              if (jm1.eq.0) go to 212
                  do 211 i=1,m
  210             b(i) = q(i,jp1)+q(i,jm1)
  211             continue
              go to 218
  212         if (nperod.eq.0) go to 215
                  do 214 i=1,m
  213             b(i) = q(i,jp1)
  214             continue
              go to 218
  215             do 217 i=1,m
  216             b(i) = q(i,jp1)+q(i,n)
  217             continue
  218         call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 220 i=1,m
  219             q(i,j) = q(i,j)+b(i)
  220             continue
              if (jp2.gt.n) go to 223
                  do 222 i=1,m
  221             b(i) = q(i,j)+q(i,jp2)
  222             continue
              go to 226
  223             do 225 i=1,m
  224             b(i) = q(i,j)
  225             continue
  226         call tri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 228 i=1,m
  227             q(i,j) = q(i,j)+b(i)
  228             continue
              if (jp2.gt.n) go to 231
                  do 230 i=1,m
  229             b(i) = q(i,j)+q(i,jp2)
  230             continue
              go to 234
  231             do 233 i=1,m
  232             b(i) = q(i,j)
  233             continue
  234         call tri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 236 i=1,m
  235             q(i,jp1) = q(i,jp1)+b(i)
  236             continue
  237         continue
  238     continue
  239 if (iex2.eq.0) go to 254
      np = n2pw
          do 253 kount=1,iex2
          k = np
          np = np/2
          k3 = k-1
          jstart = np
          if (nperod.eq.3) jstart = 1+np
              do 252 j=jstart,n,k
              jm1 = j-np
              jp1 = j+np
              if (jm1.ne.0) go to 243
              if (jp1.gt.n) go to 252
              if (nperod.eq.1.or.nperod.eq.2) go to 240
              jm1 = n
              go to 246
  240             do 242 i=1,m
  241             b(i) = q(i,jp1)
  242             continue
              go to 249
  243         if (jp1.le.n) go to 246
                  do 245 i=1,m
  244             b(i) = q(i,jm1)
  245             continue
              go to 249
  246             do 248 i=1,m
  247             b(i) = q(i,jm1)+q(i,jp1)
  248             continue
  249         call tri (np,k3,1,m,mm1,ba,bb,bc,b,twocos,d,w)
                  do 251 i=1,m
  250             q(i,j) = q(i,j)+b(i)
  251             continue
  252         continue
  253     continue
  254 return
      end
c
c
c
      subroutine sorter(cor,index,n)
c     a shell sort routine to sort index and cor down
c     according to the values of cor
      dimension cor(1),index(1)
      igap=n
    5 if (igap .le. 1) return
      igap=igap/2
      imax=n-igap
   10 iex=0
      do 20 i=1,imax
         iplusg=i+igap
         if (cor(i) .ge. cor(iplusg)) goto 20
         save=cor(i)
         cor(i)=cor(iplusg)
         cor(iplusg)=save
         isave=index(i)
         index(i)=index(iplusg)
         index(iplusg)=isave
         iex=1
   20 continue
      if (iex .ne. 0) goto 10
      goto 5
      end
c
c
c
      function xofij(i,j)
          include 'prob2.pert'
      common/cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

c
      if(icall.ne.0)goto 10
      rm=(m+1)/2.0
      rn=(n+1.)/2.
      c1=cos(theta)*hy
      c2=-sin(theta)*hy
      icall=1
   10 continue
      ri=i-rm
      rj=j-rn
      xofij=c1*ri+c2*rj
      return
      end
c
c
c
c
      subroutine alop3d
c
c --- subroutine nlop3d is the main computation loop for this
c --- code.  it advances the vorticity and solves for the new
c --- streamfunction fields.
c
          include 'prob2.pert'
c
c>>>>>>>>>>>>
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c
c>>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/dynam/alpha,beta,rf
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +  topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +  prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/pert/stf(nosp),vtf(nosp)
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/steps/sumt,sumb,cost
c
        real tempb(mn0,k0),tempc(mn0,k0),tempd(mn0,k0)
c
c
c --- new streamfunction--point by point
c
      do 350 ip=1,msq
c
        do j=1,msq
        r(j)=0.0
        enddo
c
        do 345 k=1,kz
        do 345 j=1,kz
           r(j)=r(j)+amat(k,j)*stream(ip,k)
  345   continue
c
        do 340 k=1,kz
           stream(ip,k)=r(k)
  340   continue
c
  350 continue
c
c
c
c --- solve adjoint helmholtz equation level by level.
c
      do 300 k=1,kz
c
c.......impose boundary conditions on r......
c
        do 2503 ip=1,m
 2503   stream(ip,k)=0.0
c
        do 2504 ip=msq-m+1,msq
 2504   stream(ip,k)=0.0
c
c<<<<<<<<<>>>>>>>>>>>
c
        mb=0
c
      dbasin=hy*nm1
      call pwscrt(0,0.,xbasin,mm1,mb,bda,bdb,0.,dbasin,nm1,1,bdc,bdd
     +,eigval(k),stream(1,k),m,1.,ierror,r1)
c<<<<<>>>>>>>>
c
        do j=1,n
        do i=1,m
          ip=(j-1)*m+i
          do k1=1,kz
            vort(ip,k1)=vort(ip,k1)+ainv(k,k1)*stream(ip,k)
          enddo
        enddo
        enddo
c
c
  300 continue
c
c
c ---advance the adjoint vorticity
c
c.....skip nvortq if icalc=0.....
c
      if(icalc.eq.0)goto 8080
      do 100 k=1,kz
        call avortq(k)
  100 continue
 8080 continue
c
c
      return
      end
c
	function spts(a)
c
c	function to convert to mercator 
c
	rad = abs(a)*0.01745329
	sina=sin(rad)
	spts=(7915.704*alog10(tan(0.785398+rad*0.5))
     *	-sina*(23.268932+0.0525*sina*sina))*sign(1.,a)/60.
	return
	end
c
c
c
      function yofij(i,j)
          include 'prob2.pert'
      common/cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

c
      if(icall.ne.0)goto 10
      rm=(m+1)/2.0
      rn=(n+1.)/2.
      c2=cos(theta)*hy
      c1=sin(theta)*hy
      icall=1
   10 continue
      ri=i-rm
      rj=j-rn
      yofij=c1*ri+c2*rj
      return
      end
c
c
c
      subroutine arajac(psi,q,r,m,n)
      dimension psi(m,n),q(m,n),r(m,n)
c
c --- subroutine to compute the arakawa jacobian j(psi,q)=r
c
c --- interior
c
      m1=m-1
      n1=n-1
      do 15 i=1,m
        do 10 j=2,n1
	  im1=i-1
	  ip1=i+1
	  if(i.eq.1)im1=m-1
	  if(i.eq.m)ip1=2
          r(i,j)=((psi(im1,j)-psi(i,j-1))*q(im1,j-1)+
     /      (psi(im1,j-1)+psi(im1,j)-psi(ip1,j-1)-psi(ip1,j))*q(i,j-1)+
     /      (psi(i,j-1)-psi(ip1,j))*q(ip1,j-1)+
     /      (psi(im1,j+1)+psi(i,j+1)-psi(im1,j-1)-psi(i,j-1))*q(im1,j)+
     /      (psi(i,j-1)+psi(ip1,j-1)-psi(i,j+1)-psi(ip1,j+1))*q(ip1,j)+
     /      (psi(i,j+1)-psi(im1,j))*q(im1,j+1)+
     /      (psi(ip1,j)+psi(ip1,j+1)-psi(im1,j)-psi(im1,j+1))*q(i,j+1)+
     /      (psi(ip1,j)-psi(i,j+1))*q(ip1,j+1))/12.
   10   continue
   15 continue
c
c --- left edge
c
      do 20 i=1,m
	im1=i-1
	ip1=i+1
	if(i.eq.1)im1=m-1
	if(i.eq.m)ip1=2
        r(i,1)=((psi(im1,2)+psi(i,2)-2.*psi(i,1))*q(im1,1)+
     /     (2.*(psi(im1,1)-psi(ip1,1)))*q(i,1)+
     /     (2.*psi(i,1)-psi(i,2)-psi(ip1,2))*q(ip1,1)+
     /     (psi(i,2)-psi(im1,1))*q(im1,2)+
     /     (psi(ip1,1)+psi(ip1,2)-psi(im1,1)-psi(im1,2))*q(i,2)+
     /     (psi(ip1,1)-psi(i,2))*q(ip1,2))/12.
   20 continue
c
c
c --- bottom
c
c     do 30 i=2,n1
c       r(m,i)=((psi(m-1,i-1)+psi(m-1,i)-2.*psi(m,i))*q(m,i-1)+
c    /     (2.*(psi(m,i-1)-psi(m,i+1)))*q(m,i)+
c    /     (2.*psi(m,i)-psi(m-1,i)-psi(m-1,i+1))*q(m,i+1)+
c    /     (psi(m-1,i)-psi(m,i-1))*q(m-1,i-1)+
c    /     (psi(m,i+1)+psi(m-1,i+1)-psi(m,i-1)-psi(m-1,i-1))*q(m-1,i)+
c    /     (psi(m,i+1)-psi(m-1,i))*q(m-1,i+1))/12.
c  30 continue
c
c --- right edge
c
      do 40 i=1,m
       im1=i-1
       ip1=i+1
       if(i.eq.1)im1=m-1
       if(i.eq.m)ip1=2
        r(i,n)=((psi(ip1,n1)+psi(i,n1)-2.*psi(i,n))*q(ip1,n)+
     /     (2.*(psi(ip1,n)-psi(im1,n)))*q(i,n)+
     /     (2.*psi(i,n)-psi(i,n1)-psi(im1,n1))*q(im1,n)+
     /     (psi(i,n1)-psi(ip1,n))*q(ip1,n1)+
     /     (psi(im1,n)+psi(im1,n1)-psi(ip1,n)-psi(ip1,n1))*q(i,n1)+
     /     (psi(im1,n)-psi(i,n1))*q(im1,n1))/12.
   40 continue
c
c --- top
c
c     do 50 i=2,n1
c       r(1,i)=((psi(2,i+1)+psi(2,i)-2.*psi(1,i))*q(1,i+1)+
c    /     (2.*(psi(1,i+1)-psi(1,i-1)))*q(1,i)+
c    /     (2.*psi(1,i)-psi(2,i)-psi(2,i-1))*q(1,i-1)+
c    /     (psi(2,i)-psi(1,i+1))*q(2,i+1)+
c    /     (psi(1,i-1)+psi(2,i-1)-psi(1,i+1)-psi(2,i+1))*q(2,i)+
c    /     (psi(1,i-1)-psi(2,i))*q(2,i-1))/12.
c  50 continue
c
	return
c --- corners
c
      r(1,1)=((2.*(psi(1,2)-psi(2,1)))*q(1,1)+
     /   (2.*psi(1,1)-psi(2,2)-psi(1,2))*q(2,1)+
     /   (psi(2,1)-psi(1,2))*q(2,2)+
     /   (psi(2,1)+psi(2,2)-2.*psi(1,1))*q(1,2))/12.
c
      r(m,1)=((2.*(psi(m-1,1)-psi(m,2)))*q(m,1)+
     /   (2.*psi(m,1)-psi(m-1,2)-psi(m-1,1))*q(m,2)+
     /   (psi(m,2)-psi(m-1,1))*q(m-1,2)+
     /   (psi(m,2)+psi(m-1,2)-2.*psi(m,1))*q(m-1,1))/12.
c
      r(m,n)=((2.*(psi(m,n1)-psi(m-1,n)))*q(m,n)+
     /   (2.*psi(m,n)-psi(m-1,n1)-psi(m,n1))*q(m-1,n)+
     /   (psi(m-1,n)-psi(m,n1))*q(m-1,n1)+
     /   (psi(m-1,n)+psi(m-1,n1)-2.*psi(m,n))*q(m,n1))/12.
c
      r(1,n)=((2.*(psi(2,n)-psi(1,n1)))*q(1,n)+
     /   (2.*psi(1,n)-psi(2,n1)-psi(2,n))*q(1,n1)+
     /   (psi(1,n1)-psi(2,n))*q(2,n1)+
     /   (psi(1,n1)+psi(2,n1)-2.*psi(1,n))*q(2,n))/12.
c
      return
      end
c
c
	subroutine bscal
c
        include 'prob2.pert'
c
	common/egy/peint,akeint
	common/l2n/al2s,al2v
	common/egl2/pl,akl,alls,allv
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
	common/swoi/iter
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	call bener
	call bl2n
c
        max=(tmax-tstart)/dt
	if(icalc.ne.0.and.icalc.ne.max)goto 10
	pl=0.
	akl=0.
	alls=0.
	allv=0.
 10     continue
c
	write(89,1)icalc,peint,akeint,al2s,al2v
 1      format(i5,1p,4e10.3)
	if(icalc.eq.0.or.icalc.eq.max)return
	denom=(dt*(peint+pl))
	if(denom.gt.0.)grp=2.*(peint-pl)/denom
	denom=(dt*(akeint+akl))
	if(denom.gt.0.)grk=2.*(akeint-akl)/denom
	denom=(dt*(al2s+alls))
	if(denom.gt.0.)gr2s=2.*(al2s-alls)/denom
	denom=(dt*(al2v+allv))
	if(denom.gt.0.)gr2v=2.*(al2v-allv)/denom
	write(88,1)icalc,grp,grk,gr2s,gr2v
	pl=peint
	akl=akeint
	alls=al2s
	allv=al2v
c
	return
	end
c
c
c
      subroutine filter(zz,m,n,nord)
c
c --- zz = field to be filtered
c ---  m = no of points along a side of the array
c --- nord = order of the filter
c
c --- this subroutine filters a field of data. it works in combination
c --- with subroutine shpiro. filter picks the column or row of data,
c --- calls shpiro, and then adds a factor to the data (as calculated
c --- by shpiro.
c
          include 'prob2.pert'
      dimension zz(1),g(mx0+10)
      if(m.gt.mx0+10.or.n.gt.mx0+10)stop'filerr'
      mm1=m-1
      iodev=(nord+1)/2 -nord/2
      fac=-1.+2.*float(iodev)
      fac=fac/2.**(2*nord)
      do 100 j=1,n
        irow=(j-1)*m
        do 10 i=1,m
          g(i)=zz(irow+i)
   10   continue
c
        call shprox(g,m,nord)
c
        do 50 i=1,m
          zz(irow+i)=zz(irow+i)+fac*g(i)
   50   continue
  100 continue
c
      do 200 i=1,m
        do 110 j=1,n
          g(j)=zz((j-1)*m+i)
  110   continue
c
        call shpiro(g,n,nord)
c
        do 150 j=1,n
          irow=(j-1)*m+i
          zz(irow)=zz(irow)+fac*g(j)
  150   continue
  200 continue
      return
      end
c
      subroutine invmtx (a,na,v,nv,n,d,ip,ier)
      integer na,nv,n,ip(1),ier
      real a(na,n),v(nv,n),d
      data iexmax/75/
  115 format(28h0*matrix singular in invmtx*)
  116 format(34h0*determinant too large in invmtx*)
      ier = ierinv(n,na,nv)
      if (ier .ne. 0) return
      do 102 j=1,n
         ip(j) = 0
         do 101 i=1,n
            v(i,j) = a(i,j)
  101    continue
  102 continue
      d = 1.
      iex = 0
      do 110 m=1,n
         vmax = 0.
         do 104 j=1,n
            if (ip(j) .ne. 0) go to 104
            do 103 i=1,n
               if (ip(i) .ne. 0) go to 103
               vh = abs(v(i,j))
               if (vmax .gt. vh) go to 103
               vmax = vh
               k = i
               l = j
  103       continue
  104    continue
         ip(l) = k
         npm = n+m
         ip(npm) = l
         d = d*v(k,l)
  105    if (abs(d) .le. 1.0) go to 106
         d = d*0.1
         iex = iex+1
         go to 105
  106    continue
         pvt = v(k,l)
         if (m .eq. 1) pvtmx = abs(pvt)
         if (abs(pvt/float(m))+pvtmx .eq. pvtmx) go to 113
         v(k,l) = 1.
         do 107 j=1,n
            hold = v(k,j)
            v(k,j) = v(l,j)
            v(l,j) = hold/pvt
  107    continue
         do 109 i=1,n
            if (i .eq. l) go to 109
            hold = v(i,l)
            v(i,l) = 0.
            do 108 j=1,n
               v(i,j) = v(i,j)-v(l,j)*hold
  108       continue
  109    continue
  110 continue
      m = n+n+1
      do 112 j=1,n
         m = m-1
         l = ip(m)
         k = ip(l)
         if (k .eq. l) go to 112
         d = -d
         do 111 i=1,n
            hold = v(i,l)
            v(i,l) = v(i,k)
            v(i,k) = hold
  111    continue
  112 continue
      if (iex .gt. iexmax) go to 114
      d = d*10.**iex
      return
  113 ier = 33
      print 115
      return
  114 ier = 1
      d = float(iex)
      print 116
      return
      end
c
c
      subroutine rhs(s,prv)
c
c --- this subroutine calculates the right-hand-side (rhs) of
c --- the finite-element vorticity equation
c
          include 'prob2.pert'
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      dimension s(1),prv(1)
c
      fac=-36*dt/(hy**2)
      do 10 ip=1,msq
        s(ip)=fac*s(ip)
   10 continue
c
      do 40 ip=1,msq
        fq=prv(ip)
        prv(ip)=s(ip)
        s(ip)=1.5*s(ip)-0.5*fq
   40 continue
      return
      end
c
      subroutine rhs2(s,prv)
c
c --- this subroutine calculates the right-hand-side (rhs) of
c --- the finite-element vorticity equation
c
          include 'prob2.pert'
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      dimension s(1),prv(1)
c
      fac=-36*dt/(hy**2)
      do 10 ip=1,msq
        s(ip)=fac*s(ip)
   10 continue
c
      do 40 ip=1,msq
        prv(ip)=s(ip)
   40 continue
      return
      end
c
c
      subroutine arelvor(psi,dx2psi,dy2psi,m,n,hy)
c
c --- this routine calculates the fourth order laplacian of psi
c --- at the boundaries and in the interior field
c
c
       dimension psi(1),dx2psi(1),dy2psi(1)
c
       fac=1./12./hy**2
        ip=0
      do 100 j=1,n
      do 110 i=1,m
      ip=ip+1
      if(i.eq.1.or.i.eq.m) go to 50
      if(i.eq.2.or.i.eq.m-1) go to 60
      dx2psi(ip)=fac*(-psi(ip-2)+16.*psi(ip-1)-30.*psi(ip)
     +                 +16.*psi(ip+1)-psi(ip+2))
      go to 70
c
  50  continue
      nn=1
      if(i.eq.m)nn=-1
      dx2psi(ip)=fac*(45.*psi(ip)-10.*psi(ip+nn*5)+61.*psi(ip+nn*4)
     +                -156.*psi(ip+3*nn)+214.*psi(ip+2*nn)
     +                -154.*psi(ip+nn))
       go to 70
 60   continue
      nn=1
      if(i.eq.m-1)nn=-1
      dx2psi(ip)=fac*(-15.*psi(ip)+psi(ip+4*nn)-6.*psi(ip+nn*3)
     +                 +14.*psi(ip+2*nn)-4.*psi(ip+nn)
     +                 +10.*psi(ip-nn))
c
  70  continue
      if(j.eq.1.or.j.eq.n) go to 80
      if(j.eq.2.or.j.eq.n-1) go to 90
      dy2psi(ip)=fac*(-psi(ip-2*m)+16.*psi(ip-m)-30.*psi(ip)
     +                +16.*psi(ip+m)-psi(ip+2*m))
      go to 110
c
  80  continue
      nn=1
      if(j.eq.n)nn=-1
      dy2psi(ip)=fac*(45*psi(ip)-10.*psi(ip+5*nn*m)+61.*psi(ip+4*nn*m)
     +                -156.*psi(ip+3*nn*m)+214.*psi(ip+2*m*nn)
     +                -154.*psi(ip+nn*m))
      go to 110
 90   continue
      nn=1
      if(j.eq.n-1)nn=-1
      dy2psi(ip)=fac*(-15.*psi(ip)+psi(ip+4*nn*m)-6.*psi(ip+3*m*nn)
     +                +14.*psi(ip+2*nn*m)-4.*psi(ip+nn*m)
     +                +10.*psi(ip-nn*m))
 110   continue
 100  continue
      return
      end
c
c
      subroutine getfld(anapsi,anazeta,anatop,time)
c
c
          include 'prob2.pert'
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     + ifsbl,ifeva
      common/field/psi(mn0,k0),zeta(mn0,k0),tlast,day1,tkeep,
     +  curpsi(mn0,k0),curzet(mn0,k0),psiprv(mn0,k0),zetprv(mn0,k0),
     +  topprv(mn0),topcur(mn0),topdens(mn0)
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +     amat(k0,k0)
      common/user/pram(16)
      dimension anapsi(m20,k0),anazeta(m20,k0),anatop(m20)
      dimension therm(m20,k0),dxxpsi(m20),dyypsi(m20)
      ifvort=nint(pram(2))
      iftherm=nint(pram(3))
c
        read(33,1234)mbound,nbound,day1
 1234   format(2i6,1p,1e10.3)
         print 5,mbound,nbound,day1
    5    format( ' mbound,nbound,day1 = ',2i10,e16.5)
        day1=day1+0.03
        nsq=mbound*nbound
        do 10 k=1,kz
          read(33,1235)(anapsi(i,k),i=1,nsq)
   10   continue
        if(ifvort.eq.0) then
          do 20 k=1,kz
            read(33,1235)(anazeta(i,k),i=1,nsq)
   20     continue
        endif
       if(iftop.ne.0) read(33,1235)(anatop(i),i=1,nsq)
 1235  format(1p,8e10.3)
       np=nbound
       mp=mbound 
      write(iout,111)mbound,nbound,time/day1
  111 format(2x,'infld called for m,n,time=',2i3,f5.2)
c
       if(iftherm.eq.1) call thermvt(anapsi,therm,kz,nsq)
       do 500 k=1,kz
       if(ifvort.eq.1) then
       call rel2nd(anapsi(1,k),dxxpsi,dyypsi,mp,np,hy)
c       call arelvor(anapsi(1,k),dxxpsi,dyypsi,mp,np,hy)
       call filter(dxxpsi,mp,np,4)
       call filter(dyypsi,mp,np,4)
       do 400 ip=1,nsq
 400   anazeta(ip,k)=dxxpsi(ip)+dyypsi(ip)
       endif  
       do 600 ip=1,nsq
       anazeta(ip,k)=anazeta(ip,k)+iftherm*therm(ip,k)
       if(k.eq.1.and.iftop.ne.0) anazeta(ip,k)=anazeta(ip,k)+
     1    iftherm*anatop(ip)/hz(1)
600    continue
500    continue
c     if(day1.lt.0.05) rewind 33
      if(day1.lt.0.05) day1=9999.9
c
      return
      end


      function lnblk(string, nchar)
c
cc  this function returns the character position
cc  of the last non-blank character in the string.
cc  characters are numbered starting with 1 at the
cc  left of the first word and are packed ten word
cc  with blank fill. nchar tells lnblk where to start
cc  the search for the last non blank character
c
        character string*255
        do 10 i=nchar,1,-1
        if(string(i:i).ne.' ')goto 11
10      continue
        lnblk=0
        return
11      lnblk=i
      return
      end
c
        subroutine thermvt(psi,tvort,kz,msq)
          include 'prob2.pert'
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
        dimension tvort(mn0,k0),psi(mn0,k0),vtm(k0,k0)
c
       do 5 k=1,kz
       do 6 kk=1,kz
       vtm(k,kk)=0.0
 6    continue
      do 7 ip=1,msq
      tvort(ip,k)=0.
 7    continue
 5    continue
c
c
      do 20 i=1,kz
        ip=0
        do 19 j=1,kz
          tmp=0.0
          do 18 k=1,kz
            tmp=tmp+amat(i,k)*eigval(k)*ainv(k,j)
   18     continue
          vtm(i,j)=tmp
   19   continue
   20 continue
c
c --- now calculate the thermal vorticity
c
       do 500 k=1,kz
        do 170 kk=1,kz
          tmp=vtm(k,kk)
          do 165 ip=1,msq
            tvort(ip,k)=tvort(ip,k)+tmp*psi(ip,kk)
  165     continue
  170   continue
 500    continue
c
c
        return
        end
c
	function gasdev(idum)
	data iset/0/
	if (iset.eq.0) then
1		v1=2.*ran1(idum)-1.
		v2=2.*ran1(idum)-1.
		r=v1**2+v2**2
		if(r.ge.1.)go to 1
		fac=sqrt(-2.*log(r)/r)
		gset=v1*fac
		gasdev=v2*fac
		iset=1
	else
		gasdev=gset
		iset=0
	endif
	return
	end
c
	function ran1(idum)
	dimension r(97)
	parameter (m1=259200,ia1=7141,ic1=54773,rm1=1./m1)
	parameter (m2=134456,ia2=8121,ic2=28411,rm2=1./m2)
	parameter (m3=243000,ia3=4561,ic3=51349)
	data iff /0/
	if (idum.lt.0.or.iff.eq.0) then
		iff=1
		ix1=mod(ic1-idum,m1)
		ix1=mod(ia1*ix1+ic1,m1)
		ix2=mod(ix1,m2)
		ix1=mod(ia1*ix1+ic1,m1)
		ix3=mod(ix1,m3)
		do 10 j=1,97
			ix1=mod(ia1*ix1+ic1,m1)
			ix2=mod(ia2*ix2+ic2,m2)
			r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 10             continue
		idum=1
	endif
	ix1=mod(ia1*ix1+ic1,m1)
	ix2=mod(ia2*ix2+ic2,m2)
	ix3=mod(ia3*ix3+ic3,m3)
	j=1+(97*ix3)/m3
	if(j.gt.97.or.j.lt.1)pause
	ran1=r(j)
	r(j)=(float(ix1)+float(ix2)*rm2)*rm1
	return
	end
c
	subroutine gregy
c
        include 'prob2.pert'
c
	common/grate/gre(mn0),grip(mn0),tezm(mn0),aipzm(mn0)
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
c
	real pe(mn0,k01),ke(mn0,k0)
     &    ,tez(mn0),pez(mn0),kez(mn0),aipz(mn0)
c
	do 6 ip=1,mn0
	pez(ip)=0.0
	kez(ip)=0.0
	tez(ip)=0.0
 6      aipz(mn0)=0.0
c
c.......calc. pe......
c
	do 1 j=2,n-1
	do 1 i=2,m-1
	ip=(j-1)*m+i
	if(iftop.ne.0)pe(ip,1)=0.5*scsig*sigz(1)*(topden(ip)**2)
	if(ifbot.ne.0)pe(ip,kzp1)=0.5*scsig*sigz(kzp1)*(botden(ip)**2)
	if(iftop.eq.0)pe(ip,1)=0.0
	if(ifbot.eq.0)pe(ip,kzp1)=0.0
c
	do 2 k=2,kz
	depth=0.5*(hz(k-1)+hz(k))
	psiz=(stream(ip,k-1)-stream(ip,k))/depth
c	psibz=(bstm(ip,k-1)-bstm(ip,k))/depth
 2      pe(ip,k)=0.5*scsig*sigz(k)*(psiz**2)
 1      continue
c
c.......calc. ke......
c
	do 3 j=2,n-1
	do 3 i=2,m-1
	ip=(j-1)*m+i
	do 4 k=1,kz
	u=0.5*(stream(ip-m,k)-stream(ip+m,k))/hy
	v=0.5*(stream(ip+1,k)-stream(ip-1,k))/hy
c	ub=0.5*(bstm(ip-m,k)-bstm(ip+m,k))/hy
c	vb=0.5*(bstm(ip+1,k)-bstm(ip-1,k))/hy
 4      ke(ip,k)=0.5*(u**2+v**2)
 3      continue
c
c.......calc. integrals of pe and ke.......
c
	do 11 j=2,n-1
	do 11 i=2,m-1
	do 11 k=1,kz
	ip=(j-1)*m+i
 11     pez(ip)=pez(ip)+0.5*hz(k)*(pe(ip,k)+pe(ip,k+1))
c
	do 17 j=2,n-1
	do 17 i=2,m-1
	do 17 k=1,kz-1
	ip=(j-1)*m+i
	depth=0.5*(hz(k)+hz(k+1))
 17     kez(ip)=kez(ip)+0.5*depth*(kez(k)+kez(k+1))
c
	do 18 ip=1,mn0
 18     tez(ip)=pez(ip)+kez(ip)
c
	do 20 j=2,n-1
	do 20 i=2,m-1
	do 20 k=1,kz
	ip=(j-1)*m+i
 20     aipz(ip)=aipz(ip)+stream(ip,k)*vort(ip,k)
c
c.........calc growth rates.......
c
	if(icalc.ne.1)then
	do 40 ip=1,mn0
	denome=dt*(tez(ip)+tezm(ip))
	denomi=dt*(aipz(ip)+aipzm(ip))
	if(denome.ne.0)gre(ip)=gre(ip)+2.*(tez(ip)-tezm(ip))/denome
 40	if(denomi.ne.0)grip(ip)=grip(ip)+2.*(aipz(ip)-aipzm(ip))/denomi
	endif
c
c.....swop arrays.....
c
	do 50 ip=1,mn0
	tezm(ip)=tez(ip)
 50     aipzm(ip)=aipz(ip)
c
	return
	end
c
c
c
      subroutine tridp (kstart,kstop,ising,m,mm1,a,b,c,y,twocos,d,w)
      dimension        a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     1                 twocos(1)  ,d(1)       ,w(1)
c
c uro is a machine dependent constant.  it is the smallest floating
c point number such that 1 + uro is greater than 1 in the
c arithmetic of this computer.
c
      uro = 1./2.**20
  101     do 106 k=kstart,kstop
          x = -twocos(k)
          d(1) = c(1)/(b(1)+x)
          w(1) = a(1)/(b(1)+x)
          y(1) = y(1)/(b(1)+x)
          bm = b(m)
          z = c(m)
              do 103 i=2,mm1
              den = b(i)+x-a(i)*d(i-1)
              d(i) = c(i)/den
              w(i) = -a(i)*w(i-1)/den
              y(i) = (y(i)-a(i)*y(i-1))/den
              y(m) = y(m)-z*y(i-1)
              bm = bm-z*w(i-1)
  102         z = -z*d(i-1)
  103         continue
          d(mm1) = d(mm1)+w(mm1)
          z = a(m)+z
          den = bm+x-z*d(mm1)
          y(m) = y(m)-z*y(m-1)
          y(m) = y(m)/den
          y(mm1) = y(mm1)-d(mm1)*y(m)
              do 105 ip=2,mm1
              i = m-ip
  104         y(i) = y(i)-d(i)*y(i+1)-w(i)*y(m)
  105         continue
  106     continue
  107 if (ising.eq.1) return
      d(1) = c(1)/(b(1)-2.)
      w(1) = a(1)/(b(1)-2.)
      y(1) = y(1)/(b(1)-2.)
      bm = b(m)
      z = c(m)
          do 109 i=2,mm1
          den = b(i)-2.-a(i)*d(i-1)
          d(i) = c(i)/den
          w(i) = -a(i)*w(i-1)/den
          y(i) = (y(i)-a(i)*y(i-1))/den
          y(m) = y(m)-z*y(i-1)
          bm = bm-z*w(i-1)
  108     z = -z*d(i-1)
  109     continue
      d(mm1) = d(mm1)+w(mm1)
      z = a(m)+z
      den = bm-2.-z*d(mm1)
      y(m) = y(m)-z*y(m-1)
      if (abs(den).gt.20.*(m**3+m**2*abs(a(m)))*uro) go to 111
  110 y(m) = 0.
      go to 112
  111 y(m) = y(m)/den
  112 y(mm1) = y(mm1)-d(mm1)*y(m)
          do 114 ip=2,mm1
          i = m-ip
  113     y(i) = y(i)-d(i)*y(i+1)-w(i)*y(m)
  114     continue
      return
      end
      subroutine atrdi (n,a,b,c,y,x,ks,work)
      dimension       a(n)       ,b(n)       ,c(n)       ,y(n)       ,
     1                x(n)       ,work(1)
  103 format(46h0array dimension n not greater than 2 in atrdi)
  104 format('0hks not equal to 0 or 1 in call to atrdi')
      if (n .gt. 2) go to 101
      print 103
      return
  101 continue
      if (ks.eq.0 .or. ks.eq.1) go to 102
      print 104
      return
  102 continue
      call atrdslv (n,a,b,c,y,x,ks,work(1),work(n))
      return
      end
c
      subroutine atrdslv (n,a,b,c,y,x,ks,alpha,gamma)
      dimension       a(n)       ,b(n)       ,c(n)       ,y(n)       ,
     1                x(n)       ,alpha(1)   ,gamma(1)
      nm1 = n-1
      if (ks .eq. 1) go to 102
      alpha(1) = 1./b(1)
      gamma(1) = c(1)*alpha(1)
      do 101 i=2,nm1
         alpha(i) = 1./(b(i)-a(i)*gamma(i-1))
         gamma(i) = c(i)*alpha(i)
  101 continue
  102 continue
c
      do 104 j=nm1,1,-1
         jb = n-j
         x(jb+1)=x(jb+1)-gamma(jb)*x(jb)
  104 continue
c
      y(n)=x(n)/(b(n)-a(n)*gamma(nm1))
      x(nm1)=x(nm1)-a(n)*x(n)/(b(n)-a(n)*gamma(nm1))
c
      do 103 i=nm1,2,-1
         y(i)=alpha(i)*x(i)
         x(i-1)=x(i-1)-a(i)*x(i)*alpha(i)
  103 continue
c
      y(1) = x(1)*alpha(1)
c
      return
	end
c
c
c
c
      subroutine amsolve
c
c --- this subroutine solves the adams-bashforth finite-element
c --- equations for the interior of the field
c
c --- array r (entering) = right-hand side of eq's
c --- array r (leaving) = solution
c
          include 'prob2.pert'
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/cnn/n,nm1,nm2,np1
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva

      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      dimension aux(mx0)
      dimension va(mx0),vb(mx0)
c
      do i=1,mx0
      aux(i)=0.0
      enddo
c
      do 45 j=1,n
	va(j)=van(j)
	vb(j)=vbn(j)
   45 continue
c
      do 80 i=mm1,2,-1
c
        i1=m+i
c
        do 70 j=nm2,1,-1
          ip=i1+(j-1)*m
          aux(j)=r(ip)
   70   continue
c
        do 60 j=nm2,2,-1
          jj=nm1-j
          aux(jj+1)=aux(jj+1)-aux(jj)/vb(jj)
          aux(jj)=aux(jj)/vb(jj)
   60   continue
        aux(nm2)=aux(nm2)/vb(nm2)
c
        do 50 j=nm2,2,-1
          aux(j-1)=aux(j-1)-va(j)*aux(j)
          r1(i1+(j-1)*m)=aux(j)
   50   continue
c
        r1(i1)=aux(1)
c
   80 continue
c
      do 5 j=1,m
	va(j)=vam(j)
	vb(j)=vbm(j)
    5 continue
c
      do 40 j=nm1,2,-1
c
        i1=(j-1)*m+2
c
        do 30 i=mm2,1,-1
          aux(i)=r1(i1+i-1)
   30   continue
c
        do 20 i=mm2,2,-1
          ii=mm1-i
          aux(ii+1)=aux(ii+1)-aux(ii)/vb(ii)
          aux(ii)=aux(ii)/vb(ii)
   20   continue
        aux(mm2)=aux(mm2)/vb(mm2)
c
        do 10 i=mm2,2,-1
          aux(i-1)=aux(i-1)-va(i)*aux(i)
          r(i1+i-1)=aux(i)
   10   continue
        r(i1)=aux(1)
c
   40 continue
c
      return
      end
c
c
c
c
      subroutine nvagenrl(s,zp,klevel)
c
          include 'prob2.pert'
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/bsolv/ab(mx0),bb(mx0),cb(mx0),workb(mnwork)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
      dimension s(1),zp(1)
      dimension aac(5),fx(mx0),fy(mx0),inout(mx0)
c
      real r1t(mn0)
c
      do i=1,mn0
      r1t(i)=0.0
      enddo
c
      do i=1,5
      aac(i)=0.0
      enddo
c
      do 120 iside=1,4
        ib=icorn(iside)
        i1=ib+nast(iside)
        i2=ib+mast(iside)
        if(s(i2).gt.s(i1)) go to 110
        aac(iside)=aac(iside)+0.25*r(ib)
        r(i1)=r(i1)-0.5*r(ib)
        r(i2)=r(i2)-0.5*r(ib)
        go to 120
  110   continue
        r1t(ib)=r1t(ib)+r(ib)
        zp(ib)=zp(ib)-r(ib)
c        zp(ib)=0.0
c        r(ib)=0.0
  120 continue
c
      do 100 iside=1,4
c
c.......set up coefficient arrays for atrdi.......
c
        i1=icorn(iside)
        i2=i1+mast(iside)
        i3=i2+mast(iside)
c
c --- determine inflow points
c
	ist=mast(iside)
	if(ist.eq.1.or.ist.eq.-1)mn2=mm2
	if(ist.eq.m.or.ist.eq.-m)mn2=nm2
        do 10 ii=1,mn2
          inout(ii)=0
          if(s(i3).gt.s(i1)) inout(ii)=1
          i1=i2
          i2=i3
          i3=i3+mast(iside)
   10   continue
c
c......set up fx array for input to atrdi......
c
        ib=icorn(iside)
        do 60 iz=1,mn2
          ib=ib+mast(iside)
          fx(iz)=r(ib)
   60   continue
c
        ib=icorn(iside)+mast(iside)
        if(inout(1).eq.1) go to 15
        bb(1)=7.
        cb(1)=2.
        go to 20
   15   continue
        bb(1)=1.
        cb(1)=0.
c        fx(1)=0.
   20   continue
	if(ist.eq.1.or.ist.eq.-1)mn3=m-3
	if(ist.eq.m.or.ist.eq.-m)mn3=n-3
        do 40 iz=2,mn3
          ib=ib+mast(iside)
          if(inout(iz).eq.1) go to 35
          bb(iz)=8.
          ab(iz)=2.
          cb(iz)=2.
          go to 40
   35     continue
          bb(iz)=1.
          ab(iz)=0.
          cb(iz)=0.
c          fx(iz)=0.
   40   continue
        ib=ib+mast(iside)
        if(inout(mn2).eq.1) go to 45
        bb(mn2)=7.
        ab(mn2)=2.
        go to 50
   45   continue
        bb(mn2)=1.
        ab(mn2)=0.
c        fx(mn2)=0.
   50   continue
c
c --- call the ncar tri-diagonal matrix inverter
c
        call atrdi(mn2,ab,bb,cb,fy,fx,0,workb)
c
        ib=icorn(iside)+mast(iside)
        if(inout(1).eq.1) go to 150
        i5=ib+nast(iside)
        i6=i5+mast(iside)
        r(ib)=fy(1)
        r(i5)=r(i5)-4.*fy(1)
        r(i6)=r(i6)-fy(1)
        aac(iside)=aac(iside)-0.5*fy(1)
        go to 200
  150   continue
        r1t(ib)=r1t(ib)+fy(1)
        zp(ib)=zp(ib)-fy(1)
c        zp(ib)=0.0
        r(ib)=0.0
  200   continue
	if(ist.eq.1.or.ist.eq.-1)mn3=m-3
	if(ist.eq.m.or.ist.eq.-m)mn3=n-3
        do 400 iz=2,mn3
          ib=ib+mast(iside)
          if(inout(iz).eq.1) go to 350
          i5=ib+nast(iside)
          i4=i5-mast(iside)
          i6=i5+mast(iside)
          r(ib)=fy(iz)
          r(i5)=r(i5)-4.*fy(iz)
          r(i4)=r(i4)-fy(iz)
          r(i6)=r(i6)-fy(iz)
          go to 400
  350     continue
          r1t(ib)=r1t(ib)+fy(iz)
          zp(ib)=zp(ib)-fy(iz)
c          zp(ib)=0.0
          r(ib)=0.0
  400   continue
        ib=ib+mast(iside)
        if(inout(mn2).eq.1) go to 450
        i5=ib+nast(iside)
        i4=i5-mast(iside)
        r(ib)=fy(mn2)
        r(i5)=r(i5)-4.*fy(mn2)
        r(i4)=r(i4)-fy(mn2)
        aac(iside+1)=aac(iside+1)-0.5*fy(mn2)
        go to 500
  450   continue
        r1t(ib)=r1t(ib)+fy(mn2)
        zp(ib)=zp(ib)-fy(mn2)
c        zp(ib)=0.0
        r(ib)=0.0
  500   continue
c
  100 continue
c
      aac(1)=aac(1)+aac(5)
c
      do 5 iside=1,4
        ib=icorn(iside)
        ip=ncorn(iside)
        r(ib)=aac(iside)
        r(ip)=r(ip)-aac(iside)
    5 continue
c
      return
      end
c
c
c
c
      subroutine apois (iflag,nperod,n,mperod,m,a,b,c,idimy,y,w,zz)
      external         atrid       ,atridp
      dimension        w(1)       ,b(1)       ,a(1)       ,c(1)       ,
     1                 y(idimy,1),zz(1)
      iwdim1 = 5.25*n+1
      iwdim2 = iwdim1+m
      iwdim3 = iwdim2+m
      iwdim4 = iwdim3+m
      iwdim5 = iwdim4+m
          do 106 i=1,m
          a(i) = -a(i)
          c(i) = -c(i)
          w(iwdim5+i-1) = -b(i)+2.
  106     continue
      if (mperod.eq.0) go to 103
      call apoisgn (nperod,n,iflag,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),zz,atrid)
      go to 104
  103 call apoisgn (nperod,n,iflag,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),zz,atridp)
  104     do 102 i=1,m
          a(i) = -a(i)
          c(i) = -c(i)
  102     continue
       return
      end
c
      subroutine apoisgn (nperod,n,iflag,m,ba,bb,bc,idimq,q,twocos,b,t,
     1                   d,w,zz,atri)
	external atri
      dimension        ba(1)      ,bb(1)      ,bc(1)      ,q(idimq,1) ,
     1                 twocos(1)  ,b(1)       ,t(1)       ,d(1)       ,
     2                 w(1),   zz(1)
     3                 ,amat(56,56)
c
c      print *,'m=',m,' n=',n
c      print *,'ba=',(ba(i),i=1,m)
c      print *,'bb=',(bb(i),i=1,m)
c      print *,'bc=',(bc(i),i=1,m)
c
      if (iflag.ne.0) go to 101
      call poinit (nperod,n,iex2,iex3,iex5,n2pw,n2p3p,n2p3p5,ks3,ks5,
     1             twocos)
  101 mm1 = m-1
c
c.....set up some more coeficients.......
c
c
	l2=iex2
	if(iex3+iex5.eq.0.and.nperod.eq.1)l2=l2-1
	if(l2.le.0)l2=0
c
	l3=iex3
	if(iex5.eq.0.and.nperod.eq.1)l3=l3-1
	if(l3.le.0)l3=0
c
	l5=iex5
	if(nperod.eq.1)l5=l5-1
	if(l5.le.0)l5=0
c
	npl2=2**l2
	npl2l3=npl2*(3**l3)
	npl2l3l5=npl2l3*(5**l5)
c
	np=1
c
      if (iex2.eq.0) go to 254
c
c	print *,'iex2=',iex2,' n2pw=',n2pw
c
          do 253 kount=1,iex2
c
	  k=2*np
	  k3=k-1
c	  print *,'kount=',kount,' np=',np,' k=',k,' k3=',k3
          jstart = np
          if (nperod.eq.3) jstart = 1+np
              do 252 j=jstart,n,k
c
	do 2000 i=1,m
 2000   b(i)=q(i,j)
c
         call atri (np,k3,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
              jm1 = j-np
              jp1 = j+np
c
c	print *,'At start j=',j,' jm1=',jm1,' jp1=',jp1
c
c	print *,'j=',j,' jm1=',jm1,' jp1=',jp1
c
              if (jm1.ne.0) go to 243
              if (jp1.gt.n) go to 252
              if (nperod.eq.1.or.nperod.eq.2) go to 240
              jm1 = n
              go to 246
  240             do 242 i=1,m
                  q(i,jp1)=q(i,jp1)+b(i)
  242             continue
c	      print *,'j=',j,' jp1=',jp1
	      goto 252
  243         if (jp1.le.n) go to 246
c	      print *,'j=',j,' jm1=',jm1
                  do 245 i=1,m
                  q(i,jm1)=q(i,jm1)+b(i)
  245             continue
              go to 252
  246             do 248 i=1,m
	          q(i,jm1)=q(i,jm1)+b(i)
                  q(i,jp1)=q(i,jp1)+b(i)
c	print *,'j=',j,' i=',i,' jm1=',jm1,' jp1=',jp1,
c    1           ' qh=',qh(i,j),' q2=',q2(i,jp1)
  248             continue
c	      print *,'j=',j,' jm1=',jm1,' jp1=',jp1
c
 252    continue
c
	np=2*np
c
 253    continue
c
  254 continue
c
c<<<<<<<<<Compute k2 and k3 from sections 7,6,5,3>>>>>>>>
c
c
c.......compute integer counters.......
c
	k1=0
	k2=0
	k3=0
	k4=0
	k5=0
	k6=0
	k7=0
	k8=0
c
      npp = 1
      if (iex2.eq.0) go to 1180
      l = iex2
      if (iex3+iex5.ne.0) go to 1080
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 8888
 1080    do 1170 kount=1,l
          k = npp
          npp = 2*npp
          k2 = npp-1
          k3 = npp
 1170     k4 = 2*npp-1
c
 1180 if (iex3.eq.0) go to 1340
      l = iex3
      if (iex5.ne.0) go to 1190
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 8888
 1190 k2 = npp-1
          do 1330 kount=1,l
          k = npp
          npp = 3*npp
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
 1330     k4 = k2+npp
c
 1340 l = iex5
      if (nperod.eq.1) l = l-1
      if (l.le.0) go to 8888
      k2 = (n2pw+n2p3p)/2-1
          do 1550 kount=1,l
          k = npp
          npp = 5*npp
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
 1550     k4 = k2+npp
c
 8888     continue
c
      if (iex5.eq.0) go to 2090
      npp = n2p3p5
      k8 = ks5
      if (nperod.eq.1) k3 = k3+npp/5
          do 2080 kount=1,iex5
          k = npp
          npp = npp/5
          k4 = k3-1
          k3 = k3-npp
          k1 = k8+1
          k2 = k8+2*npp
          k5 = k2+1
          k6 = k2+4*npp
          k7 = k6+1
 2080     k8 = k6+2*npp
 2090     continue
c
      if (iex3.eq.0) go to 2390
      npp = n2p3p
      k2 = ks3
      if (nperod.eq.1.and.iex5.eq.0) k3 = k3+npp/3
          do 2380 kount=1,iex3
          k = npp
          npp = npp/3
          k4 = k3-1
          k3 = k3-npp
          k1 = k2+1
 2380     k2 = k2+2*npp
c
 2390     continue
c
c
      if (iex3.eq.0) go to 239
c
	np=n2pw
c
          do 238 kount=1,iex3
c
          k = 3*np
c
	k4=k3+np-1
	k1=k2-2*np+1
c
          jstart = np
          if (nperod.eq.3) jstart = np+1
              do 237 j=jstart,n,k
              jm1 = j-np
              jp1 = j+np
              jp2 = jp1+np
c
               do 236 i=1,m
 236           b(i)=q(i,jp1)
c
           call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 228 i=1,m
 228    q(i,j)=q(i,j)+b(i)
c
	if(jp2.gt.n)goto 234
	do 230 i=1,m
 230    q(i,jp2)=q(i,jp2)+b(i)
c
 234    do 2300 i=1,m
 2300   b(i)=q(i,j)
c
           call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 220 i=1,m
 220    q(i,j)=q(i,j)+b(i)
c
	if(jp2.gt.n)goto 223
c
	do 222 i=1,m
 222    q(i,jp2)=q(i,jp2)+b(i)
c
 223    continue
c
        do 2200 i=1,m
 2200   b(i)=q(i,j)
c
           call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	if(jm1.eq.0)goto 212
c
	do 211 i=1,m
	q(i,jp1)=q(i,jp1)+b(i)
 211    q(i,jm1)=q(i,jm1)+b(i)
	goto 218
c
 212    if(nperod.eq.0)goto 215
	do 214 i=1,m
 214    q(i,jp1)=q(i,jp1)+b(i)
	goto 218
 215    do 217 i=1,m
	q(i,jp1)=q(i,jp1)+b(i)
 217    q(i,n)=q(i,n)+b(i)
c
 218    continue
c
 237    continue
c
	k3=k3+np
	k2=k2-2*np
	np=3*np
c
 238    continue
c
 239    continue
c
c
c.......compute integer counters.......
c
	k1=0
	k2=0
	k3=0
	k4=0
	k5=0
	k6=0
	k7=0
	k8=0
c
      npp = 1
      if (iex2.eq.0) go to 1181
      l = iex2
      if (iex3+iex5.ne.0) go to 1081
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 8889
 1081    do 1171 kount=1,l
          k = npp
          npp = 2*npp
          k2 = npp-1
          k3 = npp
 1171     k4 = 2*npp-1
c
 1181 if (iex3.eq.0) go to 1341
      l = iex3
      if (iex5.ne.0) go to 1191
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 8889
 1191 k2 = npp-1
          do 1331 kount=1,l
          k = npp
          npp = 3*npp
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
 1331     k4 = k2+npp
c
 1341 l = iex5
      if (nperod.eq.1) l = l-1
      if (l.le.0) go to 8889
      k2 = (n2pw+n2p3p)/2-1
          do 1551 kount=1,l
          k = npp
          npp = 5*npp
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
 1551     k4 = k2+npp
c
 8889     continue
c
      if (iex5.eq.0) go to 2091
      npp = n2p3p5
      k8 = ks5
      if (nperod.eq.1) k3 = k3+npp/5
          do 2081 kount=1,iex5
          k = npp
          npp = npp/5
          k4 = k3-1
          k3 = k3-npp
          k1 = k8+1
          k2 = k8+2*npp
          k5 = k2+1
          k6 = k2+4*npp
          k7 = k6+1
 2081     k8 = k6+2*npp
 2091     continue
c
c
      if (iex5.eq.0) go to 209
c
	np=n2p3p
c
          do 208 kount=1,iex5
c
	k=5*np
c
          k4 = k3+np-1
	  kp8=k8-8*np
          k1 = kp8+1
          k2 = kp8+2*np
          k5 = k2+1
          k6 = k2+4*np
          k7 = k6+1
c
          jstart = np
          if (nperod.eq.3) jstart = 1+np
              do 207 j=jstart,n,k
              jm1 = j-np
              jp1 = j+np
              jp2 = jp1+np
              jp3 = jp2+np
              jp4 = jp3+np
c
	do 206 i=1,m
 206    b(i)=q(i,jp3)
c
           call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 203 i=1,m
 203    q(i,jp2)=q(i,jp2)+b(i)
c
	if(jp4.gt.n)goto 201
	do 200 i=1,m
 200	q(i,jp4)=q(i,jp4)+b(i)
c
 201    continue
c
	do 198 i=1,m
 198    b(i)=q(i,jp1)
c
              call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 1960 i=1,m
	q(i,j)=q(i,j)+b(i)
	q(i,jp2)=q(i,jp2)+b(i)
 1960   b(i)=q(i,j)
c
              call atri (k7,k8,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 195 i=1,m
	q(i,jp2)=q(i,jp2)+b(i)
	q(i,j)=q(i,j)+b(i)
 195    b(i)=q(i,jp2)
c
              call atri (k5,k6,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 192 i=1,m
	q(i,jp2)=q(i,jp2)+b(i)
	q(i,j)=q(i,j)+b(i)
 192    b(i)=q(i,jp2)
c
           call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 186 i=1,m
	q(i,j)=q(i,j)+b(i)
 186    q(i,jp2)=2.*b(i)+q(i,jp2)
c
	if(jp4.gt.n)goto 187
c
	do 189 i=1,m
 189    q(i,jp4)=q(i,jp4)+b(i)
c
 187    do 182 i=1,m
 182    b(i)=q(i,jp2)
c
              call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 1820 i=1,m
	q(i,jp1)=q(i,jp1)+b(i)
	q(i,jp3)=q(i,jp3)+b(i)
 1820   b(i)=q(i,j)
c
           call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	if(jm1.ne.0)goto 177
	if(nperod.eq.0)goto 174
	do 173 i=1,m
 173    q(i,jp1)=q(i,jp1)+b(i)
	goto 180
c
 174    do 176 i=1,m
	q(i,jp1)=q(i,jp1)+b(i)
 176    q(i,n)=q(i,n)+b(i)
	goto 180
 177    do 179 i=1,m
	q(i,jp1)=q(i,jp1)+b(i)
 179    q(i,jm1)=q(i,jm1)+b(i)
c
 180    continue
c
 207	continue
c
	k3=k3+np
	k8=kp8
	np=5*np
c
 208	continue
c
 209    continue
c
c
c.......compute integer counters.......
c
	k1=0
	k2=0
	k3=0
	k4=0
	k5=0
	k6=0
	k7=0
	k8=0
c
      npp = 1
      if (iex2.eq.0) go to 1182
      l = iex2
      if (iex3+iex5.ne.0) go to 1082
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 8889
 1082    do 1172 kount=1,l
          k = npp
          npp = 2*npp
          k2 = npp-1
          k3 = npp
 1172     k4 = 2*npp-1
c
 1182 if (iex3.eq.0) go to 1342
      l = iex3
      if (iex5.ne.0) go to 1192
      if (nperod.eq.1) l = l-1
      if (l.eq.0) go to 8887
 1192 k2 = npp-1
          do 1332 kount=1,l
          k = npp
          npp = 3*npp
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
 1332     k4 = k2+npp
c
 1342 l = iex5
      if (nperod.eq.1) l = l-1
      if (l.le.0) go to 8887
      k2 = (n2pw+n2p3p)/2-1
          do 1552 kount=1,l
          k = npp
          npp = 5*npp
          k1 = k2+1
          k2 = k2+k
          k3 = k2+1
 1552     k4 = k2+npp
c
 8887     continue
c
	np=npl2l3l5
c
      if (nperod.eq.1.or.nperod.eq.2) go to 171
      if (nperod.eq.0) go to 161
          do 170 i=1,m
          b(i) = q(i,1)
  170     continue
      call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
          do 168 i=1,m
          q(i,n) = q(i,n)+2.*b(i)
  168     continue
c
  161 continue
c
	do 166 i=1,m
 166    b(i)=q(i,n)
c
	if(nperod.eq.0)goto 1610
c
      call atri (k4+1,k4+2*np-1,2,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 160 i=1,m
	q(i,n)=q(i,n)+4.*b(i)
 160    b(i)=q(i,n)
c
      call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 158 i=1,m
 158    q(i,1)=2.*b(i)+q(i,1)
c
	goto 171
c
 1610   call atri (k4+1,k4+np-1,2,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 163 i=1,m
 163    q(i,n)=q(i,n)+2.*b(i)
c
 171    continue
c
c
	if(l2.ne.0.and.l3.ne.0.and.l5.ne.0)then
c
      l = iex5
      if (nperod.eq.1) l = l-1
c
	np=npl2l3l5
c
	npp=npl2l3
	k2=(n2pw+n2p3p)/2-1
	do 2 kount=1,l5
	k2=k2+npp
 2      npp=5*npp
c
          do 155 kount=1,l
c
	k=np/5
c
          k3 = k2+1
          k4 = k2+np
	k1=k2-k+1
c
          jstart = np
          if (nperod.eq.3) jstart = 1
              do 154 j=jstart,n,np
              if (j.ne.1) go to 136
              jm1 = j+k
              jm2 = jm1+k
              jm3 = jm2+k
              jm4 = jm3+k
              go to 138
  136         jm1 = j-k
              jm2 = jm1-k
              jm3 = jm2-k
              jm4 = jm3-k
              if (j.ne.n) go to 138
              if (nperod.eq.0) go to 137
              jp1 = jm1
              jp2 = jm2
              jp3 = jm3
              jp4 = jm4
              go to 139
  137         jp1 = k
              jp2 = jp1+k
              jp3 = jp2+k
              jp4 = jp3+k
              go to 139
  138         jp1 = j+k
              jp2 = jp1+k
              jp3 = jp2+k
              jp4 = jp3+k
c
 139     continue
c
	do 153 i=1,m
 153    b(i)=5.*q(i,j)
c
              call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 151 i=1,m
	q(i,jm1)=q(i,jm1)+2.*b(i)
	q(i,jm3)=q(i,jm3)+b(i)
	q(i,jp1)=q(i,jp1)+2.*b(i)
 151    q(i,jp3)=q(i,jp3)+b(i)
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 149 i=1,m
	temp=q(i,j)
	t(i)=-1.*b(i)
	q(i,j)=4.*b(i)+temp
	q(i,jm2)=q(i,jm2)+3.*b(i)
	q(i,jm4)=q(i,jm4)+b(i)
	q(i,jp2)=q(i,jp2)+3.*b(i)
	q(i,jp4)=q(i,jp4)+b(i)
 149    b(i)=temp
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 147 i=1,m
	q(i,jm1)=q(i,jm1)+b(i)
 147    q(i,jp1)=q(i,jp1)+b(i)
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 145 i=1,m
	q(i,j)=q(i,j)+2.*b(i)
	q(i,jm2)=q(i,jm2)+b(i)
	q(i,jp2)=q(i,jp2)+b(i)
 145    b(i)=t(i)+b(i)
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 143 i=1,m
	q(i,jm1)=q(i,jm1)+3.*b(i)
	q(i,jm3)=q(i,jm3)+b(i)
	q(i,jp1)=q(i,jp1)+3.*b(i)
 143    q(i,jp3)=q(i,jp3)+b(i)
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 141 i=1,m
	q(i,j)=q(i,j)+6.*b(i)
	q(i,jm2)=q(i,jm2)+4.*b(i)
	q(i,jm4)=q(i,jm4)+b(i)
	q(i,jp2)=q(i,jp2)+4.*b(i)
 141    q(i,jp4)=q(i,jp4)+b(i)
 154    continue
c
	k2=k2-k
	np=np/5
c	
 155    continue
c
	endif
c
	if(l2.ne.0.and.l3.ne.0)then
c
      if (iex3.eq.0) go to 134
c
	np=npl2l3
c
      l = iex3
      if (iex5.ne.0) go to 119
      if (nperod.eq.1) l = l-1
c
 119	npp=npl2
	k2=npp-1
	do 1 kount=1,l3
	k2=k2+npp
 1      npp=3*npp
c
c
          do 133 kount=1,l
c
	k=np/3
c
	k1=k2-k+1
          k3 = k2+1
          k4 = k2+np
          jstart = np
          if (nperod.eq.3) jstart = 1
              do 132 j=jstart,n,np
              if (j.ne.1) go to 120
              jm1 = j+k
              jm2 = jm1+k
              go to 122
  120         jm1 = j-k
              jm2 = jm1-k
              if (j.ne.n) go to 122
              if (nperod.eq.0) go to 121
              jp1 = jm1
              jp2 = jm2
              go to 123
  121         jp1 = k
              jp2 = jp1+k
              go to 123
  122         jp1 = j+k
              jp2 = jp1+k
c
 123     continue
c
 	do 131 i=1,m
 131    b(i)=3.*q(i,j)
c
              call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 129 i=1,m
	t(i)=b(i)
 129    b(i)=q(i,j)
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 127 i=1,m
	t(i)=t(i)+b(i)
	q(i,jm1)=q(i,jm1)+t(i)
	q(i,jp1)=q(i,jp1)+t(i)
 127    b(i)=t(i)
c
              call atri (k1,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 125 i=1,m
	q(i,j)=q(i,j)+2.*b(i)
	q(i,jm2)=q(i,jm2)+b(i)
 125    q(i,jp2)=q(i,jp2)+b(i)
c
 132    continue
c
	k2=k2-k
	np=np/3
c
 133    continue
c
	np=np/3
c
 134    continue
c
	endif
c
	if(l2.ne.0)then
c
c
      if (iex2.eq.0) go to 118
c
	np=npl2
c
      l = iex2
      if (iex3+iex5.ne.0) go to 108
      if (nperod.eq.1) l = l-1
  108     do 117 kount=1,l
c
	k=np/2
c
          k2 = np-1
          k3 = np
          k4 = 2*np-1
c
          jstart = np
          if (nperod.eq.3) jstart = 1
              do 116 j=jstart,n,np
              jm1 = j-k
              jp1 = j+k
              if (j.eq.1) jm1 = jp1
              if (j.ne.n) go to 109
              jp1 = jm1
              if (nperod.eq.0) jp1 = k
  109             do 115 i=1,m
                  b(i) = 2.*q(i,j)
  115             t(i)=q(i,j)
c
              call atri (k3,k4,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
                  do 113 i=1,m
                  t(i) = t(i)+b(i)
		  q(i,j)=t(i)
                  b(i) = t(i)
  113             continue
c
              call atri (k,k2,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 111 i=1,m
	q(i,jm1)=q(i,jm1)+b(i)
 111    q(i,jp1)=q(i,jp1)+b(i)
c
  116      continue
c
	np=np/2
c
  117     continue
c
  118    continue
c
	endif
c
	do 106 j=1,n
	do 105 i=1,m
 105    b(i)=q(i,j)
c
          call atri (1,1,1,m,mm1,ba,bb,bc,b,twocos,d,w,zz)
c
	do 103 i=1,m
 103    q(i,j)=-b(i)
c
 106    continue
c
	return
	end
c
      subroutine atrid (kstart,kstop,ising,m,mm1,a,b,c,y,twocos,d,w,zz)
      dimension        a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     1                 twocos(1)  ,d(1)       ,w(1), zz(1)
c
c uro is a machine dependent constant.  it is the smallest floating
c point number such that 1 + uro is greater than 1 in the
c arithmetic of this computer.
c
      uro = 1./2.**20
c
      if (ising.ne.1)then
c
      d(1) = c(1)/(b(1)-2.)
          do 114 i=2,mm1
          z = b(i)-2.-a(i)*d(i-1)
          d(i) = c(i)/z
          zz(i)=z
          y(i) = y(i)-d(i-1)*y(i-1)
  114     continue
      y(m)=y(m)-d(mm1)*y(mm1)
      z = b(m)-2.-a(m)*d(m-1)
      if (abs(z).gt.(30.*m*abs(a(m)))*uro)then
      y(mm1)=y(mm1)-a(m)*y(m)/z
      endif
      do 109 i=mm1-1,1,-1
  109 y(i) = y(i)-a(i+1)*y(i+1)/zz(i+1)
c
      if (abs(z).gt.(30.*m*abs(a(m)))*uro)goto 111
      y(m)=0
      goto 112
 111  y(m)=y(m)/z
 112  continue
c
      y(1)=y(1)/(b(1)-2.)
c
      do 120 i=2,mm1
 120  y(i)=y(i)/zz(i)
c
	endif
c
          do 106 k=kstart,kstop
          x = -twocos(k)
          d(1) = c(1)/(b(1)+x)
              do 105 i=2,m
              z = b(i)+x-a(i)*d(i-1)
              d(i) = c(i)/z
              zz(i)=z
  105         y(i) =y(i) -d(i-1)*y(i-1)
              do 103 i=mm1,1,-1
  103         y(i) = y(i)-a(i+1)*y(i+1)/zz(i+1)
c
	y(1)=y(1)/(b(1)+x)
	do 130 i=2,m
 130    y(i)=y(i)/zz(i)
c
  106     continue
	return
	end
c
      subroutine atridp (kstart,kstop,ising,m,mm1,a,b,c,y,twocos,d,w,zz)
      dimension        a(1)       ,b(1)       ,c(1)       ,y(1)       ,
     1                 twocos(1)  ,d(1)       ,w(1),  zz(1)
c
c uro is a machine dependent constant.  it is the smallest floating
c point number such that 1 + uro is greater than 1 in the
c arithmetic of this computer.
c
      uro = 1./2.**20
c
      if (ising.ne.1)then
      d(1) = c(1)/(b(1)-2.)
      w(1) = a(1)/(b(1)-2.)
      bm = b(m)
      z = c(m)
          do 114 i=2,mm1
          den = b(i)-2.-a(i)*d(i-1)
          d(i) = c(i)/den
          w(i) = -a(i)*w(i-1)/den
          y(i) =y(i)-d(i-1)*y(i-1)
          bm = bm-z*w(i-1)
	  zz(i-1)=z
          z = -z*d(i-1)
  114     continue
	dmm1=d(mm1)+w(mm1)
	  y(m)=y(m)-dmm1*y(mm1)
          do 115 i=1,mm1-1
 115      y(m)=y(m)-w(i)*y(i)
          d(mm1) = d(mm1)+w(mm1)
      z = a(m)+z
      den = bm-2.-z*d(mm1)
      if (abs(den).gt.20.*(m**3+m**2*abs(a(m)))*uro) go to 111
      y(m) = 0.
      go to 112
  111 y(m) = y(m)/den
  112 continue
c
      y(mm1)=-z*y(m)+y(mm1)
c
          do 109 i=m-2,m-mm1,-1
          den=b(i+1)-2.-a(i+1)*d(i)
          y(i) = y(i)-a(i+1)*y(i+1)/den-zz(i)*y(m)
  109     continue
c
          y(1)=y(1)/(b(1)-2.)
c
          do 117 i=2,mm1
          den=b(i)-2.-a(i)*d(i-1)
 117      y(i)=y(i)/den
c
	endif
c
          do 106 k=kstart,kstop
          x = -twocos(k)
          d(1) = c(1)/(b(1)+x)
          w(1) = a(1)/(b(1)+x)
          bm = b(m)
          z = c(m)
              do 105 i=2,mm1
              den = b(i)+x-a(i)*d(i-1)
              d(i) = c(i)/den
              w(i) = -a(i)*w(i-1)/den
              y(i) =y(i) -d(i-1)*y(i-1)
              bm = bm-z*w(i-1)
	      zz(i-1)=z
              z = -z*d(i-1)
  105         continue
	dmm1=d(mm1)+w(mm1)
	  y(m)=y(m)-dmm1*y(mm1)
          do 118 i=m-mm1,m-2
  118     y(m)=y(m)-w(i)*y(i)
c
          d(mm1) = d(mm1)+w(mm1)
          z = a(m)+z
          den = bm+x-z*d(mm1)
          y(m) = y(m)/den
          y(mm1) = y(mm1)-z*y(m)
c
              do 103 i=m-2,m-mm1,-1
              den = b(i+1)+x-a(i+1)*d(i)
              y(i) = y(i)-a(i+1)*y(i+1)/den-zz(i)*y(m)
  103         continue
c
	  y(1)=y(1)/(b(1)+x)
	  do 120 i=2,mm1
          den = b(i)+x-a(i)*d(i-1)
  120     y(i)=y(i)/den
  106     continue
c
	return
	end
c
c
	subroutine delm1(ain,aout)
c
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
	real a(mn0,k0),b(mn0,k0)
	real c(mn0,k0),d(mn0,k0)
	real tempa(mn0,k0),tempb(mn0,k0)
        real ain(mn0,k0),aout(mn0,k0)
c
c......copy ain into a
c
	do 600 k=1,kz
	do 600 j=1,n
	do 600 i=1,m
	ip=(j-1)*m+i
	a(ip,k)=ain(ip,k)
 600    continue
c
c......copy a into b........
c
       do k=1,kz
       do j=1,n
       do i=1,m
       ip=(j-1)*m+i
       b(ip,k)=a(ip,k)
       enddo
       enddo
       enddo
c
c......set the b.c.s on b....
c
      do k=1,kz
       do ip=1,m
        b(ip,k)=0.0
       enddo
       do ip=msq-m+1,msq
        b(ip,k)=0.0
       enddo
      enddo
c
c......project onto the vertical normal modes......
c
       do k=1,kz
        do ip=1,msq
         c(ip,k)=0.0
         do j=1,kz
         c(ip,k)=c(ip,k)+ainv(k,j)*b(ip,j)
         enddo
        enddo
       enddo
c
c.....invert the Helmholtz eqn.....
c
       do k=1,kz
c
        mb=0
c
      dbasin=hy*nm1
      call pwscrt(0,0.,xbasin,mm1,mb,bda,bdb,0.,dbasin,nm1,1,bdc,bdd
     +,eigval(k),c(1,k),m,1.,ierror,r1)
c
       enddo
c
c......project onto vertical modes.....
c
       do k=1,kz
        do ip=1,msq
          d(ip,k)=0.0
          do j=1,kz
          d(ip,k)=d(ip,k)+amat(k,j)*c(ip,j)
          enddo
        enddo
       enddo
c
c......set up the output array....
c
	do 602 k=1,kz
	do 602 j=1,n
	do 602 i=1,m
	ipp=ipp+1
	ip=(j-1)*m+i
 602    aout(ip,k)=d(ip,k)
c
c
	return
	end
c
	subroutine adelm1(ain,aout)
c
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
	real a(mn0,k0),b(mn0,k0)
	real c(mn0,k0),d(mn0,k0)
	real tempa(mn0,k0),tempb(mn0,k0)
	real ain(mn0,k0),aout(mn0,k0)
c
c......set up the input array....
c
	do 602 k=1,kz
	do 602 j=1,n
	do 602 i=1,m
	ip=(j-1)*m+i
 602    d(ip,k)=aout(ip,k)
c
c......clear c array.......
c
        do k=1,kz
        do ip=1,msq
        c(ip,k)=0.0
        enddo
        enddo
c
c......project onto vertical modes.....
c
       do k=1,kz
        do ip=1,msq
          do j=1,kz
          c(ip,j)=c(ip,j)+amat(k,j)*d(ip,k)
          enddo
        enddo
       enddo
c
c.....invert the Helmholtz eqn.....
c
       do k=1,kz
c
c#######: Note!!! tests show that pwscrt if self adjoint
c         so pwscrt if used instead of apwscrt which does
c         not work properly at the moment. In this case we
c         need to set c to zero at the boundaries.
c
c......set the b.c.s on c....
c
       do ip=1,m
        c(ip,k)=0.0
       enddo
       do ip=msq-m+1,msq
        c(ip,k)=0.0
       enddo
c
c
	mb=0
c
      dbasin=hy*nm1
c
      call pwscrt(0,0.,xbasin,mm1,mb,bda,bdb,0.,
     +  dbasin,nm1,1,bdc,bdd
     +  ,eigval(k),c(1,k),m,1.,ierror,r1)
c
       enddo
c
c.....clear the b array....
c
       do k=1,kz
       do ip=1,msq
       b(ip,k)=0.0
       enddo
       enddo
c
c......project onto the vertical normal modes......
c
       do k=1,kz
        do ip=1,msq
         do j=1,kz
         b(ip,j)=b(ip,j)+ainv(k,j)*c(ip,k)
         enddo
        enddo
       enddo
c
c.....clear the a array....
c
       do k=1,kz
       do ip=1,msq
       a(ip,k)=0.0
       enddo
       enddo
c
c......copy b into a........
c
       do k=1,kz
       do j=1,n
       do i=1,m
       ip=(j-1)*m+i
       a(ip,k)=b(ip,k)
       enddo
       enddo
       enddo
c
c......copy a into ain
c
	do 600 k=1,kz
	do 600 j=1,n
	do 600 i=1,m
	ip=(j-1)*m+i
	ain(ip,k)=a(ip,k)
 600    continue
c
	return
	end
c
      subroutine rel2nd(psi,dx2psi,dy2psi,m,n,hy)
c
c --- this routine calculates the fourth order laplacian of psi
c --- at the boundaries and in the interior field
c
c
       dimension psi(1),dx2psi(1),dy2psi(1)
c
       fac=1./(hy**2)
c
       do j=1,n
       do i=2,m-1
       ip=(j-1)*m+i
       dx2psi(ip)=fac*(psi(ip+1)-2.*psi(ip)+psi(ip-1))
       enddo
       enddo
c
       do j=1,n
       ip=(j-1)*m+1
       ip1=(j-1)*m+m
       dx2psi(ip)=fac*(psi(ip+1)-psi(ip))
       dx2psi(ip1)=fac*(psi(ip1-1)-psi(ip1))
       enddo
c
       do j=2,n-1
       do i=1,m
       ip=(j-1)*m+i
       dy2psi(ip)=fac*(psi(ip+m)-2.*psi(ip)+psi(ip-m))
       enddo
       enddo
c
       do i=1,m
       ip=i
       ip1=(n-1)*m+i
       dy2psi(ip)=fac*(psi(ip+m)-psi(ip))
       dy2psi(ip1)=fac*(psi(ip1-m)-psi(ip1))
       enddo
c
      return
      end
c
      subroutine adrel2nd(psi,dx2,dy2,m,n,hy)
c
c --- this routine calculates the fourth order laplacian of psi
c --- at the boundaries and in the interior field
c
c
       dimension psi(1),dx2(1),dy2(1)
c
       fac=1./(hy**2)
c
       do j=1,n
       do i=1,m
       ip=(j-1)*m+i
       psi(ip)=0.0
       enddo
       enddo
c
       do j=1,n
       do i=2,m-1
       ip=(j-1)*m+i
c
       psi(ip+1)=psi(ip+1)+dx2(ip)
       psi(ip)=psi(ip)-2.*dx2(ip)
       psi(ip-1)=psi(ip-1)+dx2(ip)
c
       enddo
       enddo
c
       do j=1,n
       ip=(j-1)*m+1
       ip1=(j-1)*m+m
c
       psi(ip+1)=psi(ip+1)+dx2(ip)
       psi(ip)=psi(ip)-dx2(ip)
       psi(ip1-1)=psi(ip1-1)+dx2(ip1)
       psi(ip1)=psi(ip1)-dx2(ip1)
c
       enddo
c
       do j=2,n-1
       do i=1,m
       ip=(j-1)*m+i
c
       psi(ip+m)=psi(ip+m)+dy2(ip)
       psi(ip)=psi(ip)-2.*dy2(ip)
       psi(ip-m)=psi(ip-m)+dy2(ip)
c
       enddo
       enddo
c
       do i=1,m
       ip=i
       ip1=(n-1)*m+i
c
       psi(ip+m)=psi(ip+m)+dy2(ip)
       psi(ip)=psi(ip)-dy2(ip)
       psi(ip1-m)=psi(ip1-m)+dy2(ip1)
       psi(ip1)=psi(ip1)-dy2(ip1)
c
       enddo
c
       do j=1,n
       do i=1,m
       ip=(j-1)*m+i
       psi(ip)=fac*psi(ip)
       enddo
       enddo
c
      return
      end
c
	subroutine calvt(a,b)
c
c........This is a variant of delsq which computes the 3-D
c        Laplacian by projecting onto the vertical modes
c        then forming the Helmholtz equations.
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
	real a(mn0,k0),b(mn0,k0)
	real tempa(mn0,k0),tempb(mn0,k0)
        real dx2(m20),dy2(m20)
c
c......copy a into tempa
c
	ipp=0
	do 600 k=1,kz
	do 600 j=1,n
	do 600 i=1,m
	ipp=ipp+1
	ip=(j-1)*m+i
	tempa(ip,k)=a(ip,k)
 600    continue
c
c......compute the 3-D Laplacian by projecting onto the
c      vertical normal modes.
c
	do 252 k=1,kz
c
        do 250 ip=1,msq
          r(ip)=0.0
          do 210 k1=1,kz
            r(ip)=r(ip)+ainv(k,k1)*tempa(ip,k1)
  210     continue
  250   continue
c
          call rel2nd(r,dx2,dy2,m,n,hy)
c
          do ip=1,msq
          r(ip)=dx2(ip)+dy2(ip)+eigval(k)*r(ip)
          enddo
c
c.......use tempb to hold r......
c
	do 601 ip=1,msq
 601    tempb(ip,k)=r(ip)
c
 252    continue
c
c
c --- new b--point by point
c
      do 350 ip=1,msq
        do 340 k=1,kz
           r(k)=tempb(ip,k)
           tempb(ip,k)=0.
  340   continue
c
        do 345 k=1,kz
        do 345 j=1,kz
           tempb(ip,k)=tempb(ip,k)+amat(k,j)*r(j)
  345   continue
  350 continue
c
c........copy tempb into b......
c
	ipp=0
	do 602 k=1,kz
	do 602 j=1,n
	do 602 i=1,m
	ipp=ipp+1
	ip=(j-1)*m+i
 602    b(ip,k)=tempb(ip,k)
c
c
	return
	end
c
	subroutine acalvt(a,b)
c
c........This is a variant of delsq which computes the 3-D
c        Laplacian by projecting onto the vertical modes
c        then forming the Helmholtz equations.
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
	real a(mn0,k0),b(mn0,k0)
	real tempa(mn0,k0),tempb(mn0,k0)
        real dx2(m20),dy2(m20)
c
c
c........copy tempb into b......
c
	ipp=0
	do 602 k=1,kz
	do 602 j=1,n
	do 602 i=1,m
	ipp=ipp+1
	ip=(j-1)*m+i
        tempa(ip,k)=0.0
 602    tempb(ip,k)=b(ip,k)
c
      do 350 ip=1,msq
c
        do j=1,msq
        r(j)=0.0
        enddo
c
        do 345 k=1,kz
        do 345 j=1,kz
           r(j)=r(j)+amat(k,j)*tempb(ip,k)
  345   continue
c
        do 340 k=1,kz
           tempb(ip,k)=r(k)
  340   continue
c
  350 continue
c
c
c......compute the 3-D Laplacian by projecting onto the
c      vertical normal modes.
c
	do 252 k=1,kz
c
c.......use tempb to hold r......
c
	do 601 ip=1,msq
 601    r(ip)=tempb(ip,k)
c
          do ip=1,msq
          dx2(ip)=r(ip)
          dy2(ip)=r(ip)
          r(ip)=eigval(k)*r(ip)
          enddo
c
          call adrel2nd(r1,dx2,dy2,m,n,hy)
c
          do ip=1,msq
          r(ip)=r(ip)+r1(ip)
          enddo
c
        do 250 ip=1,msq
          do 210 k1=1,kz
            tempa(ip,k1)=tempa(ip,k1)+ainv(k,k1)*r(ip)
  210     continue
  250   continue
c
c
 252    continue
c
c......copy a into tempa
c
	ipp=0
	do 600 k=1,kz
	do 600 j=1,n
	do 600 i=1,m
	ipp=ipp+1
	ip=(j-1)*m+i
	a(ip,k)=tempa(ip,k)
 600    continue
c
	return
	end
	subroutine delsq1(ain,aout)
c
c........This is a variant of delsq which computes the 3-D
c        Laplacian by projecting onto the vertical modes
c        then forming the Helmholtz equations.
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
	real tempa(mn0,k0),tempb(mn0,k0)
        real dx2(m20),dy2(m20)
        real ain(mn0,k0),aout(mn0,k0)
c
c......copy ain into tempa
c
c
	do 600 k=1,kz
	do 600 j=1,n
	do 600 i=1,m
	ip=(j-1)*m+i
	tempa(ip,k)=ain(ip,k)
 600    continue
c
c......compute the 3-D Laplacian by projecting onto the
c      vertical normal modes.
c
	do 252 k=1,kz
c
        do 250 ip=1,msq
          r(ip)=0.0
          do 210 k1=1,kz
            r(ip)=r(ip)+ainv(k,k1)*tempa(ip,k1)
  210     continue
  250   continue
c
          call rel2nd(r,dx2,dy2,m,n,hy)
c
          do ip=1,msq
          r(ip)=dx2(ip)+dy2(ip)+eigval(k)*r(ip)
          enddo
c
c.......use tempb to hold r......
c
	do 601 ip=1,msq
 601    tempb(ip,k)=r(ip)
c
 252    continue
c
c
c --- new b--point by point
c
      do 350 ip=1,msq
        do 340 k=1,kz
           r(k)=tempb(ip,k)
           tempb(ip,k)=0.
  340   continue
c
        do 345 k=1,kz
        do 345 j=1,kz
           tempb(ip,k)=tempb(ip,k)+amat(k,j)*r(j)
  345   continue
  350 continue
c
c........copy tempb into aout......
c
	do 602 k=1,kz
	do 602 j=1,n
	do 602 i=1,m
	ip=(j-1)*m+i
 602    aout(ip,k)=tempb(ip,k)
c
c
	return
	end
	subroutine adelsq1(ain,aout)
c
c........This is a variant of delsq which computes the 3-D
c        Laplacian by projecting onto the vertical modes
c        then forming the Helmholtz equations.
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
	real tempa(mn0,k0),tempb(mn0,k0)
        real dx2(m20),dy2(m20)
        real ain(mn0,k0),aout(mn0,k0)
c
c
c........copy aout into tempb......
c
	do 602 k=1,kz
	do 602 j=1,n
	do 602 i=1,m
	ip=(j-1)*m+i
        tempa(ip,k)=0.0
 602    tempb(ip,k)=aout(ip,k)
c
      do 350 ip=1,msq
c
        do j=1,msq
        r(j)=0.0
        enddo
c
        do 345 k=1,kz
        do 345 j=1,kz
           r(j)=r(j)+amat(k,j)*tempb(ip,k)
  345   continue
c
        do 340 k=1,kz
           tempb(ip,k)=r(k)
  340   continue
c
  350 continue
c
c
c......compute the 3-D Laplacian by projecting onto the
c      vertical normal modes.
c
	do 252 k=1,kz
c
c.......use tempb to hold r......
c
	do 601 ip=1,msq
 601    r(ip)=tempb(ip,k)
c
          do ip=1,msq
          dx2(ip)=r(ip)
          dy2(ip)=r(ip)
          r(ip)=eigval(k)*r(ip)
          enddo
c
          call adrel2nd(r1,dx2,dy2,m,n,hy)
c
          do ip=1,msq
          r(ip)=r(ip)+r1(ip)
          enddo
c
        do 250 ip=1,msq
          do 210 k1=1,kz
            tempa(ip,k1)=tempa(ip,k1)+ainv(k,k1)*r(ip)
  210     continue
  250   continue
c
c
 252    continue
c
c......copy tempa into ain
c
	do 600 k=1,kz
	do 600 j=1,n
	do 600 i=1,m
	ip=(j-1)*m+i
	ain(ip,k)=tempa(ip,k)
 600    continue
c
	return
	end
c
c
      subroutine nvortq(klevel)
c
c --- subroutine nvortq performs the calulation of the vorticity field
c
c --- in the finite-element approximation, the order of the calculation
c --- of the vorticity is the following: interior points, regular
c --- boundary points, and corner points. afterwards, filter the field.
c
c --- array vort - vorticity
c --- arrays r,r1 - scratch arrays
c --- variable bet represents the beta effect
c
c --- fq is an arakawa jacobian which represents a finite-difference
c --- approximation. fq is anti-symmetric ith its arguments and
c --- conserves vorticity, energy, and enstrophy when integrated
c --- over a closed domain.
c
          include 'prob2.pert'
c
c>>>>>>>>>>>>
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c
c>>>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
	common /cnn/n,nm1,nm2,np1
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +  topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +  prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/philtr/iltord,iltfrq,iltcnt
      common/count/icnt,itcnt,ibcnt 
      common/dynam/alpha,beta,rf 
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +       amat(k0,k0)
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/pert/stf(nosp),vtf(nosp)
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/steps/sumt,sumb,cost
c
        common/amiv/amx(m0-1,m0-1),amy(n0,n0)
c
	real tempr(mn0),r3(mn0)
c
c
c --- computation of rhs of finite-element vorticity eq.
c
c>>>>>>>>>>>>>>>>>>>>>...
c
c --- first the arakawa jacobians
c
      call arajac(bstm(1,klevel),vort(1,klevel),r,m,n)
c
	call arajac(stream(1,klevel),bvrt(1,klevel),r1,m,n)
c
c --- normalize the jacobian by alpha
c
      do 20 ip=1,msq
        r(ip)=alpha*(r(ip)+r1(ip))
   20 continue
c
c>>>>>>>>>>>>>>
c
c
c --- include beta effect
c
      if(beta.eq.0)goto 40
      call arajac(stream(1,klevel),pvort,r2,m,n)
      do 30 ip=1,msq
        r(ip)=r(ip)+r2(ip)
   30 continue
   40 continue
c
c --- apply bottom friction. if no bottom friction wanted or if
c --- not on bottom level, skip this section.
c
      if(klevel.ne.kz)go to 65
      if(rf.eq.0.0)go to 65
c
c.....northern and southern boundaries only...
c
      cf=(rf/36.0)*(hy**2)
c
	do ip=1,msq
          r3(ip)=0.0
	enddo
c
c.....southern boundary......
c
      do ip=1,m-1
          il=ip+1
          ir=ip-1
	  iu=ip+m
          iul=iu-1
          iur=iu+1
          if(ip.eq.1)ir=m
          if(ip.eq.1)iul=2*m
          r3(ip)=cf*(vort(iul,klevel) + vort(iur,klevel) +
     +          2*(vort(il,klevel) + vort(ir,klevel)) +
     +          4*vort(iu,klevel) + 8*vort(ip,klevel))
       enddo
c
c.....northern boundary......
c
      do ip=(n-1)*m+1,n*m-1
          il=ip+1
          ir=ip-1
          iu=ip-m
          iul=iu-1
          iur=iu+1
          if(ip.eq.(n-1)*m+1)ir=m*n
          if(ip.eq.(n-1)*m+1)iul=m*(n-1)
          r3(ip)=cf*(vort(iul,klevel) + vort(iur,klevel) +
     +          2*(vort(il,klevel) + vort(ir,klevel)) +
     +          4*vort(iu,klevel) + 8*vort(ip,klevel))
       enddo

c
c --- now interior points.
c
      do 50 irow=2,n-1
        ip1=(irow-1)*m+1
        ip2=(irow*m)-1
        do 49 ip=ip1,ip2
          ipm=ip-1
          if(ip.eq.ip1)ipm=ip2+1
          ipmp=ipm+m
          ipmm=ipm-m
          r3(ip)=cf*(vort(ipmm,klevel) + vort(ip-m+1,klevel) +
     +          vort(ipmp,klevel) + vort(ip+m+1,klevel) + 
     +          4*(vort(ip-m,klevel) + vort(ip+1,klevel) +
     +          vort(ip+m,klevel) + vort(ipm,klevel)) + 
     +          16*vort(ip,klevel))
   49   continue
   50 continue
c
c.....impose periodicity.
c
        do j=1,n
        ip1=(j-1)*m+m
        ip2=(j-1)*m+1
        r3(ip1)=r3(ip2)
        enddo
c
        do ip=1,msq
          r(ip)=r(ip)+r3(ip)
        enddo
c
   65 continue
c
c
c --- now convert this term to the proper form for msolve
c --- and genrl to solve the adams-bashforth problem.
c
      if(icalc.gt.1)call rhs(r,prvjac(1,klevel))
c
      if(icalc.eq.1)call rhs2(r,prvjac(1,klevel))
c
c --- if this is the initialization step, return.
c
      if(icalc.eq.0) go to 700
c
c
c>>>>>>>>>> multiply r by the inverse matrices amx and amy.......
c
        do 2020 j=1,n
        do 2020 i=1,m-1
        iq=(j-1)*m+i
        sum=0.0
        do 2021 il=1,m-1
        ip=(j-1)*m+il
 2021   sum=sum+amx(i,il)*r(ip)
 2020   tempr(iq)=sum
c
        do 2022 j=1,n
        do 2022 i=1,m-1
        iq=(j-1)*m+i
        sum=0.0
        do 2023 il=1,n
        ip=(il-1)*m+i
 2023   sum=sum+tempr(ip)*amy(il,j)
 2022   r(iq)=sum
c
c>>>>>>>>>  UPDATE LAST COLUMN OF R >>>>>>>>>>>>>
c
        do 4452 j=1,n
        ip1=(j-1)*m+m
        ip2=(j-1)*m+1
 4452   r(ip1)=r(ip2)
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c --- update the vorticity field
c
      do 70 ip=1,msq
        vort(ip,klevel)=vort(ip,klevel)+r(ip)
   70 continue
c
c
c
c --- definitions of the filter variables:
c
c --- icnt = counter (for frequency of application)
c --- iltcnt = frequency with which filter is applied
c --- iltfrq = no of times filter is applied per time step
c --- iltord = order of the filter
c
c      icnt=icnt+1
c      if(icnt.ne.iltcnt)go to 120
c      if(iltord.eq.0) go to 120
c      icnt=0
      do 110 i=1,iltfrq
        call filter(vort(1,klevel),m,n,iltord)
  110 continue
c  120 continue
c
c --- come here (skipping the update of vorticity) if icalc .eq. 0
c
  700 continue
c
c
c --- if iftop=0 then do not want top density
c
      if(iftop.eq.0)goto 8500
c
c --- section for prognostic computation of surface
c --- temperature anomaly.
c
c --- this calculation is done with the top layer only.
c
      if(klevel.ne.1)goto 8500
c
c --- compute jacobian of top density
c
      call arajac(bstm(1,1),topden,r,m,n)
c
      call arajac(stream(1,1),btopd,r1,m,n)
c
      do ip=1,msq
         r(ip)=r(ip)+r1(ip)
      enddo
c
      do 710 ip=1,msq
        r(ip)=alpha*r(ip)
  710 continue
c
c
      if(icalc.gt.1)call rhs(r,prvtop)
c
      if(icalc.eq.1)call rhs2(r,prvtop)
c
c --- if initialization, then finished with top density
c
      if(icalc.eq.0)goto 8500
c
c --- now convert this term to the proper form for msolve
c --- and genrl to solve the adams-bashforth problem.
c
      call msolve
c
c      call infld(r1,1,tim,7,1)
c
	do ip=1,msq
           r1(ip)=0.0
	enddo
c
      call ngenrl(bstm(1,1),topden,1)
c
c --- now advance the top density
c
      do 820 ip=1,msq
        topden(ip)=topden(ip)+r(ip)
  820 continue
c
c      itcnt=itcnt+1
c      if(itcnt.ne.iltcnt)goto 851
c      if(iltord.eq.0)goto 851
c      itcnt=0
      do 830 i=1,iltfrq
        call filter(topden,m,n,iltord)
  830 continue
c  851 continue
c
c --- finished with top density calculations
c
 8500 continue
c
c
      return
      end
c
      subroutine avortq(klevel)
c
c --- subroutine avortq performs the calulation of the vorticity field
c
c --- in the finite-element approximation, the order of the calculation
c --- of the vorticity is the following: interior points, regular
c --- boundary points, and corner points. afterwards, filter the field.
c
c --- array vort - vorticity
c --- arrays r,r1 - scratch arrays
c --- variable bet represents the beta effect
c
c --- fq is an arakawa jacobian which represents a finite-difference
c --- approximation. fq is anti-symmetric ith its arguments and
c --- conserves vorticity, energy, and enstrophy when integrated
c --- over a closed domain.
c
          include 'prob2.pert'
c
c>>>>>>>>>>>>
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c
c>>>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
	common /cnn/n,nm1,nm2,np1
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +  topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +  prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +  xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/philtr/iltord,iltfrq,iltcnt
      common/count/icnt,itcnt,ibcnt 
      common/dynam/alpha,beta,rf 
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +       amat(k0,k0)
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/pert/stf(nosp),vtf(nosp)
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/steps/sumt,sumb,cost
c
        common/amiv/amx(m0-1,m0-1),amy(n0,n0)
c
	real tempr(mn0)
c
c
c......clear stream array.....
c
      do i=1,msq
         stream(i,klevel)=0.0
      enddo
c
c --- if iftop=0 then do not want top density
c
      if(iftop.eq.0)goto 8500
c
c --- section for prognostic computation of surface
c --- temperature anomaly.
c
c --- this calculation is done with the top layer only.
c
      if(klevel.ne.1)goto 8500
c
c --- if initialization, then finished with top density
c
      if(icalc.eq.0)goto 850
c
      do 830 i=1,iltfrq
        call afiltr(topden,m,n,iltord)
  830 continue
c
      do 820 ip=1,msq
        r(ip)=topden(ip)
  820 continue
c
      call nvagenrl(bstm(1,1),topden,1)
c
      call amsolve
c
 850  continue
c
      if(icalc.gt.1)call arhs(r,prvtop)
c
      if(icalc.eq.1)call arhs2(r,prvtop)
c
      do 710 ip=1,msq
        r(ip)=alpha*r(ip)
  710 continue
c
      do ip=1,msq
         r1(ip)=r(ip)
      enddo
c
      call adjacv(bstm(1,1),topden,r,m,n)
c
      call adjacs(stream(1,1),btopd,r1,m,n)
c
c
c --- finished with top density calculations
c
 8500 continue
c
c
c --- if this is the initialization step, return.
c
      if(icalc.eq.0) go to 700
c
c
c --- definitions of the filter variables:
c
c --- icnt = counter (for frequency of application)
c --- iltcnt = frequency with which filter is applied
c --- iltfrq = no of times filter is applied per time step
c --- iltord = order of the filter
c
c      icnt=icnt+1
c      if(icnt.ne.iltcnt)go to 120
c      if(iltord.eq.0) go to 120
c      icnt=0
      do 110 i=1,iltfrq
        call afiltr(vort(1,klevel),m,n,iltord)
  110 continue
c  120 continue
c
c
      do 70 ip=1,msq
        r(ip)=vort(ip,klevel)
   70 continue
c
c
c>>>>>>>>>> multiply r by the inverse matrices amx and amy.......
c
	do ip=1,mn0
	 tempr(ip)=0.0
	enddo
c
        do 2022 j=1,n
        do 2022 i=1,m-1
        iq=(j-1)*m+i
        sum=r(iq)
	r(iq)=0.0
        do 2023 il=1,n
        ip=(il-1)*m+i
 2023   tempr(ip)=tempr(ip)+amy(il,j)*sum
 2022   continue
c
        do 2020 j=1,n
        do 2020 i=1,m-1
        iq=(j-1)*m+i
        sum=tempr(iq)
	tempr(iq)=0.0
        do 2021 il=1,m-1
        ip=(j-1)*m+il
 2021   r(ip)=r(ip)+amx(i,il)*sum
 2020   continue
c
c
c>>>>>>>>>  UPDATE LAST COLUMN OF R >>>>>>>>>>>>>
c
        do 4452 j=1,n
        ip1=(j-1)*m+m
        ip2=(j-1)*m+1
 4452   r(ip1)=r(ip2)
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c
  700 continue
c
c
c --- now convert this term to the proper form for msolve
c --- and genrl to solve the adams-bashforth problem.
c
      if(icalc.gt.1)call arhs(r,prvjac(1,klevel))
c
      if(icalc.eq.1)call arhs2(r,prvjac(1,klevel))
c
c
c --- include beta effect
c
      if(beta.eq.0)goto 40
c
      do 30 ip=1,msq
        r2(ip)=r(ip)
   30 continue
c
      call adjacs(stream(1,klevel),pvort,r2,m,n)
c
   40 continue
c
c
c --- normalize the jacobian by alpha
c
      do 20 ip=1,msq
        r1(ip)=alpha*r(ip)
        r(ip)=alpha*r(ip)
   20 continue
c
c
c --- first the arakawa jacobians
c
      call adjacv(bstm(1,klevel),vort(1,klevel),r,m,n)
c
      call adjacs(stream(1,klevel),bvrt(1,klevel),r1,m,n)
c
c
      return
      end
c
c
	subroutine l2norm
c
        include 'prob2.pert'
c
	common/egy/peint,akeint
	common/l2n/al2s,al2v
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
c
	real sx(n0,k0),vx(n0,k0),sz(k0),vz(k0)
c
c.......integrate norms over whole domain.....
c
	sums=0.0
	sumv=0.0
	do k=1,kz
	do j=1,n
	do i=1,m
	ip=(j-1)*m+i
	sums=sums+stream(ip,k)**2*hz(k)
        sumv=sumv+vort(ip,k)**2*hz(k)
        enddo
        enddo
        enddo
c
	al2s=sums
	al2v=sumv
c
	return
	end
c
	subroutine energy
c
        include 'prob2.pert'
c
	common/egy/peint,akeint
	common/l2n/al2s,al2v
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
c
	total=0.0
	do 1000 k=1,kz
	do 1000 j=1,n
	do 1000 i=1,m-1
	ip=(j-1)*m+i
c1000   total=total+stream(ip,k)*vort(ip,k)*hz(k)
 1000   total=total+stream(ip,k)*vort(ip,k)*hz(k)
c
	akeint=-0.5*total
        peint=0.0
c
	return
	end
c
c
      subroutine set3d(kz,ifdiff)
c
c --- subroutine set3d sets the stratification and the
c --- associated matrices and constants
c
c
          include 'prob2.pert'
c
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/bsolv/ab(mx0),bb(mx0),cb(mx0),workb(mnwork)
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
c
	common/adm/tiam(k0,k0),taml(k0,k0)
c
        common/comm1/vmat(k0,k0),vinv(k0,k0)
c
      dimension hz1(k0),strat(k0,k0)
      dimension e1(k0),e2(k0,k0),iwork(2*k0),e3(k0,k0)
c
c --- ifdiff=1/0 - finite diff/collocation
c
      if(ifdiff.eq.0)go to 50
c
c --- calculate the related depths hz1
c
      ierr=0
      do 10 k=2,kz
        hz1(k)=0.5*(hz(k-1)+hz(k))
   10 continue
c
c --- step 1 - set the stratification matrix strat
c
      do 15 i=1,kz
      do 15 j=1,kz
        strat(i,j)=0.
   15 continue
      kzm1=kz-1
      do 20 i=1,kzm1
        a1=scsig*sigz(i+1)/(hz(i)*hz1(i+1))
        strat(i,i)=-a1
        strat(i,i+1)=a1
   20 continue
      do 25 i=2,kz
        a1=scsig*sigz(i)/(hz(i)*hz1(i))
        strat(i,i-1)=a1
        strat(i,i)=strat(i,i)-a1
   25 continue
c
c --- step 2 - get eigenvalues (eigval) and eigenvectors (amat)
c
      call eigrg1(k0,kz,strat,kz,eigval,e1,amat,e2,workb,ierr)
c
c --- check the results
c
      if(ierr.ne.0) go to 99
      do 30 i=1,kz
        if(e1(i).ne.0) go to 99
   30 continue
      xkz=float(kz)
      xkz=sqrt(xkz)
c
c --- step 3 - normalize the eigenvectors
c
      do 45 i=1,kz
        a1=0.
        do 35 i1=1,kz
          a1=a1+amat(i1,i)**2
   35   continue
        do 40 j=1,kz
          amat(j,i)=xkz*amat(j,i)/sqrt(a1)
   40   continue
   45 continue
c
c --- come directly here if using the collocation method
c
   50 continue
c
      call sort(eigval,amat,kz)
c
c......copy amat into vmat and normalise vmat so that
c      transpose(vmat) D vmat = I where D is the matrix of
c      level thickness.
c
       do j=1,kz
       do i=1,kz
       vmat(i,j)=amat(i,j)
       enddo
       enddo
c
       do j=1,kz
          sum=0.0
          do i=1,kz
          sum=sum+hz(i)*vmat(i,j)*vmat(i,j)
          enddo
          do i=1,kz
          vmat(i,j)=vmat(i,j)/sqrt(sum)
          enddo
       enddo
c
c --- step 4 - obtain inverse of amat(ainv).  since amat will be
c --- used later, a dummy array will be used in the inversion
c --- to prevent amat from being destroyed
c
      do 55 i=1,kz
      do 55 j=1,kz
        e2(i,j)=amat(i,j)
   55 continue
c
      call invmtx(e2,k0,ainv,k0,kz,det,iwork,ierr)
c
c.....calc. inverse of vmat......
c
	do 8000 i=1,kz
	do 8000 j=1,kz
 8000   e3(i,j)=vmat(i,j)
c
	call invmtx(e3,k0,vinv,k0,kz,det1,iwork,ierr)
c
c......check vmat and vinv.....
c
      print *,' trans(vmat) D vmat'
c
      do j=1,kz
         do i=1,kz
            sum=0.0
            do k=1,kz
            sum=sum+vmat(k,i)*hz(k)*vmat(k,j)
            enddo
            print *,'i=',i,' j=',j,' sum=',sum
         enddo
      enddo
c
      print *,'vinv*vmat'
c
      do j=1,kz
         do i=1,kz
            sum=0.0
            do k=1,kz
            sum=sum+vinv(i,k)*vmat(k,j)
            enddo
            print *,'i=',i,' j=',j,' sum=',sum
         enddo
      enddo
c
      if(ierr.ne.0) go to 990
      if(ifdiff.ne.0) then 
           write(4,*)' finite difference method   eigenvalues:'
           write(4,60)(eigval(i),i=1,kz)
      endif
      if(ifdiff.eq.0) then 
           write(4,*)'collocation in depth.  eigenvalues:'
           write(4,60)(eigval(i),i=1,kz)
      endif
    
   60 format(5(1x,1pe11.4,/))


      do 70 i=1,kz
        write(4,75) i
        write(4,60)(amat(j,i),j=1,kz)
   70 continue
   75 format(/'eigenvector(',i1,') :')
      return
c
c --- come here if there is a problem getting the eigenvalues
c --- or eigenvectors
c
   99 continue
      write(4,100)
  100 format(/'something wrong in eigrg1--program halted')
      write(4,101)ierr
  101 format(' ierr = ',i4)
      do 105 iii=1,kz
        write(4,110)(e2(iii,jjj),jjj=1,kz)
  105 continue
  110 format(4(2x,d13.6))
      write(4,115)
  115 format(///' e1')
      write(4,110)(e1(iii),iii=1,kz)
      stop
c
c --- come here if there is a problem in the matrix inverter
c
  990 continue
      write(4,991)
  991 format(/'something wrong in invmtx--program halted')
      stop
      end
c
      subroutine trdslv (n,a,b,c,y,x,ks,alpha,gamma)
      dimension       a(n)       ,b(n)       ,c(n)       ,y(n)       ,
     1                x(n)       ,alpha(1)   ,gamma(1)
      nm1 = n-1
      if (ks .eq. 1) go to 102
      alpha(1) = 1./b(1)
      gamma(1) = c(1)*alpha(1)
      do 101 i=2,nm1
         alpha(i) = 1./(b(i)-a(i)*gamma(i-1))
         gamma(i) = c(i)*alpha(i)
  101 continue
  102 continue
      x(1) = y(1)*alpha(1)
      do 103 i=2,nm1
         x(i) = (y(i)-a(i)*x(i-1))*alpha(i)
  103 continue
      x(n) = (y(n)-a(n)*x(nm1))/(b(n)-a(n)*gamma(nm1))
      do 104 j=1,nm1
         jb = n-j
         x(jb) = x(jb)-gamma(jb)*x(jb+1)
  104 continue
      return
	end
c
c
	subroutine matveca(npnt,right,alhs)
c
c
	include 'prob2.pert'
c
c
      character*4 restid,runid
      character*70 idstrg
c
c
        common/newbg/amms(mn0,k0),ammv(mn0,k0)
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/grate/gre(mn0),grip(mn0),tezm(mn0),aipzm(mn0)
c
c@@@@@@@@@THE FOLLOWING COMMON BLOCK IS ONLY IN CLINIC!!!!
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/egl2/pl,akl,alls,allv
c
	common/mmde/jmode
c
	common/nrest/irst
c
	common/adrst/ityrst
c
	common/cnorm/jnorm
c
	common/jctrl/ktim,icstart,iczero
c
	common/gcord/ir1,ir2,jr1,jr2
c
        common/nmmt/nmit,nmiter,idone
c
        common/comm1/vmat(k0,k0),vinv(k0,k0)
c
c>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/count/icnt,itcnt,ibcnt
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),topo(mn0),
     2     topd(mn0),botd(mn0)
c
	common/atop/radt(mn0),s2t(mn0),st(mn0),fyt(mx0,4)
c
	common/abot/radb(mn0),rad1b(mn0),rad2b(mn0),s2b(mn0),sb(mn0)
     1       ,fyb(mx0,4)
c
c
c
	common/adm/tiam(k0,k0),taml(k0,k0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
c
	common/grad/gomeg(mn0,k0),gphi(mn0,k0),ggam(mn0),
     1     gdel(mn0),dm(mn0,k0),dp(mn0,k0),dg(mn0),dd(mn0)
     2    ,cgf
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/steps/sumt,sumb,cost
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/ctest/ctol
c
	common/swoi/iter
c
	common/toff/tadj0
c
	common/istar/icall
c
        common/ccnt/callcnt
c
        common/iwfq/ibgw,iclw
c
	common/bdata/scrp(n0,k0,2),vcrp(n0,k0,2)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/lan1/maxt,ifts
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c
C************
	 real alhs(npnt),right(npnt),usea(npnt),usea1(npnt)
         real ammt(mn0,k0)
C************
c
c
c
C********INITIALISE stream and vort**************
c
        callcnt=callcnt+1
c
        print *,'callcnt=',callcnt
c
	icalc=0
c
c.......copy right into vort......
c
        do k=1,kz
        do ip=1,msq
        stream(ip,k)=0.0
        vort(ip,k)=0.0
        enddo
        enddo
c
c.......copy right into of interior of vort......
c
	ipp=0
	do 175 k=1,kz
	do 175 j=1,n
	do 175 i=1,m-1
	ip=(j-1)*m+i
	ipp=ipp+1
 175    vort(ip,k)=right(ipp)
c
c........impose periodicity......
c
        do 176 k=1,kz
        do 176 j=1,n
        ip1=j*m
        ip2=(j-1)*m+1
 176    vort(ip1,k)=vort(ip2,k)
c
c	print *,'ipp=',ipp,' npnt=',npnt
	if(ipp.ne.npnt)then
	print *,'mismatch between ipp and npnt in matveca!!!'
	stop
	endif
c
c........compute the stream......
c
        call delm1(vort,stream)
c
c
        call tlmmod
c
c.......set up the output vector.....
c
	ipp=0
	do k=1,kz
	do j=1,n
	do i=1,m-1
	ip=(j-1)*m+i
	ipp=ipp+1
        alhs(ipp)=vort(ip,k)
        enddo
        enddo
        enddo
c
c	print *,'ipp=',ipp,' npnt=',npnt
	if(ipp.ne.npnt)then
	print *,'mismatch between ipp and npnt in matveca!!!'
	stop
	endif
c
	return
        end
c
	subroutine matvecm(npnt,a,b)
c
c
           include 'prob2.pert'
c
           parameter(nptb=m0*n0*k0)
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
        common/comm1/vmat(k0,k0),vinv(k0,k0)
c
	real a(npnt),b(npnt),c(nptb),d(nptb)
	real tempa(mn0,k0),tempb(mn0,k0)
c
c........copy a into tempa.......
c
        do k=1,kz
        do i=1,msq
        tempa(i,k)=0.0
        enddo
        enddo
c
        ipp=0
        do k=1,kz
        do j=1,n
        do i=1,m
        ipp=ipp+1
        ip=(j-1)*m+i
        tempa(ip,k)=a(ipp)
        enddo
        enddo
        enddo
c
        if(ipp.ne.npnt)then
          print *,'ipp and npnt dont match in matvecm!!!'
          stop
        endif
c
        call egyop(tempa,tempb)
c
        ipp=0
        do k=1,kz
        do j=1,n
        do i=1,m
        ipp=ipp+1
        ip=(j-1)*m+i
        b(ipp)=tempb(ip,k)
        enddo
        enddo
        enddo
c
	return
	end
c
	subroutine solvem(npnt,right,alhs)
c
c
           include 'prob2.pert'
c
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
        common/comm1/vmat(k0,k0),vinv(k0,k0)
c
	 real alhs(npnt),right(npnt)
         real a(mn0,k0),b(mn0,k0),c(mn0,k0),d(mn0,k0)
c
         real dx2(mn0),dy2(mn0)
c
c.......copy right into a and divide by hz....
c
        do k=1,kz
        do i=1,msq
        a(i,k)=0.0
        enddo
        enddo
c
        ipp=0
        do k=1,kz
        do j=1,n
        do i=1,m
        ipp=ipp+1
        ip=(j-1)*m+i
        a(ip,k)=right(ipp)/hz(k)
        enddo
        enddo
        enddo
c
        if(ipp.ne.npnt)then
        print *,'ipp and npnt do not match in solvem!!'
        stop
        endif
c
c.......multiply by vinv........
c
        do k=1,kz
           do ip=1,msq
           b(ip,k)=0.
           do j=1,kz
           b(ip,k)=b(ip,k)+vinv(k,j)*a(ip,j)
           enddo
           enddo
       enddo
c
c.......compute 3-D Laplacian......
c
	do 252 k=1,kz
c
          do ip=1,msq
          r(ip)=b(ip,k)
          enddo
c
          call rel2nd(r,dx2,dy2,m,n,hy)
c
          do ip=1,msq
          c(ip,k)=dx2(ip)+dy2(ip)+eigval(k)*r(ip)
          enddo
c
 252    continue
c
c.......multiply by transpose of vinv........
c
        do k=1,kz
        do ip=1,msq
        a(ip,k)=0.0
        enddo
        enddo
c
        do k=1,kz
           do ip=1,msq
           do j=1,kz
           a(ip,j)=a(ip,j)+vinv(k,j)*c(ip,k)
           enddo
           enddo
       enddo
c
c.......copy a into alhs and divide by layer depths.....
c
        ipp=0
        do k=1,kz
        do j=1,n
        do i=1,m
        ipp=ipp+1
        ip=(j-1)*m+i
        alhs(ipp)=-2.*a(ip,k)/hz(k)
        enddo
        enddo
        enddo
c
        if(ipp.ne.npnt)then
          print *,'ipp and npnt dont match in solvem!!!'
          stop
        endif
c
	return
	end
c
	subroutine egyop(ain,aout)
c
c......This operator defines the total energy norm which is
c      given by sum over all grid boxes of -0.5*vort*stream*hz
c
c
c
           include 'prob2.pert'
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/cnn/n,nm1,nm2,np1
c
        common/comm1/vmat(k0,k0),vinv(k0,k0)
c
	real a(mn0,k0),b(mn0,k0)
	real c(mn0,k0),d(mn0,k0)
	real tempa(mn0,k0),tempb(mn0,k0)
        real ain(mn0,k0),aout(mn0,k0)
c
c......copy ain into a
c
	do k=1,kz
	do j=1,n
	do i=1,m
	ip=(j-1)*m+i
	a(ip,k)=ain(ip,k)
        enddo
        enddo
        enddo
c
c......copy interior of a into interior of b and multiply
c      by layer depths........
c
       do k=1,kz
       do j=2,n-1
       do i=2,m-1
       ip=(j-1)*m+i
       b(ip,k)=a(ip,k)*hz(k)
       enddo
       enddo
       enddo
c
c......set the b.c.s on b....
c
      do k=1,kz
       do j=1,n
       ip1=(j-1)*m+1
       ip2=j*m
       b(ip1,k)=0.
       b(ip2,k)=0.
       enddo
       do i=2,m-1
       ip1=i
       ip2=(n-1)*m+i
       b(ip1,k)=0.0
       b(ip2,k)=0.0
       enddo
      enddo
c
c......project onto the vertical normal modes using transpose of vmat......
c
       do k=1,kz
        do ip=1,msq
         c(ip,k)=0.0
         do j=1,kz
         c(ip,k)=c(ip,k)+vmat(j,k)*b(ip,j)
         enddo
        enddo
       enddo
c
c
c.....invert the Helmholtz eqn.....
c
       do k=1,kz
c
	mb=1
c
      dbasin=hy*nm1
      call pwscrt(0,0.,xbasin,mm1,mb,bda,bdb,0.,dbasin,nm1,1,bdc,bdd
     +,eigval(k),c(1,k),m,1.,ierror,r1)
c
       enddo
c
c......project onto the vertical normal modes using vmat......
c
       do k=1,kz
       do ip=1,msq
       b(ip,k)=0.0
       enddo
       enddo
c
       do k=1,kz
        do ip=1,msq
         do j=1,kz
         b(ip,j)=b(ip,j)+vmat(j,k)*c(ip,k)
         enddo
        enddo
       enddo
c
c......copy interior of b into interior of a and multiply
c      by layer depths........
c
       do k=1,kz
       do ip=1,msq
       a(ip,k)=0.0
       enddo
       enddo
c
       do k=1,kz
       do j=2,n-1
       do i=2,m-1
       ip=(j-1)*m+i
       a(ip,k)=b(ip,k)*hz(k)
       enddo
       enddo
       enddo
c
c......copy a into aout
c
	do k=1,kz
	do j=1,n
	do i=1,m
	ip=(j-1)*m+i
	aout(ip,k)=-0.5*a(ip,k)
        enddo
        enddo
        enddo
c
c
c
	return
	end
c
      function lenstr(string,ls)
c  this routine finds actual length of string by eliminating trailing blanks
      character*1 string(ls)
      common/pltmsg/msglvl
      if(msglvl.lt.-10)print 500,ls,(string(ii),ii=1,ls)
500   format(1x,'***lenstr entered with string length= ',i5,' string= ',
     1 80a1)
      do 100 i=ls,1,-1
      is=i
      if(string(i).ne.char(32)) goto 200
100   continue
      is=0
200   lenstr=is
      return
      end
c
      subroutine sgeco(a,lda,n,ipvt,rcond,z)                                    
c***begin prologue  sgeco                                                       
c***date written   780814   (yymmdd)                                            
c***revision date  820801   (yymmdd)                                            
c***category no.  d2a1                                                          
c***keywords  condition,factor,linear algebra,linpack,matrix                    
c***author  moler, c. b., (u. of new mexico)                                    
c***purpose  factors a real matrix by gaussian elimination and estimates        
c            the condition number of the matrix.                                
c***description                                                                 
c                                                                               
c     sgeco factors a real matrix by gaussian elimination                       
c     and estimates the condition of the matrix.                                
c                                                                               
c     if  rcond  is not needed, sgefa is slightly faster.                       
c     to solve  a*x = b , follow sgeco by sgesl.                                
c     to compute  inverse(a)*c , follow sgeco by sgesl.                         
c     to compute  determinant(a) , follow sgeco by sgedi.                       
c     to compute  inverse(a) , follow sgeco by sgedi.                           
c                                                                               
c     on entry                                                                  
c                                                                               
c        a       real(lda, n)                                                   
c                the matrix to be factored.                                     
c                                                                               
c        lda     integer                                                        
c                the leading dimension of the array  a .                        
c                                                                               
c        n       integer                                                        
c                the order of the matrix  a .                                   
c                                                                               
c     on return                                                                 
c                                                                               
c        a       an upper triangular matrix and the multipliers                 
c                which were used to obtain it.                                  
c                the factorization can be written  a = l*u , where              
c                l  is a product of permutation and unit lower                  
c                triangular matrices and  u  is upper triangular.               
c                                                                               
c        ipvt    integer(n)                                                     
c                an integer vector of pivot indices.                            
c                                                                               
c        rcond   real                                                           
c                an estimate of the reciprocal condition of  a .                
c                for the system  a*x = b , relative perturbations               
c                in  a  and  b  of size  epsilon  may cause                     
c                relative perturbations in  x  of size  epsilon/rcond .         
c                if  rcond  is so small that the logical expression             
c                           1.0 + rcond .eq. 1.0                                
c                is true, then  a  may be singular to working                   
c                precision.  in particular,  rcond  is zero  if                 
c                exact singularity is detected or the estimate                  
c                underflows.                                                    
c                                                                               
c        z       real(n)                                                        
c                a work vector whose contents are usually unimportant.          
c                if  a  is close to a singular matrix, then  z  is              
c                an approximate null vector in the sense that                   
c                norm(a*z) = rcond*norm(a)*norm(z) .                            
c                                                                               
c     linpack.  this version dated 08/14/78 .                                   
c     cleve moler, university of new mexico, argonne national lab.              
c                                                                               
c     subroutines and functions                                                 
c                                                                               
c     linpack sgefa                                                             
c     blas saxpy,sdot,sscal,sasum                                               
c     fortran abs,amax1,sign                                                    
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,            
c                 *linpack users  guide*, siam, 1979.                           
c***routines called  sasum,saxpy,sdot,sgefa,sscal                               
c***end prologue  sgeco                                                         
      integer lda,n,ipvt(1)                                                     
      real a(lda,1),z(1)                                                        
      real rcond                                                                
c                                                                               
      real sdot,ek,t,wk,wkm                                                     
      real anorm,s,sasum,sm,ynorm                                               
      integer info,j,k,kb,kp1,l                                                 
c                                                                               
c     compute 1-norm of a                                                       
c                                                                               
c***first executable statement  sgeco                                           
      anorm = 0.0e0                                                             
      do 10 j = 1, n                                                            
         anorm = amax1(anorm,sasum(n,a(1,j),1))                                 
   10 continue                                                                  
c                                                                               
c     factor                                                                    
c                                                                               
      call sgefa(a,lda,n,ipvt,info)                                             
c                                                                               
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .                      
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .          
c     trans(a)  is the transpose of a .  the components of  e  are              
c     chosen to cause maximum local growth in the elements of w  where          
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid            
c     overflow.                                                                 
c                                                                               
c     solve trans(u)*w = e                                                      
c                                                                               
      ek = 1.0e0                                                                
      do 20 j = 1, n                                                            
         z(j) = 0.0e0                                                           
   20 continue                                                                  
      do 100 k = 1, n                                                           
         if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))                               
         if (abs(ek-z(k)) .le. abs(a(k,k))) go to 30                            
            s = abs(a(k,k))/abs(ek-z(k))                                        
            call sscal(n,s,z,1)                                                 
            ek = s*ek                                                           
   30    continue                                                               
         wk = ek - z(k)                                                         
         wkm = -ek - z(k)                                                       
         s = abs(wk)                                                            
         sm = abs(wkm)                                                          
         if (a(k,k) .eq. 0.0e0) go to 40                                        
            wk = wk/a(k,k)                                                      
            wkm = wkm/a(k,k)                                                    
         go to 50                                                               
   40    continue                                                               
            wk = 1.0e0                                                          
            wkm = 1.0e0                                                         
   50    continue                                                               
         kp1 = k + 1                                                            
         if (kp1 .gt. n) go to 90                                               
            do 60 j = kp1, n                                                    
               sm = sm + abs(z(j)+wkm*a(k,j))                                   
               z(j) = z(j) + wk*a(k,j)                                          
               s = s + abs(z(j))                                                
   60       continue                                                            
            if (s .ge. sm) go to 80                                             
               t = wkm - wk                                                     
               wk = wkm                                                         
               do 70 j = kp1, n                                                 
                  z(j) = z(j) + t*a(k,j)                                        
   70          continue                                                         
   80       continue                                                            
   90    continue                                                               
         z(k) = wk                                                              
  100 continue                                                                  
      s = 1.0e0/sasum(n,z,1)                                                    
      call sscal(n,s,z,1)                                                       
c                                                                               
c     solve trans(l)*y = w                                                      
c                                                                               
      do 120 kb = 1, n                                                          
         k = n + 1 - kb                                                         
         if (k .lt. n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)              
         if (abs(z(k)) .le. 1.0e0) go to 110                                    
            s = 1.0e0/abs(z(k))                                                 
            call sscal(n,s,z,1)                                                 
  110    continue                                                               
         l = ipvt(k)                                                            
         t = z(l)                                                               
         z(l) = z(k)                                                            
         z(k) = t                                                               
  120 continue                                                                  
      s = 1.0e0/sasum(n,z,1)                                                    
      call sscal(n,s,z,1)                                                       
c                                                                               
      ynorm = 1.0e0                                                             
c                                                                               
c     solve l*v = y                                                             
c                                                                               
      do 140 k = 1, n                                                           
         l = ipvt(k)                                                            
         t = z(l)                                                               
         z(l) = z(k)                                                            
         z(k) = t                                                               
         if (k .lt. n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)                    
         if (abs(z(k)) .le. 1.0e0) go to 130                                    
            s = 1.0e0/abs(z(k))                                                 
            call sscal(n,s,z,1)                                                 
            ynorm = s*ynorm                                                     
  130    continue                                                               
  140 continue                                                                  
      s = 1.0e0/sasum(n,z,1)                                                    
      call sscal(n,s,z,1)                                                       
      ynorm = s*ynorm                                                           
c                                                                               
c     solve  u*z = v                                                            
c                                                                               
      do 160 kb = 1, n                                                          
         k = n + 1 - kb                                                         
         if (abs(z(k)) .le. abs(a(k,k))) go to 150                              
            s = abs(a(k,k))/abs(z(k))                                           
            call sscal(n,s,z,1)                                                 
            ynorm = s*ynorm                                                     
  150    continue                                                               
         if (a(k,k) .ne. 0.0e0) z(k) = z(k)/a(k,k)                              
         if (a(k,k) .eq. 0.0e0) z(k) = 1.0e0                                    
         t = -z(k)                                                              
         call saxpy(k-1,t,a(1,k),1,z(1),1)                                      
  160 continue                                                                  
c     make znorm = 1.0                                                          
      s = 1.0e0/sasum(n,z,1)                                                    
      call sscal(n,s,z,1)                                                       
      ynorm = s*ynorm                                                           
c                                                                               
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm                                 
      if (anorm .eq. 0.0e0) rcond = 0.0e0                                       
      return                                                                    
      end                                                                       

      subroutine sgedi(a,lda,n,ipvt,det,work,job)                               
c***begin prologue  sgedi                                                       
c***date written   780814   (yymmdd)                                            
c***revision date  820801   (yymmdd)                                            
c***category no.  d2a1,d3a1                                                     
c***keywords  determinant,factor,inverse,linear algebra,linpack,matrix          
c***author  moler, c. b., (u. of new mexico)                                    
c***purpose  computes the determinant and inverse of a matrix                   
c            using the factors computed by sgeco or sgefa.                      
c***description                                                                 
c                                                                               
c     sgedi computes the determinant and inverse of a matrix                    
c     using the factors computed by sgeco or sgefa.                             
c                                                                               
c     on entry                                                                  
c                                                                               
c        a       real(lda, n)                                                   
c                the output from sgeco or sgefa.                                
c                                                                               
c        lda     integer                                                        
c                the leading dimension of the array  a .                        
c                                                                               
c        n       integer                                                        
c                the order of the matrix  a .                                   
c                                                                               
c        ipvt    integer(n)                                                     
c                the pivot vector from sgeco or sgefa.                          
c                                                                               
c        work    real(n)                                                        
c                work vector.  contents destroyed.                              
c                                                                               
c        job     integer                                                        
c                = 11   both determinant and inverse.                           
c                = 01   inverse only.                                           
c                = 10   determinant only.                                       
c                                                                               
c     on return                                                                 
c                                                                               
c        a       inverse of original matrix if requested.                       
c                otherwise unchanged.                                           
c                                                                               
c        det     real(2)                                                        
c                determinant of original matrix if requested.                   
c                otherwise not referenced.                                      
c                determinant = det(1) * 10.0**det(2)                            
c                with  1.0 .le. abs(det(1)) .lt. 10.0                           
c                or  det(1) .eq. 0.0 .                                          
c                                                                               
c     error condition                                                           
c                                                                               
c        a division by zero will occur if the input factor contains             
c        a zero on the diagonal and the inverse is requested.                   
c        it will not occur if the subroutines are called correctly              
c        and if sgeco has set rcond .gt. 0.0 or sgefa has set                   
c        info .eq. 0 .                                                          
c                                                                               
c     linpack.  this version dated 08/14/78 .                                   
c     cleve moler, university of new mexico, argonne national lab.              
c                                                                               
c     subroutines and functions                                                 
c                                                                               
c     blas saxpy,sscal,sswap                                                    
c     fortran abs,mod                                                           
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,            
c                 *linpack users  guide*, siam, 1979.                           
c***routines called  saxpy,sscal,sswap                                          
c***end prologue  sgedi                                                         
      integer lda,n,ipvt(1),job                                                 
      real a(lda,1),det(2),work(1)                                              
c                                                                               
      real t                                                                    
      real ten                                                                  
      integer i,j,k,kb,kp1,l,nm1                                                
c                                                                               
c     compute determinant                                                       
c                                                                               
c***first executable statement  sgedi                                           
      if (job/10 .eq. 0) go to 70                                               
         det(1) = 1.0e0                                                         
         det(2) = 0.0e0                                                         
         ten = 10.0e0                                                           
         do 50 i = 1, n                                                         
            if (ipvt(i) .ne. i) det(1) = -det(1)                                
            det(1) = a(i,i)*det(1)                                              
c        ...exit                                                                
            if (det(1) .eq. 0.0e0) go to 60                                     
   10       if (abs(det(1)) .ge. 1.0e0) go to 20                                
               det(1) = ten*det(1)                                              
               det(2) = det(2) - 1.0e0                                          
            go to 10                                                            
   20       continue                                                            
   30       if (abs(det(1)) .lt. ten) go to 40                                  
               det(1) = det(1)/ten                                              
               det(2) = det(2) + 1.0e0                                          
            go to 30                                                            
   40       continue                                                            
   50    continue                                                               
   60    continue                                                               
   70 continue                                                                  
c                                                                               
c     compute inverse(u)                                                        
c                                                                               
      if (mod(job,10) .eq. 0) go to 150                                         
         do 100 k = 1, n                                                        
            a(k,k) = 1.0e0/a(k,k)                                               
            t = -a(k,k)                                                         
            call sscal(k-1,t,a(1,k),1)                                          
            kp1 = k + 1                                                         
            if (n .lt. kp1) go to 90                                            
            do 80 j = kp1, n                                                    
               t = a(k,j)                                                       
               a(k,j) = 0.0e0                                                   
               call saxpy(k,t,a(1,k),1,a(1,j),1)                                
   80       continue                                                            
   90       continue                                                            
  100    continue                                                               
c                                                                               
c        form inverse(u)*inverse(l)                                             
c                                                                               
         nm1 = n - 1                                                            
         if (nm1 .lt. 1) go to 140                                              
         do 130 kb = 1, nm1                                                     
            k = n - kb                                                          
            kp1 = k + 1                                                         
            do 110 i = kp1, n                                                   
               work(i) = a(i,k)                                                 
               a(i,k) = 0.0e0                                                   
  110       continue                                                            
            do 120 j = kp1, n                                                   
               t = work(j)                                                      
               call saxpy(n,t,a(1,j),1,a(1,k),1)                                
  120       continue                                                            
            l = ipvt(k)                                                         
            if (l .ne. k) call sswap(n,a(1,k),1,a(1,l),1)                       
  130    continue                                                               
  140    continue                                                               
  150 continue                                                                  
      return                                                                    
      end                                                                       
      subroutine sgefa(a,lda,n,ipvt,info)
c***begin prologue  sgefa
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  factor,linear algebra,linpack,matrix
c***author  moler, c. b., (u. of new mexico)
c***purpose  factors a real matrix by gaussian elimination.
c***description
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u , where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  isamax,saxpy,sscal
c***end prologue  sgefa
      integer lda,n,ipvt(1),info
      real a(lda,1)
c
      real t
      integer isamax,j,k,kp1,l,nm1
c
c     gaussian elimination with partial pivoting
c
c***first executable statement  sgefa
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end
c
	subroutine cjet(us0,stream,vort)
c
c.......creates the vorticity field for a zonally periodic zonal
c       get with Gaussian cross-section in U.......
c
           include "para.in"
c
	parameter(im=m0,jm=n0,km=k0)
c
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
c
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
c
      common/dynam/alpha,beta,rf
c
	real vort(im,jm,km),stream(im,jm,km),gamma(jm)
c
	real uprof(jm)
c
	j0=jm/2+1
c
	pi=4.*atan(1.0)
c
	read(8,*)u0
	read(8,*)ghor
c
	us0=u0
c
	ghor=ghor/dhor
	dy=hy
	u0=u0/v0
c
	do k=1,km
	  do j=1,jm
	    do i=1,im
	      stream(i,j,k)=0.0
	      vort(i,j,k)=0.0
            enddo
          enddo
        enddo
c
	do k=1,km
	if(k.eq.1)then
	do j=1,jm
	y=(float(j)-j0)*dy
	fac=-1.*y*y/(ghor*ghor)
	uprof(j)=u0*exp(fac)
	gamma(j)=beta-alpha*2.*u0*exp(fac)*(-2.*fac-1.)/(ghor*ghor)
	do i=1,im
	stream(i,j,k)=-0.5*ghor*sqrt(pi)*u0*erf(y/ghor)
  	vort(i,j,k)=2.*y*u0*exp(fac)/(ghor*ghor)
        enddo
        enddo
        endif
        enddo
c
	write(23)uprof
	write(23)gamma
c
	return
	end
c
      subroutine eigrg1 (nm,n,a,iflag,wr,wi,z,work1,work2,ierr)
      real a(nm,n),wr(n),wi(n),work1(nm,n),work2(1),z(nm,1) 
      do 50 i = 1,n
      do 50 j = 1,n
      work1(i,j) = a(i,j)
   50 continue
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
CXXX
      call balanc(nm,n,work1,kwrk1,kwrk2,work2)
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
CXXX
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
      kwrk3=work2(n+3)
CXXX
      call elmhes(nm,n,kwrk1,kwrk2,work1,kwrk3)
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
      work2(n+3)=kwrk3
CXXX
      if (iflag .ne. 0) go to 100
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
CXXX
      call hqr(nm,n,kwrk1,kwrk2,work1,wr,wi,ierr)
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
CXXX
      if (ierr .ne. 0) go to 400
      return
  100 if (iflag .ne. n) go to 200
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
      kwrk3=work2(n+3)
CXXX
      call eltran(nm,n,kwrk1,kwrk2,work1,kwrk3,z)
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
      work2(n+3)=kwrk3
CXXX
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
CXXX
      call hqr2(nm,n,kwrk1,kwrk2,work1,wr,wi,z,ierr)
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
CXXX
      if (ierr .ne. 0) go to 400
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
CXXX
      call balbak(nm,n,kwrk1,kwrk2,work2,n,z) 
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
CXXX
      return
  200 do 300 i = 1,n
         do 300 j = 1,n
      z(i,j) = work1(i,j)
  300 continue
CXXX
      kwrk1=work2(n+1)
      kwrk2=work2(n+2)
CXXX
      call hqr(nm,n,kwrk1,kwrk2,z,wr,wi,ierr)
CXXX
      work2(n+1)=kwrk1
      work2(n+2)=kwrk2
CXXX
      if (ierr .eq. 0) return
  400 ierr = iabs(ierr) + 32
      print 1001
 1001 format(' eigrg1 failed to compute all eigenvalues')
      return
      end
c
      subroutine pwscrt (intl,a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,
     1                   bdd,elmbda,f,idimf,pertrb,ierror,w)
      dimension        bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     1                 f(idimf,n) ,w(1)
c     1                 f(idimf,1) ,w(1)
      ierror = 0
      if (a.ge.b) ierror = 1
      if (mbdcnd.lt.0.or.mbdcnd.gt.4) ierror = 2
      if (c.ge.d) ierror = 3
      if (ncheck(n).gt.1.or.n.le.2) ierror = 4
      if (nbdcnd.lt.0.or.nbdcnd.gt.4) ierror = 5
      if (ierror.ne.0) return
  101 nperod = nbdcnd
      nby2 = 0
      deltax = (b-a)/m
      twdelx = 2./deltax
      delxsq = 1./deltax**2
      deltay = (d-c)/n
      twdely = 2./deltay
      delysq = 1./deltay**2
      np = nbdcnd+1
      np1 = n+1
      mp = mbdcnd+1
      mp1 = m+1
      nstart = 1
      nstop = n
      nskip = 1
      go to (105,102,103,104,105),np
  102 nstart = 2
      go to 105
  103 nstart = 2
  104 nstop = np1
      nskip = 2
  105 nunk = nstop-nstart+1
      mstart = 1
      mstop = m
      mskip = 1
      go to (122,106,107,111,112),mp
  106 mstart = 2
      go to 108
  107 mstart = 2
      mstop = mp1
      mskip = 2
  108     do 110 j=nstart,nstop
  109     f(2,j) = f(2,j)-f(1,j)*delxsq
  110     continue
      go to 115
  111 mstop = mp1
      mskip = 2
  112     do 114 j=nstart,nstop
  113     f(1,j) = f(1,j)+bda(j)*twdelx
  114     continue
  115 go to (116,119),mskip
  116     do 118 j=nstart,nstop
  117     f(m,j) = f(m,j)-f(mp1,j)*delxsq
  118     continue
      go to 122
  119     do 121 j=nstart,nstop
  120     f(mp1,j) = f(mp1,j)-bdb(j)*twdelx
  121     continue
  122 munk = mstop-mstart+1
      go to (136,123,123,126,126),np
  123     do 125 i=mstart,mstop
  124     f(i,2) = f(i,2)-f(i,1)*delysq
  125     continue
      go to 129
  126     do 128 i=mstart,mstop
  127     f(i,1) = f(i,1)+bdc(i)*twdely
  128     continue
  129 go to (130,133),nskip
  130     do 132 i=mstart,mstop
  131     f(i,n) = f(i,n)-f(i,np1)*delysq
  132     continue
      go to 136
  133     do 135 i=mstart,mstop
  134     f(i,np1) = f(i,np1)-bdd(i)*twdely
  135     continue
  136 delysq = deltay*deltay
          do 139 i=mstart,mstop
              do 138 j=nstart,nstop
  137         f(i,j) = f(i,j)*delysq
  138         continue
  139     continue
  140 id1 = 5.25*nunk+5*munk
      id2 = id1+munk
      id3 = id2+munk
      id4 = id3+munk
      s = delysq*delxsq
      st2 = 2.*s
          do 142 i=1,munk
          w(id1+i) = s
          w(id2+i) = -st2+elmbda*delysq
  141     w(id3+i) = s
  142     continue
      go to (146,146,143,144,145),mp
  143 w(id2) = st2
      go to 146
  144 w(id2) = st2
  145 w(id3+1) = st2
  146 continue
      if (nbdcnd.ne.4) go to 151
      irev = 1
      nperod = 2
      nby2 = n/2
  147     do 150 j=1,nby2
              do 149 i=mstart,mstop
              a1 = f(i,j)
              f(i,j) = f(i,n+1-j)
  148         f(i,n+1-j) = a1
  149         continue
  150     continue
      go to (164,168),irev
  151 continue
      if (elmbda) 164,153,152
  152 ierror = 6
      go to 164
  153 if ((nbdcnd.eq.0.or.nbdcnd.eq.3).and.(mbdcnd.eq.0.or.mbdcnd.eq.3))
     1    go to 154
      go to 164
  154 a1 = 1.
      a2 = 1.
      if (nbdcnd.eq.3) a2 = 2.
      if (mbdcnd.eq.3) a1 = 2.
      s1 = 0.
      msp1 = mstart+1
      mstm1 = mstop-1
      nsp1 = nstart+1
      nstm1 = nstop-1
          do 158 j=nsp1,nstm1
          s = 0.
              do 156 i=msp1,mstm1
  155         s = s+f(i,j)
  156         continue
  157     s1 = s1+s*a1+f(mstart,j)+f(mstop,j)
  158     continue
      s1 = a2*s1
      s = 0.
          do 160 i=msp1,mstm1
  159     s = s+f(i,nstart)+f(i,nstop)
  160     continue
      s1 = s1+s*a1+f(mstart,nstart)+f(mstart,nstop)+f(mstop,nstart)+
     1     f(mstop,nstop)
      s = (2.+(nunk-2)*a2)*(2.+(munk-2)*a1)
      pertrb = s1/s
          do 163 j=nstart,nstop
              do 162 i=mstart,mstop
  161         f(i,j) = f(i,j)-pertrb
  162         continue
  163     continue
      pertrb = pertrb/delysq
  164 call pois (intl,nperod,nunk,mbdcnd,munk,w(id1+1),w(id2+1),
     1           w(id3+1),idimf,f(mstart,nstart),w)
      if (nby2.eq.0) go to 165
      irev = 2
      go to 147
  165 continue
      if (nbdcnd.ne.0) go to 168
          do 167 i=mstart,mstop
  166     f(i,np1) = f(i,1)
  167     continue
  168 if (mbdcnd.ne.0) go to 171
          do 170 j=nstart,nstop
  169     f(mp1,j) = f(1,j)
  170     continue
      if (nbdcnd.eq.0) f(mp1,np1) = f(1,np1)
  171 continue
      return
      end
c
      subroutine apwscrt (intl,a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,
     1                   bdd,elmbda,f,idimf,pertrb,ierror,w,zz)
      dimension        bda(1)     ,bdb(1)     ,bdc(1)     ,bdd(1)     ,
     1                 f(idimf,1) ,w(1),zz(1)
      ierror = 0
      if (a.ge.b) ierror = 1
      if (mbdcnd.lt.0.or.mbdcnd.gt.4) ierror = 2
      if (c.ge.d) ierror = 3
      if (ncheck(n).gt.1.or.n.le.2) ierror = 4
      if (nbdcnd.lt.0.or.nbdcnd.gt.4) ierror = 5
      if (ierror.ne.0) return
  101 nperod = nbdcnd
      nby2 = 0
      deltax = (b-a)/m
      twdelx = 2./deltax
      delxsq = 1./deltax**2
      deltay = (d-c)/n
      twdely = 2./deltay
      delysq = 1./deltay**2
      np = nbdcnd+1
      np1 = n+1
      mp = mbdcnd+1
      mp1 = m+1
      nstart = 1
      nstop = n
      nskip = 1
      go to (105,102,103,104,105),np
  102 nstart = 2
      go to 105
  103 nstart = 2
  104 nstop = np1
      nskip = 2
  105 nunk = nstop-nstart+1
      mstart = 1
      mstop = m
      mskip = 1
      go to (122,106,107,111,122),mp
  106 mstart = 2
      go to 108
  107 mstart = 2
      mstop = mp1
      mskip = 2
  108     continue
      go to 122
  111 mstop = mp1
      mskip = 2
  122 munk = mstop-mstart+1
c
      id1 = 5.25*nunk+5*munk
      id2 = id1+munk
      id3 = id2+munk
      id4 = id3+munk
c     s = delysq*delxsq
      s = deltay*deltay*delxsq
      st2 = 2.*s
          do 142 i=1,munk
          w(id1+i) = s
          w(id2+i) = -st2+elmbda*deltay*deltay
  141     w(id3+i) = s
  142     continue
      go to (146,146,143,144,145),mp
  143 w(id2) = st2
      go to 146
  144 w(id2) = st2
  145 w(id3+1) = st2
  146 continue
c
	if(nbdcnd.eq.4)goto 147
c
      call apois (intl,nperod,nunk,mbdcnd,munk,w(id1+1),w(id2+1),
     1           w(id3+1),idimf,f(mstart,nstart),w,zz)
      if (nby2.eq.0) go to 165
      irev = 2
      go to 1470
  165 continue
c
      if (elmbda) 164,153,152
  152 ierror = 6
      go to 164
  153 if ((nbdcnd.eq.0.or.nbdcnd.eq.3).and.(mbdcnd.eq.0.or.mbdcnd.eq.3))
     1    go to 154
      go to 164
  154 a1 = 1.
      a2 = 1.
      if (nbdcnd.eq.3) a2 = 2.
      if (mbdcnd.eq.3) a1 = 2.
c
      msp1 = mstart+1
      mstm1 = mstop-1
      nsp1 = nstart+1
      nstm1 = nstop-1
c
	pertrb=0.0
	do 163 j=nstart,nstop
	do 163 i=mstart,mstop
 163    pertrb=pertrb-f(i,j)
c
      s = (2.+(nunk-2)*a2)*(2.+(munk-2)*a1)
c
	s1p=pertrb/s
	s=a1*s1p
	s1=s1p*a2
	s2=a1*s1
c
	f(mstart,nstart)=f(mstart,nstart)+s1p
	f(mstart,nstop)=f(mstart,nstop)+s1p
	f(mstop,nstart)=f(mstop,nstart)+s1p
	f(mstop,nstop)=f(mstop,nstop)+s1p
c
	do 160 i=msp1,mstm1
	f(i,nstart)=f(i,nstart)+s
	f(i,nstop)=f(i,nstop)+s
 160    continue
c
	do 158 j=nsp1,nstm1
	f(mstart,j)=f(mstart,j)+s1
	f(mstop,j)=f(mstop,j)+s1
	do 1600 i=msp1,mstm1
 1600   f(i,j)=f(i,j)+s2
 158    continue
c
 164    continue
c
        adelysq = deltay*deltay
          do 139 i=mstart,mstop
              do 138 j=nstart,nstop
              f(i,j) = f(i,j)*adelysq
  138         continue
  139     continue
c
	goto 151
c
 147  irev = 1
      nperod = 2
      nby2 = n/2
 1470    do 150 j=1,nby2
              do 149 i=mstart,mstop
              a1 = f(i,n+1-j)
              f(i,n+1-j) = f(i,j)
              f(i,j) = a1
  149         continue
  150     continue
      go to (168,164),irev
  151 continue
c
      go to (1220,1060,1060,115,115),mp
 1060     do 110 j=nstart,nstop
c         f(1,j) = -f(2,j)*delxsq
c         f(1,j) = -f(2,j)*delxsq+f(1,j)
 110      continue
  115 go to (116,1220),mskip
  116     do 118 j=nstart,nstop
c         f(mp1,j) = -f(m,j)*delxsq
c         f(mp1,j) = -f(m,j)*delxsq+f(mp1,j)
  118     continue
 1220 go to (1360,1230,1230,129,129),np
 1230     do 125 i=mstart,mstop
c         f(i,1) = -f(i,2)*delysq
c         f(i,1) = -f(i,2)*delysq+f(i,1)
  125     continue
  129 go to (130,1360),nskip
  130     do 132 i=mstart,mstop
c 131     f(i,np1) = -f(i,n)*delysq
c 131     f(i,np1) = -f(i,n)*delysq+f(i,np1)
  132     continue
 1360     continue
c
      if (mbdcnd.ne.0) go to 171
          do 170 j=nstart,nstop
          f(mp1,j) =f(1,j)
  170     continue
      if (nbdcnd.eq.0) f(mp1,np1) =f(1,np1)
  171 continue
c
      if (nbdcnd.ne.0) go to 168
          do 167 i=mstart,mstop
          f(i,np1) =f(i,1)
  167     continue
  168     continue
c
c
	return
	end
c
      subroutine shprox(g,m,nord)
c
c --- g = portion of a larger array to be filtered
c --- m = length of this array
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter it using subroutine filter
c
c      dimension g(1)
      dimension g(m)
      mm1=m-1
      mm2=m-2
      do 100 kord=1,nord
      g2=g(2)
        do 10 i=1,mm1
          g(i)=g(i+1)-g(i)
   10   continue
	g(m)=g2-g(m)
	g12=g(m-1)
        do 20 ii=1,mm1
          i=m+1-ii
          g(i)=g(i)-g(i-1)
   20   continue
	g(1)=g(1)-g12
  100 continue
      return
      end
c
      subroutine ashprox(g,m,nord)
c
c --- g = portion of a larger array to be filtered
c --- m = length of this array
c --- nord = order of the filter
c
c --- this subroutine calculates the factor to be added to
c --- a field in order to filter it using subroutine filter
c
c      dimension g(1)
      dimension g(m)
      mm1=m-1
      mm2=m-2
      do 100 kord=1,nord
c     g1=g(1)
c<<<<<<<<>>>>>>
	g2=g(2)
c<<<<<<>>>>>>>>>
        do 10 i=1,mm1
          g(i)=g(i)-g(i+1)
   10   continue
c	g(m-1)=g(m-1)-g1
c<<<<<<<>>>>>>
	g(m)=g(m)-g2
c<<<<<<<>>>>>>
c	g12=g(m)
c<<<<<<<<>>>>>>
	g12=g(m-1)
c<<<<<<>>>>>>>>
        do 20 ii=1,mm1
          i=m+1-ii
          g(i)=g(i-1)-g(i)
   20   continue
c	g(2)=g(2)+g12
c	g(1)=-1.*g(1)
c<<<<<>>>>>>>>
	g(1)=g12-g(1)
c<<<<<>>>>>>>
  100 continue
      return
      end
c
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,its,low,mp2,enm2,ierr
      real h(nm,n),wr(n),wi(n)
      real p,q,r,s,t,w,x,y,zz,norm,machep
      logical notlas
      machep = 2.**(-47)
      ierr = 0
      norm = 0.0
      k = 1
      do 50 i = 1, n
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0
   50 continue
      en = igh
      t = 0.0
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0) s = norm
         if (abs(h(l,l-1)) .le. machep * s) go to 100
   80 continue
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (its .eq. 30) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
      t = t + x
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75 * s
      y = x
      w = -0.4375 * s * s
  130 its = its + 1
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. machep * abs(p)
     1    * (abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))) go to 150
  140 continue
  150 mp2 = m + 2
      do 160 i = mp2, en
         h(i,i-2) = 0.0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0
  160 continue
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1) 
         r = 0.0
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         do 210 j = k, en
            p = h(k,j) + q * h(k+1,j)
            if (.not. notlas) go to 200
            p = p + r * h(k+2,j)
            h(k+2,j) = h(k+2,j) - p * zz
  200       h(k+1,j) = h(k+1,j) - p * y
            h(k,j) = h(k,j) - p * x
  210    continue
         j = min0(en,k+3)
         do 230 i = l, j
            p = x * h(i,k) + y * h(i,k+1)
            if (.not. notlas) go to 220
            p = p + zz * h(i,k+2)
            h(i,k+2) = h(i,k+2) - p * r
  220       h(i,k+1) = h(i,k+1) - p * q
            h(i,k) = h(i,k) - p
  230    continue
  260 continue
      go to 70
  270 wr(en) = x + t
      wi(en) = 0.0
      en = na
      go to 60
  280 p = (y - x) / 2.0
      q = p * p + w
      zz = sqrt(abs(q))
      x = x + t
      if (q .lt. 0.0) go to 320
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0) wr(en) = x - w / zz
      wi(na) = 0.0
      wi(en) = 0.0
      go to 330
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
 1000 ierr = en
 1001 return
      end
c
	subroutine tlmmod
c
c
	include 'prob2.pert'
c
      character*4 restid,runid
      character*70 idstrg
c
        common/newbg/amms(mn0,k0),ammv(mn0,k0)
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/grate/gre(mn0),grip(mn0),tezm(mn0),aipzm(mn0)
c
c@@@@@@@@@THE FOLLOWING COMMON BLOCK IS ONLY IN CLINIC!!!!
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/egl2/pl,akl,alls,allv
c
	common/mmde/jmode
c
	common/nrest/irst
c
	common/adrst/ityrst
c
	common/cnorm/jnorm
c
	common/jctrl/ktim,icstart,iczero
c
	common/gcord/ir1,ir2,jr1,jr2
c
        common/nmmt/nmit,nmiter,idone
c
c>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/count/icnt,itcnt,ibcnt
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),topo(mn0),
     2     topd(mn0),botd(mn0)
c
	common/atop/radt(mn0),s2t(mn0),st(mn0),fyt(mx0,4)
c
	common/abot/radb(mn0),rad1b(mn0),rad2b(mn0),s2b(mn0),sb(mn0)
     1       ,fyb(mx0,4)
c
c
c
	common/adm/tiam(k0,k0),taml(k0,k0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
c
	common/grad/gomeg(mn0,k0),gphi(mn0,k0),ggam(mn0),
     1     gdel(mn0),dm(mn0,k0),dp(mn0,k0),dg(mn0),dd(mn0)
     2    ,cgf
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/steps/sumt,sumb,cost
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/ctest/ctol
c
	common/swoi/iter
c
	common/toff/tadj0
c
	common/istar/icall
c
        common/ccnt/callcnt
c
        common/iwfq/ibgw,iclw
c
	common/bdata/scrp(n0,k0,2),vcrp(n0,k0,2)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/lan1/maxt,ifts
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
       real savstm1(mn0,k0),savvrt1(mn0,k0),savbot1(mn0),savtop1(mn0)
       real savstm2(mn0,k0),savvrt2(mn0,k0),savbot2(mn0),savtop2(mn0)
c
c
c.......clear prvjac.......
c
        do k=1,kz
        do i=1,msq
        prvjac(i,k)=0.0
        enddo
        enddo
c
        tim=0.
c
c.......write initial fields.....
c
        if(nmit.ne.0)then
          
            write(23)stream
            write(23)vort
        endif
c
c --- top of main computation loop
c
      do 9999 itime=1,maxt
c
      tim=tim+dt
c
      icalc=itime
c
	if(mod(itime-1,8).eq.0)then
c	write(53)stream
c	write(53)vort
	endif
c
c
c      if(mod(icalc,10).eq.0)print *,'icalc=',icalc,' in qg-model'
c
c
c --- call nlop3d to perform calculations
c
      call nlop3d
c
c
c.....calc. scalar diagnostics every timestep during last iteration.....
c
c
        if(nmit.eq.nmiter.and.kz.gt.1)call scalar
        if(nmit.eq.nmiter.and.kz.eq.1)call scalar
c
c
c
c --- return to top of main loop
c
 9999 continue
c
c.......write final fields.....
c
        if(nmit.ne.0)then
            write(23)stream
            write(23)vort
        endif
c
	return
        end
c
	subroutine rdback(iskip,isign,savstm,savvrt,savbot,savtop)
c
c
	include 'prob2.pert'
c
      character*4 restid,runid
      character*70 idstrg
c
        common/newbg/amms(mn0,k0),ammv(mn0,k0)
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
	common/grate/gre(mn0),grip(mn0),tezm(mn0),aipzm(mn0)
c
c@@@@@@@@@THE FOLLOWING COMMON BLOCK IS ONLY IN CLINIC!!!!
c
	common/bgsvn/btopd(mn0),bbotd(mn0)
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/egl2/pl,akl,alls,allv
c
	common/mmde/jmode
c
	common/nrest/irst
c
	common/adrst/ityrst
c
	common/cnorm/jnorm
c
	common/jctrl/ktim,icstart,iczero
c
	common/gcord/ir1,ir2,jr1,jr2
c
        common/nmmt/nmit,nmiter,idone
c
c>>>>>>>>>>
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/count/icnt,itcnt,ibcnt
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),topo(mn0),
     2     topd(mn0),botd(mn0)
c
	common/atop/radt(mn0),s2t(mn0),st(mn0),fyt(mx0,4)
c
	common/abot/radb(mn0),rad1b(mn0),rad2b(mn0),s2b(mn0),sb(mn0)
     1       ,fyb(mx0,4)
c
c
c
	common/adm/tiam(k0,k0),taml(k0,k0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
c
	common/grad/gomeg(mn0,k0),gphi(mn0,k0),ggam(mn0),
     1     gdel(mn0),dm(mn0,k0),dp(mn0,k0),dg(mn0),dd(mn0)
     2    ,cgf
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/steps/sumt,sumb,cost
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/ctest/ctol
c
	common/swoi/iter
c
	common/toff/tadj0
c
	common/istar/icall
c
        common/ccnt/callcnt
c
        common/iwfq/ibgw,iclw
c
	common/bdata/scrp(n0,k0,2),vcrp(n0,k0,2)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	common/lan1/maxt,ifts
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
	real bigs(mnbig,k0),bigv(mnbig,k0)
     &       ,bigt(mnbig),bigb(mnbig)
c
        real savstm(mn0,k0),savvrt(mn0,k0),savtop(mn0),savbot(mn0)
c
c
c
c.......read in background field......
c
	if(alpha.lt.1.e-20)goto 9876
c
	if(ktim.eq.0.and.icalc.gt.1)goto 9876
c
c.....read in qg fields at the required timestep.....
c
c     iadj=2 : bvrt
c     iadj=3 : btopd
c     iadj=4 : bbotd
c     iadj=5 : bstm
c
c	print *,'iskip in clinic=',iskip
c
	read(99,rec=iskip)itc,iadj,(bigv(ip,1),ip=1,mnbig)
c 	print *,'itc=',itc,' iadj=',iadj,' icalc=',icalc
c	write(6,*)'itc=',itc,'  iadj=',iadj
c
c.......check that the timestep of the qg fields is the same as that of
c	adjoint model
c
	if(itc-icstart-iclw+1.ne.icalc+isign*ibgw.and.
     &       icalc.ne.1.and.icalc.ne.maxt)then
         print *,'mismatch in timesteps in rdback: itc, icstart, icalc',
     &     itc,icstart,icalc
         stop
        endif
c
c.......check that the first field read is vorticity.....
c
	if(iadj.ne.2)then
          print *,'2 wrong variable read in rdback: iadj=',iadj
          stop
        endif
C>>>>>>>
c	if(iftop.eq.0)goto 5011
C>>>>>>>
	iskip=iskip+1
	read(99,rec=iskip)itc,iadj,bigt
	if(iadj.ne.3)then
          print *,'btopd mismatch in rdback: iadj=',iadj
          stop
        endif
 5011	do 5033 k=2,kz
	iskip=iskip+1
 5033   read(99,rec=iskip)itc,iadj,(bigv(ip,k),ip=1,mnbig)
C>>>>>>>
c       if(ifbot.eq.0)goto 5022
C>>>>>>>
	iskip=iskip+1
	read(99,rec=iskip)itc,iadj,bigb
c	print *,'itc=',itc,' iadj=',iadj
	if(iadj.ne.4)then
          print *,'bbotd mismatch in rdback: iadj=',iadj
          stop
        endif
 5022   continue
	do 8033 k=1,kz
	iskip=iskip+1
	read(99,rec=iskip)itc,iadj,(bigs(ip,k),ip=1,mnbig)
 8033   continue
c
c.......check to see if last array read is bstm...
c
	if(iadj.ne.5)then
          print *,'bstm not found in rdback: iadj=',iadj
          stop
        endif
c
c
c
c........copy bigs, bigv, bigt and bigb into small arrays....
c
	ip=0
	nock=0
	do 8766 j=jr1,jr2
	do 8766 i=ir1,ir2
	ip=ip+1
	nock=nock+1
	ipb=(j-1)*mbig+i
	savtop(ip)=bigt(ipb)
	savbot(ip)=bigb(ipb)
	do 8766 k=1,kz
        savstm(ip,k)=bigs(ipb,k)
 8766	savvrt(ip,k)=bigv(ipb,k)
c
	if(nock.ne.mn0)then
	print *,'wrong domain size in clinic!!!'
	stop
	endif
c
c        print *,'done background read in rdback, icalc=',icalc
c
 9876   continue
c
c
	return
        end
c
	subroutine delsq(npnt,right,nptb,alhs)
c
c
c
c
c --- this is a package of routines for a finite-element solution of
c --- the quasi-geostrophic equations,as described in "a baroclinic
c --- quasi-geostrophic open ocean model" (1981) by r.n. miller,
c --- a.r. robinson, and d.b. haidvogel of harvard.
c
c --- bottom topography (for slopes less than o(r0)) and top vertical
c --- velocity have been implemented (10/15/82).
c
c --- a driving package is required to provide initial and boundary
c --- values of the streamfunction, vorticity, and top and bottom
c --- density variations; as well as the topography. this package
c --- is called infld.
c
c --- a variety of options are available:
c
c        ifdif = 1/0     finite difference/collocation in depth
c        ifpert = 1/0    calculate perturbation fields (yes/no)
c        ifrst = 1/0     restart at time rdt (yes/no)
c        iftop = 1/0     use top density information (yes/no)
c        ifbot = 1/0     use bottom density information (yes/no)
c        iftvv = 1/0     use top vertical velocity (yes/no)
c        ifrel = 1/0     use bottom relief (yes/no)
c        ifsbl = 1/0     use surface boundary layer
c        ifeva = 1/0     execute eva preprocessor
c
c
           include 'prob2.pert'
c
	parameter(m0m1=m0-1,mmat=m0m1*n0)
c
      character*4 restid,runid
      character*70 idstrg
c
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
	common/solad/uan(n0),ubn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/count/icnt,itcnt,ibcnt
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),topo(mn0),
     2     topd(mn0),botd(mn0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
c
	common/grad/gomeg(mn0,k0),gphi(mn0,k0),ggam(mn0),
     1     gdel(mn0),dm(mn0,k0),dp(mn0,k0),dg(mn0),dd(mn0)
     2    ,cgf
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/steps/sumt,sumb,cost
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/ctest/ctol
c
	common/swoi/iter
c
	common/toff/tadj0
c
	common/istar/icall
c
	common/perst/psf(mn0,k0),pvr(mn0,k0)
c
	common/egl2/pl,akl,alls,allv
c
	common/amiv/amx(m0-1,m0-1),amy(n0,n0)
c
	common/lan1/maxt,ifts
c
c
c
c
c
      dimension therm(m20,k0),dxxpsi(m20),dyypsi(m20)
C************
	 real alhs(nptb),right(npnt),atmp(mn0,k0)
C************
c
c
      equivalence (nprnt, nplot)
c
c.......copy right into of interior of vort......
c
        do k=1,kz
        do i=1,msq
        atmp(i,k)=0.0
        enddo
        enddo
c
	ipp=0
	do 175 k=1,kz
	do 175 j=2,n-1
	do 175 i=2,m-1
	ip=(j-1)*m+i
	ipp=ipp+1
 175    atmp(ip,k)=right(ipp)
c
	if(ipp.ne.npnt)then
	print *,'mismatch between ipp and npnt!!!'
	stop
	endif
c
c.......take delsq of input........
c
       call thermvt(atmp,therm,kz,msq)
       do 5005 k=1,kz
       call rel2nd(atmp(1,k),dxxpsi,dyypsi,m,n,hy)
       call filter(dxxpsi,m,n,4)
       call filter(dyypsi,m,n,4)
       do 4005 ip=1,msq
 4005  atmp(ip,k)=(dxxpsi(ip)+dyypsi(ip)+therm(ip,k))
 5005   continue
c
	ipp=0
	do k=1,kz
	do j=1,n
	do i=1,m
	ip=(j-1)*m+i
	ipp=ipp+1
        alhs(ipp)=atmp(ip,k)
	enddo
	enddo
	enddo
c
	if(ipp.ne.nptb)then
	print *,'mismatch 2 between ipp and npnt!!!'
	stop
	endif
c
	return
	end
	subroutine adelsq(npnt,right,alhs)
c
c
c
c
c --- this is a package of routines for a finite-element solution of
c --- the quasi-geostrophic equations,as described in "a baroclinic
c --- quasi-geostrophic open ocean model" (1981) by r.n. miller,
c --- a.r. robinson, and d.b. haidvogel of harvard.
c
c --- bottom topography (for slopes less than o(r0)) and top vertical
c --- velocity have been implemented (10/15/82).
c
c --- a driving package is required to provide initial and boundary
c --- values of the streamfunction, vorticity, and top and bottom
c --- density variations; as well as the topography. this package
c --- is called infld.
c
c --- a variety of options are available:
c
c        ifdif = 1/0     finite difference/collocation in depth
c        ifpert = 1/0    calculate perturbation fields (yes/no)
c        ifrst = 1/0     restart at time rdt (yes/no)
c        iftop = 1/0     use top density information (yes/no)
c        ifbot = 1/0     use bottom density information (yes/no)
c        iftvv = 1/0     use top vertical velocity (yes/no)
c        ifrel = 1/0     use bottom relief (yes/no)
c        ifsbl = 1/0     use surface boundary layer
c        ifeva = 1/0     execute eva preprocessor
c
c
           include 'prob2.pert'
c
	parameter(m0m1=m0-1,mmat=m0m1*n0)
c
      character*4 restid,runid
      character*70 idstrg
c
c
	common/bgsv/bstm(mn0,k0),bvrt(mn0,k0)
c
      common/scrat/r(mn0),r1(mn0),r2(mn0)
      common/fields/vort(mn0,k0),stream(mn0,k0),prvjac(mn0,k0),
     +     topden(mn0),botden(mn0),relief(mn0),prvtop(mn0),
     +     prvbot(mn0),topfct(k0),botfct(k0),pvort(mn0)
      common/parm/kz,m,mm1,mm2,mp1,msq,tim,tmax,tstart,hy,rea,dt,
     +   xbasin,yy(mx0),icalc,kzp1,kzm1,iftop,ifbot,iftvv,ifrel,theta,
     +   ifsbl,ifeva
      common/solver/vam(m0),vbm(m0),van(n0),vbn(n0)
	common/solad/uan(n0),ubn(n0)
      common/edges/icorn(4),ncorn(4),mast(4),nast(4)
      common/philtr/iltord,iltfrq,iltcnt
      common/aphil/iadn,iadf,iadc,itn,itf,itz
      common/count/icnt,itcnt,ibcnt
      common/threed/hz(k0),scsig,sigz(k01),eigval(k0),ainv(k0,k0),
     +      amat(k0,k0)
      common/vert/buoysur,tsurf,xc(k01),salsur,eta(k01)
      common/user/pram(16)
      common/dynam/alpha,beta,rf
      common/unit/iin,iout,idff,idif,iorr,iowr,iplt,ifndd
      common/cnn/n,nm1,nm2,np1
      common/extpar/rlat0,rlng0,t0,v0,dhor,ht,time0,
     +                r0,f0star
      common/inform/ifdiff,ifpert,titl(20,2),runid(6),restid(9),
     +                                                     idstrg
      common/control/ddt(4),adt(4),nalev(4),ialev(4,k02),iprnt(k02),
     !      lprnt(k0)
c
	common/qgm/strm1(mn0,k0),strm2(mn0,k0),vort1(mn0,k0),
     1     vort2(mn0,k0),topo(mn0),
     2     topd(mn0),botd(mn0)
c
c
	common/adjv/omega1(mn0,k0),omega2(mn0,k0),phi1(mn0,k0),
     1     phi2(mn0,k0),prvjr(mn0,k0),apvt(mn0),
     2     delta1(mn0),delta2(mn0),gamma1(mn0),
     3     gamma2(mn0),prvjd(mn0),prvjb(mn0),icmax,rad(mn0,k0),
     4     zhat(mn0,k0),s2(mn0,k0),fys(mx0,4,k0),ss(mn0,k0)
c
c
	common/obs/obstm(nosp),obvort(nosp),iobtim(nosp),
     1     ipobs(nosp),kobs(nosp),ipobv(nosp),kobv(nosp)
c
c
	common/grad/gomeg(mn0,k0),gphi(mn0,k0),ggam(mn0),
     1     gdel(mn0),dm(mn0,k0),dp(mn0,k0),dg(mn0),dd(mn0)
     2    ,cgf
c
        common/qgbc/inflow(ibp),ipos(ibp)
c
	common/pert/stf(nosp),vtf(nosp)
c
	common/cint/st0(mn0,k0),vt0(mn0,k0),top0(mn0),bot0(mn0)
c
        common/ovar/sigs,sigv,ndats,ndatv
c
	common/steps/sumt,sumb,cost
c
c
	common/pass/ipass,icg,nsk,iflag,idis,itest,itol,numi
c
	common/ctest/ctol
c
	common/swoi/iter
c
	common/toff/tadj0
c
	common/istar/icall
c
	common/perst/psf(mn0,k0),pvr(mn0,k0)
c
	common/egl2/pl,akl,alls,allv
c
	common/amiv/amx(m0-1,m0-1),amy(n0,n0)
c
	common/lan1/maxt,ifts
c
      dimension therm(m20,k0),dxxpsi(m20),dyypsi(m20)
C************
	 real alhs(npnt),right(npnt),atmp(mn0,k0)
C************
c
c
      equivalence (nprnt, nplot)
c
	ipp=0
	do k=1,kz
	do j=1,n
	do i=1,m
	ip=(j-1)*m+i
	ipp=ipp+1
        atmp(ip,k)=alhs(ipp)
	enddo
	enddo
	enddo
c
	if(ipp.ne.npnt)then
	print *,'mismatch 2 between ipp and npnt!!!'
	stop
	endif
c
c.......take adelsq of input........
c
       do 5005 k=1,kz
       do ip=1,msq
       dxxpsi(ip)=atmp(ip,k)
       dyypsi(ip)=atmp(ip,k)
       therm(ip,k)=atmp(ip,k)
       enddo
       call afiltr(dxxpsi,m,n,4)
       call afiltr(dyypsi,m,n,4)
       call adrel2nd(atmp(1,k),dxxpsi,dyypsi,m,n,hy)
 5005   continue
c
       call adtherm(atmp,therm,kz,msq)
c
	ipp=0
	do 175 k=1,kz
	do 175 j=1,n
	do 175 i=1,m
	ip=(j-1)*m+i
	ipp=ipp+1
 175    right(ipp)=atmp(ip,k)
c
	if(ipp.ne.npnt)then
	print *,'mismatch between ipp and npnt!!!'
	stop
	endif
c
	return
	end
c
