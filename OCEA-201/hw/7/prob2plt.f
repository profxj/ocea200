         program srnplt
c
         parameter(nvmax=50)
c
c         parameter (m0=41,n0=25,mn0=m0*n0,m20=mn0,mx0=m0+n0,k0=1,
         parameter (m0=41,n0=25,mn0=m0*n0,m20=mn0,mx0=m0+n0,k0=2,
     +           k01=k0+1,k02=k0+2,mnwork=6*n0+8*m0
     +           ,nosp=1,ibp=2*m0+2*(n0-2))
c
         parameter(mm0=m0)
c
         parameter(nval=7)
c
	 real ritzr(nvmax),ritzi(nvmax)
c
         real si(mm0,n0,k0,nvmax),vi(mm0,n0,k0,nvmax)
c
         real sf(mm0,n0,k0,nvmax),vf(mm0,n0,k0,nvmax)
c
	 real top(mm0,n0,nvmax)
c
	 real bstm(mm0,n0,k0),bvrt(mm0,n0,k0)
c
         real ap(mm0,n0),asve(mm0,n0),cval(1)
c
	 real uprof(n0),gamma(n0)
c
	 real xarr(n0),yarr(n0),xlarr(2),ylarr(2)
c
         real cvalp(nval),cvalm(nval)
c
	real spv(2),cfill(3),grylev(3),cfill1(3),grylev1(3)
c
         complex cfreq
c
	character ccrap*2,critz*10,char1*4,char2*3
c
	character camuq*10,camue*10,camup*10,cnum*1
	character cper*10,cef*10
c
      COMMON/CONPAR/ISPEC,IOFFPP,SPVALL,ILEGG,ILABB,NHII,NDECCN,NLBLL,
     .              LSCAL,LDASH
      DATA ISPEC,IOFFPP,SPVALL,ILEGG,ILABB,NHII,NDECCN,NLBLL,
     .   LSCAL,LDASH
     1  /   1 ,     0,    -999., 0,    0,    0,    1,     3, 0, 0/
c
c    growth time in days.
c
        tau=1.0
c
        ispec=1
        nhii=-1
        ilegg=1
c
	fact1=1.5/float(k0)
c
	fact=3./4.
        xlen=4./fact
        ylen=1.5/fact
c        ylen=float(n0)*xlen/float(m0)
c
        call newdev('prob2.ps',7)
	call psinit(.true.)
        call factor(fact1*fact)
c
         open(unit=10,file='ftnmplt.dat',form='unformatted')
c
	read(10)uprof
	read(10)gamma
c
	read(10)nev,u0,beta
c
	print *,'nev=',nev
c
	if(nev.gt.nvmax)then
	  print *,'nev exceeds nvmax!!!'
	  stop
	endif
c
	read(10)bstm
	read(10)bvrt
c
	 do nv=1,nev
	   read(10)ritzr(nv),ritzi(nv),rcrap
           read(10)(((si(i,j,k,nv),i=1,mm0),j=1,n0),k=1,k0)
           read(10)(((vi(i,j,k,nv),i=1,mm0),j=1,n0),k=1,k0)
           read(10)(((sf(i,j,k,nv),i=1,mm0),j=1,n0),k=1,k0)
           read(10)(((vf(i,j,k,nv),i=1,mm0),j=1,n0),k=1,k0)
           print *,'read s, nv=',nv
	 enddo
c
	xstrt=1./fact
	ystrt=10.5/fact
c
        call plot(xstrt,ystrt,-3)
c
         nmode=0
c
         do iplt=1,nev,2
c
         nmode=nmode+1
c
	 nplt=0
c
	 do klv=1,k0
c
	 do nfld=1,5
c
	 nplt=nplt+1
c
         amp=200.0
         do j=1,n0
         do i=1,mm0
            if(nfld.eq.1)ap(i,j)=bstm(i,j,klv)
            if(nfld.eq.2)ap(i,j)=amp*si(i,j,klv,iplt)
            if(nfld.eq.3)ap(i,j)=bstm(i,j,klv)+amp*si(i,j,klv,iplt)
            if(nfld.eq.4)ap(i,j)=amp*vi(i,j,klv,iplt)
            if(nfld.eq.5)ap(i,j)=bvrt(i,j,klv)+amp*vi(i,j,klv,iplt)
         enddo
         enddo
c
	amuq=sqrt(ritzr(iplt)**2+ritzi(iplt)**2)
        cfreq=clog(cmplx(ritzr(iplt),ritzi(iplt)))/cmplx(tau,0.0)
c
        efold=1./real(cfreq)
        eper=1./aimag(cfreq)
c
	write(camuq,'(1p,e10.3)')amuq
	write(camue,'(1p,e10.3)')amuq
	write(camup,'(1p,e10.3)')amuq
	write(cper,'(1p,e10.3)')eper
	write(cef,'(1p,e10.3)')efold
c
	print *,'iplt=',iplt
	print *,'amuq=',amuq
	print *,'amue=',amue
	print *,'amup=',amup
c
c........determine the largest values of ap....
c
         big=-1.e20
         do j=1,n0
         do i=1,mm0
         if(abs(ap(i,j)).gt.big)big=abs(ap(i,j))
         enddo
         enddo
c
c........determine the contour intervals.....
c
         dc=big/float(nval)
         do i=1,nval
         cvalp(i)=float(i)*dc
         enddo
         do i=nval,1,-1
         cvalm(i)=-1.*float(i)*dc
         enddo
c
        cfill(1)=-1.e20
        cfill(2)=0.0
        cfill(3)=1.e20
        grylev(1)=0.75
        grylev(2)=0.75
        grylev(3)=1.0
c
        cfill1(1)=cvalm(nval)
        cfill1(2)=cvalm(1)
        cfill1(3)=1.e20
        grylev1(1)=0.75
        grylev1(2)=0.75
        grylev1(3)=1.0
c
        call confill(ap,mm0,mm0,n0,xlen,ylen,
     &    cfill1,grylev1,3,ioffpp,spvall)
c
        call conrec(ap,mm0,mm0,n0,xlen,ylen,cvalp,nval)
        call conrec(ap,mm0,mm0,n0,xlen,ylen,cvalm,nval)
c
        call border(xlen,ylen,1111,1111,1,1,1,1)
c
	ch1=0.15
	if(nfld.eq.1.and.klv.eq.1)then
	write(critz,'(1p,e10.3)')u0
	call keksym(0.1*xlen,1.05*ylen,ch1,2HU=,0.,2,0)
	call keksymc(999.,999.,ch1,critz,0.,10,0)
	write(critz,'(1p,e10.3)')beta
	call grksym(0.6*xlen,1.05*ylen,ch1,26,0.,1,0)
	call keksym(999.,999.,ch1,1H=,0.,1,0)
	call keksymc(999.,999.,ch1,critz,0.,10,0)
	endif
	if(nfld.eq.5.and.klv.eq.1)then
	dyl=0.4
	call keksym(0.5*xlen,-1.,ch1,7HPeriod=,0.,7,0)
	call keksymc(999.,999.,ch1,cper,0.,10,0)
	call keksym(0.5*xlen,-1.-dyl,ch1,15He-folding time=,0.,15,0)
	call keksymc(999.,999.,ch1,cef,0.,10,0)
	endif
c
        if(klv.eq.1.and.nfld.eq.1)then
        call keksym(xlen+0.25,1.4*ylen,0.3,13HNormal mode #,0.,13,1)
	write(cnum,'(i1)')nmode
	call keksymc(999.,999.,0.3,cnum,0.,1,0)
	endif
c
        if(klv.eq.1.and.nfld.eq.1)
     &   call keksym(0.5*xlen,1.2*ylen,0.2,7HLayer 1,0.,7,1)
        if(klv.eq.2.and.nfld.eq.1)
     &   call keksym(0.5*xlen,1.2*ylen,0.2,7HLayer 2,0.,7,1)
c
	if(klv.eq.1)then
	ch2=0.2
	if(nfld.eq.1)then
          call grksym(-0.8,0.5*ylen,ch2,23,0.,1,0)
          call keksym(0.02*xlen,0.8*ylen,ch2,1HA,0.,1,0)
        endif
	if(nfld.gt.1)then
	  if(nfld.eq.2)then
            call grksym(-0.8,0.5*ylen,ch2,47,0.,1,0)
            call keksym(0.02*xlen,0.8*ylen,ch2,1HB,0.,1,0)
          endif
	  if(nfld.eq.3)then
            call grksym(-0.8,0.5*ylen,ch2,23,0.,1,0)
            call keksym(999.,999.,ch2,1H+,0.,1,0)
            call grksym(999.,999.,ch2,47,0.,1,0)
            call keksym(0.02*xlen,0.8*ylen,ch2,1HC,0.,1,0)
          endif
	  if(nfld.eq.4)then
             call grksym(-0.8,0.5*ylen,ch2,30,0.,1,0)
             call keksym(0.02*xlen,0.8*ylen,ch2,1HD,0.,1,0)
          endif
	  if(nfld.eq.5)then
            call grksym(-0.8,0.5*ylen,ch2,16,0.,1,0)
            call keksym(999.,999.,ch2,1H+,0.,1,0)
            call grksym(999.,999.,ch2,30,0.,1,0)
            call keksym(0.02*xlen,0.8*ylen,ch2,1HE,0.,1,0)
          endif
	endif
	endif
	if(klv.eq.2)then
	if(nfld.eq.1)then
          call keksym(0.02*xlen,0.8*ylen,ch2,1HF,0.,1,0)
        endif
	if(nfld.gt.1)then
	  if(nfld.eq.2)then
            call keksym(0.02*xlen,0.8*ylen,ch2,1HG,0.,1,0)
          endif
	  if(nfld.eq.3)then
            call keksym(0.02*xlen,0.8*ylen,ch2,1HH,0.,1,0)
          endif
	  if(nfld.eq.4)then
             call keksym(0.02*xlen,0.8*ylen,ch2,1HI,0.,1,0)
          endif
	  if(nfld.eq.5)then
            call keksym(0.02*xlen,0.8*ylen,ch2,1HJ,0.,1,0)
          endif
	endif
	endif
c
         if(mod(nplt,5).ne.0)call plot(0.,-1.*(ylen+0.5),-3)
         if(mod(nplt,5).eq.0)call plot(xlen+0.5,4.*(ylen+0.5),-3)
c
         if(nplt.eq.5*k0)then
               call chopit(0.,0.)
               call plot(xstrt,ystrt,-3)
         endif
c
         enddo
c
         enddo
c
         enddo
c
         call plot(-1.*xstrt+2.,-1.*ystrt+5.,-3)
c
	 xlen=5.
	 ylen=5.
c
        call border(xlen,ylen,-1110,1111,25,1,10,1)
c
	ymax=5.
	v0=0.4
	do i=1,n0
	  xarr(i)=float(i-1)*xlen/float(n0-1)
	  yarr(i)=v0*uprof(i)*ylen/ymax+0.5*ylen
	enddo
c
	ch3=0.12
	do i=1,11
	  x=-0.3
	  y=float(i-1)*ylen/10.
	  write(char1,'(f4.1)')-ymax/2.+float(i-1)*ymax/10.
	  call keksymc(x,y-0.5*ch3,ch3,char1,0.,4,2)
	enddo
c
	call dshcrv(xarr,yarr,n0,8,0.02)
c
c
        ymax=100.
        do i=1,n0
          yarr(i)=gamma(i)*ylen/ymax+0.5*ylen
        enddo
c
        print *,'gamma=',gamma
        print *,'yarr=',yarr
c
        do i=1,11
          x=xlen+0.2
          y=float(i-1)*ylen/10.
          int=-ymax/2
          int=int+(i-1)*ymax/10
          write(char2,'(i3)')int
          call keksymc(x,y-0.5*ch3,ch3,char2,0.,3,0)
        enddo
c
        call drwcrv(xarr,yarr,n0,0.01,.false.)
c
        call keksym(0.5*xlen,-0.4,1.5*ch3,3Hy=0,0.,3,1)
c
        call keksym(-1.4,0.5*ylen,1.5*ch3,4HU(y),0.,4,0)
c
        call grksym(xlen+0.9,0.5*ylen,1.5*ch3,27,0.,1,0)
        call keksym(999.,999.,1.5*ch3,3H(y),0.,3,0)
c
        ch4=1.5*ch3
        xlarr(1)=0.8*xlen
        xlarr(2)=xlarr(1)+0.5
        ylarr(1)=-1.
        ylarr(2)=ylarr(1)
        call drwcrv(xlarr,ylarr,2,0.01,.false.)
        call grksym(xlarr(2)+0.1,ylarr(1)-0.5*ch4,ch4,27,0.,1,0)
c
        xlarr(1)=0.8*xlen
        xlarr(2)=xlarr(1)+0.5
        ylarr(1)=-1.5
        ylarr(2)=ylarr(1)
        call dshcrv(xlarr,ylarr,2,8,0.02)
        call keksym(xlarr(2)+0.1,ylarr(1)-0.5*ch4,ch4,1HU,0.,1,0)
c
c
        call plotnd
c
         stop
         end
