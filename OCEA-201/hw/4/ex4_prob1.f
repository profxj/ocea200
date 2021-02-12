        PROGRAM INDPAC
C
C
        PARAMETER(IM=51)
        PARAMETER(JM=51)
        PARAMETER(MMI=im-1)
        PARAMETER(MMJ=jm-1)
        CHARACTER*85 CDUM
        CHARACTER*4 CFIL,CRAP,CCCC
C
       REAL GDH1X(IM,JM),GDH2X(IM,JM),GDH1Y(IM,JM),
     Z GDH2Y(IM,JM),US1(IM,JM),US2(IM,JM),VS1(IM,JM),VS2(IM,JM),
     1    ECOM(JM)
        DIMENSION U1(IM,JM,3),U2(IM,JM,3),V1(IM,JM,3),V2(IM,JM,3),
     X H1(IM,JM,3),H2(IM,JM,3)
     2  ,IHM(IM,JM),IUM(IM,JM),IVM(IM,JM),WXS(IM,JM),WYS(IM,JM)
     6  ,U1MEAN(MMI,MMJ),RUS1(IM,JM),RVS1(IM,JM),V1MEAN(MMI,MMJ)
     7    ,TAUX(84,30),TAUY(84,30),IDAY(12),DUM(IM,JM)
        DIMENSION TPH(IM,JM),POTE(7),AKINE(7),ALAT1(7),ALAT2(7),
     1   ALONG1(7),ALONG2(7),IX1(7),IX2(7),IY1(7),IY2(7),ITM(IM,JM),
     2   STRKX(6),STRKY(6),TRGT(6),TRACKW(6),RMSTK(6),RMSREG(7),
     3   DAY(12),DOBS(IM,JM),TISO(6),IOBS(IM,JM),NOTAU(400),
     4   BY(IM,JM),VECW(IM,IM),DIFCORD(IM,JM),RHO(IM,JM),
     5   NOBS(7),NNOBS(6),DIFH(IM,JM),DINM(12),BX(IM,JM)
     6    ,WEIGHT(IM,JM),VARREG(7),VARTRK(6),IPOS(400),JPOS(400)
     7     ,ERR(IM,JM),ERRREG(7),ERRTK(6)
     8   ,TH1(IM,JM),TH2(IM,JM),TU1(IM,JM),TU2(IM,JM),TV1(IM,JM),
     9   TV2(IM,JM),DEX1(IM),STUA(IM),STUB(IM),STVA(IM),STVB(IM)
     9   ,DEY1(IM),IUP(IM,JM)
c
	real dex2(im),dey2(im)
C
         DIMENSION HKM(IM,JM),HK(IM,JM),HKP(IM,JM),UKM(IM,JM),
     1    UK(IM,JM),UKP(IM,JM),VKM(IM,JM),VK(IM,JM),VKP(IM,JM)
     2   ,AMEAN(IM,JM),WXT(169,61),WYT(169,61),TIX(94,58),
     3    TIY(94,58),SUPOT(6),SUKIN(6)
c
         DIMENSION HKM2(IM,JM),HK2(IM,JM),HKP2(IM,JM),UKM2(IM,JM),
     1    UK2(IM,JM),UKP2(IM,JM),VKM2(IM,JM),VK2(IM,JM),VKP2(IM,JM)
c
	real h1p(im,jm),u1p(im,jm),v1p(im,jm)
	real h2p(im,jm),u2p(im,jm),v2p(im,jm)
c
	real wx0(im,jm),wy0(im,jm)
C
C
	common/vars/beta,dx,dy,g12,g23,f0,
     &         fh(jm),fv(jm),wx(im,jm),wy(im,jm)
c
	common/vars1/je,w0
c
c.......open plot file......
c
	call newdev('prob1.ps',8)
	call psinit(.true.)
c
	open(unit=15,file='ex4_prob1.input',status='old')
c
c	open(unit=8,file='plottest.out',form='unformatted')
c	open(unit=27,file='restart.dat',form='unformatted')
c
c.......read input parameters.......
c
	read(15,*)itstop
c	read(15,*)iftype
c	read(15,*)phi0
c	read(15,*)ipert
	read(15,*)ipltm
	read(15,*)w0
	read(15,*)asubh
c	read(15,*)ifnonl
	read(15,*)r1
c
	itstop=12*itstop
	ipltm=12*ipltm
	iftype=2
	ifnonl=1
	phi0=40.
	ipert=2
	w0=10.*w0
C
        IMAX=IM
        JMAX=JM
        IMJM=IM*JM
        IMIM=IM*IM
        IR=IMAX-2
C
C
       KTAU=0
C.......ISTART=START UP TIMESTEP FOR NEXT RUN...
        ISTART=16*30*12
C.......JJSTOP=STOP TIMESTEP FOR THIS RUN.....
        JJSTOP=16*30*itstop
C......IDSTR=START UP DAY....
        IDSTR=15
C......MTHSTR=START UP MONTH.......
        MTHSTR=1
C........IYR=START UP YEAR........
        IYR=77
        MONTH1=MTHSTR
        IYEAR1=IYR
        IMTH=MTHSTR
        JMTH=MTHSTR
        ITAU=KTAU
        ISTOP=0
          LNUP=0
        CFIL='FILE'
        JSUP=0
C
C.......DEFINE MEAN ISOTHERM DEPTH IN METRES....
C
        AMLD=200
C
C......CLEAR NVAL= NO. OF TIMESTEPS USED TO DETERMINE MEAN DEPTH OF
C      LAYER AT EACH GRIDPOINT...
C
       NVAL=0
C
       DO 22 I=1,IMAX
       DO 22 J=1,JMAX
        DO 23 KK=1,3
        V1(I,J,KK)=0.
        V2(I,J,KK)=0.
        U1(I,J,KK)=0.
        U2(I,J,KK)=0.
        H1(I,J,KK)=0.
        H2(I,J,KK)=0.
 23      CONTINUE
 22      CONTINUE
       KM=1
       K=2
       KP=3
       DO 24 I=1,IM
       DO 24 J=1,JM
        AMEAN(I,J)=0.0
        WEIGHT(I,J)=0.0
        ERR(I,J)=0.0
        DIFH(I,J)=0.0
        TPH(I,J)=0.0
        ITM(I,J)=0
        DOBS(I,J)=0.0
        IOBS(I,J)=0
        IUP(I,J)=0
        BX(I,J)=0.0
        BY(I,J)=0.0
        DIFCORD(I,J)=0.0
        RHO(I,J)=1.0
       US1(I,J)=0
      US2(I,J)=0
      VS1(I,J)=0
      VS2(I,J)=0
      GDH1X(I,J)=0.
      GDH2X(I,J)=0.
      GDH1Y(I,J)=0.
      GDH2Y(I,J)=0.
        WX(I,J)=0.0
        WY(I,J)=0.0
      WXS(I,J)=0.0
        WYS(I,J)=0.0
        DUM(I,J)=0.0
 24    CONTINUE
        CALL SETRAR(WXT,10309,0.0)
        CALL SETRAR(WYT,10309,0.0)
        CALL SETRAR(TIX,5452,0.0)
        CALL SETRAR(TIY,5452,0.0)
C        CALL SETRAR(NOTAU,400,0)
C        CALL SETRAR(IPOS,400,1)
C        CALL SETRAR(JPOS,400,1)
C
        DO 759 J=1,30
        DO 759 I=1,84
        TAUX(I,J)=0.0
 759    TAUY(I,J)=0.0
      ITAPE=100
 98    FORMAT(I5)
 97    FORMAT(1X,' KTAPE=',I5,' KTAU=',I5)
C*******************************
C......SET UP G12 SUCH THAT C=280 CM/S WITH ANY HH1.......
C
       G12=7.84
c
c	g12=2.
C*******************************
       G23=7.84/2.
c	g23=2.
c
        SCV=1.0E3
        SCH=1.0E2
        DX=1.E7/0.9
        DY=1.E7/0.9
C*******DEEPEN LAYER IF NECESSARY********
       HH1=10000.
       HH2=10000.
       BETA=2.E-13
       DO 765 I=1,IM
       DEX1(I)=HH1
       DEY1(I)=HH1
       STUA(I)=0.0
       STUB(I)=0.0
       STVA(I)=0.0
       STVB(I)=0.0
 765   CONTINUE
       TDX=2.*DX
       TDY=2.*DY
       DX2=DX*DX
       DY2=DY*DY
        RDX=1.0/DX
        RDY=1.0/DY
        RDX2=1.0/DX2
        RDY2=1.0/DY2
        DXDY=DX*DY
        RDXDY=1.0/DXDY
c        DAMP=1.0E4
        DAMP=0.
      DD=0.0
       IMAXM=IMAX-1
       JMAXM=JMAX-1
       ISMTH=31
       IPR=365
      DVT=r1
      DVB=g12*dvt/g23
       DT=86400./16.
       TDT=2.*DT
C.....DEFINE VECTOR WORK SPACE
        CALL SETRAR(VECW,IMIM,0.0)
        DO 341 I=1,IM
        DO 341 II=1,IM
  341   VECW(II,I)=(I-II)
C
C******************************
        JE=(jm-1)/2+1
        IE=(im-1)/2+1
C********************************
c
c......f-plane...
c
       if(iftype.eq.1)then
c
	f0=2.*7.292e-5*sin(phi0*2.*4.*atan(1.0)/360.)
c
       DO J=1,JMAX
        FH(J)=f0
        FV(J)=f0
       ENDDO
c
       endif
c
c......mid-latitude beta-plane...
c
       if(iftype.eq.2)then
c
	f0=2.*7.292e-5*sin(phi0*2.*4.*atan(1.0)/360.)
c
       DO 15 J=1,JMAX
        FH(J)=f0+(J-JE-0.5)*DY*BETA
        FV(J)=f0+(J-JE)*DY*BETA
 15     CONTINUE
c
       endif
c
c......equatorial beta-plane...
c
       if(iftype.eq.3)then
c
       DO J=1,JMAX
        FH(J)=(J-JE-0.5)*DY*BETA
        FV(J)=(J-JE)*DY*BETA
       ENDDO
c
       endif
c
	do j=2,jmaxm
c	  ecom(j)=2.e8
	  ecom(j)=asubh*1.e4
	enddo
c
C********READ IN GEOMETRY*********
c
	do j=1,jm
	do i=1,im
	  ihm(i,j)=1
	enddo
	enddo
c
	do j=1,jm
	  ihm(1,j)=0
	  ihm(im,j)=0
	enddo
c
	do i=1,im
	  ihm(i,1)=0
	  ihm(i,jm)=0
	enddo
C
        DO 43 J=1,JM
       DO 43 I=1,IMAXM
 43     IUM(I,J)=IHM(I,J)*IHM(I+1,J)
       DO 46 J=1,JMAXM
       DO 46 I=1,IM
 46     IVM(I,J)=IHM(I,J)*IHM(I,J+1)
C
c
c########################################
c
c	skip read!!!
c
c########################################
c
c.......initialise the model with an initial H1 field....
c
	if(ipert.eq.1)then
	do j=je-5,je+5	
	do i=ie-5,ie+5	
	  h1(i,j,1)=100.0 
	  h1(i,j,2)=100.0 
	enddo
	enddo
	endif
c
	if(ipert.eq.2)then
	  kwramp=16*30*24
	  do j=1,jm
	   do i=1,im
	     wx0(i,j)=w0*sin(2.*4.*atan(1.0)*float(j-je)/float(2*jm))
c	     wx0(i,j)=-0.1
	     wy0(i,j)=0.
	   enddo
	  enddo
	endif
c
	print *,'wx0=',(wx0(3,j),j=1,jm)
c
C
 500    CONTINUE
c
       KTAU=KTAU+1
c
c	print *,'ktau=',ktau
c
       TDT=2.*DT
       IF(KTAU.EQ.1)TDT=DT
c
c......ramp up the strength of the wind......
c
	if(ipert.eq.2)then
c
	wfac=1.0
	if(ktau.le.kwramp)wfac=float(ktau)/float(kwramp)
	do j=1,jm
	  do i=1,im
	    wx(i,j)=wfac*wx0(i,j)
	    wy(i,j)=wfac*wy0(i,j)
	  enddo
	enddo
c
	endif
C
c
c#######################################
c
c       skip the wind read and interpolation
c
c######################################
c
c	if(ktau.le.16*30*2)then
c	do j=27,37
c	do i=154,164
c	  wx(i,j)=0.1
c	  wy(i,j)=0.
c	enddo
c	enddo
c	endif
c
c	if(ktau.gt.16*30*2)then
c	do j=27,37
c	do i=154,164
c	  wx(i,j)=0.
c	  wy(i,j)=0.
c	enddo
c	enddo
c	endif
c
C
C     *****
C     IF YOU CHANGE FORCING AMPLITUDE, CHANGE THE SCALING
C
      SCV=1.0E3
      SCH=1.0E2
       DO 20 J=2,JMAXM
       DO 20 I=2,IMAXM
       GDH2X(I,J)=-G23*(H2(I+1,J,K)-H2(I,J,K))/DX
       GDH2Y(I,J)=-G23*(H2(I,J+1,K)-H2(I,J,K))/DY
 20     CONTINUE
c
       DO 21 J=2,JMAXM
       DO 21 I=2,IMAXM
       GDH1X(I,J)=GDH2X(I,J)-G12*(H1(I+1,J,K)-H1(I,J,K))/DX
       GDH1Y(I,J)=GDH2Y(I,J)-G12*(H1(I,J+1,K)-H1(I,J,K))/DY
 21     CONTINUE
       DO 128 J=2,JMAXM
      DO 28 I=2,IMAXM
      H1AX=HH1-(H1(I,J,K)+H1(I+1,J,K))/2.
      H2AX=HH2+(H1(I,J,K)+H1(I+1,J,K))/2.-(H2(I,J,K)+H2(I+1,J,K))/2.
      H1AY=HH1-(H1(I,J,K)+H1(I,J+1,K))/2.
      H2AY=HH2+(H1(I,J,K)+H1(I,J+1,K))/2.-(H2(I,J,K)+H2(I,J+1,K))/2.
      US1(I,J)=U1(I,J,K)/H1AX
      VS1(I,J)=V1(I,J,K)/H1AY
      US2(I,J)=U2(I,J,K)/H2AX
      VS2(I,J)=V2(I,J,K)/H2AY
 399     FORMAT(1X,2I5,8E12.3)
 28      CONTINUE
 128    CONTINUE
C
C
C........DEFINE TEMPORARY ARRAYS....
C
         DO 6789 J=1,JM
         DO 6789 I=1,IM
         HKM(I,J)=H1(I,J,KM)
         HK(I,J)=H1(I,J,K)
         HKP(I,J)=H1(I,J,KP)
         UKM(I,J)=U1(I,J,KM)
         UK(I,J)=U1(I,J,K)
         UKP(I,J)=U1(I,J,KP)
         VKM(I,J)=V1(I,J,KM)
         VK(I,J)=V1(I,J,K)
         VKP(I,J)=V1(I,J,KP)
         HKM2(I,J)=H2(I,J,KM)
         HK2(I,J)=H2(I,J,K)
         HKP2(I,J)=H2(I,J,KP)
         UKM2(I,J)=U2(I,J,KM)
         UK2(I,J)=U2(I,J,K)
         UKP2(I,J)=U2(I,J,KP)
         VKM2(I,J)=V2(I,J,KM)
         VK2(I,J)=V2(I,J,K)
         VKP2(I,J)=V2(I,J,KP)
 6789    CONTINUE
C
C
       DO 30 J=2,JMAXM
CDIR$ IVDEP
       DO 30 I=2,IMAXM
      DEX1(I)=HH1-(HK(I,J)+HK(I+1,J))/2.
      DEX2(I)=HH2-(HK2(I,J)+HK2(I+1,J)-HK(I,J)-HK(I+1,J))/2.
      STUA(I)=DVT*(UKM(I,J) - UKM2(I,J))
      STUB(I)=DVB*UKM2(I,J)
       UKP(I,J)=UKM(I,J)+TDT*IUM(I,J)*(
     1 0.25*( FV(J)*(VK(I,J)+VK(I+1,J))
     1  +FV(J-1)*(VK(I,J-1)+VK(I+1,J-1))) -GDH1X(I,J)*DEX1(I)
     2 + WX(I,J) - STUA(I)
     3 +ECOM(J)*((UKM(I+1,J)-2.*UKM(I,J)+UKM(I-1,J))/DX2+
     5  ((UKM(I,J+1)-UKM(I,J))/DY*IUM(I,J)*IUM(I,J+1)
     5  +(UKM(I,J-1)-UKM(I,J))/DY*IUM(I,J)*IUM(I,J-1))/DY)
     9         )
       UKP2(I,J)=UKM2(I,J)+TDT*IUM(I,J)*(
     1   0.25* ( FV(J)*(VK2(I,J)+VK2(I+1,J))
     1  +FV(J-1)*(VK2(I,J-1)+VK2(I+1,J-1)))-GDH2X(I,J)*DEX2(I)
     3  +ECOM(J)*((UKM2(I+1,J)-2.*UKM2(I,J)+UKM2(I-1,J))/DX2+
     4((UKM2(I,J+1)-UKM2(I,J))/DY*IUM(I,J)*IUM(I,J+1)
     4 +(UKM2(I,J-1)-UKM2(I,J))/DY*IUM(I,J)*IUM(I,J-1))/DY)
     9    + STUA(I) - STUB(I)
     9         )
30    CONTINUE
C
c.......switch on nonlinear terms if required....
c
	if(ifnonl.ne.0)then
c
	do j=2,jmaxm
	do i=2,imaxm
c
       UKP(I,J)=UKP(I,J)+TDT*IUM(I,J)*(
     5  -(((UK(I+1,J)+UK(I,J))*(US1(I+1,J)+US1(I,J))-
     6  (UK(I-1,J)+UK(I,J))*(US1(I-1,J)+US1(I,J)))/DX
     7  +((VK(I+1,J)+VK(I,J))*(US1(I,J+1)+US1(I,J))-
     8   (VK(I+1,J-1)+VK(I,J-1))*(US1(I,J)+US1(I,J-1)))/DY)/4.
     9         )
       UKP2(I,J)=UKP2(I,J)+TDT*IUM(I,J)*(
     5  -(((UK2(I+1,J)+UK2(I,J))*(US2(I+1,J)+US2(I,J))-
     6  (UK2(I-1,J)+UK2(I,J))*(US2(I-1,J)+US2(I,J)))/DX
     7  +((VK2(I+1,J)+VK2(I,J))*(US2(I,J+1)+US2(I,J))-
     8  (VK2(I+1,J-1)+VK2(I,J-1))*(US2(I,J-1)+US2(I,J)))/DY)/4.
     9    + STUA(I) - STUB(I)
     9         )
c
	enddo
	enddo
c
	endif
c
       DO 31 J=2,JMAXM
CDIR$ IVDEP
       DO 31 I=2,IMAXM
      DEY1(I)=HH1-(HK(I,J)+HK(I,J+1))/2.
      DEY2(I)=HH2-(HK2(I,J)+HK2(I,J+1)-HK(I,J)-HK(I,J+1))/2.
      STVA(I)=DVT*(VKM(I,J) - VKM2(I,J))
      STVB(I)=DVB*VKM2(I,J)
      VKP(I,J)=VKM(I,J)+TDT*IVM(I,J)*(
     1 -0.25*(FH(J+1)*(UK(I,J+1)+UK(I-1,J+1))
     1  +FH(J)*(UK(I,J)+UK(I-1,J)))-GDH1Y(I,J)*DEY1(I)
     2  +WY(I,J) - STVA(I)
     3 +ECOM(J)*(((VKM(I+1,J)-VKM(I,J))/DX*IVM(I,J)*IVM(I+1,J)+
     3 (VKM(I-1,J)-VKM(I,J))/DX*IVM(I,J)*IVM(I-1,J))/DX+
     4         (VKM(I,J+1)-2.*VKM(I,J)+VKM(I,J-1))/DY2)
     9         )
       VKP2(I,J)=VKM2(I,J)+TDT*IVM(I,J)*(
     1   -0.25*(FH(J+1)*(UK2(I,J+1)+UK2(I-1,J+1))
     1  +FH(J)*(UK2(I,J)+UK2(I-1,J)))-GDH2Y(I,J)*DEY2(I)
     3 +ECOM(J)*(((VKM2(I+1,J)-VKM2(I,J))/DX*IVM(I,J)*IVM(I+1,J)+
     3  (VKM2(I-1,J)-VKM2(I,J))/DX*IVM(I,J)*IVM(I-1,J))/DX+
     4         (VKM2(I,J+1)-2.*VKM2(I,J)+VKM2(I,J-1))/DY2)
     9     + STVA(I) - STVB(I)
     9         )
31     CONTINUE
c
c......switch on nonlinear terms if required....
c
	if(ifnonl.ne.0)then
c
	do j=2,jmaxm
	do i=2,imaxm
c
      VKP(I,J)=VKP(I,J)+TDT*IVM(I,J)*(
     5  -(((UK(I,J)+UK(I,J+1))*(VS1(I,J)+VS1(I+1,J))-
     6   (UK(I-1,J)+UK(I-1,J+1))*(VS1(I-1,J)+VS1(I,J)))/DX
     7  +((VK(I,J+1)+VK(I,J))*(VS1(I,J+1)+VS1(I,J))-
     8  (VK(I,J-1)+VK(I,J))*(VS1(I,J-1)+VS1(I,J)))/DY)/4.
     9         )
       VKP2(I,J)=VKP2(I,J)+TDT*IVM(I,J)*(
     5  -(((UK2(I,J)+UK2(I,J+1))*(VS2(I+1,J)+VS2(I,J))-
     6  (UK2(I-1,J)+UK2(I-1,J+1))*(VS2(I-1,J)+VS2(I,J)))/DX
     7  +((VK2(I,J+1)+VK2(I,J))*(VS2(I,J+1)+VS2(I,J))-
     8  (VK2(I,J-1)+VK2(I,J))*(VS2(I,J-1)+VS2(I,J)))/DY)/4.
     9     + STVA(I) - STVB(I)
     9         )
c
	enddo
	enddo
c
	endif
C
       DO 32 J=2,JMAXM
       DO 32 I=2,IMAXM
c
       h1ch=((UK(I,J)-UK(I-1,J))/DX+
     1 (VK(I,J)-VK(I,J-1))/DY) - DD*HKM(I,J)
c     2   +DAMP*(((HKM(I+1,J)-HKM(I,J))*IHM(I,J)*IHM(I+1,J)
c     3   +(HKM(I-1,J)-HKM(I,J))*IHM(I,J)*IHM(I-1,J))*RDX2
c     4   +((HKM(I,J+1)-HKM(I,J))*IHM(I,J)*IHM(I,J+1)
c     5   +(HKM(I,J-1)-HKM(I,J))*IHM(I,J-1)*IHM(I,J))*RDY2)
c
       HKP(I,J)=HKM(I,J)+TDT*IHM(I,J)*h1ch
c
       HKP2(I,J)=HKM2(I,J)+TDT*IHM(I,J)*(
     1      h1ch+((UK2(I,J)-UK2(I-1,J))/DX+
     1 (VK2(I,J)-VK2(I,J-1))/DY) + DD*HKM(I,J)
c     2   +DAMP*(((HKM2(I+1,J)-HKM2(I,J))*IHM(I,J)*IHM(I+1,J)
c     3   +(HKM2(I-1,J)-HKM2(I,J))*IHM(I,J)*IHM(I-1,J))*RDX2
c     4   +((HKM2(I,J+1)-HKM2(I,J))*IHM(I,J)*IHM(I,J+1)
c     5   +(HKM2(I,J-1)-HKM2(I,J))*IHM(I,J-1)*IHM(I,J))*RDY2)
     6          )
  32    CONTINUE
C
C......RESET H1, U1 AND V1 ARRAYS...........
C
       DO 6790 J=1,JM
       DO 6790 I=1,IM
       H1(I,J,KP)=HKP(I,J)
       U1(I,J,KP)=UKP(I,J)
       V1(I,J,KP)=VKP(I,J)
       H2(I,J,KP)=HKP2(I,J)
       U2(I,J,KP)=UKP2(I,J)
       V2(I,J,KP)=VKP2(I,J)
 6790  CONTINUE
C
       IF(KTAU/ISMTH*ISMTH.NE.KTAU) GO TO 700
       DO 600 J=1,JMAX
       DO 600 I=1,IMAX
       U1(I,J,KP)=.5*(U1(I,J,KP)+U1(I,J,K))
       U1(I,J,K )=.5*(U1(I,J,KM)+U1(I,J,K))
       V1(I,J,KP)=.5*(V1(I,J,KP)+V1(I,J,K))
       V1(I,J,K )=.5*(V1(I,J,KM)+V1(I,J,K))
       U2(I,J,KP)=.5*(U2(I,J,KP)+U2(I,J,K))
       U2(I,J,K )=.5*(U2(I,J,KM)+U2(I,J,K))
       V2(I,J,KP)=.5*(V2(I,J,KP)+V2(I,J,K))
       V2(I,J,K )=.5*(V2(I,J,KM)+V2(I,J,K))
       H1(I,J,KP)=.5*(H1(I,J,KP)+H1(I,J,K))
       H1(I,J,K )=.5*(H1(I,J,KM)+H1(I,J,K))
       H2(I,J,KP)=.5*(H2(I,J,KP)+H2(I,J,K))
       H2(I,J,K )=.5*(H2(I,J,KM)+H2(I,J,K))
 600    CONTINUE
 700    CONTINUE
C
C.....PREVENT THERMOCLINE FROM SURFACING........
C
c        DO 9191 J=1,JM
c        DO 9191 I=1,IM
c        H1(I,J,KP)=AMIN1(H1(I,J,KP),0.9*HH1)
c 9191   H2(I,J,KP)=AMIN1(H2(I,J,KP),0.9*HH2)
c
	do j=1,jm
	do i=1,im
	  if(h1(i,j,kp).gt.0.95*hh1)h1(i,j,kp)=0.95*hh1
	  if(h2(i,j,kp).gt.0.95*hh2)h2(i,j,kp)=0.95*hh2
	enddo
	enddo
c
c.......Add extra condition along western boundary
c
c	do j=1,jm
c	  do i=1,5
c	  if(h1(i,j,kp).gt.0.75*hh1)h1(i,j,kp)=0.75*hh1
c	  if(h2(i,j,kp).gt.0.75*hh2)h2(i,j,kp)=0.75*hh2
c	  enddo
c	enddo
C
C
c
	if(mod(ktau,16*30*12).eq.0)then
c
	    bigh1=-1.e20
	    bigh2=-1.e20
	    bigu1=-1.e20
	    bigu2=-1.e20
	    do j=1,jm
	      do i=1,im
	         if(abs(h1(i,j,kp)).gt.bigh1)bigh1=abs(h1(i,j,kp))
	         if(abs(h2(i,j,kp)).gt.bigh2)bigh2=abs(h2(i,j,kp))
	         if(abs(us1(i,j)).gt.bigu1)bigu1=abs(us1(i,j))
	         if(abs(us2(i,j)).gt.bigu2)bigu2=abs(us2(i,j))
	      enddo
	    enddo
c
	    print *,'ktau=',ktau,' bigh1,bigh2,bigu1,bigu2=',
     &       bigh1,bigh2,bigu1,bigu2
	endif
c
c
c##############################
c
c.......write h1 for plotting........
c
c##############################
c
	if(mod(ktau,16*30*ipltm).eq.0)then
c
c	  write(8)ihm
c	  write(8)((h1(i,j,k),i=1,im),j=1,jm)
c
	do j=1,jm
	  do i=1,im
	    h1p(i,j)=h1(i,j,k)
	    u1p(i,j)=us1(i,j)
	    v1p(i,j)=vs1(i,j)
	  enddo
	enddo
c
	do j=1,jm
	  do i=1,im
	    h2p(i,j)=h2(i,j,k)
	    u2p(i,j)=us2(i,j)
	    v2p(i,j)=vs2(i,j)
	  enddo
	enddo
c
	  call resplt(ktau,hh1,hh2,h1p,u1p,v1p,
     &        h2p,u2p,v2p,ihm,ium,ivm)
c
	print *,'done resplt!'
c
	endif
C
C
        TDT=2.0*DT
       KS=KM
       KM=K
       K=KP
       KP=KS
 987    FORMAT(1X,13E8.2)
  323 CONTINUE
C
C
       IF(KTAU.LT.JJSTOP) GO TO 500
 398       CONTINUE
C
c......close plot file.....
c
	call plotnd
c
C
 8889      FORMAT(1P,8E10.3)
 8890      FORMAT(I5)
 8892     FORMAT(4I5)
 8893     FORMAT(8I10)
C
      STOP
      END
      SUBROUTINE GPRINT(ARR,IMP,JMP,SC)
       REAL ARR(169,61)
      REAL DUM(168)
       IMPM=IMP-1
       JMPM=JMP-1
       DO 20 JJ=1,JMP
      J=JMP+1-JJ
      DO 10 I=1,IMP
   10 DUM(I)=ARR(I,J)/SC
C     WRITE(6,3252)J,(DUM(I),I=1,6),(DUM(I),I=10,60,5),(DUM(I),I=61,63)
      WRITE(6,3252)J,(DUM(I),I=1,46,2)
 20    CONTINUE
      WRITE(6,528)
 3252    FORMAT(1X,I2,23F5.2)
 528    FORMAT(//)
      RETURN
       END
        SUBROUTINE TRKVAR(A,ITM,IOBS,IUA,VAR)
C
        REAL A(169,61)
C
        INTEGER IUA(169,61),ITM(169,61),IOBS(169,61)
C
        IM=169
        JM=61
        VAR=0.0
        N=0
        SUMSQ=0.0
C
        DO 100 J=1,JM
        DO 100 I=1,IM
        SUMSQ=SUMSQ+A(I,J)*A(I,J)*IUA(I,J)*IOBS(I,J)*ITM(I,J)
 100    N=N+ITM(I,J)*IUA(I,J)*IOBS(I,J)
        IF(N.EQ.0)GOTO 110
        AN=N
        SUMSQ=SUMSQ/AN
 110    VAR=SQRT(SUMSQ)
        RETURN
        END
        SUBROUTINE REGVAR(A,I1,I2,J1,J2,IUA,IOBS,VAR)
C
        REAL A(169,61)
C
        INTEGER IUA(169,61),IOBS(169,61)
C
        VAR=0.0
        SUMSQ=0.0
        N=0
C
        DO 100 J=J1,J2
        DO 100 I=I1,I2
        SUMSQ=SUMSQ+A(I,J)*A(I,J)*IUA(I,J)*IOBS(I,J)
 100    N=N+IUA(I,J)*IOBS(I,J)
        IF(N.EQ.0)GOTO 110
        AN=N
        SUMSQ=SUMSQ/AN
 110    VAR=SQRT(SUMSQ)
        RETURN
        END
        SUBROUTINE CALRMS(A,I1,I2,J1,J2,IUA,IOBS,N,RMS)
C
        REAL A(169,61)
        INTEGER IUA(169,61),IOBS(169,61)
C
        RMS=0.0
        N=0
        DO 100 J=J1,J2
        DO 100 I=I1,I2
        RMS=RMS+A(I,J)*A(I,J)*IUA(I,J)*IOBS(I,J)
 100    N=N+IUA(I,J)*IOBS(I,J)
        IF(N.EQ.0)GOTO 110
        AN=N
        RMS=SQRT(RMS/AN)
 110     CONTINUE
        RETURN
        END
        SUBROUTINE RMSTRK(A,ITM,IOBS,IUA,N,RMS)
C
        REAL A(169,61)
C
        INTEGER IUA(169,61),ITM(169,61),IOBS(169,61)
C
        IM=169
        JM=61
        RMS=0.0
        N=0
        DO 100 J=1,JM
        DO 100 I=1,IM
        RMS=RMS+A(I,J)*A(I,J)*ITM(I,J)*IOBS(I,J)*IUA(I,J)
 100    N=N+ITM(I,J)*IOBS(I,J)*IUA(I,J)
        IF(N.EQ.0)GOTO 110
        AN=N
        RMS=SQRT(RMS/AN)
 110      CONTINUE
        RETURN
        END
       SUBROUTINE SETRAR(A,LEN,VALUE)
       REAL A(LEN)
       DO 320 I=1,LEN
  320  A(I)=VALUE
       RETURN
       END
        SUBROUTINE STRESS(WX,WY,IUM,IVM,IM,JM,IDATE)
C
        REAL WX(IM,JM),WY(IM,JM),DUM(169,61)
        INTEGER IUM(IM,JM),IVM(IM,JM)
C
C.....DUMMY ARRAY.....
C
        DO 100 J=1,JM
        DO 100 I=1,IM
 100     DUM(I,J)=0.0
         ITEST=IDATE
C
C
C.......SET UP CORRECT WIND DATA INPUT CHANNEL...
C
        IF(IDATE.EQ.820110)NIN=3
        GOTO 998
 999    NP=NIN
        IF(NP.EQ.3)NIN=90
        IF(NP.GT.80)NIN=NIN+1
        IF(NP.EQ.96)STOP 7777
 998    CONTINUE
         READ(NIN,310,END=999)IDATE,IDTYP,
     1  ((DUM(I,J),I=2+(J/60),168),J=60,2,-1)
            READ(NIN,300,END=999)IDATE,IDTYP,
     2                      ((DUM(I,J),I=2,168),J=60,2,-1)
            READ(NIN,300,END=999)IDATE,IDTYP,
     2                      ((WX(I,J),I=2,168),J=60,2,-1)
            READ(NIN,300,END=999)IDATE,IDTYP,
     2                      ((WY(I,J),I=2,168),J=60,2,-1)
 300    FORMAT(I6,I1,35X,13F6.1/(20F6.1))
 310      FORMAT(I6,I1,41X,12F6.1/(20F6.1))
        WRITE(6,800)IDATE
 800   FORMAT(1X,'WIND READ AT   ',I7)
C
C
C.......CONVERT PSEUDO WIND STRESS TO REAL WIND STRESS IN
C       NEWTONS PER SQUARE METRE USING THE DRAG EQUATION....
C
        RHO=1.25
        DO 766 J=1,JM
        DO 766 I=1,IM
        AMODV=(WX(I,J)*WX(I,J)+WY(I,J)*WY(I,J))**0.25
        CD=1.0E-3*(0.61+0.063*AMODV)
        IF(AMODV.LE.7.8)CD=1.1E-3
        WX(I,J)=CD*RHO*WX(I,J)*IUM(I,J)
 766    WY(I,J)=CD*RHO*WY(I,J)*IVM(I,J)
C
C.......CONVERT WIND UNITS TO DYNES PER SQUARE CENTIMETRE....
        DO 767 J=1,JM
        DO 767 I=1,IM
        WX(I,J)=WX(I,J)*10.0
 767    WY(I,J)=WY(I,J)*10.0
C
        RETURN
        END
C
        SUBROUTINE TOTWIND(WX,WY,IM,JM,MONTH,IYEAR,KTAU)
C
        REAL WX(IM,JM),WY(IM,JM),TAUX(84,30),TAUY(84,30),
     1  TIX(94,58),TIY(94,58),WXT(169,61),WYT(169,61)
C
C
C.......READ IN PACIFIC WIND AND CONVERT TO STRESS
C
        NIN=3
        IF(KTAU.GE.39664)NIN=10
        READ(NIN,7614)MONTH,IYEAR,TAUX,TAUY
        WRITE(6,757)MONTH,IYEAR
 7614   FORMAT(2I5,14F5.1/(16F5.1))
 757    FORMAT(1X,'WIND READ ON MONTH=',I3,' YEAR=19',I2)
        CALL WIND(TAUX,TAUY,WXT,WYT,169,61)
C
C.......READ IN INDIAN WIND AND CONVERT TO STRESS....
C
        READ(8,7614)INDM,INDY,TIX,TIY
        IF(INDM.NE.MONTH)PRINT *,'WIND MONTHS DONT MATCH'
        IF((INDY-1900).NE.IYEAR)PRINT *,'WIND YEARS DONT MATCH'
        IF(INDM.NE.MONTH)STOP 8888
        IF((INDY-1900).NE.IYEAR)STOP 9999
C
        DO 680 J=1,58
        DO 680 I=1,94
        IF(TIX(I,J).GT.980.0)TIX(I,J)=0.0
 680    IF(TIY(I,J).GT.980.0)TIY(I,J)=0.0
C
        RHO=1.25
        DO 681 J=1,58
        DO 681 I=1,93
        AMODV=(TIX(I,J)*TIX(I,J)+TIY(I,J)*TIY(I,J))**0.25
        CD=1.0E-3*(0.61+0.063*AMODV)
        IF(AMODV.LE.7.8)CD=1.1E-3
        WX(I,J)=CD*RHO*TIX(I,J)
 681    WY(I,J)=CD*RHO*TIY(I,J)
C
C........CONVERT INDIAN WINDS TO DYNES/CM2....
C
        DO 684 J=1,58
        DO 684 I=1,93
        WX(I,J)=WX(I,J)*10.0
  684   WY(I,J)=WY(I,J)*10.0
C
C........MAP PACIFIC WIND STRESS INTO WX AND WY.....
C
        DO 682 J=1,61
        DO 682 I=2,169
        WX(I+92,J+2)=WXT(I,J)
 682    WY(I+92,J+2)=WYT(I,J)
C
C.......LINEARLY INTERPOLATE BETWEEN INDIAN AND PACIFIC WINDS
C       ALONG I=93
C
        DO 683 J=1,JM
        WX(93,J)=0.5*(WX(92,J)+WX(94,J))
 683    WY(93,J)=0.5*(WY(92,J)+WY(94,J))
        RETURN
        END
C
        SUBROUTINE WIND(TAUX,TAUY,WX,WY,IM,JM)
C
        REAL TAUX(84,30),TAUY(84,30),WX(IM,JM),WY(IM,JM)
C
        IMAX=IM
        JMAX=JM
C.......SET TAUX AND TAUY TO ZERO OVER LAND....
        DO 758 J=1,30
        DO 758 I=1,84
        IF(TAUX(I,J).GT.998.0)TAUX(I,J)=0.0
  758   IF(TAUY(I,J).GT.998.0)TAUY(I,J)=0.0
C
C.......PUT PSEUDO WIND STRESS ONTO A 2 DEGREE GRID....
C
        DO 760 J=1,30
        DO 760 I=1,84
        WX(2*I,2*J)=TAUX(I,J)
 760    WY(2*I,2*J)=TAUY(I,J)
C.......SET WX AND WY TO ZERO AT BOUNDARIES...
        DO 761 J=1,JM
        WX(1,J)=0.0
        WX(IM,J)=0.0
        WY(1,J)=0.0
 761    WY(IM,J)=0.0
C
        DO 762 I=1,IM
        WX(I,1)=0.0
        WX(I,JM)=0.0
        WY(I,1)=0.0
  762   WY(I,JM)=0.0
C......LINEARLY INTERPOLATE PSEUDO WIND STRESS ONTO A I DEGREE GRID..
C
        JMAXM2=JMAX-2
        IMAXM2=IMAX-2
        JMAXM=JMAX-1
        IMAXM=IMAX-1
        DO 763 J=2,JMAXM,2
        DO 763 I=3,IMAXM2,2
        WX(I,J)=0.5*(WX(I+1,J)+WX(I-1,J))
 763    WY(I,J)=0.5*(WY(I+1,J)+WY(I-1,J))
C
        DO 764 J=3,JMAXM2,2
        DO 764 I=2,IMAXM,2
        WX(I,J)=0.5*(WX(I,J+1)+WX(I,J-1))
 764    WY(I,J)=0.5*(WY(I,J+1)+WY(I,J-1))
C
        DO 765 J=3,JMAXM2,2
        DO 765 I=3,IMAXM2,2
        WX(I,J)=0.25*(WX(I-1,J-1)+WX(I-1,J+1)+WX(I+1,J-1)+WX(I+1,J+1))
 765    WY(I,J)=0.25*(WY(I-1,J-1)+WY(I-1,J+1)+WY(I+1,J-1)+WY(I+1,J+1))
C
C.......CONVERT PSEUDO WIND STRESS TO REAL WIND STRESS IN
C       NEWTONS PER SQUARE METRE USING THE DRAG EQUATION....
C
        RHO=1.25
        DO 766 J=1,JM
        DO 766 I=1,IM
        AMODV=(WX(I,J)*WX(I,J)+WY(I,J)*WY(I,J))**0.25
        CD=1.0E-3*(0.61+0.063*AMODV)
        IF(AMODV.LE.7.8)CD=1.1E-3
        WX(I,J)=CD*RHO*WX(I,J)
 766    WY(I,J)=CD*RHO*WY(I,J)
C
C.......CONVERT WIND UNITS TO DYNES PER SQUARE CENTIMETRE....
        DO 767 J=1,JM
        DO 767 I=1,IM
        WX(I,J)=WX(I,J)*10.0
 767    WY(I,J)=WY(I,J)*10.0
C
        RETURN
        END
        SUBROUTINE ENERGY(H1,US1,VS1,HBAR,G12
     1                   ,POTE,AKINE,IHM,IUM,IVM,IM,JM)
C
        REAL H1(IM,JM),US1(IM,JM),VS1(IM,JM),POTEGY(63),AKEGY(63)
     1       ,POTE(7),AKINE(7),ALAT1(7),ALAT2(7),ALONG1(7)
     2      ,ALONG2(7)
        INTEGER IHM(IM,JM),IUM(IM,JM),IVM(IM,JM)
     1          ,IX1(7),IX2(7),JY1(7),JY2(7)
C
        DATA ALAT1/5.,5.,-5.,-5.,-17.,-17.,-24./
        DATA ALAT2/20.,20.,5.,5.,-5.,-5.,-17./
        DATA ALONG1/144.,176.,150.,235.,154.,202.,165./
        DATA ALONG2/160.,220.,205.,280.,198.,260.,275./
C........CALC. INDICES FOR TRAPEZOIDAL SUMMATIONS.....
C
        DO 800 LMN=1,7
        IX1(LMN)=ALONG1(LMN)-122.5
        IX2(LMN)=ALONG2(LMN)-122.5
        JY1(LMN)=ALAT1(LMN)+31.5
        JY2(LMN)=ALAT2(LMN)+31.5
 800    CONTINUE
C
        DO 111 II=1,7
C
        DO 101 J=1,JM
        POTEGY(J)=0.0
 101    AKEGY(J)=0.0
        I1=IX1(II)
        I2=IX2(II)
        J1=JY1(II)
        J2=JY2(II)
C
        I1P1=I1+1
        I2M1=I2-1
        J1P1=J1+1
        J2M1=J2-1
C........CALC. PE AND KE ALONG LINES OF CONSTANT LAT...
C
        DO 10 J=J1,J2
        SUMETA2=0.0
        SUMVEL2=0.0
        DO 11 I=I1P1,I2M1
        SUMETA2=SUMETA2+H1(I,J)*H1(I,J)*IHM(I,J)
        SUMVEL2=SUMVEL2+US1(I,J)*US1(I,J)*IUM(I,J)
     1         +VS1(I,J)*VS1(I,J)*IVM(I,J)
  11       CONTINUE
        POTEGY(J)=0.25*G12*(H1(I1,J)*H1(I1,J)*IHM(I1,J)
     1           +2.*SUMETA2+H1(I2,J)*H1(I2,J)*IHM(I2,J))
        AKEGY(J)=0.25*HBAR*(US1(I1,J)*US1(I1,J)*IUM(I1,J)
     1          +VS1(I1,J)*VS1(I1,J)*IVM(I1,J)+2.*SUMVEL2
     2          +US1(I2,J)*US1(I2,J)*IUM(I2,J)
     3          +VS1(I2,J)*VS1(I2,J)*IVM(I2,J))
 10     CONTINUE
C
C......CALC. PE AND KE OVER MODEL DOMAINS.....
C
        SUMPE=0.0
        SUMKE=0.0
        DO 12 J=J1P1,J2M1
        SUMPE=SUMPE+POTEGY(J)
        SUMKE=SUMKE+AKEGY(J)
 12     CONTINUE
        POTE(II)=0.5*(POTEGY(J1)+2.*SUMPE+POTEGY(J2))
        AKINE(II)=0.5*(AKEGY(J1)+2.*SUMKE+AKEGY(J2))
 111    CONTINUE
        RETURN
        END
C
        SUBROUTINE SPINUP(H1,US1,VS1,HBAR,G12
     1                   ,POTE,AKINE,IHM,IUM,IVM,IM,JM)
C
        REAL H1(IM,JM),US1(IM,JM),VS1(IM,JM),POTEGY(63),AKEGY(63)
     1       ,POTE(6),AKINE(6),ALAT1(6),ALAT2(6),ALONG1(6)
     2      ,ALONG2(6)
        INTEGER IHM(IM,JM),IUM(IM,JM),IVM(IM,JM)
     1          ,IX1(6),IX2(6),JY1(6),JY2(6)
C
        DATA ALAT1/1.,28.,38.,1.,28.,38./
        DATA ALAT2/27.,37.,63.,27.,37.,63./
        DATA ALONG1/1.,1.,1.,93.,93.,93./
        DATA ALONG2/94.,94.,94.,261.,261.,261./
C........CALC. INDICES FOR TRAPEZOIDAL SUMMATIONS.....
C
        DO 800 LMN=1,6
        IX1(LMN)=ALONG1(LMN)
        IX2(LMN)=ALONG2(LMN)
        JY1(LMN)=ALAT1(LMN)
        JY2(LMN)=ALAT2(LMN)
 800    CONTINUE
C
        DO 111 II=1,6
C
        DO 101 J=1,JM
        POTEGY(J)=0.0
 101    AKEGY(J)=0.0
        I1=IX1(II)
        I2=IX2(II)
        J1=JY1(II)
        J2=JY2(II)
C
        I1P1=I1+1
        I2M1=I2-1
        J1P1=J1+1
        J2M1=J2-1
C........CALC. PE AND KE ALONG LINES OF CONSTANT LAT...
C
        DO 10 J=J1,J2
        SUMETA2=0.0
        SUMVEL2=0.0
        DO 11 I=I1P1,I2M1
        SUMETA2=SUMETA2+H1(I,J)*H1(I,J)*IHM(I,J)
        SUMVEL2=SUMVEL2+US1(I,J)*US1(I,J)*IUM(I,J)
     1         +VS1(I,J)*VS1(I,J)*IVM(I,J)
  11       CONTINUE
        POTEGY(J)=0.25*G12*(H1(I1,J)*H1(I1,J)*IHM(I1,J)
     1           +2.*SUMETA2+H1(I2,J)*H1(I2,J)*IHM(I2,J))
        AKEGY(J)=0.25*HBAR*(US1(I1,J)*US1(I1,J)*IUM(I1,J)
     1          +VS1(I1,J)*VS1(I1,J)*IVM(I1,J)+2.*SUMVEL2
     2          +US1(I2,J)*US1(I2,J)*IUM(I2,J)
     3          +VS1(I2,J)*VS1(I2,J)*IVM(I2,J))
 10     CONTINUE
C
C......CALC. PE AND KE OVER MODEL DOMAINS.....
C
        SUMPE=0.0
        SUMKE=0.0
        DO 12 J=J1P1,J2M1
        SUMPE=SUMPE+POTEGY(J)
        SUMKE=SUMKE+AKEGY(J)
 12     CONTINUE
        POTE(II)=0.5*(POTEGY(J1)+2.*SUMPE+POTEGY(J2))
        AKINE(II)=0.5*(AKEGY(J1)+2.*SUMKE+AKEGY(J2))
 111    CONTINUE
        RETURN
        END
c
	subroutine resplt(ktau,hh1,hh2,h1,u1,v1,h2,u2,v2,
     &             ihm,ium,ivm)
c
         parameter(im=51,jm=51,imm=im-2,jmm=jm-2)
c
	 parameter(nval=5)
	 parameter(nval1=30)
c
	real h1(im,jm),hp(imm,jmm)
c
	real u1(im,jm),up(im,jm)
c
	real v1(im,jm),vp(im,jm)
c
	real h2(im,jm),h2p(imm,jmm)
c
	real u2(im,jm),u2p(imm,jmm)
c
	real v2(im,jm),v2p(imm,jmm)
c
	real speed(im,jm)
c
	real q1(im,jm),q2(im,jm),qgeos(im,jm)
     &    ,psib(im,jm)
c
	real psibt(im,jm)
c
	real spv(2),cfill(3),grylev(3),cval(nval),cval1(1)
c
	real conp(nval)
	real conm(nval)
	real conp1(2*nval)
	real conp2(nval1)
	real conm2(nval1)
c
	real sver1(im,jm),sver2(im,jm),wcurl(im,jm),sver(im,jm)
c
	real sum1(im,jm),sum2(im,jm),sumc(im,jm)
c
	real xlrr(4),ylrr(4)
c
	integer ihm(im,jm)
c
	integer ium(im,jm)
c
	integer ivm(im,jm)
c
	character carrow*10,char1*3,char2*7,char4*7,char3*17
c
	common/vars/beta,dxm,dym,g12,g23,f0,
     &         fh(jm),fv(jm),wx(im,jm),wy(im,jm)
c
	common/vars1/je,w0
c
      COMMON/CONPAR/ISPEC,IOFFPP,SPVALL,ILEGG,ILABB,NHII,NDECCN,NLBLL,
     .              LSCAL,LDASH,HGTLAB
c
      DATA ISPEC,IOFFPP,SPVALL,ILEGG,ILABB,NHII,NDECCN,NLBLL,
     .   LSCAL,LDASH,HGTLAB
     1  /   1 ,     1,  -99999., 1,    1,   -1,    1,     2, 1, 0, 0/
c
	char1='t= '
	char4=' months'
	ptime=float(ktau)/(16.*30.)
	write(char2,'(f7.1)')ptime
	char3=char1//char2//char4
c
	ch1=0.25
	ch2=0.2
	ch3=0.2
c
	pscal=1.e12
c
c.......compute the potential vorticity in each layer.....
c
	do j=1,jm
	   do i=1,im
	     q1(i,j)=fv(j)/hh1
	     q2(i,j)=fv(j)/hh2
	     qgeos(i,j)=fh(j)-f0
	   enddo
	enddo
c
	do j=1,jm-1
	 do i=1,im-1
	   hfl=0.25*(h1(i,j)+h1(i+1,j)+h1(i,j+1)+h1(i+1,j+1))
	   zeta1=(v1(i+1,j)-v1(i,j))/dxm-(u1(i,j+1)-u1(i,j))/dym
	   q1(i,j)=(zeta1+fv(j))/(hh1-hfl)
	 enddo
	enddo
c
	do j=1,jm-1
	 do i=1,im-1
	   hfl1=0.25*(h1(i,j)+h1(i+1,j)+h1(i,j+1)+h1(i+1,j+1))
	   hfl2=0.25*(h2(i,j)+h2(i+1,j)+h2(i,j+1)+h2(i+1,j+1))
	   zeta2=(v2(i+1,j)-v2(i,j))/dxm-(u2(i,j+1)-u2(i,j))/dym
	   q2(i,j)=(zeta2+fv(j))/(hh2+hfl1-hfl2)
	 enddo
	enddo
c
c.......compute the sverdrup transport and wind stress curl....
c
	do j=1,jm-1
	 do i=1,im-1
	    hfl1=hh1-0.25*(h1(i,j)+h1(i+1,j)+h1(i,j+1)+h1(i+1,j+1))
	    hfl2=hh2+hfl1
     &           -0.25*(h2(i,j)+h2(i+1,j)+h2(i,j+1)+h2(i+1,j+1))
	    vfl1=0.5*(v1(i,j)+v1(i+1,j))
	    vfl2=0.5*(v2(i,j)+v2(i+1,j))
	    wcurl(i,j)=(wy(i+1,j)-wy(i,j))/dxm-(wx(i,j+1)-wx(i,j))/dym
	    sver1(i,j)=hfl1*vfl1
	    sver2(i,j)=hfl2*vfl2
	    sver(i,j)=sver1(i,j)+sver2(i,j)
	 enddo
	enddo
c
c.......zonally integrate the Sverdrup transport and wind stress curl....
c
	do j=1,jm-1
	   do i=im,1,-1
	      sum1(i,j)=0.0
	      sum2(i,j)=0.0
	      sumc(i,j)=0.0
	      do ii=im,i+1,-1
	         sum1(i,j)=sum1(i,j)+sver1(ii,j)*dxm/pscal
	         sum2(i,j)=sum2(i,j)+sver2(ii,j)*dxm/pscal
	         sumc(i,j)=sumc(i,j)+wcurl(ii,j)*dxm/(beta*pscal)
	      enddo
	   enddo
	enddo
c
	do j=1,jm
	  do i=1,im
	    sver1(i,j)=-1.*sum1(i,j)
	    sver2(i,j)=-1.*sum2(i,j)
	    wcurl(i,j)=-1.*sumc(i,j)
	    sver(i,j)=sver1(i,j)+sver2(i,j)
	  enddo
	enddo
c
c.......compute the geostrophic contours......
c
	fhat=f0*f0*(hh1+hh2)/(g12*hh1*hh2)
	cr=beta/fhat
c
	print *,'cr=',cr
c
	do j=1,jm
	  do i=1,im
c	    psi1=g12*(hh1-h1(i,j))+g23*(hh2-h2(i,j))
c	    psi2=g23*(hh2-h2(i,j))
	    psi1=g12*h1(i,j)+g23*h2(i,j)
	    psi2=g23*h2(i,j)
	    psi1=-1.*psi1/f0
	    psi2=-1.*psi2/f0
	    psib(i,j)=(hh1*psi1+hh2*psi2)/(hh1+hh2)
	    psibt(i,j)=w0*cos(2.*4.*atan(1.0)*float(j-je)/float(2*jm))
     &       *f0*float(im-i)*dxm/(beta*(hh1+hh2))
	    qgeos(i,j)=fh(j)-f0+fhat*psib(i,j)
	  enddo
	enddo
c
	cval1(1)=0.5
c
	do i=1,nval
	   cval(i)=float(i-1)*0.2
	enddo
c
        xlen=5.
        ylen=float(jm)*xlen/float(im)
c
	call factor(3./5.)
c
        call plot(1.,9.,-3)
c
	chh1=hh1/100.
	chh2=hh2/100.
c
c......convert h to m......
c
	do j=1,jm
	  do i=1,im
	    h1(i,j)=h1(i,j)/100.
	    h2(i,j)=h2(i,j)/100.
	  enddo
	enddo
c
c......find the largest value of h and compute the
c      contour interval
c
	bigh1=-1.e20
	bigh2=-1.e20
	do j=1,jm
	  do i=1,im
	   if(abs(h1(i,j)).gt.bigh1)bigh1=abs(h1(i,j))
	   if(abs(h2(i,j)).gt.bigh2)bigh2=abs(h2(i,j))
	  enddo
	enddo
c
	print *,'bigh1=',bigh1
	print *,'bigh2=',bigh2
c
c	dc=abs(chh1-bigh)/float(nval)
	bigh=bigh1
	dc=bigh/float(nval)
	do i=1,nval
c	  conp(i)=float(i)*dc+chh1
c	  conm(i)=-1.*abs(chh1-bigh)+float(i-1)*dc+chh1
	  conp(i)=float(i)*dc
	  conm(i)=-1.*bigh+float(i-1)*dc
	enddo
c
c.......calc perturbation h and set special values....
c
	do j=1,jm
	 do i=1,im
	  if(ihm(i,j).eq.0)h1(i,j)=spvall
	 enddo
	enddo
c
	do j=2,jm-1
	   jj=j-1
	   do i=2,im-1
           ii=i-1
	   hp(ii,jj)=h1(i,j)
	   enddo
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
c        call confill(hp,imm,imm,jmm,xlen,ylen,
c     &    cfill,grylev,3,ioffpp,spvall)
c
        call confill(h1,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
c
c        call conrec(hp,imm,imm,jmm,xlen,ylen,conp,nval)
c        call conrec(hp,imm,imm,jmm,xlen,ylen,conm,nval)
c
        call conrec(h1,im,im,jm,xlen,ylen,conp,nval)
        call conrec(h1,im,im,jm,xlen,ylen,conm,nval)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HA,0.,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1Hh,0.,1,1)
	call subber(1H1,1,ch1,0.)
c
	call keksymc(xlen+0.5,1.2*ylen,1.2*ch1,char3,0.,17,1)
c
	call coast(xlen,ylen)
c
c
	call plot(xlen+1.,0.,-3)
c
	bigh=bigh2
	dc=bigh/float(nval)
	do i=1,nval
c	  conp(i)=float(i)*dc+chh1
c	  conm(i)=-1.*abs(chh1-bigh)+float(i-1)*dc+chh1
	  conp(i)=float(i)*dc
	  conm(i)=-1.*bigh+float(i-1)*dc
	enddo
c
c.......calc perturbation h and set special values....
c
	do j=1,jm
	 do i=1,im
	  if(ihm(i,j).eq.0)h2(i,j)=spvall
	 enddo
	enddo
c
	do j=2,jm-1
	   jj=j-1
	   do i=2,im-1
           ii=i-1
	   hp(ii,jj)=h2(i,j)
	   enddo
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
c        call confill(hp,imm,imm,jmm,xlen,ylen,
c     &    cfill,grylev,3,ioffpp,spvall)
c
        call confill(h2,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
c
c        call conrec(hp,imm,imm,jmm,xlen,ylen,conp,nval)
c        call conrec(hp,imm,imm,jmm,xlen,ylen,conm,nval)
c
        call conrec(h2,im,im,jm,xlen,ylen,conp,nval)
        call conrec(h2,im,im,jm,xlen,ylen,conm,nval)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HB,0.,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1Hh,0.,1,1)
	call subber(1H2,1,ch1,0.)
c
c
	call coast(xlen,ylen)
c
	call plot(-1.*(xlen+1.),-1.*(ylen+1.0),-3)
c
	do j=1,jm
	 do i=1,im
	   speed(i,j)=sqrt(u1(i,j)*u1(i,j)+v1(i,j)*v1(i,j))/100.
	   up(i,j)=0.0
	   vp(i,j)=0.0
	 enddo
	enddo
c
	do j=2,jm-1
	   do i=2,im-1
	   up(i,j)=0.5*(u1(i,j)+u1(i-1,j))/100.
	   vp(i,j)=0.5*(v1(i,j)+v1(i-1,j))/100.
	   enddo
	enddo
c
        call border(xlen,ylen,0000,1111,1,1,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1HV,0.,1,1)
	call subber(1H1,1,ch1,0.)
c
	call keksym(-0.5,1.05*ylen,ch3,1HC,0.,1,1)
c
	bigwn=-1.e20
	do j=1,jm
	  do i=1,im
c	   spd=sqrt(up(i,j)*up(i,j)+vp(i,j)*vp(i,j))/100.
	   spd=sqrt(u1(i,j)*u1(i,j)+v1(i,j)*v1(i,j))/100.
	   if(spd.gt.bigwn)bigwn=spd
	  enddo
	enddo
c
	print *,'bigwn=',bigwn
c
	bigh=bigwn
	dc=bigh/float(nval)
	do i=1,nval
	  conp(i)=float(i)*dc
	enddo
c
	call setlw(0.008)
        call conrec(speed,im,im,jm,xlen,ylen,conp,nval)
	call setlw(0.0)
c
        dx=xlen/float(im-1)
        dy=ylen/float(jm-1)
        scale=0.3
        do j=2,jm-1
        y=float(j-1)*dy
        do i=2,im-1
        x=float(i-1)*dx
        spd=sqrt(up(i,j)*up(i,j)+vp(i,j)*vp(i,j))/100.
        alen=scale*spd
	x1=x
	y1=y
	if(bigwn.gt.0.0)then
        x1=x+up(i,j)*scale/bigwn
        y1=y+vp(i,j)*scale/bigwn
	endif
c
c	if(spd.gt.0.)then
c           x1=x+up(i,j)*scale
c           y1=y+vp(i,j)*scale
c	endif
c
c        if(spd.ge.0.3*bigwn)call arrow(x,y,x1,y1,0.08,20.,0)
c        if(spd.ge.1.e-6)call arrow(x,y,x1,y1,0.08,20.,0)
c        if(spd.ge.0.03*bigwn.and.spd.gt.0.0)then
        if(spd.gt.0.0)then
        ang=atan2((y1-y),(x1-x))
	pi=4.*atan(1.0)
        ang=0.5*pi-ang
        ang=360.*ang/(2.*pi)
c        call sldlin(x,y,x1,y1,0.01)
        call farohed(x1,y1,ang,0.08,20.,0,.true.)
        endif
        enddo
        enddo
c
c        call sldlin(0.,-0.3,scale,-0.3,0.01)
c        call farohed(scale,-0.3,90.,0.08,20.,0,.true.)
c
c	write(carrow,'(1p,1e10.3)')bigwn
c
c	call keksymc(scale+0.1,-0.3,ch2,carrow,0.,10,0)
c
	call coast(xlen,ylen)
c
	call plot((xlen+1.),0.,-3)
c
	do j=1,jm
	 do i=1,im
	   speed(i,j)=sqrt(u2(i,j)*u2(i,j)+v2(i,j)*v2(i,j))/100.
	   up(i,j)=0.0
	   vp(i,j)=0.0
	 enddo
	enddo
c
	do j=2,jm-1
	   do i=2,im-1
	   up(i,j)=0.5*(u2(i,j)+u2(i-1,j))/100.
	   vp(i,j)=0.5*(v2(i,j)+v2(i-1,j))/100.
	   enddo
	enddo
c
        call border(xlen,ylen,0000,1111,1,1,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1HV,0.,1,1)
	call subber(1H2,1,ch1,0.)
c
	call keksym(-0.5,1.05*ylen,ch3,1HD,0.,1,1)
c
	bigwn=-1.e20
	do j=1,jm
	  do i=1,im
c	   spd=sqrt(up(i,j)*up(i,j)+vp(i,j)*vp(i,j))/100.
	   spd=sqrt(u2(i,j)*u2(i,j)+v2(i,j)*v2(i,j))/100.
	   if(spd.gt.bigwn)bigwn=spd
	  enddo
	enddo
c
	print *,'bigwn=',bigwn
c
	bigh=bigwn
	dc=bigh/float(nval)
	do i=1,nval
	  conp(i)=float(i)*dc
	enddo
c
	call setlw(0.008)
        call conrec(speed,im,im,jm,xlen,ylen,conp,nval)
	call setlw(0.0)
c
        dx=xlen/float(im-1)
        dy=ylen/float(jm-1)
        scale=0.3
        do j=2,jm-1
        y=float(j-1)*dy
        do i=2,im-1
        x=float(i-1)*dx
        spd=sqrt(up(i,j)*up(i,j)+vp(i,j)*vp(i,j))/100.
        alen=scale*spd
	x1=x
	y1=y
	if(bigwn.gt.0.0)then
        x1=x+up(i,j)*scale/bigwn
        y1=y+vp(i,j)*scale/bigwn
	endif
c
c
c	if(spd.gt.0.)then
c           x1=x+up(i,j)*scale/spd
c           y1=y+vp(i,j)*scale/spd
c	endif
c
c        if(spd.ge.0.3*bigwn)call arrow(x,y,x1,y1,0.08,20.,0)
c        if(spd.ge.1.e-6)call arrow(x,y,x1,y1,0.08,20.,0)
c        if(spd.ge.0.1*bigwn.and.spd.gt.0.0)then
        if(spd.ge.0.0)then
        ang=atan2((y1-y),(x1-x))
	pi=4.*atan(1.0)
        ang=0.5*pi-ang
        ang=360.*ang/(2.*pi)
c        call sldlin(x,y,x1,y1,0.01)
        call farohed(x1,y1,ang,0.08,20.,0,.true.)
        endif
        enddo
        enddo
c
c        call sldlin(0.,-0.3,scale,-0.3,0.01)
c        call farohed(scale,-0.3,90.,0.08,20.,0,.true.)
c
c	write(carrow,'(1p,1e10.3)')bigwn
c
c	call keksymc(scale+0.1,-0.3,ch2,carrow,0.,10,0)
c
	call coast(xlen,ylen)
c
	call chopit(0.,0.)
c
	lscal=0
C
        call plot(1.,9.,-3)
c
c......find the largest value of q1 and q2 and compute the
c      contour interval
c
	bigq1=-1.e20
	bigq2=-1.e20
	bigpb=-1.e20
	smlq1=1.e20
	smlq2=1.e20
	smlpb=1.e20
	bigqg=-1.e20
	smlqg=1.e20
	do j=1,jm
	  do i=1,im
	   if(abs(psib(i,j)).gt.bigpb)bigpb=abs(psib(i,j))
	   if(abs(psib(i,j)).lt.smlpb)smlpb=abs(psib(i,j))
	  enddo
	enddo
	do j=1,jm
	  do i=1,im
	   if(abs(q1(i,j)).gt.bigq1)bigq1=abs(q1(i,j))
	   if(abs(q2(i,j)).gt.bigq2)bigq2=abs(q2(i,j))
	   if(abs(qgeos(i,j)).gt.bigqg)bigqg=abs(qgeos(i,j))
	   if(abs(q1(i,j)).lt.smlq1)smlq1=abs(q1(i,j))
	   if(abs(q2(i,j)).lt.smlq2)smlq2=abs(q2(i,j))
	   if(abs(qgeos(i,j)).lt.smlqg)smlqg=abs(qgeos(i,j))
	  enddo
	enddo
c
	print *,'bigq1=',bigq1
	print *,'bigq2=',bigq2
	print *,'bigqg=',bigqg
c
	bigh=bigq1-smlq1
	dc=bigh/float(2*nval)
	do i=1,2*nval
	  conp1(i)=smlq1+float(i-1)*dc
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
c        call confill(q1,im,im,jm,xlen,ylen,
c     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(q1,im,im-1,jm,xlen,ylen,conp1,2*nval)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HE,0.,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1Hq,0.,1,1)
	call subber(1H1,1,ch1,0.)
c
	call keksymc(xlen+0.5,1.2*ylen,1.2*ch1,char3,0.,17,1)
c
	call coast(xlen,ylen)
c
	call plot(xlen+1.,0.,-3)
c
	bigh=bigq2-smlq2
	dc=bigh/float(2*nval)
	do i=1,2*nval
	  conp1(i)=smlq2+float(i-1)*dc
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
c        call confill(q2,im,im,jm,xlen,ylen,
c     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(q2,im,im-1,jm,xlen,ylen,conp1,2*nval)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HF,0.,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1Hq,0.,1,1)
	call subber(1H2,1,ch1,0.)
c
	call coast(xlen,ylen)
c
	call plot(-1.*(xlen+1.),-1.*(ylen+1.5),-3)
c
c	bigh=bigqg-smlqg
c	dc=bigh/float(2*nval)
c	do i=1,2*nval
c	  conp1(i)=smlqg+float(i-1)*dc
c	enddo
c
	bigh=bigqg
	dc=bigh/float(nval)
	do i=1,nval
	  conp(i)=float(i)*dc
	  conm(i)=-1.*bigh+float(i-1)*dc
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
        call confill(qgeos,im,im-1,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(qgeos,im,im-1,jm,xlen,ylen,conp,nval)
        call conrec(qgeos,im,im-1,jm,xlen,ylen,conm,nval)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HG,0.,1,1)
c
	call keksym(0.5*xlen,1.05*ylen,ch1,1Hq,0.,1,1)
	call subber(1h2,1,ch1,0.)
	call setfnt(29)
	call keksym(0.5*xlen,1.05*ylen+ch1,ch1,331,0.,-999,1)
	call setfnt(20)
c
	call coast(xlen,ylen)
c
	call plot((xlen+1.),0.,-3)
c
c	bigh=bigpb-smlpb
c	dc=bigh/float(2*nval)
c	do i=1,2*nval
c	  conp1(i)=smlpb+float(i-1)*dc
c	enddo
c
	bigh=bigpb
	dc=bigh/float(nval)
	do i=1,nval
	  conp(i)=float(i)*dc
	  conm(i)=-1.*bigh+float(i-1)*dc
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
        call confill(psib,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
c        call conrec(psib,im,im,jm,xlen,ylen,conp1,2*nval)
        call conrec(psib,im,im,jm,xlen,ylen,conp,nval)
        call conrec(psib,im,im,jm,xlen,ylen,conm,nval)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HH,0.,1,1)
c
	call grksym(0.5*xlen,1.05*ylen,ch1,47,0.,1,1)
	call subber(1hB,1,ch1,0.)
c
	call coast(xlen,ylen)
c
	call chopit(0.,0.)
C
        call plot(1.,9.,-3)
c
c......find the largest value of sver1 and sver2, etc and compute the
c      contour interval
c
	bigs1=-1.e20
	bigs2=-1.e20
	bigs=-1.e20
	bigwc=-1.e20
	do j=1,jm
	  do i=1,im
	   if(abs(sver1(i,j)).gt.bigs1)bigs1=abs(sver1(i,j))
	   if(abs(sver2(i,j)).gt.bigs2)bigs2=abs(sver2(i,j))
	   if(abs(sver(i,j)).gt.bigs)bigs=abs(sver(i,j))
	   if(abs(wcurl(i,j)).gt.bigwc)bigwc=abs(wcurl(i,j))
	  enddo
	enddo
c
	lscal=1
c
	bigwc=bigwc
c
	print *,'bigs1=',bigs1
	print *,'bigs2=',bigs2
	print *,'bigs=',bigs
	print *,'bigwc=',bigwc
c
c	bigh=bigs1
	bigh=2.*bigwc
c	dc=bigh/float(nval1)
	dc=0.4*w0/0.5
	do i=1,nval1
	  conp2(i)=float(i)*dc
c	  conm2(i)=-1.*bigh+float(i-1)*dc
	  conm2(i)=-1.*float(nval1)*dc+float(i-1)*dc
	enddo
c
        cfill(1)=-1.e20
        cfill(2)=conp2(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
        call confill(sver1,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(sver1,im,im,jm,xlen,ylen,conp2,nval1)
        call conrec(sver1,im,im,jm,xlen,ylen,conm2,nval1)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HI,0.,1,1)
c
c	call keksym(0.5*xlen,1.05*ylen+0.3,0.7*ch1,
c     &         28HZonally Integrated Northward,0.,28,1)
c	call keksym(0.5*xlen,1.05*ylen,0.7*ch1,
c     &       30HTransport from East in Layer 1,0.,30,1)
c
	call grksym(0.5*xlen,1.05*ylen,ch1,46,0.,1,1)
	call subber(1H1,1,ch1,0.)
c
	call keksymc(xlen+0.5,1.2*ylen,1.2*ch1,char3,0.,17,1)
c
	call coast(xlen,ylen)
c
	call plot(xlen+1.,0.,-3)
c
        cfill(1)=-1.e20
        cfill(2)=conp2(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
        call confill(sver2,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(sver2,im,im,jm,xlen,ylen,conp2,nval1)
        call conrec(sver2,im,im,jm,xlen,ylen,conm2,nval1)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HJ,0.,1,1)
c
c
c	call keksym(0.5*xlen,1.05*ylen+0.3,0.7*ch1,
c     &         28HZonally Integrated Northward,0.,28,1)
c	call keksym(0.5*xlen,1.05*ylen,0.7*ch1,
c     &       30HTransport from East in Layer 2,0.,30,1)
c
	call grksym(0.5*xlen,1.05*ylen,ch1,46,0.,1,1)
	call subber(1H2,1,ch1,0.)
c
	call coast(xlen,ylen)
C
        call plot(-1.*(xlen+1.),-1.*(ylen+1.5),-3)
c
        cfill(1)=-1.e20
        cfill(2)=conp2(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
        call confill(sver,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(sver,im,im,jm,xlen,ylen,conp2,nval1)
        call conrec(sver,im,im,jm,xlen,ylen,conm2,nval1)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HK,0.,1,1)
c
c
c	call keksym(0.5*xlen,1.05*ylen+0.3,0.7*ch1,
c     &         28HZonally Integrated Northward,0.,28,1)
c	call keksym(0.5*xlen,1.05*ylen,0.7*ch1,
c     &       26HTransport from East: Total,0.,26,1)
c
	call grksym(0.5*xlen,1.05*ylen,ch1,46,0.,1,0)
	call subber(1H1,1,ch1,0.)
	call keksym(999.,999.,ch1,1H+,0.,1,0)
	call grksym(999.,999.,ch1,46,0.,1,0)
	call subber(1H2,1,ch1,0.)
c
	call coast(xlen,ylen)
c
	call plot(xlen+1.,0.,-3)
c
        cfill(1)=-1.e20
        cfill(2)=conp2(1)
        cfill(3)=1.e20
        grylev(1)=1.0
        grylev(2)=1.0
        grylev(3)=0.75
c
        call confill(wcurl,im,im,jm,xlen,ylen,
     &    cfill,grylev,3,ioffpp,spvall)
c
        call conrec(wcurl,im,im,jm,xlen,ylen,conp2,nval1)
        call conrec(wcurl,im,im,jm,xlen,ylen,conm2,nval1)
c
        call border(xlen,ylen,0000,1111,1,1,2,1)
c
	call keksym(-0.5,1.05*ylen,ch3,1HL,0.,1,1)
c
c
c	call keksym(0.5*xlen,1.05*ylen+0.3,0.7*ch1,
c     &         18HZonally Integrated,0.,18,1)
c	call keksym(0.5*xlen,1.05*ylen,0.7*ch1,
c     &       16HWind Stress Curl,0.,16,1)
c
	call grksym(0.5*xlen,1.05*ylen,ch1,47,0.,1,1)
	call subber(1HS,1,ch1,0.)
c
	call coast(xlen,ylen)
c
	call chopit(0.,0.)
c
	lscal=0
c
c
	print *,'interior ',u1(25,25),u2(25,25)
	print *,'SE ',u1(42,6),u2(42,6)
c
	return
        end
c
	subroutine coast(xlen,ylen)
c
        PARAMETER(IM=51)
        PARAMETER(JM=51)
c
	real xlrr(4),ylrr(4)
c
c       west coast
c
         dlx=xlen/(float(im-1))
         dly=ylen/float(jm-1)
	 xlrr(1)=0.0
	 ylrr(1)=0.0
	 xlrr(2)=0.5*dlx
	 ylrr(2)=ylrr(1)
	 xlrr(3)=xlrr(2)
	 ylrr(3)=ylen
	 xlrr(4)=0.
	 ylrr(4)=ylrr(3)
         call filrgn(xlrr,ylrr,4,0.)
c
c       east coast
c
	 xlrr(1)=xlen-0.5*dlx
	 ylrr(1)=0.0
	 xlrr(2)=xlen
	 ylrr(2)=ylrr(1)
	 xlrr(3)=xlrr(2)
	 ylrr(3)=ylen
	 xlrr(4)=xlen-0.5*dlx
	 ylrr(4)=ylrr(3)
         call filrgn(xlrr,ylrr,4,0.)
c
	 xlrr(1)=0.0
	 ylrr(1)=0.0
	 xlrr(2)=xlen
	 ylrr(2)=ylrr(1)
	 xlrr(3)=xlrr(2)
	 ylrr(3)=0.5*dly
	 xlrr(4)=0.
	 ylrr(4)=ylrr(3)
         call filrgn(xlrr,ylrr,4,0.)
c
	 xlrr(1)=0.0
	 ylrr(1)=ylen-0.5*dly
	 xlrr(2)=xlen
	 ylrr(2)=ylrr(1)
	 xlrr(3)=xlrr(2)
	 ylrr(3)=ylen
	 xlrr(4)=0.
	 ylrr(4)=ylrr(3)
         call filrgn(xlrr,ylrr,4,0.)
c
	return
	end
c
