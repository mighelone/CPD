C  THIS IS THE CPDCP-NLG MODEL
C
C  THIS MODEL WAS DEVELOPED BY SANDIA NATIONAL LABORATORIES UNDER 
C  FWP 0709 FOR THE DEPARTMENT OF ENERGY'S PITTSBURGH ENERGY
C  TECHNOLOGY CENTER AND THE DOE DIVISION OF ENGINEERING AND GEOSCIENCES
C  THROUGH THE OFFICE OF BASIC ENERGY SCIENCES;
C  AND BY THE UNIVERSITY OF UTAH THROUGH FUNDING FROM 
C  THE ADVANCED COMBUSTION ENGINEERING RESEARCH CENTER (ACERC), WHICH 
C  IS PRINCIPALLY SPONSORED BY THE NATIONAL SCIENCE FOUNDATION, THE
C  STATE OF UTAH, AND BY A CONSORTIUM OF INDUSTRIAL COMPANIES.
C  THE CODE WILL NOT BE FORMALLY LICENSED.  NEITHER THE U.S. OR THE 
C  DOE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, 
C  OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, 
C  COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR 
C  PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD INFRINGE PRIVATELY 
C  OWNED RIGHTS.

C
C  The CPD model is intended to solve pyrolysis rate equations
C  based on a percolative bond-breaking scheme.  This version includes the 
C  flash distillation program to distinguish between tar and metaplast.
C  This program also includes a crosslinking scheme.  (January, 1991)
C
C
C  Most recent modifications to this model include (a) the nitrogen release
C  model of Perry, (b) the model of Genetti to break the light gas into 
C  species based on a correlation, and (C) slight modification to mdel
C  to account for the mass associated with c0.  
C  These modifications were made at BYU by Dominic Genetti in his 
C  M.S. Thesis work (1999) and Steve Perry in his Ph.D. work (1999).
C
C
      program CPDCP
C
C  This version is coupled with a solver for the particle energy and 
C   momentum equations.
C       UNITS ARE G,K,CM,S,CAL
C  PARAMETER INITIALIZATION
C
C   A BLOWING CORRECTION TO THE HEAT TRANSFER MODEL IS EMPLOYED.
C
C   Merrick heat capacity correlations are used
C   This version uses the heat capacity at the initial temperature as recommended
C   by Maloney, D. J., R. Sampath and J. W. Zondlo, Combustion and Flame, 116, 94-104 (1999).
C
C   Water is included through Omegaw
C
      IMPLICIT real*4 (A-H,O-Z)
      save
C  Input parameters are found in PERKIN
C  Outut parameters are found in PERKOT
      DIMENSION XT(50),TGC(50),zv(50),vpz(50)
      CHARACTER*80 A,vname,tgasin,output,infile,outgas,outnit,hrfile
      REAL KG,NU,KMJ1,KMJ2
      character*80 perkin,perkot,perkot2,perkot3
      character*1 ans
      dimension y(4),yp(4),ypp(4),ypred(4),tim(50),tem(50),gasf(9,10)
      dimension yygas(5),ffgas(5)
      real*4 l0,l,ma,kp,ft(35),mt(35),mgas,mtot
      real*4 mwchar,machar,mwcharold,mw1,mdel,mb,trateold,volc
      real*4 sp,smin,sigmdel,cscale,svar,cHR,cvisc,cext,sratio,sc
      integer*2 lib
C  intar = .true. calculates tar molecular weight distribution
C                 in subroutine perkp.
      logical intar,idiff,imolw,ipred,inside
C  ipred = .true. when on the predictor step      
      common/init/l0,c0,g0,ma,rba,finf,sig,siginv,nmax,pstar
      common/rate/ab,eb0,ebsig,ac,ec0,ecsig,ag,eg0,egsig,rg
      common/tarmw/press
      common/timer/tms,time
      common/spheat/yelem(5)
      common/nitmod/fnca0,fnca,an,en0,ensig
      common/lgmod/yf,xoc,yhc
      data xt/50*0./,tgc/50*0./,zv/50*0./,vpz/50*0./
      data tim/50*0./,tem/50*0./,ft/35*0./,mt/35*0./
      data y/4*0./,yp/4*0./,ypp/4*0./,ypred/4*0./
      data nmax/10/,nt/1/,intar/.false./,idiff/.false./
      data ftar0,fgas0,fchar0/0,0,1./,ip/30/,zero/0.0/
      data small/1.e-7/ftarold/0.0/
C  y1 = l       labile bridges
C  y2 = del     ends
C  y3 = C       char links
C  y4 = fnca    mass fraction of nitrogen in site
      data g0,g1,g2/0.,0.,0./
      DATA PI/3.14159/,DELHV/0./,PR/0.7/,delhw/-540./,cpw/1.0/
      data TMS,Xmm,V,FRAC,BLOW/4*0.0,1./
	  xwb=0.0
	  rtot=0.0
	  tratem=0.0
	  trateold=0.0
	  trcheck=0
	  volc=0.0
	  tpmax=0.0
C
C
      IX = 1
      IV = 1
      X = 0.0

C
C   input parameters
C
      print*,' Enter main input file'
C      read(*,1)infile
      read(*,*)infile
      print*,infile
      open(unit=21,name=infile,type='unknown')
      print*,' Enter input particle velocity file'
C      read (21,1)vname
      read (21,*)vname
      open(unit=25,name=vname,type='unknown')
	  read(25,*)nv
C      do 710 i=1,20
C         read(25,1)a
C         write(*,1)a
C         if(a(1:1).eq.' ')go to 711
C710    continue
C711    CONTINUE
      do 1001 i=1,nv
         read(25,*)zv(i),vpz(i)
C	 write(*,*)zv(i),vpz(i)
         zv(i) = zv(i)/10.
         write(*,*)zv(i),vpz(i)
1001  continue
C      print*,' Error: Too many velocity data points'
C1112  continue
C      nv = i-1
C      Enter input gas temperature file
      read (21,*)tgasin
C      read (21,1)tgasin
      open(unit=22,name=tgasin,type='unknown')
C   read in gas temperature profile
	  read(22,*)nx
C      do 10 i=1,20
C         read(22,1)a
C1        format(a)
C         if(a(1:1).eq.' ')go to 11
C10    continue
C11    CONTINUE
      do 1111 i=1,nx
         read(22,*)xt(i),tgc(i)
         xt(i) = xt(i)/10.
	 write(*,*)xt(i),tgc(i)
1111  continue
C         print*,'too many data points for gas temperature'
C         stop
C1002  continue
C      nx = i-1
      print*,' Enter name for output file'
      read (21,*)output
      write(*,*)output
      if(output(1:1).eq.' ')output(1:11)='ubcp.out'
      open(unit=1,name=output,type='unknown')
            print*,' Enter name for nitrogen release output file'
      read (21,*)outnit
      write(*,*)outnit
      open(unit=23,name=outnit,type='unknown')
      read (21,*)outgas
      write(*,*)outgas
      open(unit=27,name=outgas,type='unknown')
C
      read(21,*)TIMAX
C      read(21,*)TG0
C      read(21,*)VG0  !cm/s
       TG0=tgc(1)
       VG0=vpz(1)
       write(*,*)Tg0,vg0
      read(21,*)RHOP !G/CM**3
      read(21,*)DP !CM
C      read(21,*)Tau0 !adjustable to match particle temperatures (or mass release)
C	  Rhop=rhop*tau0 !this is one way to implement the correction
      read(21,*)Swell !fraction final diameter increase
	  swell=swell-1.0
	  if((swell+1.0).gt.0.5) then
		sratio=swell+1.0 !use experimental swelling ratio
	  else
	    swell=0.0
	    sratio=0.0 !triggers use of swelling correlation by Randy Shurtz
	  endif
	  sc=sratio
      read(21,*)DELHV !CAL/G (- MEANS ENDOTHERMIC)
      read(21,*)OMEGAW
      read(21,*)OMEGAA 
      read(21,*)EMIS 
      read(21,*)TWALL 
      read(21,*)dt,dtmax,iprint
C   coal input parameters
C  read in measured p0 from NMR data
      read(21,*)p0
      print*,'p0',p0
C  estimate c0
      read(21,*)c0
C  read in sigma+1 (coordination number) from NMR data
      read(21,*)sigp1
C  read in real MW/Cluster from NMR data (includes side chains)
      read(21,*)mw1
C  read in real mdel from NMR data
      read(21,*)mdel         
       print*,'p0,c0,sigp1,mw1,mdel',p0,c0,sigp1,mw1,mdel
C  kinetic parameters
      read(21,*)ab
      read(21,*)eb0
      read(21,*)ebsig
      read(21,*)ac
      read(21,*)ec0
      read(21,*)ag
      read(21,*)eg0
      read(21,*)egsig
      read(21,*)acr
      read(21,*)ecr
C  kinetic parameters for light gas nitrogen release from char:
C  A. release by reaction with free radicals (fast)
      read(21,*)arad
      read(21,*)erad
      read(21,*)fstable
C  B. release by high T decomposition (slow)
      read(21,*)an
      read(21,*)en0
      read(21,*)ensig
C  pressure in atmospheres
      read(21,*)press
C  read in DAF composition of the coal (CHNOS)
      do 12 i=1,5
        read(21,*)yelem(i)
12    continue
	  read(21,*)vastm
	  read(21,*)cpswitch
	  read(21,*)hrfile
C      READ(21,*)THTR !These temperatures are for a Sandia drop-tube (CDL);
C      READ(21,*)TTUB !setting them to the the wall temperature makes the view factors irrelavent
       THTR=TWALL
       TTUB=TWALL
	  open(unit=29,name=hrfile,type='unknown')
	  write(29,*) "Pressure (atm) ",press	  
	  write(29,*) "Particle Diameter (microns) ",(DP*1.e4)
      write(29,*) "Moisture Fraction ", OMEGAW
	  write(29,*) "Dry Ash Fraction ",OMEGAA
      write(29,*) "CHNOS daf Percentages ",
     & yelem(1),yelem(2),yelem(3),yelem(4),yelem(5)
	  write(29,*)"ASTM daf Volatiles Percent",vastm
C    Predict NMR parameters from Genetti's correlation if not provided, convert composition from percentages to fractions
      if((p0.lt.0.0).or.(mdel.lt.0.0).or.
     & (mw1.lt.0.0).or.(sigp1.lt.0.0).or.(c0.lt.0.0)) then
        mdel=(421.957)+(-8.64692)*yelem(1)
     & +(0.0463894)*yelem(1)*yelem(1)+(-8.47272)*yelem(2)
     & +(1.18173)*yelem(2)*yelem(2)+(1.15366)*yelem(4)
     & +(-0.0434024)*yelem(4)*yelem(4)
     & +(0.556772)*vastm+(-0.00654575)*vastm*vastm
        mw1=(1301.41)+(16.3879)*yelem(1)+(-0.187493)*yelem(1)*yelem(1)
     & +(-454.773)*yelem(2)+(51.7109)*yelem(2)*yelem(2)
     & +(-10.072)*yelem(4)+(0.0760827)*yelem(4)*yelem(4)
     & +(1.36022)*vastm+(-0.0313561)*vastm*vastm
        p0=(0.489809)+(-0.00981566)*yelem(1)
     & +(0.000133046)*yelem(1)*yelem(1)+(0.155483)*yelem(2)
     & +(-0.0243873)*yelem(2)*yelem(2)+(0.00705248)*yelem(4)
     & +(0.000219163)*yelem(4)*yelem(4)+(-0.0110498)*vastm
     & +(0.000100939)*vastm*vastm
        sigp1=(-52.1054)+(1.63872)*yelem(1)
     & +(-0.0107548)*yelem(1)*yelem(1)+(-1.23688)*yelem(2)
     & +(0.0931937)*yelem(2)*yelem(2)+(-0.165673)*yelem(4)
     & +(0.00409556)*yelem(4)*yelem(4)+(0.00926097)*vastm
     & +(-8.26717E-05)*vastm*vastm
        if (yelem(1).gt.85.9) then
          c0=min((0.1183*yelem(1)-10.16),0.36)
        elseif (yelem(4).gt.12.5) then
          c0=min((0.014*yelem(4)-0.175),0.15)
	else
	  c0=0.0
	endif
      endif
C	calculate terms in swelling model as developed by Randy Shurtz
      write(29,*)"c0 ",c0
      write(29,*)"p0 ",p0
      write(29,*)"(sigma+1) ",sigp1
      write(29,*)"MW ",mw1
      write(29,*)"Mdel ",mdel
	  smin=((1.-0.01*vastm)*(1.-omegaa)+omegaa)**(0.3333333333)
	  sigmdel=sigp1/mdel
	  if((sigmdel.le.0.01831).or.(sigmdel.ge.0.30084)) then
          svar=0.0
      elseif (sigmdel.le.0.2067) then
          svar=1.6852*sigmdel-0.030856
      else
          svar=-3.3722*sigmdel+1.0145
      endif
	  write(29,*)"(sigma+1)/mdel ",sigmdel
	  write(29,*)"smin ",smin
	  write(29,*)"svar ",svar
	  if ((sigmdel.le.0.1062).or.(sigmdel.ge.0.2541)) then
		cHR=0.0
	  else
	    cHR=-191.293*sigmdel*sigmdel+68.921*sigmdel-5.161
	  endif
   	  write(29,*)"cHR ",cHR
		  cvisc=1.0
		  cext=1.0
	  cscale=0.0
      if(press.le.1.0) then
		  sp=1.0
	  else
		  cvisc=7.77
		  cext=3.47
		  cscale=38.89*sigp1-167.10
		  if(sigp1.lt.4.297) then
			cscale = 0.0
		  endif
		  sp=1.0+cscale*((log(press))**cvisc)/(press**cext)
	  endif
	  write(29,*)"cvisc ",cvisc
	  write(29,*)"cext ",cext
	  write(29,*)"cscale ",cscale
	  write(29,*)"sp (swelling pressure factor) ",sp
      yelem(1)=yelem(1)*.01
      yelem(2)=yelem(2)*.01
      yelem(3)=yelem(3)*.01
      yelem(4)=yelem(4)*.01
      yelem(5)=yelem(5)*.01
       print*,'p0,c0,sigp1,mw1,mdel',p0,c0,sigp1,mw1,mdel
       print*,'C,H,N,O,S fracs',yelem(1),yelem(2),yelem(3),yelem(4),
     & yelem(5)
C  save initial char NMR parameters as coal NMR parameters
C  (char parameters calculated independent of those using
C  empirical correlation for mdel)
      mwchar = mw1
      mwcharold = mwchar
      machar = mwchar-sigp1*mdel
C  adjust mdel to correct for c0 (Steve Perry, May 99)
      mdel = mdel/(1.0-c0)
C  empirical correlation to allow a small portion of alpha-carbon to stay with
C  the aromatic cluster
      mdel = mdel-7
C  now calculate other chemical structure coefficients
      l0 = p0 - c0
      mb = 2.*mdel
      ma = mw1-sigp1*mdel
      sma = mw1-sigp1*(mdel+7)
      sig = sigp1-1
      finf = 1./(2.*ma/(mb*sigp1*(1-c0))+1)
      rba = mb/ma
      print*,'rba=',rba
      fnit = yelem(3)
      fnt = fnit
      fnca0 = fnit*mw1/machar
      Grav = 980.
      dp0 = dp
C   initialize variables
      y(1) = l0
      y(2) = 2.*(1.-c0-l0)
      y(3) = c0
      y(4) = fnca0
      aind0 = l0 + (1.-c0-l0)
      siginv = 1./sig
      pstar = 0.5*siginv
      ynhcn = 0.0
      yntar = 0.0
      yytar = 0.0
      yf = 0.0
      inside=.true.
C   START calculation DO LOOP
      Rg = 1.987                       !CAL/GMOLE K
      iii = 0
      fcross = 0.0
C   CALCULATE INITIAL PARTICLE VELOCITY
      CALL PROPS(TG0,TG0,RHOG,UG,KG,CPG,diffw,press)
C      VP = VG0 - DP**2*RHOP*Grav/(18*UG)
      VP = Vg0
      VG = VG0
C   START DO LOOP
      V = 0.
      Tp = tg0
      Tg = tg0
      MOVE = 1
      TIME = 0.0
C  for now, assume that the apparent density is indicative of the as 
C  received coal.
      ALPHAP = (4./3.)*PI*(DP/2.)**3*RHOP
      alfa = alphap*omegaa
      alfw = alphap*omegaw
      ALFC = ALPHAP*(1-OMEGAA-Omegaw)
      alfc0 = alfc
      ALFV = 0.
      ALPHA = Alphap
      FRAC = ALFC/ALPHAP
      W = ALFW/ALPHAP*100.
      AP = PI*DP**2
      RE = RHOG*ABS(VG-VP)*DP/UG
      NU = 2. + 0.6*RE**.5*PR**.333
      H = NU*KG/DP
      SIGMA = 1.335E-12               !CAL/S CM**2 K**4
      HEAT = 0.0
      v = 0.0
C   calculate O/C and H/C ratios for light gas model
      xoc = (yelem(4)/16)/(yelem(1)/12)
      yhc = yelem(2)/(yelem(1)/12)
C
      R = 1.987                       !CAL/GMOLE K
      dpm = dp*1.e4
      write(1,5000)dpm
5000  format('C  Coal Calculations'/
     x'C  Initial Particle Diameter = ',f6.2,' microns'/'C')
      write(6,5010)
5010  format('   Tim (ms)   Dis (mm)   Tp (K)    Tg (K)',
     x '       Vol.          Char       Blow')
      write(1,202)
202   format('C time(ms)    dis(mm)  Tp(K)      Tg(K)        % Vol.',
     x'        fchar        fcross     ftar    fmet      Trate')
      write(23,204)
204   format('C time(ms)   temp(K)    MW        Nsite    Nchar   fnchar'
     x,'     fntar     fnhcn',
     x'     fntot')
      write(27,205)
205   format('C time(ms)      fgas     fh20     fco2     fch4',
     x'     fco     fother    yh20     yco2     ych4',
     x'     yco     yother    Xgas')
        call heatcp(tp,cpc)
		IF (CPSWITCH.GT.0.0) THEN
			CPC=CPC*CPSWITCH
		ENDIF
		cpcold=cpc
        call heatap(tp,cpa)
      fntar = 0.0
C        CALL HEATWP(tp,cpw)
        cp = (alfc*cpc + alfa*cpa + alfw*cpw)/(alpha)
            WRITE(6,1000)TMS,Xmm,Tp,TG0,v,fchar,BLOW
            write(1,1000)TMS,Xmm,Tp,TG,v,fchar,fcross,
     x     ftar,fmet,trate
      DO 100 I=1,100000 !goes to continue statement 100, exits to statement 3333

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C PREDICTOR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        fvolold = fvol
        fcrossold = fcross
        fmetold = fmet
        XP = X + VP*DT
C   CALCULATE GAS VELOCITY FROM Thermocouple Measurements
        IF(XP.LE.XT(nx))THEN
          TG = (XT(IX+1)-XP)/(XT(IX+1)-XT(IX))*(TGC(IX)-TGC(IX+1)) 
     X       + TGC(IX+1)
        ELSE
          PRINT*,'REACHED END OF GAS TEMPERATURE CORRELATION'
           if(inside.eqv..false.)then
             print*,'!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!'
             print*,'O/C and H/C ratios are outside the  '
             print*,'bounds of the library coals. '
             print*,'Estimation of Light gas distribution '
             print*,'is based on library coal no.',lib
           endif
          go to 3333
        ENDIF
        if(tg.gt.4000.)then
           print*,'     >>>>   WARNING    <<<<<'
           print*,' gas temperature too high----',tg
        endif
C interpolate to get particle velocity
        IF(XP.LE.zv(nv))THEN
          vpp = (zv(iv+1)-xP)/(zv(iv+1)-zv(iv))*(vpz(iv)-vpz(iv+1)) 
     X       + vpz(iv+1)
        ELSE
          PRINT*,'REACHED END OF PARTICLE VELOCITY CORRELATION'
          go to 3333
        ENDIF
C  CALCULATE GAS PROPERTIES
        CALL PROPS(TG,Tp,RHOG,UG,KG,CPG,diffw,press)
C   INTEGRATE PARTICLE VELOCITY
C       Use this section if gas velocity .NE. particle velocity
C        VLAG = (VP-VG)
C        RE = RHOG*DP*ABS(VLAG)
C        RE = RE/UG
C        SUMP = -18.*UG*VLAG*(1.+.15*RE**.687)/(DP**2*RHOP)+Grav
C        SUMP = -18.*UG*VLAG/(DP**2*RHOP)+Grav
C        VPP = VP + SUMP*DT
          RE = 0.
C ENERGY EQUATION
C
C--   CONVECTION
        NU = 2. + .6*RE**.5*PR**.333
        B = CPG*(rtot)/(2.*PI*DP*KG)
        IF(B.GE.1.E-4)THEN
          BLOW = B/(EXP(B)-1)
        ELSE
          BLOW = 1.0
        ENDIF
        H = BLOW*NU*KG/DP
        QCONV = H*AP*(TG-Tp)
C--   mass transfer
        if(alfw.gt.0.)then
          Bw = (rtot)/(2.*PI*DP*diffw*rhog)
          IF(Bw.GE.1.E-4)THEN
            BLOwW = Bw/(EXP(Bw)-1)
          ELSE
            BLOwW = 1.0
          ENDIF
       else
          bloww = 1.0
       endif
C--    RADIATION
        Z = XP
        RAD = 2.54 ! 1" radius burner or injector probe?
        if(z.gt.0.)then
          R2 = RAD/Z
          R3 = .125*2.54/Z !1/8" radius feeder tube?
        else
          r2 = 1.e12
          r3 = 1.e12
        endif
        f12 = .5*(1.-1./sqrt(1.+r2**2)) !view factor to burner/drop tube injector
        f13 = .5*(1.-1./sqrt(1.+r3**2)) !view factor to some other tube that is the same distance; feeder tube?
        FAC = AP*SIGMA*EMIS
        QTOP = (f12-F13)*(THTR**4-Tp**4)*FAC
        QTUBE = F13*(TTUB**4-Tp**4)*FAC
        QWALL = (1.-F12)*(TWALL**4-Tp**4)*FAC
        QRAD = QWALL + QTOP + QTUBE
C--    Water Evaporation Rate
        if(alfw.gt.0.)then
          xw0 = exp(18.3036-3816.44/(Tp-46.13))/760.
          xw0 = min(xw0/press,1.) !It should be the total pressure to get the mole fraction, not just 1 atm
          xw0 = max(xw0,0.)
          rwp = bloww*2.*rhog*diffw*pi*dp*(xw0-xwb)/(1.-xw0)
        else
          rwp = 0.
        endif
C--  Coal pyrolysis rate
        call perks(y,ypp,Tp)
C  free radical light gas nitrogen release mechanism
      if((mw1-mwchar)/mw1.GT.fstable)ypp(4)=ypp(4)-y(4)*arad*
     x exp(-erad/rg/Tp)*(mwcharold-mwchar)/mwchar*
     x machar/mwchar/dt
C COMPONENT MASS CONSERVATION
        do 50 j=1,4
           ypred(j) = y(j) + dt*ypp(j)
           ypred(j) = max(ypred(j),zero)
50      continue
C  crosslinking rate
        fracr = 1.
        if(fmetold.gt.small.and.acr.gt.0.)then
           ratecr = acr*exp(-ecr/rg/tp)*fmetold*dt
           fracr = 1.-ratecr/fmetold
           fmet = fmetold-ratecr
           fcross = fcrossold+ratecr
           if(fmet.lt.0.)then
             fcross = fcrossold-ratecr+fmet+ratecr
             fmet = 0.0
             fracr = 0.
           endif
         endif
C  get product distribution from ypred array
        if(ypred(1).gt.small)intar = .true.
        call perkp(ypred,mgas,ftar,ftart,fgas,fchart,ft,mt,intar)
        intar = .false.
          gasmw = rba*ma/2.
C  flash distillation
          if(fgas.ge.1.e-5)then
           imolw = .false.
           if(mod(iii,iprint).eq.0)imolw=.true.
           ipred = .true.
           call flash(fgas,gasmw,ft,mt,fracr,ftar,fmet,
     x                tp,press,nmax,imolw,ipred)
          elseif(fgas.lt.1.e-5)then
           fmet = ftart
           ftar = 0.
          endif
          intar = .false.
          fvol = fgas+ftar
          fchar = 1.-fvol
          dvdt = (fvol-fvolold)/dt*alfc0
C  done with pyrolysis rate!
        rtotp = dvdt + rwp
	delhw= 6.7845E-08*(tp**4) - 1.1123E-04*(tp**3)
     x  + 6.8595E-02*(tp**2) - 1.8093E+01*tp + 1.1229E+03
        HEAT = dvdt*delhv + rwp*delhw
C calculate heat capacity
	fmetmaxp=max(fmetmaxp,fmet)
	if (CPSWITCH.ge.0.0) then
C        call heatcp(tp,cpc)
C        	cpc=(fmet/fmetmaxp)*cpcold+(1.0-fmet/fmetmaxp)*cpc
		cpc=cpcold
	else
		CALL HEATCP(TP,CPC)
	endif
        call heatap(tp,cpa)
        cp = (alfc*cpc + alfa*cpa + alfw*cpw)/(alpha)
C
        TRATEP = (QCONV+QRAD+HEAT)/(ALPHA*CP)
        TPRED = Tp + TRATEP*DT
C        print*,'<trate> qconv,qrad,heat,alpha,cp,tratep,tp,tg',
C     x        qconv,qrad,heat,alpha,cp,tratep,tp,tg
C COMPONENT MASS CONSERVATION
        alfwp = alfw - rwp*dt
        ALFVP = fvol*alfc0
        ALFcP = fchar*alfc0
        ALFCP = MAX(ALFCP,0.0)
        alfwp = max(alfwp,0.)
        ALPHA = (ALFCP + alfa + alfwp)
        omegaa = alfa/alpha
C  particle swelling
        l = ypred(1)
        dp = dp0*(1.+swell*(1.-l/l0))
        AP = PI*DP**2
C PARTICLE DENSITY CHANGES DURING DEVOLATILIZATION
        RHOP = ALPHA/((4./3.)*PI*(DP/2)**3)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C CORRECTOR
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        X = X + 0.5*(VPP+VP)*DT
C   Interpolate to get gas temperature
        IF(X.LE.XT(nx))THEN
          TG = (XT(IX+1)-X)/(XT(IX+1)-XT(IX))*
     X       (TGC(IX)-TGC(IX+1)) + TGC(IX+1)
        ELSE
          PRINT*,' REACHED END OF GAS TEMPERATURE CORRELATION'
          go to 3333
        ENDIF
C interpolate to get particle velocity
        IF(X.LE.zv(nv))THEN
          vp = (zv(iv+1)-x)/(zv(iv+1)-zv(iv))*(vpz(iv)-vpz(iv+1)) 
     X       + vpz(iv+1)
        ELSE
          PRINT*,'REACHED END OF PARTICLE VELOCITY CORRELATION'
          go to 3333
        ENDIF
C   CALCULATE GAS PROPERTIES
        CALL PROPS(TG,TPRED,RHOG,UG,KG,CPG,diffw,press)
C       INTEGRATE PARTICLE VELOCITY
C       Use this section if gas velocity .NE. particle velocity
C        VLAG = (VP-VG)
C        RE = rhog*DP*ABS(VLAG)
C        RE = RE/UG
C        SUM = -18.*UG*VLAG*(1.+.15*RE**.687)/(DP**2*RHOP)+
C        SUM = -18.*UG*VLAG/(DP**2*RHOP)+Grav
C        VP = VP + 0.5*(SUMP+SUM)*DT
          RE = 0
C ENERGY EQUATION
C
C--    CONVECTION
        NU = 2. + 0.6*RE**.5*PR**.333
        B = CPG*(rtotp)/(2.*PI*DP*KG)
        IF(B.GE.1.E-4)THEN
          BLOW = B/(EXP(B)-1)
        ELSE
          BLOW = 1.0
        ENDIF
        H = BLOW*NU*KG/DP
        QCONV = H*AP*(TG-TPRED)
C--   mass transfer
        if(alfw.gt.0.)then
          Bw = (rtotp)/(2.*PI*DP*diffw*rhog)
          IF(Bw.GE.1.E-4)THEN
            BLOwW = Bw/(EXP(Bw)-1)
          ELSE
            BLOwW = 1.0
          ENDIF
       else
          bloww = 1.
       endif
C--    RADIATION
        Z = X
        RAD = 2.54
        if(z.gt.0.)then
          R2 = RAD/Z
          R3 = .125*2.54/Z
        else
          r2 = 1.e12
          r3 = 1.e12
        endif
        f12 = .5*(1.-1./sqrt(1.+r2**2))
        f13 = .5*(1.-1./sqrt(1.+r3**2))
        FAC = AP*SIGMA*EMIS
        QTOP = (f12-F13)*(THTR**4-TPRED**4)*FAC
        QTUBE = F13*(TTUB**4-TPRED**4)*FAC
        QWALL = (1.-F12)*(TWALL**4-TPRED**4)*FAC
        QRAD = QWALL + QTOP + QTUBE
C--    water evaporation rate
        if(alfw.gt.0.)then
          xw0 = exp(18.3036-3816.44/(Tpred-46.13))/760.
          xw0 = min(xw0/press,1.)
          xw0 = max(xw0,0.)
          rw = bloww*2.*rhog*diffw*pi*dp*(xw0-xwb)/(1.-xw0)
        else
          rw = 0.
        endif
C--  coal pyrolysis rate
        Call perks(ypred,yp,Tpred)
C   time step control
           if(y(1).gt.5.e-3)then
              dy1 = dt*0.5*(yp(1)+ypp(1))
           else
              dy1 = dt*0.5*(yp(3)+ypp(3))
           endif
           if(abs(dy1).lt.0.001)then
              dt = dt*2.
              if(dt.lt.dtmax)print*,'at time=',time, 
     x                        '  dt changed to ',dt
           elseif(abs(dy1).gt.0.02)then
              dt = 0.01/abs(dy1)*dt
              print*,'at time=',time, 'dt changed to ',dt
           endif
           dt = min(dt,dtmax)
C  free radical light gas nitrogen release mechanism
      if((mw1-mwchar)/mw1.GT.fstable)yp(4)=yp(4)-y(4)*arad*
     x exp(-erad/rg/Tpred)*(mwcharold-mwchar)/mwchar*
     x machar/mwchar/dt
C COMPONENT MASS CONSERVATION
        do 70 j=1,4
           y(j) = y(j) + dt*0.5*(yp(j)+ypp(j))
           y(j) = max(zero,y(j))
70      continue
C  update current and old mwchar
      mwcharold = mwchar
      g = 2.*(1.-y(1)-c0)-y(2)
      mwchar = mw1-g*mdel*sigp1/2.
C  crosslinking rate
        fracr = 1.
        if(fmetold.gt.small.and.acr.gt.0.)then
           ratecr = acr*exp(-ecr/rg/tpred)*fmetold*dt
           fracr = 1.-ratecr/fmetold
           fmet = fmetold-ratecr
           fcross = fcrossold+ratecr
           if(fmet.lt.0.)then
             fcross = fcrossold-ratecr+fmet+ratecr
             fmet = 0.0
             fracr = 0.
           endif
         endif
C  get product distribution from y array
        if(y(1).gt.small)intar = .true.
        call perkp(y,mgas,ftar,ftart,fgas,fchart,ft,mt,intar)
        intar = .false.
          gasmw = rba*ma/2.
C  flash distillation
          if(fgas.ge.small)then
           imolw = .false.
           if(mod(iii,iprint).eq.0)imolw=.true.
           ipred = .false.
           call flash(fgas,gasmw,ft,mt,fracr,ftar,fmet,
     x                tpred,press,nmax,imolw,ipred)
          elseif(fgas.lt.1.e-5)then
           fmet = ftart
           ftar = 0.
          endif
          intar = .false.
          fvol = fgas+ftar
          fchar = 1.-fvol
          dvdt = (fvol-fvolold)/dt*alfc0
C  done with pyrolysis rate!
C        rtot = dvdt + rwp !this looks wrong
        rtot = dvdt + rw
       	delhw= 6.7845E-08*(tpred**4) - 1.1123E-04*(tpred**3)
     x  + 6.8595E-02*(tpred**2) - 1.8093E+01*tpred + 1.1229E+03
C        HEAT = dvdt*delhv + rwp*delhw
        HEAT = dvdt*delhv + rw*delhw
C
C calculate heat capacity
	fmetmaxc=max(fmetmaxc,fmet)
	if (CPSWITCH.ge.0.0) then
C        call heatcp(tpred,cpc)
C		cpc=(fmet/fmetmaxc)*cpcold+(1.0-fmet/fmetmaxc)*cpc
		cpc=cpcold
	else
		CALL HEATCP(TPRED,CPC)
	endif
        call heatap(tpred,cpa)
        cp = (alfcp*cpc + alfa*cpa + alfw*cpw)/(alpha)
C
        if (trate.gt.tratem) trateold=tratem
        TRATE = (QCONV+QRAD+HEAT)/(ALPHA*CP)
	TRATEM=max(trate,tratem)
        Tp = Tp + 0.5*(TRATE+TRATEP)*DT
		tpmax=max(tp,tpmax)
C        print*,'<trate> qconv,qrad,heat,alpha,cp,trate,tp,tg',
C     x        qconv,qrad,heat,alpha,cp,trate,tp,tg
        Alfw = alfw - rw*dt
        ALFV = fvol*alfc0
        ALFC = fchar*alfc0
        ALFC = MAX(ALFC,0.0)
        ALFw = MAX(ALFw,0.0)
        ALPHA = (ALFC + alfa + alfw)
        omegaa = alfa/alpha
C particle diameter changes due to swelling during devolatilization
	volc=fvol-fvolold
        if((alfw.eq.0.0).and.(fvol.ge.1.e-7)) then !check for true maximum in heating rate
	  if ((sratio.lt.0.5).and.(tratem.gt.trate)) then
	     sratio=smin+svar*sp*((5.79e4/tratem)**cHR)
	     write(29,*) "Max HR Time ",(time+dt)*1000.
             write(29,*)"HR for Swelling Correlation ",tratem
             write(29,*)"Correlated Swelling Ratio ",sratio
             swell=sratio-1.0 !use swelling correlation from Randy Shurtz
			 trateold=tratem
	  elseif ((sc.lt.0.5).and.(trateold.ne.tratem)) then
	     if (((volc).gt.1.e-6).and.(tratem.gt.trate)) then
	       sratio=smin+svar*sp*((5.79e4/tratem)**cHR)
	       write(29,*) "Max HR Time ",(time+dt)*1000.
               write(29,*)"HR for Swelling Correlation ",tratem
               write(29,*)"Correlated Swelling Ratio ",sratio
	       swell=sratio-1.0 !update swelling ratio with new maximum heating rate
		   trateold=tratem
	     endif
	   endif
	endif
        l = y(1)
        fnca = y(4)
        del2 = y(2)/2.
        aind = del2 + l
        dp = dp0*(1.+swell*(1.-l/l0))
        AP = PI*DP**2
C     
C  NITROGEN RELEASE CALCULATIONS
C
C  calculate tar nitrogen yield
C  (assumes tar is released before light gas and HCN)
      yntar = yntar + (ftar-ftarold)*fnt
      ftarold = ftar
C  nitrogen content of char and metaplast
      fnt = y(4)*machar/mwchar
C  nitrogen remaining in char
      ynchar = fchar*fnt
C  fraction of original nitrogen remaining in char
      fnchar = ynchar/fnit
C  fraction of original nitrogen released as tar
      fntar = yntar/fnit
C  fraction of original nitrogen released as light gas (diff.)
      fnhcn = 1.0 - fnchar - fntar
C  total fractional release of nitrogen
      fntot = (fnit - fnt*fchar)/fnit
      
    
C  DISTRIBUTE LIGHT GAS INTO H20, CO2, CO, CH4, & other HC's

C  yf is a CPD indicator of the fraction of total light gas
C   that has been released. The look up table on light gas 
C   composition is based on yf. 
      yf = 1-aind/aind0
       
      call lightgas(yygas,inside,lib)
      
C  calculate fraction of total mass release that is h2o, co2, ch4,
C   co, and other light gases

       do 218 ik=1,5
         ffgas(ik)=fgas*yygas(ik)
218    continue 
C
C PARTICLE DENSITY CHANGES DURING DEVOLATILIZATION
        RHOP = ALPHA/((4./3.)*PI*(DP/2)**3)
C
        FRAC = ALFC/ALPHAP
        omegawf = alfaw/alpha
        w = (alfw/alphap)*100.
        TIME = TIME + DT
        tms = time*1000.
         IF(MOD(I,iprint).EQ.0)THEN
          xmm = x*10.
          v = fvol*100.
          write(1,1000)TMS,Xmm,Tp,TG,v,fchar,fcross,ftar,fmet,trate
          write(23,1021)tms,tp,mwchar,y(4),fnt,fnchar,fntar,fnhcn,fntot
          write(27,1022)tms,fgas,ffgas(1),ffgas(2),ffgas(3),ffgas(4),
     x                ffgas(5),yygas(1),yygas(2),yygas(3),yygas(4),
     x                yygas(5),yf
            if(mod(i,iprint).eq.0)then
               WRITE(6,1000)TMS,Xmm,Tp,TG,v,fchar,BLOW
            endif
         ENDIF
1000    FORMAT(' ',F9.3,2X,F9.4,2X,F6.0,5X,F6.0,5X,1PE10.3,5X,1PE10.3
     X         ,3X,0PF6.3,3X,0PF6.2,3X,0PF6.3,3x,1pe10.3)
1021    format(' ',2(1pe10.3),0pf8.2,6(0pf10.6))
1022    format(' ',1pe10.3,0pf9.4,11(0pf9.4))
        IF(TIME.GE.timax)go to 3333
C   CHECK TO SEE IF INTERPOLATION FACTORS NEED UPDATE
        IF(X.GT.XT(IX+1))THEN
          IX = IX + 1
          IF(IX.GE.50)go to 3333
        ENDIF
        IF(X.GT.zv(Iv+1))THEN
          Iv = Iv + 1
          IF(Iv.GE.50)go to 3333
        ENDIF
100   CONTINUE
3333  continue
            v = fvol*100.
            xmm = x*10.
            WRITE(6,1000)TMS,Xmm,Tp,TG,v,fchar,BLOW
            write(1,1000)TMS,Xmm,Tp,TG,v,fchar,fcross,
     x           ftar,fmet,trate
	   write(29,*)"Final heat capacity (cal/gm-K) ",cpc
	   write(29,*)"Maximum particle temperature (K) ",tpmax
	   write(29,*)"Final Volatiles yield (%daf) ",v
	   write(29,*)"Maximum Heating Rate ", TRATEM
	   write(29,*)"Swelling Ratio Used ",(swell+1.0)
      STOP
      END
C
      SUBROUTINE PROPS(TG,TP,RHOG,UG,KG,CPG,diffw,press)
C   CALCULATION OF GAS PROPERTIES OF N2 AT ATMOSPHERIC PRESSURE; upgraded for pressure by Randy Shurtz 2010
C   FROM HOLMAN, P. 505
      PARAMETER (NDAT=10)
      REAL TG0(NDAT),UG0(NDAT),KG,KG0(NDAT)
      DATA TG0/300,400,500,600,700,800,900,1000,1100,1200/
      DATA RHO/.3412E-3/,TRHO/1000./
      DATA UG0/1.784E-4,2.198E-4,2.57E-4,2.911E-4,3.213E-4,3.484E-4,
     X3.749E-4,4.E-4,4.228E-4,4.450E-4/
C  UNITS ON KG0 ARE W/M/C, AND MUST BE DIVIDED BY 418.4 TO GET CAL/CM/S/C
      DATA KG0/.0262,.03335,.03984,.0458,.0512,.0561,.0607,.0648,.0685,
     X.07184/
C  HEAT CAPACITY DATA
      DATA CPA,CPB,CPC,CPD/6.529,.1488e-2,-0.02771e-5,0./
C   CALCULATE GAS DENSITY FROM IDEAL GAS LAW
      RHOG = RHO*TRHO*press/TG
C   T = FILM TEMPERATURE (TP+TG)/2
      T = (TP+TG)/2.
C   CALCULATE CPG FROM HEAT CAPACITY EXPRESSION IN Himmelblau
      CPG = (CPA+CPB*T+CPC*T**2)/28.
C   CALCULATE UG AND KG FROM INTERPOLATION OF TABLE VALUES FROM HOLMAN
C   FIND INTERVAL WHERE TEMPERATURE LIES.  
      IF(T.GT.1200.)THEN
         UG = UG0(NDAT)*(T/TG0(NDAT))**.67
         KG = KG0(NDAT)*(T/TG0(NDAT))**.58/418.4
         RETURN
      ENDIF
      J=1
      DO 10 I=1,9
        IF(T.GT.TG0(I))J=J+1
10    CONTINUE
      FAC = (TG0(J)-T)/(TG0(J)-TG0(J+1))
      UG = -FAC*(UG0(J)-UG0(J+1)) + UG0(J)
      KG = (-FAC*(KG0(J)-KG0(J+1)) + KG0(J))/418.4
C      p = 1.
      mw = 28.
C      diffw = 1.905e-5*t**1.67/p
	   diffw = 4.3991e-7*t**2.334/press !Randy Shurtz 2010 from Bird, Stewart and Lightfoot 2002, H2O and N2
      RETURN
      END
C  this program calculates the heat capacity of a particle from Merrick's 
C  correlations.
      Subroutine heatcp(TP,cp)
C  calculates daf coal heat capacity
C   tp          Particle Temperature (K)
C   cp          Particle heat capacity (cal/g/K)
C   u(i)        Atomic Weights of elements
C   y(i)        Mass fractions of elements (C,H,N,O,S)
C   RGAS        Gas constant (1.987 cal/gmole/K)
C   RGASJ       Gas constant (8314.3 J/kg/K)
      dimension y(5),u(5)
      common/spheat/yelem(5)
      data u/12.,1.,14.,16.,32./
      data RGAS/1.987/,RGASJ/8314.3/
C  calculate mean atomic weight, a
      a = 0.0
      do 10 i=1,5
        a = a+(yelem(i)/u(i))
10    continue
      a = 1./a
      f1 = 380./tp
      f2 = 1800./tp
      cp = rgas/a *(g1(f1)+2.*g1(f2))
      return
      end
C
      Function g1(z)
      g1=exp(z)/((exp(z)-1)/z)**2
      return
      end
      Subroutine heatap(TP,cpa)
C  calculates ash heat capacity
C   tp          Particle Temperature (K)
C   cpa         Ash heat capacity (cal/g/K)
C   RGAS        Gas constant (1.987 cal/g/K)
      cpa = (754+.586*(tp-273))/(4.184e3)
      return
      end

C
C------------------------------------------------------------------
C
      subroutine perks(y,yp,T)
      IMPLICIT real*4 (A-H,O-Z)
      save
C
C  this subroutine is the meat of the devolatilization model
C      
C  y1 = l       labile bridges
C  y2 = del     ends
C  y3 = C       char links
C  y4 = fnca    mass nitrogen per mass (aromatic) site
C  yp(i) = derivative of y(i) in time
C
C     nmax = number of terms in expansion for mol. wt. distribution
      real*4 l0,l,ma,kb,kc,kg,kn,kp,ft(35),mt(35),y(4),yp(4)
      common/init/l0,c0,g0,ma,rba,finf,sig,siginv,nmax,pstar
      common/rate/ab,eb0,ebsig,ac,ec0,ecsig,ag,eg0,egsig,rg
      common/nitmod/fnca0,fnca,an,en0,ensig
      data fx/0./,ft/35*0./,mt/35*0./
C
      l = y(1)
      del = y(2)
      C = y(3)
      fnca = y(4)
      p = l+C
      g1 = 2.*(1-p)-del
      g2 = 2.*(C-c0)
      g = g1+g2
C  calculate current activation energy using error function solution
      if(c0.lt.1.)fx = g/(1.-c0)/2.
      call inverf(fx,x)
      eg = eg0 + x*egsig
      if (fnca0.gt.0)fx = 1.-fnca/fnca0
      call inverf(fx,x)
      en = en0 + x*ensig
      if(l0.gt.0.)fx = 1.-l/l0
      call inverf(fx,x)
      eb = eb0+x*ebsig
      ec = ec0
C  calculate rate constants
      rt = rg*t
      kb = ab*exp(-eb/rt)
      rho = ac*exp(-ec/rt)
      kg = ag*exp(-eg/rt)
      kn = an*exp(-en/rt)
C  calculate rate of destruction of labile bridges
      yp(1) = -(kb)*l
C  calculate rate of formation of ends (danglers)
      yp(2) = 2.*rho*kb*l/(rho+1.) - kg*del
C  calculate rate of formation of char
      yp(3) = kb*l/(rho+1.)
C  calculate rate of high T (slow) nitrogen loss (as HCN)
      yp(4) = -kn*fnca
      return
      end
C------------------------------------------------------------------
C
      subroutine perkp(y,mgas,ftar,ftart,fgas,fchar,ft,mt,intar)
      IMPLICIT real*4 (A-H,O-Z)
      save
C
C   calculates fractions of tar, gas, and char from p, sig, l, and C
C
      real*4 l0,l,ma,mb,kp,y(3),yp(3),ft(35),mt(35),
     x    mgas,mtot
      logical intar
      common/init/l0,c0,g0,ma,rba,finf,sig,siginv,nmax,pstar
      common/rate/ab,eb0,ebsig,ac,ec0,ecsig,ag,eg0,egsig,rg
      common/tarmw/press
      data pi/3.141592/,iprint/1/
C
      l = y(1)
      del = y(2)
      C = y(3)
      p = l+C
      if(intar)then
        if(p.gt.0.9999)then
          delfac = 1.
        else
          delfac = del/(1.-p)
        endif
        a = 1.+rba*(l/p + (sig-1.)/4. *delfac)
        b = (delfac/2. - l/p)
C  find pstar
        pstar0 = pstar
        pinv = siginv+1.e-4
        if(p.ge.0.9999)then
           pstar = 0.
        elseif(p.ge.pinv)then
           do 10 i=1,25
              f = pstar*(1-pstar)**(sig-1) - p*(1-p)**(sig-1)
              fp = (1-pstar)**(sig-1)-pstar*(sig-1)*(1-pstar)**(sig-2)
              ppstar = pstar - f/fp
              err= abs(1.-ppstar/pstar)
              if(err.le.1.e-4)go to 11
              pstar = ppstar
10         continue
           print*,' warning--pstar did not converge'
           print*,' p=',p,' sig=',sig,' pstar=',pstar,
     x       '  pstar0=',pstar0
11         continue
        else
           pstar = p
        endif
C
C  check to see if pstar is in the right range
      if(pstar.lt.0. .or. (p.ne.pstar .and. pstar.ge.siginv))then
           print*,' error--pstar out of acceptable ranges!'
           print*,'     pstar=',pstar,';  HOW DID YOU DO THAT?'
           print*,'p=',p,'  sig=',sig,'  pstar0=',pstar0
           stop
        endif
C
        sfac = (sig+1.)/(sig-1.)
        fp = (pstar/p)**(sfac)
        kp = (pstar/p)**sfac*(1.-(sig+1)/2. *pstar)
C  calculate wt fraction tar, gas, and char
        ftart = 2.*(a*fp+rba*b*kp)/(2.+rba*(1.-c0)*(sig+1.))
      endif
      tarfac = 1.-ftar
      g1 = (2.*(1.0 - p) - del)
      g2 = 2.*(C-c0)
      g = g1+g2
C
      mgas = rba*ma*g*(sig+1)/4.*tarfac
      mtot = ma + rba*ma*(sig+1)/2. *(1.-c0)
      fgas = mgas/mtot
      fchar = 1.-ftar-fgas
C  calculate tar molecular weight distribution
      if(.not.intar)return
C      open(unit=25,name='perktar.out',type='unknown')
      ftsum = 0.0
      do 20 n=1,nmax
         tn = n*(sig-1.)+2
         xm = n*sig+1.
         yk = n-1.
         xm1 = xm+1.
C  gamln is the solution to the gamma function in the Sandia Math Library
         fg1 = gamln(xm1)
         if(fg1.le.1.e-10) then
            fgam = 0.
         else
           yk1 = yk+1.
           fg2 = gamln(yk1)
           xmyk = xm-yk+1.
           fg3 = gamln(xmyk)
           fgam = exp(fg1-fg2-fg3)
         endif
         bnn = (sig+1.)/(n*sig+1.)*fgam
         qn = bnn*(p**(n-1))*((1-p)**tn)/n
C  ft(n) = weight fraction of each tar bin
         ft(n) = 2.*(n*a*qn+rba*b*qn)/(2.+rba*(1.-c0)*(sig+1.))
         ftsum = ftsum + ft(n)
C  check to not divide by zero
         if(p.le.1.e-9)then
            fac = 0
         else
            fac = l/p
         endif
         tst = 1.-p
         if(tst.le.1.e-9)then
            fac1 = 0.
         else
            fac1 = del/(1.-p)
         endif
C  mt(n) = molecular weight of each tar bin
         mt(n) = n*ma+(n-1)*rba*ma*fac+tn*rba*ma/4.*fac1
         if(iprint.eq.1)then
C           print*,n,ft(n),ftsum,mt(n)
C           write(24,*)ft(n),ftsum,mt(n)
           iprint = iprint + 1
         endif
20    continue
      return
      end

      SUBROUTINE INVERF(Y,X)
      IMPLICIT real*4 (A-H,O-Z)
      save
C  this program calculates the inverse of the area under the normal curve.
C  if y=area(x), then given y, this program will calculate x.
C  A table lookup is performed.
      dimension xx(18),yy(18)
      data xx/3.4,3.2,3.,2.8,2.6,2.4,2.2,2.,1.8,1.6,1.4,
     x        1.2,1.,.8,.6,.4,.2,0./
      data yy/.9997,.9993,.9987,.9974,.9953,.9918,.9861,.9772,.9641,
     x        .9452,.9192,.8849,.8413,.7881,.7257,.6554,.5793,.5/
C
      fac = 1.
C  check to see if y is within range
C      if(y.lt.0.0228)then
C         x = -2.0
C         return
      if(y.lt.0.0003)then
         x = -3.4
         return
      elseif(y.lt.0.5)then
        yp = 1.-y
        fac = -1.
      elseif(y.gt.0.9997)then
        x = 3.5
        return
      else
        yp = y
      endif
C  search for range
      do 10 i=18,1,-1
        if(yp.le.yy(i-1))then
          x = xx(i) + (yp-yy(i))*(xx(i-1)-xx(i))/(yy(i-1)-yy(i))
          x = fac*x
          return
        endif
10    continue
      return
      end
C
      function gamln(x)
C   this is a program to calculate the ln of the gamma function,
C   taken from Abramowitz, p. 257, 6.1.41
      data pi/3.141592/
         gamln = (x-.5)*alog(x)-x+.5*alog(2.*pi)+1./(12.*x)
     x        -1./(360.*x**3)+1./(1260.*x**5)-1./(1680.*x**7)
      return
      end
C
      subroutine flash(fgas,gasmw,ft,mt,fracr,ftar,fmet,
     x                  temp,press,nmax,imolw,ipred)
      IMPLICIT real*4 (A-H,O-Z)
      save
C  flash distillation of metaplast to form liquid and tar vapor
      real*4 k(36),l(36),Ltot,mt(35),metold(35),metold2(35)
      dimension f(36),ft(35),xmw(36),z(36),pv(36),x(36),y(36),
     x          v(36),ftold(35),tarold(35),ftold2(35),tarold2(35)
      common/timer/tms,time
      logical imolw,ipred
C  ipred = .true. only on predictor step, when old values are updated      
      data small/1.e-3/,a,b,g/87058,299,0.5903/,x3,x2/.2,.3/,
     x     zero/0.0/,ip/30/,x/36*0./,y/36*0./
      data k/36*0/,l/36*0/,v/36*0/,metold/35*0./,ftold/35*0./,
     x     xmw/36*0./,tarold/35*0./,f/36*0./,z/36*0./,pv/36*0./
     x     ,fgasold/0./,metold2/35*0./,ftold2/35*0./,tarold2/35*0./
     x     ,fgasold2/0./,ftarold/0./,ftarold2/0./
C
C  xxxold2 is the value at the previous time step, while
C  xxxold  is the value from the last time the subroutine was called
C  metold(i) = mass fraction of coal contained in metaplast of mer size i
C  fracr = fraction to account for reduction of metaplast by crosslinking
C          in latest time step
C  
C  renormalize in tar fragment data to molar basis
C  f(i) = moles of mer i
C  ft(i) = mass of mer i
C  
      Ftot = 0.0
      do 10 i=1,nmax
         i1=i+1
         xmw(i1) = mt(i)
         if(ipred)then
            ftold2(i) = ftold(i)
            metold2(i) = metold(i)
            tarold2(i) = tarold(i)
            fgasold2 = fgasold
            ftarold2 = ftarold
         endif
         dif = ft(i)-ftold2(i)
         dif = max(dif,zero)
         f(i1) = (dif+metold2(i)*fracr)/mt(i)
         ftold(i) = ft(i)
         Ftot = Ftot + f(i1)
10    continue
c20    continue
      ntot = nmax + 1
      f(1) = (fgas-fgasold2)/gasmw
      f(1) = max(f(1),0.)
      fgasold = max(fgas,fgasold2)
      xmw(1) = gasmw
      Ftot = Ftot + f(1)
C  get mole fraction of components in the feed stream
C  and compute equilibrium contants k(i) from vapor pressure and
C  Raoults law expression
      sum = 0.0
      do 30 ii=1,ntot
        sum = sum + f(ii)
        pv(ii) = a*exp(-b*xmw(ii)**g/temp)
        k(ii) = pv(ii)/press
        if(k(ii).lt.0.001)k(ii) = 0.0
30    continue
      if(sum.le.1.e-8)return
      do 40 ii=1,ntot
        z(ii) = f(ii)/sum
40    continue
C  use the Rachford-Rice formulation for flash distillation
C  x = V/F, first guess
      x1 = x3
C  calculate sum (Eq. 11-24, Separation Processes, by King, 1971)
      f1 = 0.0
      do 50 ii=1,ntot
         f1 = f1 + z(ii)*(k(ii)-1)/((k(ii)-1)*(x1)+1)
50    continue
C  Secant Method for convergence
      if(x2.eq.x1)x2 = x1+.005
      do 70 iter=1,100
C  calculate sum (Eq. 11-24, Separation Processes, by King, 1971)
        f2 = 0.0
        do 60 ii=1,ntot
           f2 = f2 + z(ii)*(k(ii)-1)/((k(ii)-1)*(x2)+1)
60      continue
        if(abs(f2).le.small.or.abs(f2-f1).le.small**2)go to 100
        x3 = x2-f2*(x2-x1)/(f2-f1)
        if(x3.gt.1.0)x3 = 1.0-small**2
        if(x3.lt.0.)x3 = 0.0+small**2
        if(x3.eq.x2)then
          print*,'problem---f(V/F) not converged, but x3=x2',x3,x2
          if(x2.ge.small)then
            x3=x2-small
          else
            x3 = x2+small
          endif
        endif
        if(x2.le.1.e-5.and.x1.le.1.e-5)then
           x2 = 1.e-7
           go to 100
        endif
        if(x2.ge..9999 .and. x1.ge..9999)then
           x2 = .9999
           go to 100
        endif
        f1 = f2
        x1 = x2
C  under-relax solution (especially near the V/F=1 limit
        x2 = 0.2*x2+0.8*x3
70      continue
        print*,'Convergence not achieved in Flash distillation'
        print*,' last two guesses for V/F were',x1,x3
        do 81 ii=1,ntot
          write(30,*)z(ii),k(ii)
81      continue
        stop
100     continue
C  now calculate molecular weight distributions on a
C  light-gas free basis, wt fractions
        Vtot = Ftot*x2
        Ltot = Ftot-Vtot
        VoL = Vtot/Ltot
        sumx = 0.0
        sumy = 0.0
        xmwtot = 0.0
        ttot = 0.0
        do 200 ii=2,ntot
          i = ii-1
          l(ii) = f(ii)/(1.+k(ii)*VoL)
          v(ii) = f(ii)-l(ii)
          x(ii) = l(ii)*xmw(ii)
          y(ii) = v(ii)*xmw(ii)
          metold(i) = max(x(ii),zero)
          tarold(i) = tarold2(i)+y(ii)
          xmwtot = xmwtot+tarold(i)*xmw(ii)
          ttot = ttot+tarold(i)
C          if(imolw)write(ip,213)xmw(ii),tarold(i),metold(i)
C213       format(' ',f7.0,3(1pe10.3,3x))
          sumx = sumx + x(ii)
          sumy = sumy + y(ii)
C          write(ip,210)ii,xmw(ii),k(ii),f(ii),v(ii),l(ii)
C          write(ip,210)i,xmw(ii),k(ii),ft(i),tarold(i),metold(i)
200     continue
        if(ttot.gt.0.)xmwtot = xmwtot/ttot
        do 250 ii=2,ntot
           if(sumx.ne.0.0)x(ii) = x(ii)/sumx
           if(sumy.ne.0.0)y(ii) = y(ii)/sumy
C          write(ip,210)ii,xmw(ii),k(ii),f(ii),yv,xl
250     continue
        ftar = ftarold2+sumy
        ftarold = ftar
        fmet = sumx
        if(imolw)then
          ip = ip+1
          if(mod(ip,10).eq.0)then
C              print*,'Weight Av. Molecular Wt.=',xmwtot
C     x         ,' Gas MW =',gasmw
           endif
        endif
C210     format(' ',i2,2x,f5.0,4(1pe10.3,2x))
C211     format(' ',f5.0,4(1pe10.3,2x))
        return
        end
C  ------------------------------------------------------------------
C  --------------LIGHT GAS DISTRIBUTION MODEL------------------------
C  ------------------------------------------------------------------
        subroutine lightgas(yygas,inside,lib)
C  This program calculates the distribution of light gas species
C  based on a look up table of reference data

        IMPLICIT real*4 (A-H,O-Z)
        real*4 ind
        integer*2 lib
        logical inside
        common/lgmod/yf,xoc,yhc
        save
        dimension xx(12),yy(12),a(23),b(23),s1(12),s2(12),s3(12)
        dimension p1(12),p2(12),p3(12),n(12) 
        dimension xz(12,12),fh2o(12,12),fco2(12,12)
        dimension out(5,14),fco(12,12),fch4(12,12),ygas(4,3),yygas(5)
        
C *******************************************************************
C ****************LIGHT GAS DISTRIBUTION REFERENCE LIBRARY***********
C *******************************************************************

C  This library can be moved outside of submodel as long as
C  it is linked to the light gas sub-model

C  xz = The index used in correlation.  In main program it 
C  corresponds to variable "yf"

C  f*** = fraction of light gas that is species ***

C  the data is organized in the following order with 12 ordered
C  pairs for each species (ordered with xz)

C  each table is organized in rows in the following order 
C       1 Lower Kittaning (Chen)
C       2 Pocahontas #3 (ANL)
C       3 Upper Freeport (ANL)
C       4 Pittsburgh (Chen)
C       5 Lewis Stockton (ANL)  
C       6 Utah Blind Canyon (ANL)
C       7 Illinois #6 (ANL)
C       8 Illinois #6 (Chen)
C       9 Wyodak (ANL)
C       10 Beulah Zap (ANL)
C       11 Dietz (Chen)
C       12 PSOC 1448 (Serio)

C  reference data for xz = yf, the fractional light gas released

        data xz/0.,.04,.11,.14,.21,.27,.34,.675,.9,1.,.0,.0,
     x  .0,.161,.442,.663,.777,.874,.921,.967,1.,.0,.0,.0,
     x  .0,.022,.20,.430,.526,.64,.787,.875,.927,.955,1.,.0,
     x  .0,.04,.12,.15,.23,.29,.36,.68,.9,1.,.0,.0,
     x  .0,.018,.058,.21,.417,.572,.696,.778,.821,.883,.932,1.,
     x  .0,.052,.144,.291,.498,.639,.746,.859,.925,.949,.966,1.,
     x  .0,.063,.178,.33,.506,.612,.706,.813,.895,.94,1.,.0,
     x  .0,.04,.12,.15,.23,.29,.36,.68,.9,1.,.0,.0,
     x  .0,.061,.146,.374,.535,.622,.714,.8,.883,.931,.964,1.,
     x  .0,.034,.087,.179,.316,.472,.585,.694,.777,.872,.935,1.,
     x  .0,.04,.12,.16,.25,.31,.37,.68,.9,1.,.0,.0,
     x  .0,.02,.055,.17,.313,.434,.546,.716,.874,.935,.973,1./
     
C  fraction of light gas that is water
     
          data fh2o/.772,.772,.738,.455,.371,.304,.290,.273,.218,
     x              .218,.0,.0,
     x  .699,.632,.299,.269,.247,.249,.236,.225,.226,.0,.0,.0,
     x  .0,.0,.35,.297,.301,.299,.284,.291,.306,.297,.283,.0,
     x  .636,.636,.646,.550,.436,.320,.186,.199,.195,.195,.0,.0,
     x  1.,.983,.754,.488,.413,.385,.373,.382,.377,0.362,.367,.348,
     x  .665,.636,.604,.508,.435,.409,.383,.362,.351,.343,.342,.339,
     x  .763,.737,.698,.572,.527,.470,.438,.411,.411,.396,.378,.0,
     x  .748,.748,.637,.704,.490,.446,.348,.268,.266,.266,.0,.0,
     x  .0,.0,.385,.461,.396,.369,.344,.323,.292,.277,.266,.257,
     x  .0,.0,.197,.267,.26,.333,.361,.369,.346,.306,.285,.267,
     x  .521,.521,.55,.523,.511,.46,.414,.388,.313,0.313,.0,.0,
     x  .0,.0,.291,.335,.264,.271,.261,.211,.171,.160,.153,.149/
     
C  fraction of light gas that is carbon dioxide

          data fco2/.0,.0,.0,.174,.174,.167,.129,.102,.071,.071,.0,.0,
     x  .259,.234,.113,.086,.097,.109,.116,.118,.122,.0,.0,.0,
     x  .333,.327,.070,.052,.057,.06,.059,.062,.066,.08,0.115,.0,
     x  .194,.194,.152,.117,.116,.122,.081,.092,.065,.065,.0,.0,
     x  .0,.0,.0,.122,.103,.086,.083,.082,.085,.086,.093,.128,
     x  .332,.318,.165,.141,.120,.108,.105,.119,.120,.122,.125,.130,
     x  .229,.221,.125,.09,.07,.073,.083,.133,.132,.13,.147,.0,
     x  .111,.111,.142,.175,.149,.155,.136,.122,.133,.133,.0,.0,
     x  .98,.984,.55,.345,.317,.285,.286,.277,.273,.264,.254,.255,
     x  .993,.989,.786,.572,.519,.416,.375,.345,.335,.32,.303,.299,
     x  .363,.363,.353,.325,.321,.35,.318,.251,.249,.249,.0,.0,
     x  1.,.983,.448,.179,.104,.09,.104,.151,.166,.160,.158,.154/
     
C  fraction of light gas that is methane

          data fch4/.203,.203,.078,.160,.180,.219,.258,.294,.320,
     x               .320,.0,.0,
     x  .041,.037,.388,.389,.359,.332,.323,.307,.299,.0,.0,.0,
     x  .667,.655,.42,.454,.444,.419,.382,.353,.331,.321,.306,.0,
     x  .055,.055,.073,.088,.116,.124,.170,.15,.189,.189,.0,.0,
     x  .0,.0,.188,.195,.234,.243,.224,.21,.2,.186,.177,.167,
     x  .0,.0,.11,.155,.176,.172,.185,.173,.163,.159,.156,.151,
     x  .0,.0,.075,.136,.159,.178,.174,.157,.143,.141,.132,.0,
     x  .02,.02,.026,.042,.045,.049,.064,.1,.128,.128,.0,.0,
     x  .0,.0,.0,.029,.048,.067,.069,.072,.069,.066,.063,.061,
     x  .0,.0,.0,.0,.035,.05,.061,0.058,.057,.053,.049,.046,
     x  .01,.01,.011,.016,.011,.021,.023,.035,.06,.06,.0,.0,
     x  .0,.0,.216,.262,.362,.327,.307,.25,.203,.189,.182,.177/
     
C  fraction of light gas that is carbon monoxide
     
          data fco/.0,.0,.157,.121,.141,.112,.139,.085,.145,
     x              .145,.0,.0,
     x  .0,.0,.0,.057,.097,.109,.124,.15,.153,.0,.0,.0,
     x  .0,.0,.0,.0,.0,.024,.078,.097,.099,.104,.099,.0,
     x  .083,.083,.038,.066,.032,.168,.286,.324,.313,.313,.0,.0,
     x  .0,.0,.0,.0,.055,.091,.124,.131,.142,.171,.168,.162,
     x  .0,.0,.0,.028,.093,.129,.142,.162,.181,.191,.193,.195,
     x  .0,.0,.0,.075,.099,.122,.139,.133,.148,.167,.177,.0,
     x  .101,.101,.173,.054,.219,.247,.335,.349,.28,.280,.0,.0,
     x  .0,.0,.055,.115,.151,.168,.172,.2,.236,.264,.287,.298,
     x  .0,.0,.0,.133,.142,.150,.15,.173,.206,.265,.307,.331,
     x  .096,.096,.066,.113,.123,.13,.2,.281,.334,.334,.0,.0,
     x  .0,.0,.0,.084,.078,.115,.130,.191,.262,.294,.311,.322/
     
C       ************************************************************
C       **********End of Reference library**************************
C       ************************************************************
        
C  *********************************************************************        
C  *******DETERMINE THE APPROPRIATE TRIANGLE FOR INTERPOLATION**********
C  *********************************************************************

C define the equations of line, distance and area

      yyy(aa,xd,bb)=aa*xd+bb
      xxx(aa,yd,bb)=(yd-bb)/aa
      d(x2,x1,y2,y1)=((x2-x1)**2+(y2-y1)**2)**.5
      at(aa,bb,cc)=0.5*bb*cc*(1-((bb**2+cc**2-aa**2)/(2*bb*cc))**2)**.5

C  inititalize variables
        x = xoc
        y = yhc
        ind = yf
        
C  look up table of the reference points of the triangular mesh

        data xx/.0177734,.0203654,.0659401,.0773465,.0893623,.1077369,
     x          .1305803,.1357569,.1803479,.2093441,.2603201,.0687/
        data yy/.6717240,.5810955,.6550527,.8088697,.7575807,.8506428,
     x          .7671163,.8523190,.8499221,.7890888,.8572938,.863/
     
C  look up table for a and b of equations of all triangle sides

        data a/-34.965,1.6228,-.34612,2.3021,1.7337,1.1993,4.3774,
     x       -4.2685,.23134,5.0647,1.3746,-3.6565,.059818,16.459,
     x       1.6638,-.05375,.27897,-2.0979,.092179,1.3380,3.7559,
     x       -6.2604,-.31655/
        data b/1.2932,.54805,.67788,.63081,.54074,.65041,.36641,
     x       1.1390,.73691,.30499,.70255,1.2446,.84420,-1.3821,
     x       .54985,.85962,.73069,1.2283,.83330,.50899,.60497,
     x       1.2931,.88475/
     
C  look up table for the three sides that correspond to each triangle

        data s1/1,3,4,7,8,10,12,14,15,18,21,22/
        data s2/2,7,6,5,10,9,14,15,17,20,4,11/
        data s3/3,6,8,9,11,12,13,16,18,19,22,23/
        
C  do loop to find the appropriate triangle for interpolation
C       if(iii.eq.1)then
        m=0.0
        inside = .true.   
        DO 10 i=1,12
           z1=xxx(a(s1(i)),y,b(s1(i)))
           z2=xxx(a(s2(i)),y,b(s2(i)))
           z3=yyy(a(s3(i)),x,b(s3(i)))
          if(x.GE.z1.AND.x.LE.z2.AND.y.LE.z3)THEN
             m=i  
          go to 30
          endif
          if(i.EQ.12.AND.m.EQ.0.0)Then
             print*,'one or both ratios are out of bounds'
             inside = .false.
          endif    
10      continue
30      continue
C        print*,' triangle =',m
C       endif
C  *****************************************************************    
C  *****************TRIANGULAR INTERPOLATION************************
C  *****************************************************************
        if(inside)then
C  This interpolation scheme is taken from Zhao et al., 25th Symp.
C  on Comb. (Int.), 1994/pp. 553-560. 

C  look up table for points 1,2, and 3 for each triangle

        data p1/2,3,1,3,5,5,7,7,7,10,1,4/
        data p2/1,1,4,5,4,6,6,8,9,9,12,12/
        data p3/3,5,5,7,6,7,8,9,10,11,4,6/
        
C  calculate the length of each side
C       if(iii.eq.1)then
        ds1=d(xx(p1(i)),xx(p2(i)),yy(p1(i)),yy(p2(i)))
        ds2=d(xx(p3(i)),xx(p1(i)),yy(p3(i)),yy(p1(i)))
        ds3=d(xx(p3(i)),xx(p2(i)),yy(p3(i)),yy(p2(i)))
        ds4=d(x,xx(p2(i)),y,yy(p2(i)))
        ds5=d(xx(p1(i)),x,yy(p1(i)),y)
        ds6=d(xx(p3(i)),x,yy(p3(i)),y)
        
C       print*,' ds1,ds2,ds3 = ',ds1,ds2,ds3

C  calculate the area of each triangle used in interpolation scheme

        A1=at(ds1,ds2,ds3)
        A2=at(ds1,ds5,ds4)
        A3=at(ds5,ds2,ds6)
        
C       print*,' A1,A2,A3 = ',A1,A2,A3
        
C  calculate S and R, the weighted fraction of two of the points
C  the weighted fraction of other point will be 1-R-S

        S=A2/A1
        R=A3/A1
C       endif
C       ************************************************************
C       *************CALCULATE LIGHT GAS DISTRIBUTION***************
C       ************************************************************
 
C  n is the number of order pairs of data for coals 1-12

        data n/10,9,11,10,12,12,11,10,12,12,10,12/

C  do loop to calculate the light gas composition of each reference
C  coal (point) of triangle.  this is accomplished using linear 
C  interpolation between points in reference library
        
C  j specifies the point (1,2, or3) of the triangle selected above

        DO 35 j=1,3
          if(j.EQ.1)THEN
            lib=p1(i)
          elseif(j.EQ.2)THEN
            lib=p2(i)
          else
            lib=p3(i)
          endif
        
C  do loop to find the two xz points that surround ind

        DO 40 ii=1,12
          if(ind.GE.xz(ii,lib).AND.ind.LE.xz(ii+1,lib))THEN
          go to 45
          endif
40      continue

C       for ygas(k,j)
C               k=1;fh2o
C               k=2;fco2
C               k=3;fch4
C               k=4;fco
C               k=5;other light gases (HC'2, H2, parafins, olifins)

C  linear interpolation to find reference values as a function of ind

45      ygas(1,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fh2o(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fh2o(ii+1,lib)
        ygas(2,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fco2(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fco2(ii+1,lib)
        ygas(3,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fch4(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fch4(ii+1,lib)
        ygas(4,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fco(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fco(ii+1,lib)
     
C         print*,'lib=',lib
C         print*,'ii=',ii
C         print*,'xz=',xz(ii,lib)

35      continue
        endif
C       ************************************************************
C       *******CALCULATE GAS COMPOSITION FROM LIBRARY COALS*********
C       ************************************************************
        if(inside)then
C  yygas(k) = the fraction of light gas that is k
C  1=h20, 2=co2, 3=ch4, 4=co

        DO 50 k=1,4           
          yygas(k)=(1-R-S)*ygas(k,1)+R*ygas(k,2)+S*ygas(k,3)
C         print*,ygas(k,1),ygas(k,2),ygas(k,3)
C         print*,'xz(5,1)',xz(5,1)
50      continue
        endif
        
C       ************************************************************
C       *****Estimate composition for coals outside mesh************
C       ************************************************************
        
        if(inside.eqv..false.)then
        data out/.0,.085,.835,1.5,12,.085,.12,.835,1.5,6,.12,.155,
     x          .835,1.5,8,.155,.222,.835,1.5,9,.222,.285,.75,1.5,11,
     x          .0,.089,.75,.835,4,.0,.05,.63,.75,1,.05,.175,.63,
     x          .69,3,.066,.175,.69,.835,7,.175,.222,.0,.835,10,
     x          .0,.175,.55,.63,2,.222,1.,.0,.75,10,.285,1,0.75,1.5,13,
     x           .0,.175,.0,.55,14/
     
        DO 60 kk=1,14
            if(x.GT.out(1,kk).AND.x.LE.out(2,kk))then
               if(y.GE.out(3,kk).AND.y.LT.out(4,kk))then
            lib=out(5,kk)
            if(lib.eq.13)then
              yygas(1)=0.24
              yygas(2)=0.37
              yygas(3)=0.06
              yygas(4)=0.28
              go to 70
            endif
            if(lib.eq.14)then
              yygas(1)=0.18
              yygas(2)=0.08
              yygas(3)=0.37
              yygas(4)=0.18
              go to 70
            endif
C  do loop to find the two xz points that surround ind

        DO 63 ii=1,12
          if(ind.GE.xz(ii,lib).AND.ind.LE.xz(ii+1,lib))THEN
          go to 64
          endif
63      continue

C  linear interpolation to find reference values as a function of ind

64      yygas(1)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fh2o(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fh2o(ii+1,lib)
        yygas(2)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fco2(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fco2(ii+1,lib)
        yygas(3)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fch4(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fch4(ii+1,lib)
        yygas(4)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))
     x                  *fco(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)
     x            -xz(ii,lib)))*fco(ii+1,lib)
               print*,out(1,kk),out(2,kk),out(3,kk),out(4,kk),out(5,kk)
               go to 65
               endif
            endif
60      continue
65      print*,'Light gas distribution is based on ref. # ',lib
C       print*,out(1,6),out(2,6),out(3,6),out(4,6),out(5,6)
        print*,x,y
        endif
        
70      yygas(5)=1-yygas(1)-yygas(2)-yygas(3)-yygas(4)
C       print*,'R and S = ',R,S
C       print*,'fh2o, fco2, fch4, fco, & fother =',yygas(1),
C    x             yygas(2),yygas(3),yygas(4),yygas(5)
C       print*,inside,m,x,y,ind
        nn=ntmax-2
        return    
90      end
