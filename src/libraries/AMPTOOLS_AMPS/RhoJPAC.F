c------------------------------------       
      double complex function ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,M)
c      implicit double precision (a-h,o-z)
      double complex p1,p2,p3,p4
      double complex param 
      double precision b
      double precision csig
      common /pars/ param(4)
      common /slopes/ b(4) 
      double complex a_ro_pomex,a_ro_f2ex,a_ro_piex,a_ro_sigmaex
        
!        csig=CONV*0.5D0*PIK/((2.d0*PI)**3*32.d0*PM2*OMEGA**2)
!     factor 0.5 changed to 0.25 as both photon helicities are included in the c-s
      pim=0.13957d0
      pm =  0.93827d0
      pm2 = pm**2
      CONV=389.38d0
      pi = 4.d0*datan(1.d0)
       pik=dsqrt(xm**2/4.d0-pim**2)
! lab photon energy
       omega=(s-pm**2)/(2.d0*pm)
       csig=CONV*0.25D0*pik/((2.d0*PI)**3*32.d0*PM2*OMEGA**2)

c       test=a_ro_pomex(s,t,xm,lbd_g,lbd_t,lbd_r,M)        
c       print*,"s, t, xm", s, t, xm, p1, test

       p1 =  param(1)*dexp(b(1)*t)
       p2 =  param(2)*dexp(b(2)*t)
       p3 =  param(3)*dexp(b(3)*t)
       p4 =  param(4)*dexp(b(4)*t) 
c       cisg = 1.
c       p2 = 0.
c       p3 = 0.
c       p4 = 0.
c       p1 = 1.
      ampltot = dsqrt(csig)*(
     $      p1*a_ro_pomex(s,t,xm,lbd_g,lbd_t,lbd_r,M)
     $     +p2*a_ro_f2ex(s,t,xm,lbd_g,lbd_t,lbd_r,M)
     $     +p3*a_ro_piex(s,t,xm,lbd_g,lbd_t,lbd_r,M)
     $     +p4*a_ro_sigmaex(s,t,xm,lbd_g,lbd_t,lbd_r,M) )

c        print*,"ampltot",ampltot
c        print*,"ro_piex", a_ro_piex(s,t,xm,lbd_g,lbd_t,lbd_r,M)
c        print*,"ro_sigmaex", a_ro_sigmaex(s,t,xm,lbd_g,lbd_t,lbd_r,M)

        return
        end 

c--------------------------


      
      function xnorm(s,t,xm,xnorm_sum)
c     implicit double precision (a-h,o-z)
      double complex  ala,ampltot

c     Lukasz 
      
c       dsdtdm=dsdtdm+csig*cdabs(
c     & a*dexp(b*et)*a_ro_pomex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)+
c     & a_ro_piex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)+
c     & a_ro_f2ex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)+
c     & c*a_ro_sigmaex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)+
c     & dexp(x*et)*(y+z*empipi))**2

c      print*,"xnorm s, t, xm", s, t, xm

      xnorm = 0.d0
      xnorm_sum = 0.d0
      do lbd_g=-1,1,2
      do lbd_t=-1,1,2
      do lbd_r=-1,1,2
      do M=-1,1
      ala =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,M)
c      print*,"g, t, r, M, ala", lbd_g, lbd_t, lbd_r, M, ala
      xnorm_sum = xnorm_sum + dble(ala*dconjg(ala))
      enddo
      enddo
      enddo
      enddo

      return
      end 

      
      double complex function rho_h(k,l,lp,s,t,xm)
      implicit double precision (a-h,u-w)
      double complex xi,rho
      double complex  ala,ola,ampltot,a_ro_pomex
      xi = (0.d0,1.d0)

c      print*,"s, t, xm", s, t, xm

c      xtemp = 0.d0
c      xn = xnorm(s,t,xm,xtemp)
c      xn = xtemp/2.

c      xn = 0. 
c      do lbd_g=-1,1,2
c      do lbd_t=-1,1,2
c      do lbd_r=-1,1,2
c      do M=-1,1
c      xn = xn + dble(ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,M)
c     $  *dconjg(ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,M)))
c      enddo
c      enddo
c      enddo
c      enddo
c      xn = xn/2.       
      
c---------------- rho^0 
      if (k.eq.0) then
      rho = (0.d0,0.d0)
      do lbd_g=-1,1,2
      do lbd_t=-1,1,2
      do lbd_r=-1,1,2
      do M=-1,1
         do Mp=-1,1
         if ((Mp.eq.lp).and.(M.eq.l)) then
      ala =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,M)
c            ala = a_ro_pomex(s,t,xm,lbd_g,lbd_t,lbd_r,M)
      ola =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,Mp)
c            ola = a_ro_pomex(s,t,xm,lbd_g,lbd_t,lbd_r,Mp)
      rho = rho + ala*dconjg(ola)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c----------------------------------------      
c---------------- rho^1
      if (k.eq.1) then
      rho = (0.d0,0.d0)
      do lbd_g=-1,1,2
      do lbd_t=-1,1,2
      do lbd_r=-1,1,2
      do M=-1,1
      do Mp=-1,1
      if ((Mp.eq.lp).and.(M.eq.l)) then
      if (lbd_g.eq.1) then
      lbd_gp = -1
      else
      lbd_gp = 1
      endif
      ala =  ampltot(s,t,xm,lbd_gp,lbd_t,lbd_r,M)
      ola =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,Mp)
      rho = rho + ala*dconjg(ola)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c----------------------------------------            
c---------------- rho^2
      if (k.eq.2) then
      rho = (0.d0,0.d0)
      do lbd_g=-1,1,2
      do lbd_t=-1,1,2
      do lbd_r=-1,1,2
      do M=-1,1
      do Mp=-1,1
      if ((Mp.eq.lp).and.(M.eq.l)) then
      if (lbd_g.eq.1) then
      lbd_gp = -1
      else
      lbd_gp = 1
      endif
      ala =  ampltot(s,t,xm,lbd_gp,lbd_t,lbd_r,M)
      ola =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,Mp)
      rho = rho + xi*lbd_g*ala*dconjg(ola)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c---------------- rho^3
      if (k.eq.3) then
      rho = (0.d0,0.d0)
      do lbd_g=-1,1,2
      do lbd_t=-1,1,2
      do lbd_r=-1,1,2
      do M=-1,1
      do Mp=-1,1
      if ((Mp.eq.lp).and.(M.eq.l)) then
         ala =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,M)
         ola =  ampltot(s,t,xm,lbd_g,lbd_t,lbd_r,Mp)
      rho = rho + lbd_g*ala*dconjg(ola)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c-----------
      
      rho_h = rho 
c      /(2.*xn)

      return
      end 
      
      
      

      double complex function rho_gj(k,l,lp,s,t,xm)
      implicit double precision (a-h,o-z)
      dimension d(3,3) 
      double complex rho_h
      double complex tmp(3,3)
      
      xmn=0.93827d0
      tgj = 2.*s*(t - xm**2+(s + xm**2 - xmn**2)*(s-xmn**2)/(2.*s))
     $   /(dsqrt((s-(xmn+xm)**2)*(s-(xmn-xm)**2))*(s-xmn**2))
      beta = dsqrt((s-(xm+xmn)**2)*(s-(xm-xmn)**2))
     $     /(s+xm**2-xmn**2)
      tg = (beta-tgj)/(beta*tgj-1)

      d(1+2,1+2) = (1. + tg)/2.
      d(1+2,0+2) = -dsqrt(1.-tg**2)/dsqrt(2.d0)
      d(1+2,-1+2) = (1. - tg)/2.

      d(0+2,1+2) = dsqrt(1.-tg**2)/dsqrt(2.d0)
      d(0+2,0+2) = tg
      d(0+2,-1+2) = -dsqrt(1.-tg**2)/dsqrt(2.d0)

      d(-1+2,1+2) = (1. - tg)/2.
      d(-1+2,0+2) = dsqrt(1.-tg**2)/dsqrt(2.d0)
      d(-1+2,-1+2) =  (1. + tg)/2.
      
      do l1=-1,1
      do l2=-1,1
      tmp(l1+2,l2+2) = (0.d0,0.d0)
      do j1=-1,1
      do j2=-1,1   
      tmp(l1+2,l2+2) = tmp(l1+2,l2+2) 
     $  + d(l1+2,j1+2)*d(l2+2,j2+2)*rho_h(k,j1,j2,s,t,xm)
      enddo
      enddo
      enddo
      enddo

      rho_gj = tmp(l+2,lp+2)
      return
      end

c----------------------------------------------------------------------------      
!This is the ro photoproduction amplitude with t-channel sigma exchange and
! subsequent ro->pi^+pi^-  decay
!The amplitude is based on the model by Oh and Lee Phys.Rev. C69, 025201 (2004)
C--------------------------------------------------------
! This is the version of the amplitude which computes 4 nucleon spin combinations at once
C--------------------------------------------------------
       SUBROUTINE AP_F2(BPP,BMP,BPM,BMM)
       implicit double precision(a-h,o-z)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & APP(3),APM(3),AMP(3),AMM(3),FFR(3),THETA(3),UPS(3),XIX(3),
     & XIY(3),XIZ(3),X0(3),Xx(3),Xy(3),Xz(3),C0(3),CX(3),CY(3),CZ(3)
       COMPLEX*16 BPP(3),BMP(3),BPM(3),BMM(3),CI,BWIG
       DOUBLE PRECISION LBDFNN,LBDFROGAM,LBD1,LBD2
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG
! st. sprzeżenia i parametry formfaktorów z pracy Oh i Lee
       FRO=5.33d0
       GFRORO=5.76d0
       EMF2=1.2751d0
   
       GFROGAM=E*GFRORO/FRO
       GFNN=2.12d0
       LBDFNN=1.4d0
       LBDFROGAM=1.4d0
       GROPIPI=6.04d0

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=EFFM
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E1=(S-PM2+T)/2.d0/EMPI
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1KW=E1**2-PM**2
       P1Z=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1X=DSQRT(P1KW-P1Z**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq
       QX=qxver*Q
       QZ=qzver*Q
       xLBD_G=DBLE(LBDG)
       epsx=-xLBD_G*costhq/DSQRT(2.d0)
       epsz=-xLBD_G*sinthq/DSQRT(2.d0)

       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       RHO(1)=sinthq/DSQRT(2.D0)
       RHO(2)=costhq
       RHO(3)=-RHO(1)

       SIG1=(1.d0+EMPI2/EMF2**2)**2+(EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q)
     & /(3.d0*EMF2**2)
       SIG2=EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q
       SIG3=1.d0+EMPI2/EMF2**2+2.d0*EMPI*Q/EMF2**2
       SIG4=2.d0*EMPI*Q

       LBD1=Q+EMPI+2.d0*(EMPI/EMF2)**2*Q-(EMPI/EMF2)**4*(EMPI-Q)+
     & (EMPI-Q)/(3.d0*EMF2**2)*(-EMPI2-6.d0*EMPI*Q+(EMPI2/EMF2)**2)
       LBD2=Q+EMPI-(EMPI-Q)/EMF2**2*(EMPI2+2.d0*EMPI*Q)

       DLT1=(P2MOD**2-P1MOD**2)*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & ((P2MOD**2+P1Z*P2MOD)/E2M-(P1Z*P2MOD+P1KW)/E1M)
       DLT2=(epsx*P1X+epsz*(P1Z-P2MOD))*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
       DLT3=(E1+E2)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)-
     & (P2MOD**2-P1KW)*P1X*P2MOD/EE
       DLT4=(E1+E2)*(-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M)-
     & (epsx*P1X+epsz*(P1Z-P2MOD))*P1X*P2MOD/EE
       DLT5=-(E1+E2)/DSQRT(2.d0)*(P2MOD/E2M+P1Z/E1M)
       DLT6=(E1+E2)/DSQRT(2.d0)*P1X/E1M

       PSI1=(P2MOD**2-P1KW)*
     & (P2MOD**2/E2M-P1KW/E1M-P2MOD*P1Z*(E2-E1)/EE)
       PSI2=-(P2MOD**2/E2M+P1KW/E1M-P2MOD*P1Z*EPLUS)/3.d0
       PSI3=-((P2MOD**2-P1KW)*(-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
     & +(epsx*P1X+epsz*(P1Z-P2MOD))*((P2MOD**2+P1Z*P2MOD)/E2M-
     & (P1Z*P2MOD+P1KW)/E1M))
       PSI5Y=(P2MOD**2-P1KW)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)
       PSI6Y=(P1X*(P2MOD/E2M+P1Z/E1M)-(P1Z-P2MOD)*P1X/E1M)/3.d0
       PSI7X=(P2MOD**2-P1KW)*(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       PSI7Y=(P2MOD**2-P1KW)*(epsx*(P2MOD/E2M+P1Z/E1M)-epsz*P1X/E1M)+
     & (epsx*P1X+epsz*(P1Z-P2MOD))*(QX*(P2MOD/E2M+P1Z/E1M)-QZ*P1X/E1M)
       PSI7Z=-(P2MOD**2-P1KW)*P1X/E1M/DSQRT(2.d0)

       THETA(1)=Q*RHO(1)
       THETA(2)=2.d0*P2MOD+Q*RHO(2)
       THETA(3)=Q*RHO(3)

       UPS(1)=Q/E1M*RHO(1)
       UPS(2)=P2MOD*EPLUS+Q/E1M*RHO(2)
       UPS(3)=Q/E1M*RHO(3)

       XIX(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIX(2)=0.d0
       XIX(3)=XIX(1)
       XIY(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIY(2)=P1X/E1M
       XIY(3)=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIZ(1)=-P1X/E1M/DSQRT(2.D0)
       XIZ(2)=0.d0
       XIZ(3)=-P1X/E1M/DSQRT(2.D0)

       Ex=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       Ey=-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M
       Ez=P1X/E1M/DSQRT(2.d0)

       do I=1,3
       X0(I)=(PSI1*SIG1+PSI2*SIG2)*BE(I)+PSI3*SIG3*Q*RHO(I)+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*UPS(I)+
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)*THETA(I))*SIG4

       Xx(I)=-PSI7X*Q*RHO(I)*SIG3+((epsx*P1X+epsz*(P1Z-P2MOD))*XIX(I)
     & -Ex*THETA(I))*SIG4

       Xy(I)=(PSI5Y*SIG1+PSI6Y*SIG2)*BE(I)+PSI7Y*Q*SIG3*RHO(I)-
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIY(I)-THETA(I)*Ey)*SIG4

       Xz(I)=-PSI7Z*RHO(I)*SIG3+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIZ(I)-THETA(I)*Ez)*SIG4
       enddo

       U=(EMPI+Q-(EMPI/EMF2)**2*(EMPI-Q))**2-
     & (1.d0-(EMPI-Q)**2/EMF2**2)/3.d0*
     & (EMPI2+6.d0*EMPI*Q-(EMPI2/EMF2)**2)

       do I=1,3
       C0(I)=(E1+E2)*U*BE(I)*(1.d0-P1Z*P2MOD)/EE-LBD1*DLT1*BE(I)+
     & Q*RHO(I)*LBD2*DLT2+X0(I)
       CX(I)=-Q*RHO(I)*LBD2*DLT5+Xx(I)
       CY(I)=((E1+E2)*U*P1X*P2MOD/EE-LBD1*DLT3)*BE(I)+
     & Q*RHO(I)*LBD2*DLT4+Xy(I)
       CZ(I)=-Q*RHO(I)*LBD2*DLT6+Xz(I)
       enddo

       R=-4.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)*GFNN/PM*GROPIPI
     & *GFROGAM/EMF2/DSQRT(4.d0*PI)
       CT=DSQRT((1.D0+P1Z/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)
! 2012.02.21 - uwaga przyjmuję konwencję indeksowania spinów taką jak w amplitudach drella, z wymianą pi i sigma
       PRF2=1.d0/(T-EMF2**2)

       FF_NNF=(LBDFNN**2-EMF2**2)/(LBDFNN**2-T)
       FF_GAMROF=(LBDFROGAM**2-EMF2**2)/(LBDFROGAM**2-T)

       DO 101 M=1,3
       APP(M)=C0(M)*ST+CX(M)*CT-CY(M)*CT-CZ(M)*ST
       AMP(M)=-C0(M)*CT-CX(M)*ST-CY(M)*ST-CZ(M)*CT
       APM(M)=C0(M)*CT-CX(M)*ST+CY(M)*ST-CZ(M)*CT
       AMM(M)=C0(M)*ST-CX(M)*CT-CY(M)*CT+CZ(M)*ST
c       MM=M
       
       BPP(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*APP(M)
       BMP(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*AMP(M)
       BPM(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*APM(M)
       BMM(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*AMM(M)
101    CONTINUE
!        print*,t,EFFM,FORMA,BSLOPE,BPP(1),BPP(2)
       RETURN
       END

c-------------------------------------------------------------------------------------
! Use routines AP_F2_PP,AP_F2_MP,AP_F2_PM,AP_F2_MM below to compute amplitudes for individual spin combinations
c-------------------------------------------------------------------------------------

       SUBROUTINE AP_F2_PP(BPP, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & APP(3),THETA(3),UPS(3),XIX(3),
     & XIY(3),XIZ(3),X0(3),Xx(3),Xy(3),Xz(3),C0(3),CX(3),CY(3),CZ(3)
       COMPLEX*16 BPP(3),CI,BWIG
       DOUBLE PRECISION LBDFNN,LBDFROGAM,LBD1,LBD2
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG
! couplings and ff. params from Oh and Lee paper
       FRO=5.33d0
       GFRORO=5.76d0
       EMF2=1.2751d0
  
c       print*,"s, t, xm", S, T, zmpipi
c       print*,"stale",CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,REZ_GAM2
c     & ,E,F,PIM
 
       GFROGAM=0.30282d0*GFRORO/FRO
c       print*,"test GFROGAM",E,GFRORO,FRO
       GFNN=2.12d0
       LBDFNN=1.4d0
       LBDFROGAM=1.4d0
       GROPIPI=6.04d0

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E1=(S-PM2+T)/2.d0/EMPI
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1KW=E1**2-PM**2
       P1Z=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1X=DSQRT(P1KW-P1Z**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq
       QX=qxver*Q
       QZ=qzver*Q
       xLBD_G=DBLE(LBDG)
       epsx=-xLBD_G*costhq/DSQRT(2.d0)
       epsz=-xLBD_G*sinthq/DSQRT(2.d0)

       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       RHO(1)=sinthq/DSQRT(2.D0)
       RHO(2)=costhq
       RHO(3)=-RHO(1)

       SIG1=(1.d0+EMPI2/EMF2**2)**2+(EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q)
     & /(3.d0*EMF2**2)
       SIG2=EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q
       SIG3=1.d0+EMPI2/EMF2**2+2.d0*EMPI*Q/EMF2**2
       SIG4=2.d0*EMPI*Q

       LBD1=Q+EMPI+2.d0*(EMPI/EMF2)**2*Q-(EMPI/EMF2)**4*(EMPI-Q)+
     & (EMPI-Q)/(3.d0*EMF2**2)*(-EMPI2-6.d0*EMPI*Q+(EMPI2/EMF2)**2)
       LBD2=Q+EMPI-(EMPI-Q)/EMF2**2*(EMPI2+2.d0*EMPI*Q)

       DLT1=(P2MOD**2-P1MOD**2)*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & ((P2MOD**2+P1Z*P2MOD)/E2M-(P1Z*P2MOD+P1KW)/E1M)
       DLT2=(epsx*P1X+epsz*(P1Z-P2MOD))*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
       DLT3=(E1+E2)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)-
     & (P2MOD**2-P1KW)*P1X*P2MOD/EE
       DLT4=(E1+E2)*(-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M)-
     & (epsx*P1X+epsz*(P1Z-P2MOD))*P1X*P2MOD/EE
       DLT5=-(E1+E2)/DSQRT(2.d0)*(P2MOD/E2M+P1Z/E1M)
       DLT6=(E1+E2)/DSQRT(2.d0)*P1X/E1M

       PSI1=(P2MOD**2-P1KW)*
     & (P2MOD**2/E2M-P1KW/E1M-P2MOD*P1Z*(E2-E1)/EE)
       PSI2=-(P2MOD**2/E2M+P1KW/E1M-P2MOD*P1Z*EPLUS)/3.d0
       PSI3=-((P2MOD**2-P1KW)*(-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
     & +(epsx*P1X+epsz*(P1Z-P2MOD))*((P2MOD**2+P1Z*P2MOD)/E2M-
     & (P1Z*P2MOD+P1KW)/E1M))
       PSI5Y=(P2MOD**2-P1KW)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)
       PSI6Y=(P1X*(P2MOD/E2M+P1Z/E1M)-(P1Z-P2MOD)*P1X/E1M)/3.d0
       PSI7X=(P2MOD**2-P1KW)*(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       PSI7Y=(P2MOD**2-P1KW)*(epsx*(P2MOD/E2M+P1Z/E1M)-epsz*P1X/E1M)+
     & (epsx*P1X+epsz*(P1Z-P2MOD))*(QX*(P2MOD/E2M+P1Z/E1M)-QZ*P1X/E1M)
       PSI7Z=-(P2MOD**2-P1KW)*P1X/E1M/DSQRT(2.d0)

       THETA(1)=Q*RHO(1)
       THETA(2)=2.d0*P2MOD+Q*RHO(2)
       THETA(3)=Q*RHO(3)

       UPS(1)=Q/E1M*RHO(1)
       UPS(2)=P2MOD*EPLUS+Q/E1M*RHO(2)
       UPS(3)=Q/E1M*RHO(3)

       XIX(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIX(2)=0.d0
       XIX(3)=XIX(1)
       XIY(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIY(2)=P1X/E1M
       XIY(3)=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIZ(1)=-P1X/E1M/DSQRT(2.D0)
       XIZ(2)=0.d0
       XIZ(3)=-P1X/E1M/DSQRT(2.D0)

       Ex=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       Ey=-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M
       Ez=P1X/E1M/DSQRT(2.d0)

       do I=1,3
       X0(I)=(PSI1*SIG1+PSI2*SIG2)*BE(I)+PSI3*SIG3*Q*RHO(I)+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*UPS(I)+
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)*THETA(I))*SIG4

       Xx(I)=-PSI7X*Q*RHO(I)*SIG3+((epsx*P1X+epsz*(P1Z-P2MOD))*XIX(I)
     & -Ex*THETA(I))*SIG4

       Xy(I)=(PSI5Y*SIG1+PSI6Y*SIG2)*BE(I)+PSI7Y*Q*SIG3*RHO(I)-
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIY(I)-THETA(I)*Ey)*SIG4

       Xz(I)=-PSI7Z*RHO(I)*SIG3+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIZ(I)-THETA(I)*Ez)*SIG4
       enddo

       U=(EMPI+Q-(EMPI/EMF2)**2*(EMPI-Q))**2-
     & (1.d0-(EMPI-Q)**2/EMF2**2)/3.d0*
     & (EMPI2+6.d0*EMPI*Q-(EMPI2/EMF2)**2)

       do I=1,3
       C0(I)=(E1+E2)*U*BE(I)*(1.d0-P1Z*P2MOD)/EE-LBD1*DLT1*BE(I)+
     & Q*RHO(I)*LBD2*DLT2+X0(I)
       CX(I)=-Q*RHO(I)*LBD2*DLT5+Xx(I)
       CY(I)=((E1+E2)*U*P1X*P2MOD/EE-LBD1*DLT3)*BE(I)+
     & Q*RHO(I)*LBD2*DLT4+Xy(I)
       CZ(I)=-Q*RHO(I)*LBD2*DLT6+Xz(I)
       enddo

c       print*,"test R",EE,PIK,GFNN,PM,GROPIPI,GFROGAM,EMF2
       R=-4.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)*GFNN/PM*GROPIPI
     & *GFROGAM/EMF2/DSQRT(4.d0*PI)
       CT=DSQRT((1.D0+P1Z/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)
! 2012.02.21 - uwaga przyjmuję konwencję indeksowania spinów taką jak w amplitudach drella, z wymianą pi i sigma
       PRF2=1.d0/(T-EMF2**2)

       FF_NNF=(LBDFNN**2-EMF2**2)/(LBDFNN**2-T)
       FF_GAMROF=(LBDFROGAM**2-EMF2**2)/(LBDFROGAM**2-T)

c       print*,"BPP", BPP
       DO 102 M=1,3
       APP(M)=C0(M)*ST+CX(M)*CT-CY(M)*CT-CZ(M)*ST
       
       BPP(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*APP(M)
c       print*,"test BPP",R,BWIG,PRF2,FF_NNF,FF_GAMROF,APP(M)
c       print*,"M, BPP", M, BPP(M)
102    CONTINUE

       RETURN
       END

       SUBROUTINE AP_F2_MP(BMP, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & AMP(3),THETA(3),UPS(3),XIX(3),
     & XIY(3),XIZ(3),X0(3),Xx(3),Xy(3),Xz(3),C0(3),CX(3),CY(3),CZ(3)
       COMPLEX*16 BMP(3),CI,BWIG
       DOUBLE PRECISION LBDFNN,LBDFROGAM,LBD1,LBD2
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG
! couplings and ff. params from Oh and Lee paper
       FRO=5.33d0
       GFRORO=5.76d0
       EMF2=1.2751d0
   
       GFROGAM=0.30282d0*GFRORO/FRO
       GFNN=2.12d0
       LBDFNN=1.4d0
       LBDFROGAM=1.4d0
       GROPIPI=6.04d0

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E1=(S-PM2+T)/2.d0/EMPI
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1KW=E1**2-PM**2
       P1Z=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1X=DSQRT(P1KW-P1Z**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq
       QX=qxver*Q
       QZ=qzver*Q
       xLBD_G=DBLE(LBDG)
       epsx=-xLBD_G*costhq/DSQRT(2.d0)
       epsz=-xLBD_G*sinthq/DSQRT(2.d0)

       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       RHO(1)=sinthq/DSQRT(2.D0)
       RHO(2)=costhq
       RHO(3)=-RHO(1)

       SIG1=(1.d0+EMPI2/EMF2**2)**2+(EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q)
     & /(3.d0*EMF2**2)
       SIG2=EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q
       SIG3=1.d0+EMPI2/EMF2**2+2.d0*EMPI*Q/EMF2**2
       SIG4=2.d0*EMPI*Q

       LBD1=Q+EMPI+2.d0*(EMPI/EMF2)**2*Q-(EMPI/EMF2)**4*(EMPI-Q)+
     & (EMPI-Q)/(3.d0*EMF2**2)*(-EMPI2-6.d0*EMPI*Q+(EMPI2/EMF2)**2)
       LBD2=Q+EMPI-(EMPI-Q)/EMF2**2*(EMPI2+2.d0*EMPI*Q)

       DLT1=(P2MOD**2-P1MOD**2)*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & ((P2MOD**2+P1Z*P2MOD)/E2M-(P1Z*P2MOD+P1KW)/E1M)
       DLT2=(epsx*P1X+epsz*(P1Z-P2MOD))*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
       DLT3=(E1+E2)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)-
     & (P2MOD**2-P1KW)*P1X*P2MOD/EE
       DLT4=(E1+E2)*(-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M)-
     & (epsx*P1X+epsz*(P1Z-P2MOD))*P1X*P2MOD/EE
       DLT5=-(E1+E2)/DSQRT(2.d0)*(P2MOD/E2M+P1Z/E1M)
       DLT6=(E1+E2)/DSQRT(2.d0)*P1X/E1M

       PSI1=(P2MOD**2-P1KW)*
     & (P2MOD**2/E2M-P1KW/E1M-P2MOD*P1Z*(E2-E1)/EE)
       PSI2=-(P2MOD**2/E2M+P1KW/E1M-P2MOD*P1Z*EPLUS)/3.d0
       PSI3=-((P2MOD**2-P1KW)*(-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
     & +(epsx*P1X+epsz*(P1Z-P2MOD))*((P2MOD**2+P1Z*P2MOD)/E2M-
     & (P1Z*P2MOD+P1KW)/E1M))
       PSI5Y=(P2MOD**2-P1KW)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)
       PSI6Y=(P1X*(P2MOD/E2M+P1Z/E1M)-(P1Z-P2MOD)*P1X/E1M)/3.d0
       PSI7X=(P2MOD**2-P1KW)*(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       PSI7Y=(P2MOD**2-P1KW)*(epsx*(P2MOD/E2M+P1Z/E1M)-epsz*P1X/E1M)+
     & (epsx*P1X+epsz*(P1Z-P2MOD))*(QX*(P2MOD/E2M+P1Z/E1M)-QZ*P1X/E1M)
       PSI7Z=-(P2MOD**2-P1KW)*P1X/E1M/DSQRT(2.d0)

       THETA(1)=Q*RHO(1)
       THETA(2)=2.d0*P2MOD+Q*RHO(2)
       THETA(3)=Q*RHO(3)

       UPS(1)=Q/E1M*RHO(1)
       UPS(2)=P2MOD*EPLUS+Q/E1M*RHO(2)
       UPS(3)=Q/E1M*RHO(3)

       XIX(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIX(2)=0.d0
       XIX(3)=XIX(1)
       XIY(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIY(2)=P1X/E1M
       XIY(3)=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIZ(1)=-P1X/E1M/DSQRT(2.D0)
       XIZ(2)=0.d0
       XIZ(3)=-P1X/E1M/DSQRT(2.D0)

       Ex=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       Ey=-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M
       Ez=P1X/E1M/DSQRT(2.d0)

       do I=1,3
       X0(I)=(PSI1*SIG1+PSI2*SIG2)*BE(I)+PSI3*SIG3*Q*RHO(I)+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*UPS(I)+
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)*THETA(I))*SIG4

       Xx(I)=-PSI7X*Q*RHO(I)*SIG3+((epsx*P1X+epsz*(P1Z-P2MOD))*XIX(I)
     & -Ex*THETA(I))*SIG4

       Xy(I)=(PSI5Y*SIG1+PSI6Y*SIG2)*BE(I)+PSI7Y*Q*SIG3*RHO(I)-
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIY(I)-THETA(I)*Ey)*SIG4

       Xz(I)=-PSI7Z*RHO(I)*SIG3+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIZ(I)-THETA(I)*Ez)*SIG4
       enddo

       U=(EMPI+Q-(EMPI/EMF2)**2*(EMPI-Q))**2-
     & (1.d0-(EMPI-Q)**2/EMF2**2)/3.d0*
     & (EMPI2+6.d0*EMPI*Q-(EMPI2/EMF2)**2)

       do I=1,3
       C0(I)=(E1+E2)*U*BE(I)*(1.d0-P1Z*P2MOD)/EE-LBD1*DLT1*BE(I)+
     & Q*RHO(I)*LBD2*DLT2+X0(I)
       CX(I)=-Q*RHO(I)*LBD2*DLT5+Xx(I)
       CY(I)=((E1+E2)*U*P1X*P2MOD/EE-LBD1*DLT3)*BE(I)+
     & Q*RHO(I)*LBD2*DLT4+Xy(I)
       CZ(I)=-Q*RHO(I)*LBD2*DLT6+Xz(I)
       enddo

       R=-4.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)*GFNN/PM*GROPIPI
     & *GFROGAM/EMF2/DSQRT(4.d0*PI)
       CT=DSQRT((1.D0+P1Z/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)
! 2012.02.21 - uwaga przyjmuję konwencję indeksowania spinów taką jak w amplitudach drella, z wymianą pi i sigma
       PRF2=1.d0/(T-EMF2**2)

       FF_NNF=(LBDFNN**2-EMF2**2)/(LBDFNN**2-T)
       FF_GAMROF=(LBDFROGAM**2-EMF2**2)/(LBDFROGAM**2-T)

       DO 103 M=1,3
       AMP(M)=-C0(M)*CT-CX(M)*ST-CY(M)*ST-CZ(M)*CT
       
       BMP(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*AMP(M)
103    CONTINUE

       RETURN
       END

       SUBROUTINE AP_F2_PM(BPM, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & APM(3),THETA(3),UPS(3),XIX(3),
     & XIY(3),XIZ(3),X0(3),Xx(3),Xy(3),Xz(3),C0(3),CX(3),CY(3),CZ(3)
       COMPLEX*16 BPM(3),CI,BWIG
       DOUBLE PRECISION LBDFNN,LBDFROGAM,LBD1,LBD2
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG
! couplings and ff. params from Oh and Lee paper
       FRO=5.33d0
       GFRORO=5.76d0
       EMF2=1.2751d0
   
       GFROGAM=0.30282d0*GFRORO/FRO
       GFNN=2.12d0
       LBDFNN=1.4d0
       LBDFROGAM=1.4d0
       GROPIPI=6.04d0

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E1=(S-PM2+T)/2.d0/EMPI
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1KW=E1**2-PM**2
       P1Z=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1X=DSQRT(P1KW-P1Z**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq
       QX=qxver*Q
       QZ=qzver*Q
       xLBD_G=DBLE(LBDG)
       epsx=-xLBD_G*costhq/DSQRT(2.d0)
       epsz=-xLBD_G*sinthq/DSQRT(2.d0)

       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       RHO(1)=sinthq/DSQRT(2.D0)
       RHO(2)=costhq
       RHO(3)=-RHO(1)

       SIG1=(1.d0+EMPI2/EMF2**2)**2+(EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q)
     & /(3.d0*EMF2**2)
       SIG2=EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q
       SIG3=1.d0+EMPI2/EMF2**2+2.d0*EMPI*Q/EMF2**2
       SIG4=2.d0*EMPI*Q

       LBD1=Q+EMPI+2.d0*(EMPI/EMF2)**2*Q-(EMPI/EMF2)**4*(EMPI-Q)+
     & (EMPI-Q)/(3.d0*EMF2**2)*(-EMPI2-6.d0*EMPI*Q+(EMPI2/EMF2)**2)
       LBD2=Q+EMPI-(EMPI-Q)/EMF2**2*(EMPI2+2.d0*EMPI*Q)

       DLT1=(P2MOD**2-P1MOD**2)*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & ((P2MOD**2+P1Z*P2MOD)/E2M-(P1Z*P2MOD+P1KW)/E1M)
       DLT2=(epsx*P1X+epsz*(P1Z-P2MOD))*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
       DLT3=(E1+E2)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)-
     & (P2MOD**2-P1KW)*P1X*P2MOD/EE
       DLT4=(E1+E2)*(-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M)-
     & (epsx*P1X+epsz*(P1Z-P2MOD))*P1X*P2MOD/EE
       DLT5=-(E1+E2)/DSQRT(2.d0)*(P2MOD/E2M+P1Z/E1M)
       DLT6=(E1+E2)/DSQRT(2.d0)*P1X/E1M

       PSI1=(P2MOD**2-P1KW)*
     & (P2MOD**2/E2M-P1KW/E1M-P2MOD*P1Z*(E2-E1)/EE)
       PSI2=-(P2MOD**2/E2M+P1KW/E1M-P2MOD*P1Z*EPLUS)/3.d0
       PSI3=-((P2MOD**2-P1KW)*(-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
     & +(epsx*P1X+epsz*(P1Z-P2MOD))*((P2MOD**2+P1Z*P2MOD)/E2M-
     & (P1Z*P2MOD+P1KW)/E1M))
       PSI5Y=(P2MOD**2-P1KW)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)
       PSI6Y=(P1X*(P2MOD/E2M+P1Z/E1M)-(P1Z-P2MOD)*P1X/E1M)/3.d0
       PSI7X=(P2MOD**2-P1KW)*(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       PSI7Y=(P2MOD**2-P1KW)*(epsx*(P2MOD/E2M+P1Z/E1M)-epsz*P1X/E1M)+
     & (epsx*P1X+epsz*(P1Z-P2MOD))*(QX*(P2MOD/E2M+P1Z/E1M)-QZ*P1X/E1M)
       PSI7Z=-(P2MOD**2-P1KW)*P1X/E1M/DSQRT(2.d0)

       THETA(1)=Q*RHO(1)
       THETA(2)=2.d0*P2MOD+Q*RHO(2)
       THETA(3)=Q*RHO(3)

       UPS(1)=Q/E1M*RHO(1)
       UPS(2)=P2MOD*EPLUS+Q/E1M*RHO(2)
       UPS(3)=Q/E1M*RHO(3)

       XIX(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIX(2)=0.d0
       XIX(3)=XIX(1)
       XIY(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIY(2)=P1X/E1M
       XIY(3)=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIZ(1)=-P1X/E1M/DSQRT(2.D0)
       XIZ(2)=0.d0
       XIZ(3)=-P1X/E1M/DSQRT(2.D0)

       Ex=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       Ey=-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M
       Ez=P1X/E1M/DSQRT(2.d0)

       do I=1,3
       X0(I)=(PSI1*SIG1+PSI2*SIG2)*BE(I)+PSI3*SIG3*Q*RHO(I)+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*UPS(I)+
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)*THETA(I))*SIG4

       Xx(I)=-PSI7X*Q*RHO(I)*SIG3+((epsx*P1X+epsz*(P1Z-P2MOD))*XIX(I)
     & -Ex*THETA(I))*SIG4

       Xy(I)=(PSI5Y*SIG1+PSI6Y*SIG2)*BE(I)+PSI7Y*Q*SIG3*RHO(I)-
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIY(I)-THETA(I)*Ey)*SIG4

       Xz(I)=-PSI7Z*RHO(I)*SIG3+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIZ(I)-THETA(I)*Ez)*SIG4
       enddo

       U=(EMPI+Q-(EMPI/EMF2)**2*(EMPI-Q))**2-
     & (1.d0-(EMPI-Q)**2/EMF2**2)/3.d0*
     & (EMPI2+6.d0*EMPI*Q-(EMPI2/EMF2)**2)

       do I=1,3
       C0(I)=(E1+E2)*U*BE(I)*(1.d0-P1Z*P2MOD)/EE-LBD1*DLT1*BE(I)+
     & Q*RHO(I)*LBD2*DLT2+X0(I)
       CX(I)=-Q*RHO(I)*LBD2*DLT5+Xx(I)
       CY(I)=((E1+E2)*U*P1X*P2MOD/EE-LBD1*DLT3)*BE(I)+
     & Q*RHO(I)*LBD2*DLT4+Xy(I)
       CZ(I)=-Q*RHO(I)*LBD2*DLT6+Xz(I)
       enddo

       R=-4.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)*GFNN/PM*GROPIPI
     & *GFROGAM/EMF2/DSQRT(4.d0*PI)
       CT=DSQRT((1.D0+P1Z/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)
! 2012.02.21 - uwaga przyjmuję konwencję indeksowania spinów taką jak w amplitudach drella, z wymianą pi i sigma
       PRF2=1.d0/(T-EMF2**2)

       FF_NNF=(LBDFNN**2-EMF2**2)/(LBDFNN**2-T)
       FF_GAMROF=(LBDFROGAM**2-EMF2**2)/(LBDFROGAM**2-T)

       DO 104 M=1,3
       APM(M)=C0(M)*CT-CX(M)*ST+CY(M)*ST-CZ(M)*CT
       
       BPM(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*APM(M)
104    CONTINUE

       RETURN
       END

       SUBROUTINE AP_F2_MM(BMM, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & AMM(3),THETA(3),UPS(3),XIX(3),
     & XIY(3),XIZ(3),X0(3),Xx(3),Xy(3),Xz(3),C0(3),CX(3),CY(3),CZ(3)
       COMPLEX*16 BMM(3),CI,BWIG
       DOUBLE PRECISION LBDFNN,LBDFROGAM,LBD1,LBD2
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG
! couplings and ff. params from Oh and Lee paper
       FRO=5.33d0
       GFRORO=5.76d0
       EMF2=1.2751d0
   
       GFROGAM=0.30282d0*GFRORO/FRO
       GFNN=2.12d0
       LBDFNN=1.4d0
       LBDFROGAM=1.4d0
       GROPIPI=6.04d0

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E1=(S-PM2+T)/2.d0/EMPI
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1KW=E1**2-PM**2
       P1Z=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1X=DSQRT(P1KW-P1Z**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq
       QX=qxver*Q
       QZ=qzver*Q
       xLBD_G=DBLE(LBDG)
       epsx=-xLBD_G*costhq/DSQRT(2.d0)
       epsz=-xLBD_G*sinthq/DSQRT(2.d0)

       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       RHO(1)=sinthq/DSQRT(2.D0)
       RHO(2)=costhq
       RHO(3)=-RHO(1)

       SIG1=(1.d0+EMPI2/EMF2**2)**2+(EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q)
     & /(3.d0*EMF2**2)
       SIG2=EMPI2-(EMPI2/EMF2)**2+6.d0*EMPI*Q
       SIG3=1.d0+EMPI2/EMF2**2+2.d0*EMPI*Q/EMF2**2
       SIG4=2.d0*EMPI*Q

       LBD1=Q+EMPI+2.d0*(EMPI/EMF2)**2*Q-(EMPI/EMF2)**4*(EMPI-Q)+
     & (EMPI-Q)/(3.d0*EMF2**2)*(-EMPI2-6.d0*EMPI*Q+(EMPI2/EMF2)**2)
       LBD2=Q+EMPI-(EMPI-Q)/EMF2**2*(EMPI2+2.d0*EMPI*Q)

       DLT1=(P2MOD**2-P1MOD**2)*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & ((P2MOD**2+P1Z*P2MOD)/E2M-(P1Z*P2MOD+P1KW)/E1M)
       DLT2=(epsx*P1X+epsz*(P1Z-P2MOD))*(1.d0-P1Z*P2MOD/EE)+(E1+E2)*
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
       DLT3=(E1+E2)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)-
     & (P2MOD**2-P1KW)*P1X*P2MOD/EE
       DLT4=(E1+E2)*(-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M)-
     & (epsx*P1X+epsz*(P1Z-P2MOD))*P1X*P2MOD/EE
       DLT5=-(E1+E2)/DSQRT(2.d0)*(P2MOD/E2M+P1Z/E1M)
       DLT6=(E1+E2)/DSQRT(2.d0)*P1X/E1M

       PSI1=(P2MOD**2-P1KW)*
     & (P2MOD**2/E2M-P1KW/E1M-P2MOD*P1Z*(E2-E1)/EE)
       PSI2=-(P2MOD**2/E2M+P1KW/E1M-P2MOD*P1Z*EPLUS)/3.d0
       PSI3=-((P2MOD**2-P1KW)*(-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)
     & +(epsx*P1X+epsz*(P1Z-P2MOD))*((P2MOD**2+P1Z*P2MOD)/E2M-
     & (P1Z*P2MOD+P1KW)/E1M))
       PSI5Y=(P2MOD**2-P1KW)*(-QX*(P1Z/E1M+P2MOD/E2M)+QZ*P1X/E1M)
       PSI6Y=(P1X*(P2MOD/E2M+P1Z/E1M)-(P1Z-P2MOD)*P1X/E1M)/3.d0
       PSI7X=(P2MOD**2-P1KW)*(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       PSI7Y=(P2MOD**2-P1KW)*(epsx*(P2MOD/E2M+P1Z/E1M)-epsz*P1X/E1M)+
     & (epsx*P1X+epsz*(P1Z-P2MOD))*(QX*(P2MOD/E2M+P1Z/E1M)-QZ*P1X/E1M)
       PSI7Z=-(P2MOD**2-P1KW)*P1X/E1M/DSQRT(2.d0)

       THETA(1)=Q*RHO(1)
       THETA(2)=2.d0*P2MOD+Q*RHO(2)
       THETA(3)=Q*RHO(3)

       UPS(1)=Q/E1M*RHO(1)
       UPS(2)=P2MOD*EPLUS+Q/E1M*RHO(2)
       UPS(3)=Q/E1M*RHO(3)

       XIX(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIX(2)=0.d0
       XIX(3)=XIX(1)
       XIY(1)=(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIY(2)=P1X/E1M
       XIY(3)=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       XIZ(1)=-P1X/E1M/DSQRT(2.D0)
       XIZ(2)=0.d0
       XIZ(3)=-P1X/E1M/DSQRT(2.D0)

       Ex=-(P2MOD/E2M+P1Z/E1M)/DSQRT(2.D0)
       Ey=-epsx*(P2MOD/E2M+P1Z/E1M)+epsz*P1X/E1M
       Ez=P1X/E1M/DSQRT(2.d0)

       do I=1,3
       X0(I)=(PSI1*SIG1+PSI2*SIG2)*BE(I)+PSI3*SIG3*Q*RHO(I)+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*UPS(I)+
     & (-epsz*P2MOD/E2M+(epsx*P1X+epsz*P1Z)/E1M)*THETA(I))*SIG4

       Xx(I)=-PSI7X*Q*RHO(I)*SIG3+((epsx*P1X+epsz*(P1Z-P2MOD))*XIX(I)
     & -Ex*THETA(I))*SIG4

       Xy(I)=(PSI5Y*SIG1+PSI6Y*SIG2)*BE(I)+PSI7Y*Q*SIG3*RHO(I)-
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIY(I)-THETA(I)*Ey)*SIG4

       Xz(I)=-PSI7Z*RHO(I)*SIG3+
     & ((epsx*P1X+epsz*(P1Z-P2MOD))*XIZ(I)-THETA(I)*Ez)*SIG4
       enddo

       U=(EMPI+Q-(EMPI/EMF2)**2*(EMPI-Q))**2-
     & (1.d0-(EMPI-Q)**2/EMF2**2)/3.d0*
     & (EMPI2+6.d0*EMPI*Q-(EMPI2/EMF2)**2)

       do I=1,3
       C0(I)=(E1+E2)*U*BE(I)*(1.d0-P1Z*P2MOD)/EE-LBD1*DLT1*BE(I)+
     & Q*RHO(I)*LBD2*DLT2+X0(I)
       CX(I)=-Q*RHO(I)*LBD2*DLT5+Xx(I)
       CY(I)=((E1+E2)*U*P1X*P2MOD/EE-LBD1*DLT3)*BE(I)+
     & Q*RHO(I)*LBD2*DLT4+Xy(I)
       CZ(I)=-Q*RHO(I)*LBD2*DLT6+Xz(I)
       enddo

       R=-4.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)*GFNN/PM*GROPIPI
     & *GFROGAM/EMF2/DSQRT(4.d0*PI)
       CT=DSQRT((1.D0+P1Z/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)
! 2012.02.21 - uwaga przyjmuję konwencję indeksowania spinów taką jak w amplitudach drella, z wymianą pi i sigma
       PRF2=1.d0/(T-EMF2**2)

       FF_NNF=(LBDFNN**2-EMF2**2)/(LBDFNN**2-T)
       FF_GAMROF=(LBDFROGAM**2-EMF2**2)/(LBDFROGAM**2-T)

       DO 105 M=1,3
       AMM(M)=C0(M)*ST-CX(M)*CT-CY(M)*CT+CZ(M)*ST

       BMM(M)=R*BWIG*PRF2*FF_NNF*FF_GAMROF*AMM(M)
105    CONTINUE

       RETURN
       END

       complex*16 function a_ro_f2ex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)
       implicit double precision (f-h,o-z)
       complex*16 ampl(3),ci
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       
       ci=DCMPLX(0.d0,1.d0)

       a_ro_f2ex=dcmplx(0.d0,0.d0)

!        call this subroutine to initialise constants       
       call InitConstants(es,et,empipi)

       if (.not. (lbd_g.eq.1 .or. lbd_g.eq.-1)) then
       print*,"Invalid photon helicity"
       return
       endif

       if (.not. (lbd_t.eq.1 .or. lbd_t.eq.-1)) then
       print*,"Invalid target helicity"
       return
       endif

       if (.not. (lbd_r.eq.1 .or. lbd_r.eq.-1)) then
       print*,"Invalid recoil helicity"
       return
       endif

       if (.not. (M.ge.-1 .and. M.le.1)) then
       print*,"Invalid rho helicity"
       return
       endif

       if (empipi.lt.2.d0*pim .or. empipi.gt.(dsqrt(s)-pm)) then
c       print*,"Dipion effective mass out of range"
c       print*,"1 s, empipi:", s, empipi
       return
       endif

!        initial cm momentum
       xk=(s-pm**2)/dsqrt(4.d0*s)
!        final cm momentum
       xq=dsqrt((s-(pm-empipi)**2)*(s-(pm+empipi)**2))/dsqrt(4.d0*s)
       tmax=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)-xq)
       tmin=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)+xq)
c       if (et.lt.tmin .or. et.gt.tmax) then
c       print*,"4-momentum transfer out of range"
c       return
c       endif


       LBDG=lbd_g
       if (lbd_r.eq.1 .and. lbd_t.eq.1) then
         call AP_F2_PP(ampl, empipi)
       else if (lbd_r.eq.-1 .and. lbd_t.eq.1) then
         call AP_F2_MP(ampl, empipi)
       else if (lbd_r.eq.1 .and. lbd_t.eq.-1) then
         call AP_F2_PM(ampl, empipi)
       else
         call AP_F2_MM(ampl, empipi)
       endif

c       print*,"ampl=",ampl

       if (M.eq.-1) then
         a_ro_f2ex=ampl(3)
       else if (M.eq.0) then
         a_ro_f2ex=ampl(2) 
       else
         a_ro_f2ex=ampl(1)
       endif
       return
       end!This is the ro photoproduction amplitude with t-channel pi exchange and
! subsequent ro->pi^+pi^-  decay
!The amplitude is based on the model by Oh and Lee Phys.Rev. C69, 025201 (2004)
C--------------------------------------------------------
       SUBROUTINE AP_PI(BPP,BMP,BPM,BMM)
       implicit double precision(a-h,o-z)
       COMPLEX*16 BPP(3),BMP(3),BPM(3),BMM(3),czyn(3),CI,BWIG,Gx,Gy,Gz,
     & DEL(3),G_T,APP(3),APM(3),AMP(3),AMM(3)
       DOUBLE PRECISION LAMBDA_PI, LAMBDA_ROPIGAM,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       LAMBDA_PI=0.6d0
       LAMBDA_ROPIGAM=0.77d0
       G1=6.04d0
! GPIGARO=0.756d0 - this is the value of the coupling cinstant used by Oh and Lee
       GPIGARO=0.639d0
! this is the coupl. const. value calculated based on the masses and widths from pdg 2010 and Eq. (13) of Oh and Lee
       GPINN=DSQRT(4.d0*PI*14.3d0)

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=EFFM
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1KW=E1**2-PM**2
       P1X=DSQRT(P1KW-P1ZET**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq

       qx=qxver*Q
       qz=qzver*Q

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       Gx=-CI*qz/DSQRT(2.d0)
       Gy=DBLE(LBDG)*Q/DSQRT(2.d0)
       Gz=CI*qx/DSQRT(2.d0)
       
       DEL(1)=-(Gx-CI*Gy)/DSQRT(2.d0)
       DEL(2)=Gz
       DEL(3)=(Gx+CI*Gy)/DSQRT(2.d0)

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       R1=(qx+P1X*(E1-PM)/(E1+PM))*DSQRT(EE)
       R2=(qz-P2MOD*(E1-E2)/(E2+PM)+P1ZET*(E1-E2)/(E1+PM)+
     & P2MOD*(E1**2-PM**2+P1ZET*P2MOD)/((E2+PM)*(E1+PM)))*
     & DSQRT(EE)

       FF_PINN=(LAMBDA_PI**2-PIM**2)/(LAMBDA_PI**2-T)
       FF_ROPIGAM=(LAMBDA_ROPIGAM**2-PIM**2)/(LAMBDA_ROPIGAM**2-T)

       FFR=FF_PINN*FF_ROPIGAM
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c-------------------- 
       G_T=CI*E*FFR*GPIGARO*GPINN/(2.d0*PM*EFFM)/(T-PIM**2)
     & /DSQRT(4.d0*PI)

       DO 101 M=1,3
       APP(M)=(R1*CT-R2*ST)*G_T
       AMP(M)=(-R1*ST-R2*CT)*G_T
       APM(M)=(-R1*ST-R2*CT)*G_T
       AMM(M)=(-R1*CT+R2*ST)*G_T

       BPP(M)=-2.d0*G1*EFFM*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*APP(M)*DEL(M)
       BMP(M)=-2.d0*G1*EFFM*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*AMP(M)*DEL(M)
       BPM(M)=-2.d0*G1*EFFM*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*APM(M)*DEL(M)
       BMM(M)=-2.d0*G1*EFFM*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*AMM(M)*DEL(M)
101    CONTINUE
       RETURN
       END

       SUBROUTINE AP_PI_PP(BPP, zmpipi)
       implicit double precision(a-h,o-y)
       COMPLEX*16 BPP(3),czyn(3),CI,BWIG,Gx,Gy,Gz,
     & DEL(3),G_T,APP(3)
       DOUBLE PRECISION LAMBDA_PI, LAMBDA_ROPIGAM,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       EM=0.30282d0

       LAMBDA_PI=0.6d0
       LAMBDA_ROPIGAM=0.77d0
       G1=6.04d0
! GPIGARO=0.756d0 - this is the value of the coupling cinstant used by Oh and Lee
       GPIGARO=0.639d0
! this is the coupl. const. value calculated based on the masses and widths from pdg 2010 and Eq. (13) of Oh and Lee
       GPINN=DSQRT(4.d0*PI*14.3d0)

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1KW=E1**2-PM**2
       P1X=DSQRT(P1KW-P1ZET**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq

       qx=qxver*Q
       qz=qzver*Q

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       Gx=-CI*qz/DSQRT(2.d0)
       Gy=DBLE(LBDG)*Q/DSQRT(2.d0)
       Gz=CI*qx/DSQRT(2.d0)
       
       DEL(1)=-(Gx-CI*Gy)/DSQRT(2.d0)
       DEL(2)=Gz
       DEL(3)=(Gx+CI*Gy)/DSQRT(2.d0)

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       R1=(qx+P1X*(E1-PM)/(E1+PM))*DSQRT(EE)
       R2=(qz-P2MOD*(E1-E2)/(E2+PM)+P1ZET*(E1-E2)/(E1+PM)+
     & P2MOD*(E1**2-PM**2+P1ZET*P2MOD)/((E2+PM)*(E1+PM)))*
     & DSQRT(EE)

       FF_PINN=(LAMBDA_PI**2-PIM**2)/(LAMBDA_PI**2-T)
       FF_ROPIGAM=(LAMBDA_ROPIGAM**2-PIM**2)/(LAMBDA_ROPIGAM**2-T)
       FFR=FF_PINN*FF_ROPIGAM
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c-------------------- 
       G_T=CI*EM*FFR*GPIGARO*GPINN/(2.d0*PM*EMPI)/(T-PIM**2)
     & /DSQRT(4.d0*PI)

       DO 102 M=1,3
       APP(M)=(R1*CT-R2*ST)*G_T

       BPP(M)=-2.d0*G1*EMPI*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*APP(M)*DEL(M)
102    CONTINUE
       RETURN
       END

      SUBROUTINE AP_PI_MP(BMP, zmpipi)
       implicit double precision(a-h,o-y)
       COMPLEX*16 BMP(3),czyn(3),CI,BWIG,Gx,Gy,Gz,
     & DEL(3),G_T,AMP(3)
       DOUBLE PRECISION LAMBDA_PI, LAMBDA_ROPIGAM,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       EM=0.30282d0

       LAMBDA_PI=0.6d0
       LAMBDA_ROPIGAM=0.77d0
       G1=6.04d0
! GPIGARO=0.756d0 - this is the value of the coupling cinstant used by Oh and Lee
       GPIGARO=0.639d0
! this is the coupl. const. value calculated based on the masses and widths from pdg 2010 and Eq. (13) of Oh and Lee
       GPINN=DSQRT(4.d0*PI*14.3d0)

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1KW=E1**2-PM**2
       P1X=DSQRT(P1KW-P1ZET**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq

       qx=qxver*Q
       qz=qzver*Q

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       Gx=-CI*qz/DSQRT(2.d0)
       Gy=DBLE(LBDG)*Q/DSQRT(2.d0)
       Gz=CI*qx/DSQRT(2.d0)
       
       DEL(1)=-(Gx-CI*Gy)/DSQRT(2.d0)
       DEL(2)=Gz
       DEL(3)=(Gx+CI*Gy)/DSQRT(2.d0)

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       R1=(qx+P1X*(E1-PM)/(E1+PM))*DSQRT(EE)
       R2=(qz-P2MOD*(E1-E2)/(E2+PM)+P1ZET*(E1-E2)/(E1+PM)+
     & P2MOD*(E1**2-PM**2+P1ZET*P2MOD)/((E2+PM)*(E1+PM)))*
     & DSQRT(EE)

       FF_PINN=(LAMBDA_PI**2-PIM**2)/(LAMBDA_PI**2-T)
       FF_ROPIGAM=(LAMBDA_ROPIGAM**2-PIM**2)/(LAMBDA_ROPIGAM**2-T)
       FFR=FF_PINN*FF_ROPIGAM
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c-------------------- 
       G_T=CI*EM*FFR*GPIGARO*GPINN/(2.d0*PM*EMPI)/(T-PIM**2)
     & /DSQRT(4.d0*PI)

       DO 103 M=1,3
       AMP(M)=(-R1*ST-R2*CT)*G_T

       BMP(M)=-2.d0*G1*EMPI*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*AMP(M)*DEL(M)
103    CONTINUE
       RETURN
       END

      SUBROUTINE AP_PI_PM(BPM, zmpipi)
       implicit double precision(a-h,o-y)
       COMPLEX*16 BPM(3),czyn(3),CI,BWIG,Gx,Gy,Gz,
     & DEL(3),G_T,APM(3)
       DOUBLE PRECISION LAMBDA_PI, LAMBDA_ROPIGAM,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       EM=0.30282d0

       LAMBDA_PI=0.6d0
       LAMBDA_ROPIGAM=0.77d0
       G1=6.04d0
! GPIGARO=0.756d0 - this is the value of the coupling cinstant used by Oh and Lee
       GPIGARO=0.639d0
! this is the coupl. const. value calculated based on the masses and widths from pdg 2010 and Eq. (13) of Oh and Lee
       GPINN=DSQRT(4.d0*PI*14.3d0)

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1KW=E1**2-PM**2
       P1X=DSQRT(P1KW-P1ZET**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq

       qx=qxver*Q
       qz=qzver*Q

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       Gx=-CI*qz/DSQRT(2.d0)
       Gy=DBLE(LBDG)*Q/DSQRT(2.d0)
       Gz=CI*qx/DSQRT(2.d0)
       
       DEL(1)=-(Gx-CI*Gy)/DSQRT(2.d0)
       DEL(2)=Gz
       DEL(3)=(Gx+CI*Gy)/DSQRT(2.d0)

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       R1=(qx+P1X*(E1-PM)/(E1+PM))*DSQRT(EE)
       R2=(qz-P2MOD*(E1-E2)/(E2+PM)+P1ZET*(E1-E2)/(E1+PM)+
     & P2MOD*(E1**2-PM**2+P1ZET*P2MOD)/((E2+PM)*(E1+PM)))*
     & DSQRT(EE)

       FF_PINN=(LAMBDA_PI**2-PIM**2)/(LAMBDA_PI**2-T)
       FF_ROPIGAM=(LAMBDA_ROPIGAM**2-PIM**2)/(LAMBDA_ROPIGAM**2-T)
       FFR=FF_PINN*FF_ROPIGAM
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c-------------------- 
       G_T=CI*EM*FFR*GPIGARO*GPINN/(2.d0*PM*EMPI)/(T-PIM**2)
     & /DSQRT(4.d0*PI)

       DO 104 M=1,3
       APM(M)=(-R1*ST-R2*CT)*G_T

       BPM(M)=-2.d0*G1*EMPI*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*APM(M)*DEL(M)
104    CONTINUE
       RETURN
       END

       SUBROUTINE AP_PI_MM(BMM, zmpipi)
       implicit double precision(a-h,o-y)
       COMPLEX*16 BMM(3),czyn(3),CI,BWIG,Gx,Gy,Gz,
     & DEL(3),G_T,AMM(3)
       DOUBLE PRECISION LAMBDA_PI, LAMBDA_ROPIGAM,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       EM=0.30282d0

       LAMBDA_PI=0.6d0
       LAMBDA_ROPIGAM=0.77d0
       G1=6.04d0
! GPIGARO=0.756d0 - this is the value of the coupling cinstant used by Oh and Lee
       GPIGARO=0.639d0
! this is the coupl. const. value calculated based on the masses and widths from pdg 2010 and Eq. (13) of Oh and Lee
       GPINN=DSQRT(4.d0*PI*14.3d0)

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       P1KW=E1**2-PM**2
       P1X=DSQRT(P1KW-P1ZET**2)

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       qxver=-sinthq
       qzver=costhq

       qx=qxver*Q
       qz=qzver*Q

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1
       Gx=-CI*qz/DSQRT(2.d0)
       Gy=DBLE(LBDG)*Q/DSQRT(2.d0)
       Gz=CI*qx/DSQRT(2.d0)
       
       DEL(1)=-(Gx-CI*Gy)/DSQRT(2.d0)
       DEL(2)=Gz
       DEL(3)=(Gx+CI*Gy)/DSQRT(2.d0)

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       R1=(qx+P1X*(E1-PM)/(E1+PM))*DSQRT(EE)
       R2=(qz-P2MOD*(E1-E2)/(E2+PM)+P1ZET*(E1-E2)/(E1+PM)+
     & P2MOD*(E1**2-PM**2+P1ZET*P2MOD)/((E2+PM)*(E1+PM)))*
     & DSQRT(EE)

       FF_PINN=(LAMBDA_PI**2-PIM**2)/(LAMBDA_PI**2-T)
       FF_ROPIGAM=(LAMBDA_ROPIGAM**2-PIM**2)/(LAMBDA_ROPIGAM**2-T)
       FFR=FF_PINN*FF_ROPIGAM
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c-------------------- 
       G_T=CI*EM*FFR*GPIGARO*GPINN/(2.d0*PM*EMPI)/(T-PIM**2)
     & /DSQRT(4.d0*PI)

       DO 105 M=1,3
       AMM(M)=(-R1*CT+R2*ST)*G_T

       BMM(M)=-2.d0*G1*EMPI*DSQRT(4.d0*PI/3.d0)*PIK*BWIG*AMM(M)*DEL(M)
105    CONTINUE
       RETURN
       END

       complex*16 function a_ro_piex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)
       implicit double precision (f-h,o-z)
       complex*16 ampl(3),ci
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG

       ci=DCMPLX(0.d0,1.d0)

       a_ro_piex=dcmplx(0.d0,0.d0)

!        call this subroutine to initialise constants       
       call InitConstants(es,et,empipi)

       if (.not. (lbd_g.eq.1 .or. lbd_g.eq.-1)) then
       print*,"Invalid photon helicity"
       return
       endif

       if (.not. (lbd_t.eq.1 .or. lbd_t.eq.-1)) then
       print*,"Invalid target helicity"
       return
       endif

       if (.not. (lbd_r.eq.1 .or. lbd_r.eq.-1)) then
       print*,"Invalid recoil helicity"
       return
       endif

       if (.not. (M.ge.-1 .and. M.le.1)) then
       print*,"Invalid rho helicity"
       return
       endif

       if (empipi.lt.2.d0*pim .or. empipi.gt.(dsqrt(s)-pm)) then
c       print*,"Dipion effective mass out of range"
c       print*,"2 s, empipi:", s, empipi
       return
       endif

!        initial cm momentum
       xk=(s-pm**2)/dsqrt(4.d0*s)
!        final cm momentum
       xq=dsqrt((s-(pm-empipi)**2)*(s-(pm+empipi)**2))/dsqrt(4.d0*s)
       tmax=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)-xq)
       tmin=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)+xq)
c       if (et.lt.tmin .or. et.gt.tmax) then
c       print*,"4-momentum transfer out of range"
c       return
c       endif


       LBDG=lbd_g
       if (lbd_r.eq.1 .and. lbd_t.eq.1) then
         call AP_PI_PP(ampl, empipi)
       else if (lbd_r.eq.-1 .and. lbd_t.eq.1) then
         call AP_PI_MP(ampl, empipi)
       else if (lbd_r.eq.1 .and. lbd_t.eq.-1) then
         call AP_PI_PM(ampl, empipi)
       else
         call AP_PI_MM(ampl, empipi)
       endif

       if (M.eq.-1) then
         a_ro_piex=ampl(3)
       else if (M.eq.0) then
         a_ro_piex=ampl(2) 
       else
         a_ro_piex=ampl(1)
       endif
       return
       end




C--------------------------------------------------------      
       SUBROUTINE AMPLP(BPP,BMP,BPM,BMM)
       implicit double precision(a-h,o-z)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & APP(3),APM(3),AMP(3),AMM(3)
       COMPLEX*16 BPP(3),BMP(3),BPM(3),BMM(3),CI,BWIG,PROP_POM,
     & EMULT,FDELTA
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       S0=4.0D0

       PIM=.13957d0
       EMPI=EFFM
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       P0P=(S-PM2-EMPI2)/2.d0/EMPI
       P0=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1=P0
       E2=P0P
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+P0*P0P
       PQ=(P0P**2-P0**2-Q**2)/2.d0/Q
       PPQ=PQ+Q
       PKDQ=2.d0*PIK/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 
       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
       ALFA=DBLE(LBDG)*P2MOD*QX/Q*EPLUS
       BETA=P2MOD*EMINUS+QZET/E1M
       GAMMA=DBLE(LBDG)*(P2MOD*QZET/Q*EMINUS+Q/E1M)
       DELTA=-QX/E1M
       AWSP=(CU+PP/EE)*EMPI
       BWSP=QX*P2MOD*EMPI/EE
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))
       CE1=-ALFA/DSQRT(2.D0)
       DX=-BETA/DSQRT(2.D0)
       DYMI=-GAMMA/DSQRT(2.D0)
       DZ=-DELTA/DSQRT(2.D0)
       RHO(1)=-QX/DSQRT(2.D0)
       RHO(2)=QZET
       RHO(3)=-RHO(1)
       ER=2.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)
       C=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-C**2)
       FFR=CU/(CU-T/0.7d0)**2

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)

       BWIG=REZ_M*REZ_GAM/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

       WU=1.0D0	

c--------------------   loop over M -phi spin projections
       ALFA_P=1.08D0+0.25D0*T
       PROP_POM=-(S/S0)**(ALFA_P-1.0D0)*CDEXP(-.5D0*CI*PI*
     & (ALFA_P-1.0D0))

       DO 101 M=1,3
       ALFAP(M)=AWSP*BE(M)+RHO(M)*CE1
       BETAP(M)=RHO(M)*DX
       GAMAP(M)=BWSP*BE(M)+RHO(M)*DYMI
       DELTAP(M)=RHO(M)*DZ
       APP(M)=ER*(ALFAP(M)*ST+BETAP(M)*C-GAMAP(M)*C-DELTAP(M)*ST)
!        AMP(M)=ER*(ALFAP(M)*C-BETAP(M)*ST+GAMAP(M)*ST-DELTAP(M)*C)
!        APM(M)=ER*(-ALFAP(M)*C-BETAP(M)*ST-GAMAP(M)*ST-DELTAP(M)*C)
! 2012.02.21 - changed to adjust to nucleon spin indexing convention in pion and sigma exchange
       APM(M)=ER*(ALFAP(M)*C-BETAP(M)*ST+GAMAP(M)*ST-DELTAP(M)*C)
       AMP(M)=ER*(-ALFAP(M)*C-BETAP(M)*ST-GAMAP(M)*ST-DELTAP(M)*C)
       AMM(M)=ER*(ALFAP(M)*ST-BETAP(M)*C-GAMAP(M)*C+DELTAP(M)*ST) 
c       MM=M
       
       BPP(M)=CI*PROP_POM*FFR*BWIG*APP(M)
       BMP(M)=CI*PROP_POM*FFR*BWIG*AMP(M)
       BPM(M)=CI*PROP_POM*FFR*BWIG*APM(M)
       BMM(M)=CI*PROP_POM*FFR*BWIG*AMM(M)
       
101    CONTINUE

       RETURN
       END

       SUBROUTINE AP_POM_PP(BPP,zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & APP(3)
       COMPLEX*16 BPP(3),CI,BWIG,PROP_POM,
     & EMULT,FDELTA
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

c       print*,"s, t, xm", S, T, zmpipi
c       print*,"stale", CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,REZ_GAM2,E,F,PIM

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       S0=4.0D0

       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       P0P=(S-PM2-EMPI2)/2.d0/EMPI
       P0=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1=P0
       E2=P0P
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+P0*P0P
       PQ=(P0P**2-P0**2-Q**2)/2.d0/Q
       PPQ=PQ+Q
       PKDQ=2.d0*PIK/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 
       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
       ALFA=DBLE(LBDG)*P2MOD*QX/Q*EPLUS
       BETA=P2MOD*EMINUS+QZET/E1M
       GAMMA=DBLE(LBDG)*(P2MOD*QZET/Q*EMINUS+Q/E1M)
       DELTA=-QX/E1M
       AWSP=(CU+PP/EE)*EMPI
       BWSP=QX*P2MOD*EMPI/EE
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))
       CE1=-ALFA/DSQRT(2.D0)
       DX=-BETA/DSQRT(2.D0)
       DYMI=-GAMMA/DSQRT(2.D0)
       DZ=-DELTA/DSQRT(2.D0)
       RHO(1)=-QX/DSQRT(2.D0)
       RHO(2)=QZET
       RHO(3)=-RHO(1)
       ER=2.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)
       C=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-C**2)
       FFR=CU/(CU-T/0.7d0)**2

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)

       BWIG=REZ_M*REZ_GAM/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

       WU=1.0D0

c--------------------   loop over M -phi spin projections
       ALFA_P=1.08D0+0.25D0*T
       PROP_POM=-(S/S0)**(ALFA_P-1.0D0)*CDEXP(-.5D0*CI*PI*
     & (ALFA_P-1.0D0))

       DO 102 M=1,3
       ALFAP(M)=AWSP*BE(M)+RHO(M)*CE1
       BETAP(M)=RHO(M)*DX
       GAMAP(M)=BWSP*BE(M)+RHO(M)*DYMI
       DELTAP(M)=RHO(M)*DZ
       APP(M)=ER*(ALFAP(M)*ST+BETAP(M)*C-GAMAP(M)*C-DELTAP(M)*ST)

       BPP(M)=CI*PROP_POM*FFR*BWIG*APP(M)
102    CONTINUE

       RETURN
       END

       SUBROUTINE AP_POM_MP(BMP,zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & AMP(3)
       COMPLEX*16 BMP(3),CI,BWIG,PROP_POM,
     & EMULT,FDELTA
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       S0=4.0D0

       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       P0P=(S-PM2-EMPI2)/2.d0/EMPI
       P0=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1=P0
       E2=P0P
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+P0*P0P
       PQ=(P0P**2-P0**2-Q**2)/2.d0/Q
       PPQ=PQ+Q
       PKDQ=2.d0*PIK/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       

       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 
       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
       ALFA=DBLE(LBDG)*P2MOD*QX/Q*EPLUS
       BETA=P2MOD*EMINUS+QZET/E1M
       GAMMA=DBLE(LBDG)*(P2MOD*QZET/Q*EMINUS+Q/E1M)
       DELTA=-QX/E1M
       AWSP=(CU+PP/EE)*EMPI
       BWSP=QX*P2MOD*EMPI/EE
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))
       CE1=-ALFA/DSQRT(2.D0)
       DX=-BETA/DSQRT(2.D0)
       DYMI=-GAMMA/DSQRT(2.D0)
       DZ=-DELTA/DSQRT(2.D0)
       RHO(1)=-QX/DSQRT(2.D0)
       RHO(2)=QZET
       RHO(3)=-RHO(1)
       ER=2.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)
       C=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-C**2)
       FFR=CU/(CU-T/0.7d0)**2
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)

       BWIG=REZ_M*REZ_GAM/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

       WU=1.0D0	

c--------------------   loop over M -phi spin projections
       ALFA_P=1.08D0+0.25D0*T
       PROP_POM=-(S/S0)**(ALFA_P-1.0D0)*CDEXP(-.5D0*CI*PI*
     & (ALFA_P-1.0D0))

       DO 103 M=1,3
       ALFAP(M)=AWSP*BE(M)+RHO(M)*CE1
       BETAP(M)=RHO(M)*DX
       GAMAP(M)=BWSP*BE(M)+RHO(M)*DYMI
       DELTAP(M)=RHO(M)*DZ
       AMP(M)=ER*(-ALFAP(M)*C-BETAP(M)*ST-GAMAP(M)*ST-DELTAP(M)*C)
       
       BMP(M)=CI*PROP_POM*FFR*BWIG*AMP(M)
103    CONTINUE

       RETURN
       END

       SUBROUTINE AP_POM_PM(BPM,zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & APM(3)
       COMPLEX*16 BPM(3),CI,BWIG,PROP_POM,
     & EMULT,FDELTA
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       S0=4.0D0

       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       P0P=(S-PM2-EMPI2)/2.d0/EMPI
       P0=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1=P0
       E2=P0P
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+P0*P0P
       PQ=(P0P**2-P0**2-Q**2)/2.d0/Q
       PPQ=PQ+Q
       PKDQ=2.d0*PIK/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 
       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
       ALFA=DBLE(LBDG)*P2MOD*QX/Q*EPLUS
       BETA=P2MOD*EMINUS+QZET/E1M
       GAMMA=DBLE(LBDG)*(P2MOD*QZET/Q*EMINUS+Q/E1M)
       DELTA=-QX/E1M
       AWSP=(CU+PP/EE)*EMPI
       BWSP=QX*P2MOD*EMPI/EE
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))
       CE1=-ALFA/DSQRT(2.D0)
       DX=-BETA/DSQRT(2.D0)
       DYMI=-GAMMA/DSQRT(2.D0)
       DZ=-DELTA/DSQRT(2.D0)
       RHO(1)=-QX/DSQRT(2.D0)
       RHO(2)=QZET
       RHO(3)=-RHO(1)
       ER=2.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)
       C=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-C**2)
       FFR=CU/(CU-T/0.7d0)**2

       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)

       BWIG=REZ_M*REZ_GAM/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

       WU=1.0D0	

c--------------------   loop over M -phi spin projections
       ALFA_P=1.08D0+0.25D0*T
       PROP_POM=-(S/S0)**(ALFA_P-1.0D0)*CDEXP(-.5D0*CI*PI*
     & (ALFA_P-1.0D0))

       DO 104 M=1,3
       ALFAP(M)=AWSP*BE(M)+RHO(M)*CE1
       BETAP(M)=RHO(M)*DX
       GAMAP(M)=BWSP*BE(M)+RHO(M)*DYMI
       DELTAP(M)=RHO(M)*DZ
       APM(M)=ER*(ALFAP(M)*C-BETAP(M)*ST+GAMAP(M)*ST-DELTAP(M)*C)
       
       BPM(M)=CI*PROP_POM*FFR*BWIG*APM(M)
       
104    CONTINUE

       RETURN
       END

       SUBROUTINE AP_POM_MM(BMM,zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),RHO(3),ALFAP(3),BETAP(3),GAMAP(3),DELTAP(3),
     & AMM(3)
       COMPLEX*16 BMM(3),CI,BWIG,PROP_POM,
     & EMULT,FDELTA
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/MASAEF/EFFM,CSIG
       COMMON/TRANS/T
       COMMON/POLAR/LBDG

       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       S0=4.0D0

       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       P0P=(S-PM2-EMPI2)/2.d0/EMPI
       P0=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1=P0
       E2=P0P
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+P0*P0P
       PQ=(P0P**2-P0**2-Q**2)/2.d0/Q
       PPQ=PQ+Q
       PKDQ=2.d0*PIK/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)
       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 
       EPLUS=1.D0/E1M+1.D0/E2M
       EMINUS=1.D0/E1M-1.D0/E2M
       EE=E1M*E2M
       ALFA=DBLE(LBDG)*P2MOD*QX/Q*EPLUS
       BETA=P2MOD*EMINUS+QZET/E1M
       GAMMA=DBLE(LBDG)*(P2MOD*QZET/Q*EMINUS+Q/E1M)
       DELTA=-QX/E1M
       AWSP=(CU+PP/EE)*EMPI
       BWSP=QX*P2MOD*EMPI/EE
C                                  1 M=1  2 M=0  3 M=-1
       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))
       CE1=-ALFA/DSQRT(2.D0)
       DX=-BETA/DSQRT(2.D0)
       DYMI=-GAMMA/DSQRT(2.D0)
       DZ=-DELTA/DSQRT(2.D0)
       RHO(1)=-QX/DSQRT(2.D0)
       RHO(2)=QZET
       RHO(3)=-RHO(1)
       ER=2.D0*DSQRT(EE)*PIK*DSQRT(4.D0*PI/3.D0)
       C=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-C**2)
       FFR=CU/(CU-T/0.7d0)**2
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)

       BWIG=REZ_M*REZ_GAM/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

       WU=1.0D0	

c--------------------   loop over M -phi spin projections
       ALFA_P=1.08D0+0.25D0*T
       PROP_POM=-(S/S0)**(ALFA_P-1.0D0)*CDEXP(-.5D0*CI*PI*
     & (ALFA_P-1.0D0))

       DO 105 M=1,3
       ALFAP(M)=AWSP*BE(M)+RHO(M)*CE1
       BETAP(M)=RHO(M)*DX
       GAMAP(M)=BWSP*BE(M)+RHO(M)*DYMI
       DELTAP(M)=RHO(M)*DZ
       AMM(M)=ER*(ALFAP(M)*ST-BETAP(M)*C-GAMAP(M)*C+DELTAP(M)*ST)

       BMM(M)=CI*PROP_POM*FFR*BWIG*AMM(M)
       
105    CONTINUE

       RETURN
       END

       complex*16 function a_ro_pomex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)
       implicit double precision (f-h,o-z)
       complex*16 ampl(3),ci
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       
       ci=DCMPLX(0.d0,1.d0)

       a_ro_pomex=dcmplx(0.d0,0.d0)

!        call this subroutine to initialise constants       
       call InitConstants(es,et,empipi)

       if (.not. (lbd_g.eq.1 .or. lbd_g.eq.-1)) then
       print*,"Invalid photon helicity"
       return
       endif

       if (.not. (lbd_t.eq.1 .or. lbd_t.eq.-1)) then
       print*,"Invalid target helicity"
       return
       endif

       if (.not. (lbd_r.eq.1 .or. lbd_r.eq.-1)) then
       print*,"Invalid recoil helicity"
       return
       endif

       if (.not. (M.ge.-1 .and. M.le.1)) then
       print*,"Invalid rho helicity"
       return
       endif

       if (empipi.lt.2.d0*pim .or. empipi.gt.(dsqrt(s)-pm)) then
c       print*,"Dipion effective mass out of range"
c       print*,"3 s, empipi:", s, empipi
       return
       endif

! !        initial cm momentum
!        xk=(s-pm**2)/dsqrt(4.d0*s)
! !        final cm momentum
!        xq=dsqrt((s-(pm-empipi)**2)*(s-(pm+empipi)**2))/dsqrt(4.d0*s)
!        tmax=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)-xq)
!        tmin=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)+xq)
!        if (et.lt.tmin .or. et.gt.tmax) then
!        print*,"4-momentum transfer out of range"
!        return
!        endif


c       print*,"empipi",empipi

       LBDG=lbd_g
       if (lbd_r.eq.1 .and. lbd_t.eq.1) then
         call AP_POM_PP(ampl, empipi)
       else if (lbd_r.eq.-1 .and. lbd_t.eq.1) then
         call AP_POM_MP(ampl, empipi)
       else if (lbd_r.eq.1 .and. lbd_t.eq.-1) then
         call AP_POM_PM(ampl, empipi)
       else
         call AP_POM_MM(ampl, empipi)
       endif

c       print*,"ampl=",ampl

       if (M.eq.-1) then
         a_ro_pomex=ampl(3)
       else if (M.eq.0) then
         a_ro_pomex=ampl(2) 
       else
         a_ro_pomex=ampl(1)
       endif
       return
       end!This is the ro photoproduction amplitude with t-channel sigma exchange and
! subsequent ro->pi^+pi^-  decay
!The amplitude is based on the model by Oh and Lee Phys.Rev. C69, 025201 (2004)
C--------------------------------------------------------
       SUBROUTINE AP_SIGMA(BPP,BMP,BPM,BMM)
       implicit double precision(a-h,o-z)
       DIMENSION BE(3),APP(3),APM(3),AMP(3),AMM(3)
       COMPLEX*16 BPP(3),BMP(3),BPM(3),BMM(3),CI,BWIG
       DOUBLE PRECISION LAMBDA_SIGMA, LAMBDA_ROSIGAM,MSIGMA,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       COMMON/MASAEF/effm,csig
       COMMON/TRANS/t

       LAMBDA_SIGMA=1.d0
       LAMBDA_ROSIGAM=0.9
       MSIGMA=.5d0
       G1=6.04d0
       GSIGARO=3.0d0
       GSIGNN=DSQRT(4.d0*PI*8.d0)
       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=EFFM
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+E1*E2
       PQ=(E2**2-E1**2-Q**2)/2.d0/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2
       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       OMEGA1=1.d0+P1ZET*P2MOD/EE
       OMEGA2=P1X*P2MOD/EE

       FF_SIGMANN=(LAMBDA_SIGMA**2-MSIGMA**2)/(LAMBDA_SIGMA**2-T)
       FF_ROSIGAM=(LAMBDA_ROSIGAM**2-MSIGMA**2)/(LAMBDA_ROSIGAM**2-T)
       FFR=FF_SIGMANN*FF_ROSIGAM

! masę rezonansu zastępuję masą efektywną 2012.01.26
       G_T=E*GSIGARO*GSIGNN/EFFM/(T-MSIGMA**2)/DSQRT(4.d0*PI)
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c--------------------   petla po M -rzutach spinu fi

       DO 101 M=1,3
       APP(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(OMEGA1*ST-OMEGA2*CT)
       AMP(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(-OMEGA1*CT-OMEGA2*ST)
       APM(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(OMEGA1*CT+OMEGA2*ST)
       AMM(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(OMEGA1*ST-OMEGA2*CT)

       BPP(M)=BWIG*Q*EFFM*PIK*DSQRT(4.d0*PI/3.d0)*APP(M)*BE(M)
       BMP(M)=BWIG*Q*EFFM*PIK*DSQRT(4.d0*PI/3.d0)*AMP(M)*BE(M)
       BPM(M)=BWIG*Q*EFFM*PIK*DSQRT(4.d0*PI/3.d0)*APM(M)*BE(M)
       BMM(M)=BWIG*Q*EFFM*PIK*DSQRT(4.d0*PI/3.d0)*AMM(M)*BE(M)
101    CONTINUE
       RETURN
       END

       SUBROUTINE AP_SIGMA_PP(BPP, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),APP(3)
       COMPLEX*16 BPP(3),CI,BWIG
       DOUBLE PRECISION LAMBDA_SIGMA, LAMBDA_ROSIGAM,MSIGMA,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       COMMON/MASAEF/effm,csig
       COMMON/TRANS/t

       EM=0.30282d0

       LAMBDA_SIGMA=1.d0
       LAMBDA_ROSIGAM=0.9
       MSIGMA=.5d0
       G1=6.04d0
       GSIGARO=3.0d0
       GSIGNN=DSQRT(4.d0*PI*8.d0)
       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+E1*E2
       PQ=(E2**2-E1**2-Q**2)/2.d0/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2
       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       OMEGA1=1.d0+P1ZET*P2MOD/EE
       OMEGA2=P1X*P2MOD/EE

       FF_SIGMANN=(LAMBDA_SIGMA**2-MSIGMA**2)/(LAMBDA_SIGMA**2-T)
       FF_ROSIGAM=(LAMBDA_ROSIGAM**2-MSIGMA**2)/(LAMBDA_ROSIGAM**2-T)
       FFR=FF_SIGMANN*FF_ROSIGAM

! masę rezonansu zastępuję masą efektywną 2012.01.26
       G_T=EM*GSIGARO*GSIGNN/EMPI/(T-MSIGMA**2)/DSQRT(4.d0*PI)
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c--------------------   petla po M -rzutach spinu fi

       DO 102 M=1,3
       APP(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(OMEGA1*ST-OMEGA2*CT)

       BPP(M)=BWIG*Q*EMPI*PIK*DSQRT(4.d0*PI/3.d0)*APP(M)*BE(M)
102    CONTINUE
       RETURN
       END

       SUBROUTINE AP_SIGMA_MP(BMP, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),AMP(3)
       COMPLEX*16 BMP(3),CI,BWIG
       DOUBLE PRECISION LAMBDA_SIGMA, LAMBDA_ROSIGAM,MSIGMA,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       COMMON/MASAEF/effm,csig
       COMMON/TRANS/t

       EM=0.30282d0

       LAMBDA_SIGMA=1.d0
       LAMBDA_ROSIGAM=0.9
       MSIGMA=.5d0
       G1=6.04d0
       GSIGARO=3.0d0
       GSIGNN=DSQRT(4.d0*PI*8.d0)
       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+E1*E2
       PQ=(E2**2-E1**2-Q**2)/2.d0/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2
       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       OMEGA1=1.d0+P1ZET*P2MOD/EE
       OMEGA2=P1X*P2MOD/EE

       FF_SIGMANN=(LAMBDA_SIGMA**2-MSIGMA**2)/(LAMBDA_SIGMA**2-T)
       FF_ROSIGAM=(LAMBDA_ROSIGAM**2-MSIGMA**2)/(LAMBDA_ROSIGAM**2-T)

       FFR=FF_SIGMANN*FF_ROSIGAM

! masę rezonansu zastępuję masą efektywną 2012.01.26
       G_T=EM*GSIGARO*GSIGNN/EMPI/(T-MSIGMA**2)/DSQRT(4.d0*PI)
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c--------------------   petla po M -rzutach spinu fi

       DO 103 M=1,3
       AMP(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(-OMEGA1*CT-OMEGA2*ST)

       BMP(M)=BWIG*Q*EMPI*PIK*DSQRT(4.d0*PI/3.d0)*AMP(M)*BE(M)
103    CONTINUE
       RETURN
       END

       SUBROUTINE AP_SIGMA_PM(BPM, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),APM(3)
       COMPLEX*16 BPM(3),CI,BWIG
       DOUBLE PRECISION LAMBDA_SIGMA, LAMBDA_ROSIGAM,MSIGMA,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       COMMON/MASAEF/effm,csig
       COMMON/TRANS/t

       EM=0.30282d0

       LAMBDA_SIGMA=1.d0
       LAMBDA_ROSIGAM=0.9
       MSIGMA=.5d0
       G1=6.04d0
       GSIGARO=3.0d0
       GSIGNN=DSQRT(4.d0*PI*8.d0)
       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+E1*E2
       PQ=(E2**2-E1**2-Q**2)/2.d0/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2
       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       OMEGA1=1.d0+P1ZET*P2MOD/EE
       OMEGA2=P1X*P2MOD/EE

       FF_SIGMANN=(LAMBDA_SIGMA**2-MSIGMA**2)/(LAMBDA_SIGMA**2-T)
       FF_ROSIGAM=(LAMBDA_ROSIGAM**2-MSIGMA**2)/(LAMBDA_ROSIGAM**2-T)
       FFR=FF_SIGMANN*FF_ROSIGAM

! masę rezonansu zastępuję masą efektywną 2012.01.26
       G_T=EM*GSIGARO*GSIGNN/EMPI/(T-MSIGMA**2)/DSQRT(4.d0*PI)
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c--------------------   petla po M -rzutach spinu fi

       DO 104 M=1,3
       APM(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(OMEGA1*CT+OMEGA2*ST)

       BPM(M)=BWIG*Q*EMPI*PIK*DSQRT(4.d0*PI/3.d0)*APM(M)*BE(M)
104    CONTINUE
       RETURN
       END

       SUBROUTINE AP_SIGMA_MM(BMM, zmpipi)
       implicit double precision(a-h,o-y)
       DIMENSION BE(3),AMM(3)
       COMPLEX*16 BMM(3),CI,BWIG
       DOUBLE PRECISION LAMBDA_SIGMA, LAMBDA_ROSIGAM,MSIGMA,FFR
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       COMMON/MASAEF/effm,csig
       COMMON/TRANS/t

       EM=0.30282d0

       LAMBDA_SIGMA=1.d0
       LAMBDA_ROSIGAM=0.9
       MSIGMA=.5d0
       G1=6.04d0
       GSIGARO=3.0d0
       GSIGNN=DSQRT(4.d0*PI*8.d0)
       CI=DCMPLX(0.d0,1.d0)
       CU=DCMPLX(1.d0,0.d0)
       PIM=.13957d0
       EMPI=zmpipi
       EMPI2=EMPI**2
       PIK=DSQRT(EMPI2/4.d0-PIM**2)
       E2=(S-PM2-EMPI2)/2.d0/EMPI
       E1=(S-PM2+T)/2.d0/EMPI
       Q=(EMPI2-T)/2.d0/EMPI
       E1M=E1+PM
       E2M=E2+PM
       PP=T/2.d0-PM2+E1*E2
       PQ=(E2**2-E1**2-Q**2)/2.d0/Q
       P2MOD=DSQRT(E2**2-PM**2)
       P1MOD=DSQRT(E1**2-PM**2)
       P1ZET=(-T/2.D0+PM**2-E1*E2)/P2MOD
       QZET=-(E2**2-E1**2+Q**2)/2.D0/P2MOD
       P1KW=E1**2-PM**2
       
       P1X=DSQRT(P1KW-P1ZET**2)
       QX=(PQ*Q-P1ZET*QZET)/P1X 

       EE=E1M*E2M
C                                  1 M=1  2 M=0  3 M=-1

       costhq=(E1**2-E2**2-Q**2)/(2.d0*Q*P2MOD)
       sinthq=DSQRT(1.d0-costhq**2)

       BE(1)=DBLE(kdelta(1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG-1))
       BE(2)=-DBLE(LBDG)*sinthq/DSQRT(2.D0)
       BE(3)=DBLE(kdelta(-1,LBDG))+
     & 0.5D0*(costhq-CU)*(-1.d0)**(.5d0*DBLE(LBDG+1))

       CT=DSQRT((1.D0+P1ZET/P1MOD)/2.D0)
       ST=DSQRT(1.D0-CT**2)

       OMEGA1=1.d0+P1ZET*P2MOD/EE
       OMEGA2=P1X*P2MOD/EE

       FF_SIGMANN=(LAMBDA_SIGMA**2-MSIGMA**2)/(LAMBDA_SIGMA**2-T)
       FF_ROSIGAM=(LAMBDA_ROSIGAM**2-MSIGMA**2)/(LAMBDA_ROSIGAM**2-T)
       FFR=FF_SIGMANN*FF_ROSIGAM

! masę rezonansu zastępuję masą efektywną 2012.01.26
       G_T=EM*GSIGARO*GSIGNN/EMPI/(T-MSIGMA**2)/DSQRT(4.d0*PI)
       
       PIRO=DSQRT(REZ_M2/4.d0-PIM**2)
       GAM_MPIPI=REZ_GAM*((PIK/PIRO)**3)*2.0d0/(1.0d0+(PIK/PIRO)**2)
       BWIG=1.d0/(REZ_M2-EMPI2-CI*REZ_M*GAM_MPIPI)

c--------------------   petla po M -rzutach spinu fi

       DO 105 M=1,3
       AMM(M)=2.d0*FFR*G_T*G1*DSQRT(EE)*(OMEGA1*ST-OMEGA2*CT)

       BMM(M)=BWIG*Q*EMPI*PIK*DSQRT(4.d0*PI/3.d0)*AMM(M)*BE(M)
105    CONTINUE
       RETURN
       END

       complex*16 function 
     & a_ro_sigmaex(es,et,empipi,lbd_g,lbd_t,lbd_r,M)
       implicit double precision (f-h,o-z)
       complex*16 ampl(3),ci
       COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
       COMMON/POLAR/LBDG
       
       ci=DCMPLX(0.d0,1.d0)

       a_ro_sigmaex=dcmplx(0.d0,0.d0)

!        call this subroutine to initialise constants       
       call InitConstants(es,et,empipi)

       if (.not. (lbd_g.eq.1 .or. lbd_g.eq.-1)) then
       print*,"Invalid photon helicity"
       return
       endif

       if (.not. (lbd_t.eq.1 .or. lbd_t.eq.-1)) then
       print*,"Invalid target helicity"
       return
       endif

       if (.not. (lbd_r.eq.1 .or. lbd_r.eq.-1)) then
       print*,"Invalid recoil helicity"
       return
       endif

       if (.not. (M.ge.-1 .and. M.le.1)) then
       print*,"Invalid rho helicity"
       return
       endif

       if (empipi.lt.2.d0*pim .or. empipi.gt.(dsqrt(s)-pm)) then
c       print*,"Dipion effective mass out of range"
c       print*,"4 s, empipi:", s, empipi
       return
       endif

!        initial cm momentum
       xk=(s-pm**2)/dsqrt(4.d0*s)
!        final cm momentum
       xq=dsqrt((s-(pm-empipi)**2)*(s-(pm+empipi)**2))/dsqrt(4.d0*s)
       tmax=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)-xq)
       tmin=empipi**2-2.d0*xk*(dsqrt(xq**2+empipi**2)+xq)
c       if (et.lt.tmin .or. et.gt.tmax) then
c       print*,"4-momentum transfer out of range"
c       return
c       endif


       LBDG=lbd_g
       if (lbd_r.eq.1 .and. lbd_t.eq.1) then
         call AP_SIGMA_PP(ampl, empipi)
       else if (lbd_r.eq.-1 .and. lbd_t.eq.1) then
         call AP_SIGMA_MP(ampl, empipi)
       else if (lbd_r.eq.1 .and. lbd_t.eq.-1) then
         call AP_SIGMA_PM(ampl, empipi)
       else
         call AP_SIGMA_MM(ampl, empipi)
       endif

       if (M.eq.-1) then
         a_ro_sigmaex=ampl(3)
       else if (M.eq.0) then
         a_ro_sigmaex=ampl(2) 
       else
         a_ro_sigmaex=ampl(1)
       endif
       return
       end
       integer function kdelta(M,N)
       if (M.eq.N) then
       kdelta=1
       else
       kdelta=0
       endif 
       end

      subroutine InitConstants(es,et,empipi)
      implicit double precision (f-h,o-z)
      COMMON /STALE/ CONV,PM,PM2,DS2,S,PI,REZ_M,REZ_GAM,REZ_M2,
     & REZ_GAM2,E,F,PIM
      COMMON/MASAEF/effm,csig
      COMMON/TRANS/t
      CONV=389.38d0  !conversion GeV to mikrobarns
      PM=0.93827d0 !proton rest mass
      PM2=PM**2
      DS2=SQRT(2.D0)
      S=es
      t=et
      effm=empipi
c      print*,"5 S, t, effm:", S, t, effm
      PI=3.141592654d0
      REZ_M=0.7755D0 !rho mass
      REZ_GAM=0.1494d0 ! rho width
      REZ_M2=REZ_M**2
      REZ_GAM2=(REZ_M*REZ_GAM)**2
      E=0.30282d0  !e-m coupling
      F=6.05d0
      PIM=.13957d0
      end
