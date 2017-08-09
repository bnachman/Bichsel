c
c Plot Bichsel single collision xsec and moments
c
c plotsig setup the options  
c (all arguments in integer for kumac)  
c
c  Iptype  = 1   proton
c          = 2   pion
c          = 3   alpha
c          = 4   electron
c          = 5   kaon
c          = 6   muon   
c  MevP    = total momentum in MeV (integer) 
c  Indcomp = Xsec components 1-4; 5=total Xsec; 
c               6=Xsec integral; 7=Eloss integral; 8=2nd moment integral 
c  Levdump = Dump level   >8  to print 1 point per E point 
c
      Subroutine sigplot(Iptype, MevP, indcomp, Levdump)  
C
C Plot collision cross sections of Bichsel Calculation 
C
      Common/ INTFB / KPTYPE, betgam
      COMMON/PLTMOD/ INDP  
C
      dimension APMASS(6)
      Data APMASS/938.256, 139.578, 3727.328, 
     &           0.511004, 497.034,  105.658/
      REAL SCALOG(8)
      DATA SCALOG/1.0,1.0,4.0,9.0,4.0,1.0,2.2,9.0/
C
      KPTYPE = Iptype
      APMOM  = MevP
      INDP   = indcomp
      LEVDMP = Levdump
C
C Interface to H.Bichsel's program for the beta*gamma factor
C
      betgam = APMOM/apmass(KPTYPE) 
      return
      END  
C
C  INPUT:    XeV = energy loss in eV  
C OUTPUT: COVSIG = Collision cross-section 
C
      REAL FUNCTION COVSIG( XeV )    
      Common/ INTFB / KPTYPE, betgam
      COMMON/PLTMOD/ INDP 
      Common/ INTSP / SPINT0(1252),SPINT1(1252),SPINT2(1252)
C
      INCLUDE 'bichinc/cov.cmm'   
C  
      LOGICAL FIRST
      data     FIRST/.true./
      Parameter ( NEbin  = 1250 ) 
      Parameter ( Nrange = 25   ) 
      dimension Elow(25)
      dimension ISDEF(5)
      data ISDEF/5,1,3,4,2/
c
      character*8 pname(6) 
      data pname/'proton','pion','alpha','electron','kaon','muon'/  
C         
      COVSIG = 0 
      CALL SIGCAL( KPTYPE, betgam )  
      CALL SPECTINT
C
      If(FIRST) then  
c
c Normalisation factor = k * Na / beta**2
C k  = Rutherford X-sec formula (RPM eqn 2.2)  = 2.5496E-19 eV*cm**2 
C Na = Number of atoms/cm**3                   = 4.9938E+22 cm**-3  
c
        Cnorm = 4993.8 * 2.5496 
C
        NSKIP = NEbin/Nrange         
        do i=1,Nrange     
          IP = (i-1)*NSKIP + 1  
          if( IP.gt.NEbin ) IP = 1250  
          Elow(i) = E( IP )  
        enddo
        write(*,10) KPTYPE,pname(KPTYPE),betgam,FSG,dEdx
 10     format(1x,'Particle type   =',i3,3x,a8,
     &       /,1x,'beta*gamma      =',f9.3,
     &       /,1x,'Ncollision/cm   =',f9.2,
     &       /,1x,'dE/dx <1.346MeV =',f9.4,' MeV/cm')         
        First = .false. 
      Endif  
C
      IR  = Nrange + 1 
      Do i=1,Nrange
        if( XeV.lt.Elow(i) ) then
          IR = i 
          goto 40 
        endif
      Enddo
40    continue   
C
      Ibin = 0 
      If( IR.gt.1.and.IR.le.Nrange ) then        
        Istart = (IR-2)*Nskip + 1   
        Do i=Istart,NEbin-1 
          if(XeV.ge.E(i).and.Xev.lt.E(i+1)) then 
            Ibin = i 
            goto 50
          endif
        Enddo
      Endif
50    continue 
C 
      If( Ibin.gt.0 ) then 
       delE  = E(Ibin+1) - E(Ibin) 
       func1 = 0
       func2 = 0  
       fnorm = 1
       if(INDP.le.5) then  
         ISIND  = ISDEF(INDP)
         func1  = sig(ISIND,Ibin)  
         func2  = sig(ISIND,Ibin+1)  
       elseif(INDP.eq.6) then 
         func1  = SPINT0(Ibin)  
         func2  = SPINT0(Ibin+1)
       elseif(INDP.eq.7)   
         func1  = SPINT1(Ibin)  
         func2  = SPINT1(Ibin+1)
         fnorm  = 1.0E-6
       elseif(INDP.eq.8)   
         func1  = SPINT2(Ibin)  
         func2  = SPINT2(Ibin+1) 
         fnorm  = 1.0E-6 
       endif 
       beta2 = betgam**2/(1+betgam**2) 
       fact  = fnorm*Cnorm/beta2
       COVSIG = fact*(func1*(E(Ibin+1)-XeV)+func2*(XeV-E(Ibin)))/delE 
      Endif
C
      if( LEVDMP.gt.8 ) then 
        write(6,91) XeV,Ibin,COVSIG
91      format(1x,'COVSIG: XeV,Ibin,COVSIG =',E10.4,I6,E12.4)
      endif    
C
      Return
      End                      
c
c One time integral table 
c
      subroutine spectint
      Common/ INTSP / SPINT0(1252),SPINT1(1252),SPINT2(1252)
C
      INCLUDE 'bichinc/cov.cmm'   
c
      SPINT0(1) = 0 
      SPINT1(1) = 0 
      SPINT2(1) = 0      
      do i=2,1250
        delE = E(i) - E(i-1)  
        Eval = 0.5*(E(i)+E(i-1))
        Xsec = 0.5*(sig(5,i-1)+sig(5,i)) 
        SPINT0(i) = SPINT0(i-1) + delE*Xsec/(Eval*Eval)
        SPINT1(i) = SPINT1(i-1) + delE*Xsec/Eval
        SPINT2(i) = SPINT2(i-1) + delE*Xsec
      enddo 
c
      Return
      End
c
      include 'sigcal.for'
      include 'convolu.for'
      include 'cov3.for'
