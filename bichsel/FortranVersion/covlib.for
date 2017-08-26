C
C Basic utilitiy library extracted from the COV program 
C
C---------------------------------------------------------------------
C 
        subroutine EPRED
C
C Routine to read in heps.tab containing the tabulated result of 
C dielectric constant (epsilon). This is used for the cross section
C of small momentum transfer excitations.    
C
C
        include 'bichinc/cov.cmm'
        open (unit=14,file='bichdat/heps.tab',status='old')
                read (14,*) n2t,numt
                print*, ' EPRED, n2t,numt=',n2t,numt
        if (nume .ne. numt) print*, ' CAUTION: nume & numt differ'
        if (n2 .ne. n2t) print*, ' CAUTION: n2 & n2t differ'
        do 4 j=1,numt
                read (14,*) jt,etbl,ep(1,j),ep(2,j),rim(j)
        dfdE(j) = rim(j) * 0.0092456 * E(j)
c       if ((j/20)*20 .eq. j) print 304, j,jt,E(j),etbl,rim(j),dfdE(j)
c304            format (' EP:',2i4,2f11.2,2f12.6)
 4      continue
        END
C 
C--------------------------------------------------------------------
C 
        subroutine AERED
C
C Routine to read in MACOM.TAB containing tabulated result of 
C Generalised Oscillation Strength (GOS) calculations for longitudinal 
C excitation (K+L shells) with large momentum tansfers.   
C
        include 'bichinc/cov.cmm'
        open (unit=15,file='bichdat/macom.tab',status='old')
 
                read(15,*) n2t,numt
                print*, ' AERED, n2t,numt=',n2t,numt
        do 4 j=1,numt
                read (15,*) jt,etbl,sig(6,j)
c       if ((j/20)*20 .eq. j) print 304, j,jt,E(j),etbl,sig(6,j)
c304            format (' AE:',2i4,2f11.2,f12.6)
 4      continue
        END
C
C-----------------------------------------------------------------
C
        subroutine EMRED
C
C
C Routine to read in EMERC.TAB containing tabulated result of 
C Generalised Oscillation Strength (GOS) calculations for longitudinal 
C excitation (M shell) with large momentum tansfers.   
C
        include 'bichinc/cov.cmm'
        character*8 tit(9)
        open (unit=16,file='bichdat/emerc.tab',status='old')
 
        do 1 j=1,4
                read (16,1615) tit
 1615           format (9a8)
 1              print 1615, tit
        do 4 j=1,175
                read (16,*) jt,etbl,sig(6,j),xkmn(j)
        if (jt .lt. 20) go to 4
        if (jt .le. 31) print 304, j,jt,E(j),etbl,sig(6,j),xkmn(j)
c               if ((j/20)*20 .eq. j) print 304, j,jt,E(j),etbl,sig(6,j)
 304            format (' EM:',2i4,2f11.2,2f12.6)
 4      continue
        END
C
C---------------------------------------------------------------------
C 
        subroutine PREP
C
C Initialization routine for user parameter input & basic constants 
C Particle types   NPM = 1    proton
C                      = 2    pion
C                      = 3    alpha
C                      = 4    electron (poistron) 
C                      = 5    kaon
C
        include 'bichinc/cov.cmm'
C
        Data pi/3.14159265359/, Ry/13.6058/
C
        Dimension PMASS(5)
        Data PMASS/938.256, 139.578, 3727.328, 0.511004, 497.034/

        write (3,10) Ry
10      format(/,1x,'PREP:   Ry = ',F12.5,' eV')
C 
        write(*,'(a,$)') 
     &      'Particle type (1=P, 2=Pi, 3=Alpha, 4=e, 5=K) : '
        Read*, npm
        write(*,'(a,$)') 'Silicon thickness (microns) : '
        Read*, exth
C
C  convert microns into cm
C
        exth = exth / 1e4
C
C  set particle mass values (MeV) according to input code number  
C 
        PTM = PMASS( npm )
C
C  charge of incident particle
C
        zi = 1.
        if (npm .eq. 3) zi =2.
C
C properties of absober material (Silicon) 
C    ZA = atomic number 
C    AW = atomic weight
C   rho = density  (g/cm**3) 
C  atnu = # of atoms/cm**3 
C
        ZA = 14.
        AW = 28.086
        rho = 2.329
        atnu = 6.0222e23 * rho / Aw
C 
        write (3,601) PTM,zi
 601    format (1x,'PREP:  particle mass=',F10.3,' MeV, charge=',f3.0)
C
C Initialization kinematic parameters
C 
        call EVANS
C
        Efin = Emax
        saxk = 153540. * zi**2 * rho / (betasq*AW)
C  Saxon Eq (3a) for k
        Emk  = saxk * Emax
 
        write (3,602) Za,Aw,exth,saxk,Emk
 602    format (/,1x,'Z=',F6.2,'  A=',f9.4,3x,'t=',f9.5,'cm',
     1  '  k/Z=',f10.2,3x,'k*Emax/Z=',e12.5)
        return
        END
C
C----------------------------------------------------------------------
C 
        subroutine EVANS
C
C Initialization of kinematic parameters 
C
        include 'bichinc/cov.cmm'
C Evans, p 891, Uehling Eq (4a).  Date : 15 June 1984
 
      print*, ' Evans: particle mass (PTM) = ',ptM,' MeV'
      print*, ' Input option: 1=E_kin (MeV), 2=P_mom(MeV/c), 3=bet*gam'
      write(*,'(a,$)') 'Select option: '
      read*,    jkm
      write(*,'(a,$)') 'Input value  : '
      read*,    xxx
      go to (11,22,33), jkm
C
C Calculate the following     W = Gamma  (E_total/mass) 
C                           pKe = Particle kinetic energy 
C                            bg = beta*gamma
C
C-- input was kinetic energy
C 
 11     pkE = xxx
        W = xxx/ptM + 1.0
        bg = sqrt(W**2 - 1.)
                go to 34
c
c-- input was momentum  
c 
 22     pmom = xxx
        bg = xxx / ptM
        W = sqrt(bg**2 + 1.0)
        pkE = ptM * (W - 1.0)
        go to 35 
c
c-- input was beta*gamma 
c                
 33     bg = xxx                
        W = sqrt(bg**2 + 1.0)
        pkE = ptM * (W - 1.0)
c
c further parameters:  pmom = particle momentum (MeV)  
c                      ptE  = particle total energy (MeV)
 34     pmom = ptM * bg
 35     betasq = bg**2 / (1 + bg**2)
        beta = bg / W
        gam  = W
        ptE  = ptM * W
C
C Maximum energy transfer   Emax  (MeV) 
C Uehling, also Sternheimer & Peierls Eq.(53)
C
        telm = 2 * 0.511004
        Emax = ptM * (W**2 - 1) / (ptM/telm + telm/ptM + W)
C cannot distinguish i/o electrons
        if (npm .eq. 4) Emax = pkE / 2
C        
        print*, 'particle type = ', npm,'   Emax = ',Emax,' MeV'
        Emx  = telm * bg**2
c 
                write (3,333) bg,pmom,pkE
 333    format (/,3x,'beta*gamma=',f11.4,3x,'momentum=',f13.4,' MeV/c',
     1       3x,'E kinetic of incident particle=',f15.2,' MeV')
                write (3,334) betasq,gam,Emax,Emx
 334    format (3x,'beta**2=',f9.6,3x,'gamma=',f12.5,3x,'Emax=',
     1		2e12.4,' MeV'/)
C
C  Emax in eV
C
        Emax = 1.e6 * Emax
C
        return
        END
C
C------------------------------------------------------------------------
C 
        subroutine PREPE
C
C Definitions of energy scale (log) bin size 
C
        include 'bichinc/cov.cmm'
C
C n2 = number of bins for each factor of 2 in energy 
C
        n2   = 64
        nume = 650
        if (n2 .eq. 64) nume = 1250
        u    = log(2.) / N2
        um   = exp(u)
        eken  = 1.*log(1839. / 1.5) / u
        Emin = 1839. / 2**(eken/N2)
        E(1) = Emin
        EXS  = 1.
C
        print*, ' PREPE', n2,ken,E(1),emin
                write (3,609) N2,Emin,u,um
 609            format(/1x,'PREPE:   N2=',I4,3x,
     1  3x,'Emin=',f8.3,3X,'u=',F9.6,' e**u=',f10.6)
C
        lemx = nume + 450
        do  11 L=1,lemx
        EXS = EXS * um
        E(L+1) = E(L) * um
        if ((L/50)*50 .eq. L) print*, ' L,E=',L,E(L),exs,um
        DI(L)  = -alog(1.0 - 1.0/EXS) / u
        DE(L)  = E(L+1) - E(L)
        if (L .le. nume) H(L)   = 0.
        if (E(L) .le. Emax) leh = L
 11     continue
C
        if (leh .gt. nume) leh = nume
        Etop = E(nume) * sqrt(um)
        if (Efin .gt. Etop) Efin = Etop
                write (3,*) ' PREPE: Efin, Etop, Emax=',Efin,Etop,Emax
                print*, ' PREPE: Efin, Etop, Emax=',Efin,Etop,Emax
        RETURN
        END
C
C---------------------------------------------------------------
C 
        subroutine SPECT
C
C generate collision spectrum from ep-1,2 and ae
C
        include 'bichinc/cov.cmm'
C 
        elm = 511004.
        fac = 8. * pi * Ry**2 * (0.529177e-8)**2 / (elm * betasq)
        DEC = zi**2 * atnu * fac
 
        blg = alog ((2.*elm) * bg**2) - betasq
                write (3,307) betasq,atnu,blg
 307            format (/4X,'SPECT F.307:  beta**2=',F12.10,4X,'# of',
     2  ' atoms per cm**3=',e12.4,3x,'blg=',f9.4,/)
c               write (3,308)
c308    format (3x,'j',5x,'E/eV',5x,'df/dE ',5x,'sgg',6x,'sgh',
c    1  7x,'S 1',6x,'S 3',6x,'S 4',5x,'sum S',5x,'S(0)',3x,'dE/dx'/)
 
                jpr = 5
                jpd = 5
        do 3 L=1,5
        rM2(L)  = 0
        Tsig(L) = 0
 3      STP(L)  = 0
        S0 = 0
        S1 = 0
        avI = 0
        avI1 = 0
        bemx = betasq / Emax
        pf = pkE * 1.e6
        tmcb = 2. * elm * betasq
 
C do loop for Fano Eq 47
        do 5 j=1,nume
        if (E(j) .gt. Emax) go to 11
        if (npm .eq. 4) uef = 1 + (E(j)/(pf-E(j)))**2 + (((gam-1)/gam)
     1  * E(j)/pf)**2 - (2*gam - 1)*E(j)/(gam**2 * (pf - E(j)))
        if (npm .lt. 4) uef = 1 - E(j) * bemx     
C  uef from Uehling Eqs. 9 & 2
        if (j .eq. 1) print*, ' uef=',uef
        S0   = S0   + dfdE(j) * dE(j)
        avI  = avI  + dfdE(j) * alog(E(j)) * dE(j)
        avI1 = avI1 + dfdE(j) * E(j) * alog(E(j)) * dE(j)
        S1   = S1   + dfdE(j) * E(j) * dE(j)
        Q1 = Ry
cee  red CCS-33, 39 & 47
        if (E(j) .lt. 100.)  Q1 = 0.025**2 * Ry
        if (E(j) .lt. 11.9) Q1 = xkmn(j)**2 * Ry
        Qmin = E(j)**2 / tmcb
        sig(1,j) = 0
        if (E(j) .lt. 11.9 .and. Q1 .le. Qmin) go to 14
        sig(1,j) = E(j) * dfdE(j) * alog(Q1 / Qmin) 
 14     epbe = 1 - betasq * ep(1,j)
C  Fano Eq 47
        if (epbe .eq. 0) epbe = 1e-20
        sgg = E(j) * dfdE(j)*(-.5)*alog(epbe**2+(betasq*ep(2,j))**2)
        thet = atan (ep(2,j) * betasq / epbe)
        if (thet .lt. 0) thet = thet + pi         
C  plausible-otherwise I'd have a jump
C  Fano says [p 21]: 'arctan approaches pi for betasq*eps1 > 1'
        sgh = E(j)**2 *(betasq-ep(1,j)/(ep(1,j)**2+ep(2,j)**2))*thet
        sgh = 0.0092456 * sgh
        sig(3,j) = sgg + sgh
        sig(4,j) = 2. * sig(6,j) * uef
C the integral was over  d lnK rather than  d lnQ
c       if ((j/10)*10 .eq. j) print 327, E(j),sgg,sgh,(sig(ii,j),ii=1,4)
c327            format (1x,f11.2,1p6e11.3)
        sig(2,j) = 0
        sig(5,j) = 0
        do  27 L=1,4
                Tsig(L)  = Tsig(L)  + sig(L,j) * dE(j) / E(j)**2
                STP(L)   = STP(L)   + sig(L,j) * dE(j) / E(j)
                rM2(L)   = rM2(L)   + sig(L,j) *  dE(j)
 27             sig(5,j) = sig(5,j) + sig(L,j)
        Tsig(5) = Tsig(5) + sig(5,j) * dE(j) / E(j)**2
        STP(5)  = STP(5)  + sig(5,j) * dE(j) / E(j)  
        rM2(5)  = rM2(5)  + sig(5,j) * dE(j)
c               if (j .eq. 1) go to 28
c               if (j .ge. 320 .and. j .le. 326) go to 28
c               if ((j/10)*10 .ne. j) go to 5
c28     write (3,608) j,E(j),dfdE(j),sgg,sgh,sig(1,j),(sig(L,j),L=3,5),
c    1          S0,STP(5)
c608            format (1x,i4,f9.1,1pe11.3,0p9f9.4)
 5      continue
 11             write (3,*) '  uef=',uef
                write (3,374) Tsig,STP,rm2
 374    format (/9x,'Integ. over sig =',5F12.4/2(28x,5f12.3/))
                write (3,375) S0,avI,S1,avI1
 375    format (/9x,' S(0)=',f9.5,3x,'ln(I)=',f10.5,3x,
     1          'S(1)=',f10.3,3x,'L(1)=',f10.3/)
                write (3,*) '  following data without density effect'
                write (3,*) '  S(0)*blg=',S0*blg,'   2*L(0)=',2*avI
                write (3,*) '  S(1)*blg=',S1*blg,'   2*L(1)=',2*avI1
                print*, ' S(0)=',S0,'  L(0)=',avI
        FSG  = Tsig(5) * DEC
        dEdx = STP(5) * (dec/1.E6)
        rmf  = rM2(5) * (dec/1.E6)
                write (3,388) S0,FSG,dEdx,rmf
 388            format (/,10X,'Zeff=',F7.3,4X,'# coll/cm=',f11.3,4x,
     1          'dE/dx=',F9.4,' MeV/cm',3x,'M2=',f12.4,' keV**2/cm')

        write (3,*) ' DEC=',dec,'  # atoms/cm**3=',atnu,'  fac=',fac
        call SPTS
        RETURN
        END
C
C------------------------------------------------------------------
C 
        subroutine SPTS
        include 'bichinc/cov.cmm'
 
c               write (3,333)
c333    format (/4X,'SPTS,F.333:',/15X,'E',7x,'sig*E**2',6x,'sig',
c    1       10X,'M 0',11X,'M 1',11X,'M 2',10X,'<E>',/)
                SGM  = 0
                Stpw = 0
                SECM = 0
                jpr = 20
                ja  = 20
        do 75 j=1,nume
        if (E(j) .gt. Emax) go to 77
        nlast= j
        he2  = sig(5,j) * dec
        H(j) = he2 / E(j)**2
        SGM  = SGM + H(j)*dE(j)
        STPW = STPW + H(j) * E(j) * dE(j)
        SECM = SECM + he2 * dE(j)
        eps  = STPW / SGM
c               if (j .lt. 5) go to 11
c               if (j .eq. nume) go to 11
c               if (jpr .ne. j) go to 75
c               jpr = jpr + ja
c11             write (3,654) j,E(j),he2,H(j),SGM,STPW,SECM,eps
c654            format (1x,i6,f12.2,1p7e13.5)
 75     continue
 77             write (3,610)  nlast,SGM,STPW,SECM
 610    format (/1X,'SPTS: nlast  ',i5,' total cross section=',E15.5,
     1       3X,'dE/dx=',E15.5,', M2=',E15.5/2x,'see FSR-99 and CCS-9'/)
                write (3,*) ' final E=', Efin,' Emax=',Emax,' he2=',he2
Cee FSR-99
        sbb = 720.
        bbb = 1. - sbb*betasq / Emax
        fft = 14. * dec / 1e6
                write (3,*) ' sbb=',sbb,' eV,   bbb=',bbb,'  fft=',fft
        rm0 = bbb * ((1/Efin - 1/Emax) + 2 * (1/Efin**2 - 1/Emax**2))
     1        - betasq * alog(Emax/Efin) / Emax
                write (3,*) ' residual M0=',rm0,rm0*fft,'/cm'
                write (3,*) '  If residual M0 is large, look for error'
        call HPART(bbb,fft,sbb,stpw,secm)
        RETURN
        END
C
C-----------------------------------------------------------------------
C  
        subroutine HPART(bbb,fft,sbb,stpw,secm)
        include 'bichinc/cov.cmm'
 
        if (npm .lt. 4) rst = bbb * alog (Emax/Efin) + 
     1          sbb * (1./Efin - 1./Emax) - betasq*(1. - Efin/Emax)

        if (npm .lt. 4) go to 6 
cee FSR-143 and Uehling Eq 9
        TE = pkE * 1.e6
                print*,' ele',efin,emax,TE,gam
        rst=alog(Emax/Efin)+log(Emax)-alog(TE-Efin)  - 1/(1.-Efin/TE)
     1       + 2 + ((gam-1)/gam)**2 * (1./8. - .5 * (Efin/TE)**2)
     2       + ((2*gam - 1) / gam**2) * (alog(Emax) - alog(TE-Efin))
 
 6              write (3,*) ' residual dE/dx=',rst,rst*fft,' MeV/cm'
 
        print*,rst,fft,STPW
        tdedx = stpw/1.e6 + rst*fft
        write(3,*) ' dE/dx=',tdedx,' MeV/cm ',tdedx/2.329,' MeV cm**2/g'
                 print*, ' dE/dx=',tdedx,tdedx/rho

        secmv = secm / 1.e6
        rM2p  = Efin - 0.5 * betasq * Efin**2 / Emax
        del2  = secmv - fft * rM2p
        print *,    ' M2=',secmv,'  M2"=',fft*rM2p,' M2-M2"=',del2
        write (3,*) ' M2=',secmv,'  M2"=',fft*rM2p,' M2-M2"=',del2
        RETURN
        END
