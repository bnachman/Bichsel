        SUBROUTINE CONVLU 
C
C Main drive 
C
C Convolution program with realistic energy loss spectrum for silicon
C Link COV.FOR & convolution subroutines: COV3.FOR
C
        include 'bichinc/cov.cmm'
        logical first, newpar
        data    first/.true./ 
        data npmold/0/, bgold/0.0/
C
        newpar = .false.
        if( npm .ne. npmold ) newpar = .true. 
        if(  bg .ne. bgold  ) newpar = .true. 
C
        If( first .or. newpar ) then                  
          call PREP
          call PREPE                        
          npmold = npm 
          bgold  = bg
        Endif 
C
c
c Read in 1/epsilon  & GOS data   
C       Input files:
C       14      HEPS.TAB      in subr EPRED ( dielec.const. epsilon )
C       15      MACOM.TAB     in subr AERED ( GOS K,L shells  data  ) 
C       16      EMERC.TAB     in subr EMRED ( GOS  M  shell   data  )
c
        If( First ) then 
          call EPRED 
          call AERED
          call EMRED
        Endif 
C
c Calculate single collision spectrum 
c
        If( First .or. newpar ) call SPECT
        First = .false.
C
        call CONV
C 
        END
C
C---------------------------------------------------------------------
C 
        subroutine EPRED
C
C Routine to read in HEPS.TAB containing the tabulated result of 
C dielectric constant (epsilon). This is used for the cross section
C of transverse excitations.    
C
        include 'bichinc/cov.cmm'
C
        LUN = 15 
        open ( unit=LUN ,file='bichdat/heps.tab',status='old')
C
        read (LUN,*) n2t,numt
        if (nume .ne. numt) print*, ' CAUTION: nume & numt differ'
        if (n2 .ne. n2t) print*, ' CAUTION: n2 & n2t differ'
C
        do 4 j=1,numt
          read (LUN,*) jt,etbl,ep(1,j),ep(2,j),rim(j)
          dfdE(j) = rim(j) * 0.0092456 * E(j)
 4      continue
C
        Close( unit=LUN )
        return
        END
C 
C--------------------------------------------------------------------
C 
        subroutine AERED
C
C Routine to read in MACOM.TAB containing tabulated result of 
C Generalised Oscillation Strength (GOS) calculations for longitudinal 
C excitation function A(E).
C (Only K+L shells used eventually , M shell values replaced by EMRED).   
C
        include 'bichinc/cov.cmm'
C
        LUN = 15 
        open ( unit=LUN,file='bichdat/macom.tab', status='old' )
 
        read(LUN,*) n2t,numt
        do 4 j=1,numt
          read (LUN,*) jt,etbl,sig(6,j)
 4      continue
C
        Close( unit=LUN )
        return 
        END
C
C-----------------------------------------------------------------
C
        subroutine EMRED
C
C
C Routine to read in EMERC.TAB containing tabulated result of 
C Generalised Oscillation Strength (GOS) calculations for longitudinal 
C excitation function A(E). (M shell only to overwrite MACOM.TAB)   
C
        include 'bichinc/cov.cmm'
        character*8 tit(9)
C
        LUN = 15 
        open ( unit=LUN, file='bichdat/emerc.tab',status='old' )
C 
        do 1 j=1,4
                read (LUN,1615) tit
 1615           format (9a8)
 1      continue
        do 4 j=1,175
          read (LUN,*) jt,etbl,sig(6,j),xkmn(j)
 4      continue
C
        Close( unit=LUN )
        Return 
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
C                      = 6    muon 
C
        include 'bichinc/cov.cmm'
C
        Data pi/3.14159265359/, Ry/13.6058/
C
        Dimension PMASS(6)
        Data PMASS/938.256, 139.578, 3727.328 , 
     &            0.511004, 497.034,  105.658 /
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
        if(LEVDMP.gt.0) write (6,601) PTM,zi
 601    format (1x,'PREP : particle mass=',F10.3,
     &              ' MeV, charge=',f3.0)
C
C Initialization kinematic parameters
C 
        call EVANS
C
C  Saxon Eq (3a) for k
        Efin = Emax
        saxk = 153540. * zi**2 * rho / (betasq*AW)
C
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
C
C Evans, p 891, Uehling Eq (4a).  Date : 15 June 1984 
C
C We assume particle beta*gamma is set outside this routine 
C Calculate the following     W = Gamma  (E_total/mass) 
C                           pke = Particle kinetic energy 
c                          pmom = particle momentum (MeV)  
c                          ptE  = particle total energy (MeV)
c
        W = sqrt(bg**2 + 1.0)
        pkE = ptM * (W - 1.0) 
        pmom = ptM * bg
        betasq = bg**2 / (1 + bg**2)
        beta = bg / W
        gam  = W
        ptE  = ptM * W
C
C Maximum energy transfer   Emax  (MeV) 
C Uehling, also Sternheimer & Peierls Eq.(53)
C
C [SD Aug/90 Correct formula bug but minor effect ]
C
cc      telm = 2 * 0.511004
cc      Emax = ptM * (W**2 - 1) / (ptM/telm + telm/ptM + W)
c
        elm  = 2 * 0.511004
        Emax = ptM * (W**2 - 1) / (ptM/(2*elm) + elm/(2*ptM) + W)
C cannot distinguish i/o electrons
        if (npm .eq. 4) Emax = pkE / 2
C Emax in eV
        Emax = 1.e6 * Emax
c 
        if(LEVDMP.gt.0) write (6,333) bg,pmom
 333    format (1x,'EVANS: beta*gamma=',f11.4,
     &             '  momentum=',f13.4,' MeV/c')
C
        return
        END
C
C------------------------------------------------------------------------
C 
        subroutine PREPE
C
C Definitions of energy scale (log) bin size 
C Calculation make different treatment at different energy range: 
C Low  energy range:    Emin---THETAk ( ~1.8eV---1839eV          )  
C High energy range:  THETAk---Emax   ( 1839eV---depend on P_mom ) 
C
        include 'bichinc/cov.cmm'
C
C n2    = number of bins for each factor of 2 in energy 
C nume  = total number of energy bins in the low energy range   
C nhbin = number of high energy bins (bin szie=def low energy)   
C
        n2   = 64
        nume = 1250
        nhbin= 450  
        u    = log(2.) / N2
        um   = exp(u)
C
C [SD Aug/90] 
C This section has a problem as the result Emin=1.8 is actually due to 
C single precision rounding error while analytical calculation indicates
C that Emin should be == 1.5 
C However, since the data tables are created according to this scale
C nothing should be changed to preserve consistency.     
C
        ThetaK = 1839.
        ken  = log(ThetaK / 1.5) / u
        Emin = ThetaK / 2**(ken/N2)
        E(1) = Emin
        EXS  = 1.0 
C
C  E(L) = (Bin centre ?) energy (eV) of bin #L 
C DE(L) = linear energy bin width (eV) for bin #L 
C
        lemx = nume + Nhbin 
        do  11 L=1,lemx
          EXS = EXS * um
          E(L+1) = E(L) * um
          DI(L)  = -alog(1.0 - 1.0/EXS) / u
          DE(L)  = E(L+1)-E(L)
C
          if (L .le. nume) H(L)   = 0.
          if (E(L) .le. Emax) leh = L
 11     continue
C
        if (leh .gt. nume) leh = nume
        Etop = E(nume) * sqrt(um)
        if (Efin .gt. Etop) Efin = Etop
C
        if(LEVDMP.gt.0) Write(6,100) n2, nume, Emin, Etop, Emax 
100     format(1x,'PREPE: # energy bin for each factor of 2 =',I5,
     &       /,1x,'       # energy bins in low energy range =',I5,
     &       /,1x,'       Emin, Etop, Emax (eV) = ',3E12.4)
C
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
C 
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
c
ccAllow other particles with npm>4 to be added 
cc        if (npm .lt. 4) uef = 1 - E(j) * bemx     
        if (npm .ne. 4) uef = 1 - E(j) * bemx     
c
C  uef from Uehling Eqs. 9 & 2
        S0   = S0   + dfdE(j) * dE(j)
        avI  = avI  + dfdE(j) * alog(E(j)) * dE(j)
        avI1 = avI1 + dfdE(j) * E(j) * alog(E(j)) * dE(j)
        S1   = S1   + dfdE(j) * E(j) * dE(j)
        Q1 = Ry
c
cee  red CCS-33, 39 & 47
        if (E(j) .lt. 100.)  Q1 = 0.025**2 * Ry
        if (E(j) .lt. 11.9) Q1 = xkmn(j)**2 * Ry
        Qmin = E(j)**2 / tmcb
        sig(1,j) = 0
        if (E(j) .lt. 11.9 .and. Q1 .le. Qmin) go to 14
        sig(1,j) = E(j) * dfdE(j) * alog(Q1 / Qmin) 
 14     epbe = 1 - betasq * ep(1,j)
c
C  Fano Eq 47
        if (epbe .eq. 0) epbe = 1e-20
        sgg = E(j) * dfdE(j)*(-.5)*alog(epbe**2+(betasq*ep(2,j))**2)
        thet = atan (ep(2,j) * betasq / epbe)
        if (thet .lt. 0) thet = thet + pi         
c
C  plausible-otherwise I'd have a jump
C  Fano says [p 21]: 'arctan approaches pi for betasq*eps1 > 1'
        sgh = E(j)**2 *(betasq-ep(1,j)/(ep(1,j)**2+ep(2,j)**2))*thet
        sgh = 0.0092456 * sgh
        sig(3,j) = sgg + sgh
        sig(4,j) = 2. * sig(6,j) * uef
C the integral was over  d lnK rather than  d lnQ
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
 5      continue
 11     continue 
C
        If( LEVDMP.gt.0 ) Then 
          AIeff = exp(avI/S0) 
          Write(6,65) S0, AIeff
65        format(1x,'SPECT: Z_eff =',F10.4,'   I_eff =',F10.2, ' eV ')
        Endif
C
        If(LEVDMP.gt.1) then 
          Write(6,'(/a)') '    E(eV)      sig**E**2  1,3,4,5'
          Write(6,'(/a)') '    -----      ------------------'
          Do k=1,1250,10
            write(6,70) E(k),sig(1,k),sig(3,k),sig(4,k),sig(5,k) 
70          format(1x,5E12.4) 
          Enddo
        Endif
C
        FSG  = Tsig(5) * DEC
        dEdx = STP(5) * (dec/1.E6)
        rmf  = rM2(5) * (dec/1.E6)
        call SPTS
c
        RETURN
        END
C
C------------------------------------------------------------------
C 
        subroutine SPTS
        include 'bichinc/cov.cmm'
C
C Calculate corrected results to include high E-loss tail & 
C Print out summary results 
C
        do 75 j=1,nume
        if (E(j) .gt. Emax) go to 77
        nlast= j
        he2  = sig(5,j) * dec
        H(j) = he2 / E(j)**2
 75     continue
 77     continue
C
Cee FSR-99
C sbb is the <d1> over shells in the RPM paper as a correction factor 
C in A(E) applied to sig_u  ( 1 + d1/E + d2/E**2 )  d2 term negligible
C
        sbb = 720.
        bbb = 1. - sbb*betasq / Emax
        fft = Za * dec / 1e6
c 
c Residual cross section from Efin --> Emax 
cc bug 
cc      rm0 = bbb * ((1/Efin - 1/Emax) + 2 * (1/Efin**2 - 1/Emax**2))           
cc     1      - betasq * alog(Emax/Efin) / Emax                                 
 
        rm0 = bbb * ((1/Efin - 1/Emax) + 0.5 * (1/Efin**2 - 1/Emax**2))
     1      - betasq * alog(Emax/Efin) / Emax                                 
C
C Calculate corrected dE/dx up to Emax as the value of dedx is only
C summed to the maximum loss of Efin  
C 
        if (npm .ne. 4) rst = bbb * alog (Emax/Efin) + 
     1          sbb * (1./Efin - 1./Emax) - betasq*(1. - Efin/Emax)

        if (npm .ne. 4) go to 106 
cee FSR-143 and Uehling Eq 9
        TE = pkE * 1.e6
        rst=alog(Emax/Efin)+log(Emax)-alog(TE-Efin)  - 1/(1.-Efin/TE)
     1       + 2 + ((gam-1)/gam)**2 * (1./8. - .5 * (Efin/TE)**2)
     2       + ((2*gam - 1) / gam**2) * (alog(Emax) - alog(TE-Efin))
 
106     continue    
C
        tdedx = dedx + rst*fft
C
C 2nd moment   del2 = M2 - M2"
C
        secmv = DEC*RM2(5) / 1.e6
        rM2p  = Efin - 0.5 * betasq * Efin**2 / Emax
        del2  = secmv - fft * rM2p
C
C Print out results  
C
        If(LEVDMP.gt.0) then 
C
          write (6,'(a)') ' SPTS: summary information '
          write (6,610) FSG, rm0*fft                              
610       format(1x,'Total cross section  M0   = ',F12.4,' #coll/cm' ,          
     &         /,1x,'High E-loss residual M0   = ',F12.4,' #coll/cm' , 
     &         /,1x,'If residual M0 is large, look for error' ) 
          write(6,620) dedx, dedx/rho, tdedx, tdedx/rho 
620       format(1x,'Restricted dE/dx = ',F12.4,' MeV/cm   '   ,
     &                                    F12.4,' MeV*cm**2/g' ,  
     &         /,1x,'Corrected  dE/dx = ',F12.4,' MeV/cm   '   ,
     &                                    F12.4,' MeV*cm**2/g' )  
          write(6,630) secmv, fft*rM2p, del2
630       format(1x,'Summed 2nd moment     M2  = ',E12.4,' KeV**2/cm' ,
     &         /,1x,'Rutherford 2nd moment M2" = ',E12.4,' KeV**2/cm' ,  
     &         /,1x,'del-2  M2 - M2"           = ',E12.4,' KeV**2/cm' )
C
        Endif 
C
        RETURN
        END
