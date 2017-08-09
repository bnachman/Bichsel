C
C------------------------------------------------------------------------
C Main interface routine to initialize & perform dE/dx spectrum folding 
C------------------------------------------------------------------------
C
        subroutine CONV
C
C  Convolution subroutines for actual dE/dx straggling function 
C 
        include 'bichinc/cov.cmm'
C 
        XMC = exth * FSG
C  number of collisions in thickness exth
                jxt = log(xmc) / log(2.0) + 1
                xx = xmc / 2.0**jxt
 86     CZ0 = xx / 2.**10
CC        NU  = jxt + 10 + 1        ! Bug ? 
        NU  = jxt + 10 
        write (6,386) jxt,NU,exth,XMC,CZ0                            
 386    format (/,1x,'# folds above 1 collision  jxt =',I4,
     &          /,1x,'Total # fold loops          NU =',I4,
     &          /,1x,'Silicon thickness            t =',f10.7,' cm',      
     &          /,1x,'Average # collisions       XMC =',F11.3,
     &          /,1x,'Fold starting thickness    CZ0 =',e12.5,' cm')
C
        H0  = 0.
        CM1 = 1.0
        CM2 = 1.0
        XN  = 1.
        EX  = 1.e-15              
        MIE = 0                 
        MIF = 0
        MIH = 0
C
        call NORMAL
C
        xn = 1.
        D1 = CM1
        D2 = CM2 + CM1**2
        D3 = CM3 + 3.0*CM2*CM1 + CM1**3
        D4 = CM4 + 4.0*CM3*CM1 + 6.0*CM1**2*CM2 + CM1**4
        S  = D2 / D1
C
C    set parameters from single collision spectrum
C
        EAV  = D1
        ONRM = CM0
        STPP = 2. * EAV * saxk * ONRM / ZA
        DK   = 2. * dEdx / STPP / ZA * CMA
        DQ   = DK / Emax - 1.

        call SHRINK
 
        H0 = 1. - CZ0
        do 15 L=1,leh
 15     H(L) = H(L) * CZ0         
C   first convolution: H0 + H(E)
        CN  = CZ0
        CM1 = CZ0 * D1
        CM2 = CZ0 * D2
        thi = exth / 2**(jxt+10)
        xi  = saxk * thi * Za
        rkap = xi / Emax
C
        do 16 N1 = 1,NU
C
        thi  = 2. * thi
        xi   = 2. * xi
        rkap = 2. * rkap
        CN   = 2. * CN
C
        write (6,317) N1,MIE,MIH,Leh,CN,H0                              
 317    format (' CONVOL #',i3,' MIE,MIH,Leh=',                
     1  3i5,2x,'<coll.#> =',E12.4,'  H0=',f8.5)               
C  
        N2P = N2
        call FOLD
C  Landau-Vavilov parameters
        PB = -.4227843351 - log(rkap) - betasq
        dmmpl = xi * (pb + .225)
C
 631    N = MIH - MIE
        call OUTPUT
 16     CONTINUE
        RETURN
        END
C
C---------------------------------------
C Routine to output fold results
C---------------------------------------
C
        subroutine OUTPUT               
 
        include 'bichinc/cov.cmm'
        common/ BICHSP / OUTE(1252), OUTH(1252),
     &                   OUTASP(1252),OUTASS(1252)   
        dimension ASP(1252), ASS(1252) 
C
        X = 1.
        N = MIH - MIE
        ASP(1) = 0
        ASS(1) = 0
        hmax   = 0
        do 44 L=1,leh
          if (h(L) .lt. hmax) go to 53
          hmax = h(L)
          lmax = L          
 53       J = L + N
          ASP(L+1) = ASP(L) + H(L) * DE(J)
          ASS(L+1) = ASS(L) + H(L) * E(J) * DE(J)
 44     continue 
C
C Only output points which is not too small compared to the max amplitude
C (previous default was 1%, changed to 0.2%  SD Sep/27/96) 
C
ccc     hhun = 0.01 * h(lmax)
        hhun = 0.002 * h(lmax) 
        llow = 1                  ! default SD Apr/05
        lup  = 1250               ! 
        do 55 L=1,leh  
        if (L .lt. lmax .and. h(L) .lt. hhun) llow = L 
        if (L .gt. lmax .and. h(L) .gt. hhun) lup  = L
 55     continue                                       
C
        NPP   = (lup + llow) / 2  
        nskip = 1                  
        Mpm = (lup - llow) / 2   
        if (Mpm .gt. 70) nskip = 2     
        Mpm = Mpm + nskip       
C
C Zero the interfacing output arrays 
C
        do 57 L=1,1252
          OUTE(L) = 0.0 
          OUTH(L) = 0.0 
          OUTASP(L) = 0.0 
          OUTASS(L) = 0.0 
 57     continue        
C
C Fill the interfacing output arrays  
C
        write(*,99) llow,NPP
 99     format('llw,npp=',2i10)
        do 60 L=llow,NPP,nskip     
          J  = L + N                             
          J2 = J + Mpm                                 
          L2 = L + Mpm    
          OUTE(J)    = E(J)
          OUTH(J)    = H(L)
          OUTASP(J)  = ASP(L)
          OUTASS(J)  = ASS(L)
          OUTE(J2)   = E(J2)
          OUTH(J2)   = H(L2)
          OUTASP(J2) = ASP(L2)
          OUTASS(J2) = ASS(L2)
 60     continue      
C
        RETURN
        END
C
C--------------------------------------------------------------------
C Routine to do self-folding of dE/dx spectrum at current thickness
C-------------------------------------------------------------------- 
C
        subroutine FOLD
C On 18 July 1984, I have some questions whether ST 62 is correct
        include 'bichinc/cov.cmm'
 
        do 60 L=1,1250
        F(L) = H(L)
 60     H(L) = 0.
        F0  = H0
        LEF = leh
        MIF = MIH
C
C Magnify energy scale if spectrum too narrow 
C
        if (LEF .lt. 80) call RESET
C
C Shift starting bin number by N2 for folded spectrum = double energy 
C
        H0  = F0**2
        MIH = MIF + N2
        LEH = LEF
        if (leh .gt. 1250) leh = 1250
C
C Loop over new spectrum energy points (DELTA)  
C
        do 61 LH=1,leh
        JH = LH + MIH
C
C Loop over old energy spectrum and integrate to half way  
C   H(DELTA) = Integral [ F(DELTA-x) * F(x) * dx ]     x=0,DELTA/2 
C
        do 62 LF=1,LH
        JF = LF + MIF                        ! x energy bin
        K = JH - JF
        FLF = JH - MIF - DI(K) + 1.E-20      ! DELTA-x energy bin 
        LFF = FLF
        LE  = JF - MIE
        S   = FLF - LFF
CC        if (LFF .eq. 0) LFF = 1     ! SD 3/Sept/90
        If( LEF.eq.0 ) then 
          LEF = 1 
          S   = 0.0 
        Endif 
 62     H(LH) = H(LH) + F(LF) * ((1.0-S)*F(LFF)+S*F(LFF+1)) * DE(LE)
C
C The DELTA/2 bin should not be double counted when x2 in next stage  
C
 61     H(LH) = H(LH) - F(LH)**2 * 0.5*DE(LE)
C
C Other half of the integral is symmetric, just x2 
C
        do 63 L=1,leh
 63     H(L) = H(L) * 2.0
C
C Take into account the no-collision probability if it is still significant 
C
        if (F0 .gt. EX) call ZERO
C
C Trim out negligible tails in the spectrum & normalize 
C
 65     call SHRINK
        call NORMAL
C
        RETURN
        END
C
C------------------------------------------------------------------
C Routine to Reset/magnify energy scale in case spectrum too narrow   
C------------------------------------------------------------------
C
        subroutine RESET                
C   called from FOLD
        include 'bichinc/cov.cmm'             
C   LEF = LEH initially
 
        if (N2 .ge. 128) go to 702
 701    N2 = N2 * 2
        U = log(2.) / float(N2)
C
        do 72 LL=1,LEF                  
C   this is from top down
        L = LEF + 1 - LL
 72     F(2*L) = F(L)
        N = 2*LEF
        do 73 L=4,N,2                   
C   N is just some number
73      F(L-1) = (F(L) + F(L-2)) / 2.
        LEF = 2 * LEF + 1
        MIF = 2 * MIF
        MIE = MIF
        do 74 J=1,leh
        S     = float(J+MIE)
        E(J)  = exp(S*U) * Emin
        DE(J) = E(J)*U
        S = J
 74     DI(J) = - log(1. - exp(-S*U)) / U
 702    MIE = MIF                       
C   from 4th Line above
        do 76  J=1,leh
        S     = float(J+MIE)
        E(J)  = exp(S*U) * Emin
 76     DE(J) = E(J)*U
        RETURN
        END
C
C-------------------------------------------------------------------
C Routine to normalize just folded dE/dx spetrum to include the  
C effect of no-collision probability for extremely thin layers  
C-------------------------------------------------------------------
C
        subroutine ZERO
        include 'bichinc/cov.cmm'
 
 81     xs = 0.
        N = MIH - MIE
C
        do 82 L=1,leh
        K = L+N
 82     xs = xs + H(L)*DE(K)
        xs = (1. - F0)**2 / xs 
C
        do 83 L=1,leh
 83     H(L) = H(L) * xs
        N   = MIH - MIF
        MIH = MIF
        LEH = LEH+N
        if (leh .gt. nume) leh = nume
        LEX = LEH + 1
        do 84 LL=1,leh
        LA = LEX - LL
        K  = LA + N
 84     H(K) = H(LA)
        do 85 LB=1,N
 85     H(LB) = 0.
C
        do 86 M=1,LEF
 86     H(M) = H(M) + 2*F0*F(M)         
C   F(M) is the H of the previous convolution
C This seems to be correct (18 July 1984)
        RETURN
        END
C
C-----------------------------------------------------------------
C Routine to normalize dE/dx spectrum to probability distribution 
C-----------------------------------------------------------------
C 
        subroutine NORMAL
        include 'bichinc/cov.cmm'
 
        Y = CM1 * 2.0
        Z = CM2 * 2.0
C
        CM0 = H0
        CM1 = 0.
        CMA = 0.
        N = MIH - MIE
        do 35 L = 1,leh
        LE  = L + N
        S   = H(L) * dE(LE)
        CM0 = CM0 + S
        CM1 = CM1 + S*E(LE)
 35     cma = cma + s*E(LE)**2
        if (cm0 .ne. 0) cmq = cm1 / cm0
C
        if (cm0 - H0 .ne. 0) T = (1. - H0) / (cm0 - H0)
        CM1 = CM1*T
        CM2 = 0.
        CM3 = 0.
        CM4 = 0.
        do 31  L=1,leh
 31     H(L) = H(L) * T
C
        do 30  L=1,leh
        LE = L + N
        EC = E(LE) - CM1
        S  = H(L) * DE(LE)
        CM2 = CM2 + S*EC**2
        CM3 = CM3 + S*EC**3
 30     CM4 = CM4 + S*EC**4
        XN  = XN * CM0
        if (Y .ne. 0) Y = CM1 / Y
        if (z .ne. 0) Z = CM2 / Z
        RETURN
        END
C
C----------------------------------------------------------------
C Shift down the dE/dx spectrum points to eliminate the leading 
C low energy bins with negligible E-loss probability 
C----------------------------------------------------------------
C 
        subroutine SHRINK
        include 'bichinc/cov.cmm'
 
        S = 0.
        N = MIH - MIE
        do 40 L=1,leh
        lla = L
        K = L + N
        S = S + H(L)*DE(K)
        if (S .gt. EX) GO TO 42
 40     CONTINUE
 42     M = lla - 1
        MIH = MIH + M
        S  = 0.
        LA = LEH + 1
        do 43 K=1,leh
        L  = LA - K
        KK = L + N
        S = S + H(L) * DE(KK)
        if (S - EX)  43,43,44
 43     CONTINUE
 44     LEH = L - M
        do 45  L=1,leh
        K = L + M
 45     H(L) = H(K)
        K = LEH + 1
        do 46  L=K,nume
 46     H(L) = 0.
        RETURN
        END
