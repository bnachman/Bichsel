C   COV3.FOR;             Convolution subroutines
 
        subroutine CONV
        include 'cov.cmm'

                write (3,303)
 303            format (/1x,50(' *'))
                write (3,*) '  CONV','  t=',exth,'  fsg=',fsg
        XMC = exth * FSG
C  number of collisions in thickness exth
                jxt = log(xmc) / log(2.0) + 1
                        print*, '  jxt=',jxt
                xx = xmc / 2.0**jxt
 86     CZ0 = xx / 2.**10
        NU  = jxt + 10 + 1
                write (3,386) jxt,NU,exth,XMC,XX,CZ0
 386            format (3x,'F.386: NU=',2i3,4x,'t=',f10.7,' cm',
     1          4X,'# coll=',F11.3,4X,'xx=',F9.5,4X,'CZ0=',1pe12.5)
                write (9,*) CZ0,Emax,betasq,pkE,bg,PTM,zi
                write (9,*) lemx,NU,U,um,Emin
        do 83 k=1,lemx,6
 83             write (9,384) (E(kl),kl=k,k+5)
 384            format (1p6e13.6)
        H0  = 0.
        CM1 = 1.0
        CM2 = 1.0
        XN  = 1.
        EX  = 1.e-15              
        MIE = 0                 
        MIF = 0
        MIH = 0
        call NORMAL
        xn = 1.
        D1 = CM1
        D2 = CM2 + CM1**2
        D3 = CM3 + 3.0*CM2*CM1 + CM1**3
        D4 = CM4 + 4.0*CM3*CM1 + 6.0*CM1**2*CM2 + CM1**4
        S  = D2 / D1
                write (3,612) D1,S,um
 612    format (2x,'conv  F.612:',3X,'Initial distribution',/,
     1  '  delta 1=',1PE12.4,'  delta 2 =',E12.4,'  exp(u)=',E12.5)
C               write (3,635) (H(J),J=1,nume,5)
C    set parameters from single collision spectrum
        EAV  = D1
        ONRM = CM0
        STPP = 2. * EAV * saxk * ONRM / ZA
                print*, ' Eav',d1,cm0,stpp,cma
        DK   = 2. * dEdx / STPP / ZA * CMA
        DQ   = DK / Emax - 1.
c               write (3,615) ONRM,CMA,EAV,STPP,DK,DQ
c615            format(/,' Orig spec - M0,M2,EAV,STP,M2,D2',/1P6E15.6)

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
                write (3,401)  NU,leh,Emin,EAV,DQ,Emax
 401            format(/2x,'F 401: ',2i5,1p5e14.7)
        do 16 N1 = 1,NU
                print*, 'convol: N1',N1
        thi  = 2. * thi
        xi   = 2. * xi
        rkap = 2. * rkap
        CN   = 2. * CN
                write (3,303)
                write (3,317) N1,leh,MIE,MIH,CN,CZ0
 317    format (/,' CONVOL NUMBER=',i3,' leh=',i4,2x,'MIE,MIH=',
     1  2i5,4x, 'mean collision number=',1PE12.4,'  CZ0=',f9.6,/)
C               write (3,635) (H(J),J=1,leh,5)
C635            format(' H= ',1P10E12.5)
        N2P = N2
        call FOLD
C  Landau-Vavilov parameters
        PB = -.4227843351 - log(rkap) - betasq
        dmmpl = xi * (pb + .225)
                write (3,310) thi,xi,rkap,pb,dmmpl
 310    format (/1x,'t =',1pe12.5,2x,'xi=',e12.5,2x,'kappa=',e12.5,
     1  2x,'<lam>=',0pf12.6/3x,'Landau theory <del>-dmp=',f12.3)
 631    N = MIH - MIE
c               write (3,408) N1,leh,CN,xi,rkap,PB,PZERO   ! PZERO undefined ?
c408            format (' CONV, F.408',2i5,f9.3,1p4e11.4/)
                write (3,408) N1,leh,CN,xi,rkap,PB 
 408            format (' CONV, F.408',2i5,f9.3,1p3e11.4/)
        call OUTPUT
 16     CONTINUE
        RETURN
        END
 
        subroutine OUTPUT               
 
        include 'cov.cmm'
 
        dimension ASP(1252),ASS(1252)
        character*22 TEE 

	call TIMDAT(TEE)
		print*, TEE
        B = CM2 / (CM1**2)
        C = CM3 / sqrt(CM2**3)
        D = CM4 / (CM2**2)
        S1 = CN * D1
                        write (3,*) ' OUT',b,c,d,s1,cn
        S2 = CN * D2 / S1**2
        S3 = D3 / sqrt(D2**3 * CN)
        S4 = D4 / (D2**2 * CN) + 3.
                write (3,350) H0,cm1,B,C,D,S1,S2,S3,S4
 350            format (2x,'OUTP: zero component =',1pe16.4/30x,'mean',
     1          9x,'variance/mean**2',4x,'central3/var**1.5',
     2          20H     central4/var**2/3x,'actual values=',4e20.4,/
     3          3x,'theoret values',4e20.4,/)
                write (3,340) cm1,s1,cm1/s1
                print 340, cm1,s1,cm1/s1
 340            format (1x,'OUTPUT:',1p2e13.6,'  ratio=',0pf10.6)
        X = 1.
        N = MIH - MIE
        if (N1 .le. nu-ners) RETURN
                write (9,*) leh,N1,N2,N2P,N,H0,thi,xi,rkap
                if (N2 .eq. N2P) go to 42
                write (9,*) lemx
        do 43 k=1,lemx,8
 43             write (9,950) (E(kl),kl=k,k+7)
 950            format (8f10.3)
 42     if (N2 .ne. N2P) N2P = N2
        ASP(1) = 0
        ASS(1) = 0
        hmax   = 0
        do 44 L=1,leh
        if (h(L) .lt. hmax) go to 53
        hmax = h(L)
        lmax = L
 53     J = L + N
        ASP(L+1) = ASP(L) + H(L) * DE(J)
 44     ASS(L+1) = ASS(L) + H(L) * E(J) * DE(J)
                print 346,  lmax+N,E(lmax+N),h(lmax)
 346            format (1x,'lmax+N=',i4,'  E=',f11.2,'  h=',f12.6)
                write (3,346) lmax+N,E(lmax+N),h(lmax)
        bax = (h(lmax) - h(lmax-1)) / (h(lmax) - h(lmax+1) )
        dmp = E(lmax+N) + 0.5 * dE(lmax+N) * (bax - 1.) / (bax + 1.)
        bax1 = (h(lmax-1) - h(lmax-2)) / (h(lmax-1) - h(lmax) )
        dmp1 = E(lmax+N-1) + 0.5 * dE(lmax+N-1)*(bax1 - 1.)/(bax1 + 1.)
        write (3,*) ' bax,dmp=',bax,dmp,'  lower:',bax1,dmp1
        hhun = 0.01 * h(lmax)
        do 55 L=1,leh
        if (L .lt. lmax .and. h(L) .lt. hhun) llow = L
        if (L .gt. lmax .and. h(L) .gt. hhun) lup  = L
 55     continue
        write (3,355) llow,E(llow+N),h(llow),lup,E(lup+N),h(lup)
        print 355, llow,E(llow+N),h(llow),lup,E(lup+N),h(lup)
 355    format (' lower & upper cutoff:',2(i6,0pf9.1,1pe12.4))
                write (9,*) lmax+N,E(lmax+N),h(lmax),' lmax+N, E, h'
                write (9,*) leh,dmp,'  leh, dmp'
        do 47 kl=1,leh,6
 47             write (9,944) (H(kk),kk=kl,kl+5)
 944            format (1p6e12.5)
        NPP   = (lup + llow) / 2
        nskip = 1
        npm = (lup - llow) / 2
        if (npm .gt. 70) nskip = 2
        npm = npm + nskip
                print*, ' N,npp,npm,nskip=',N,npp,npm,nskip
                write (3,344)
 344    format ( 2(7x,'j',4x,'E/eV',8x,'phi(E)',7x,'dE/dx',9x,'M2',3x))
        do 56 L=llow,NPP,nskip
        J  = L + N
        J2 = J + npm
        L2 = L + npm
 56     write (3,356) J,E(J),H(L),ASP(L),ASS(L),J2,E(J2),H(L2),
     1  ASP(L2),ASS(L2)
 356    format (3x,i5,0pF10.1,1p3E13.4,3x,i5,0pF10.1,1p3E13.4)
        RETURN
        END
 
        subroutine FOLD
C On 18 July 1984, I have some questions whether ST 62 is correct
        include 'cov.cmm'

        do 60 L=1,1250
        F(L) = H(L)
 60     H(L) = 0.
        F0  = H0
        LEF = leh
        MIF = MIH
        if (LEF .lt. 80) call RESET
        H0  = F0**2
        MIH = MIF + N2
        LEH = LEF
        if (leh .gt. 1250) leh = 1250
        print*,     '  FOLD: MIH,lef,leh,nume=',MIH,lef,leh,nume
        write (3,*) '  FOLD: MIH,lef,leh,nume=',MIH,lef,leh,nume
        do 61 LH=1,leh
        JH = LH + MIH
        do 62 LF=1,LH
        JF = LF + MIF
        K = JH - JF
        FLF = JH - MIF - DI(K) + 1.E-20
        LFF = FLF
        LE  = JF - MIE
        S   = FLF - LFF
        if (LFF .eq. 0) LFF = 1
 62     H(LH) = H(LH) + F(LF) * ((1.0-S)*F(LFF)+S*F(LFF+1)) * DE(LE)
 61     H(LH) = H(LH) - F(LH)**2 * 0.5*DE(LE)
        do 63 L=1,leh
 63     H(L) = H(L) * 2.0
        if (F0 .gt. EX) call ZERO
 65     call SHRINK
        call NORMAL      

        RETURN
        END
 
        subroutine RESET                
C   called from FOLD
        include 'cov.cmm'             
C   LEF = LEH initially
 
        if (N2 .ge. 128) go to 702
 701    N2 = N2 * 2
        U = log(2.) / float(N2)
                write (3,375) MIE,MIH,N2
 375    format  (1x,'RESET: ',2i4,'  Coordinate change: doubling of ',
     1          'point grid',i4,' points on a factor of 2',/)
                write (3,*) ' LEH,MIH,MIE=',LEH,MIH,MIE
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
                write (3,*) ' LEH,MIH,MIE=',LEH,MIH,MIE
        RETURN
        END
 
        subroutine ZERO
        include 'cov.cmm'
 
 81     xs = 0.
        N = MIH - MIE
                write (3,*) ' zero ',mih,mie,n,leh,F0
        do 82 L=1,leh
        K = L+N
 82     xs = xs + H(L)*DE(K)
        xs = (1. - F0)**2 / xs 
                write (3,*) k,xs,H(l),dE(K)
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
                print 684, N,LEF,leh,LEX,MIE,MIF,MIH
 684            format  (' zero:',i3,2(3x,3i5))
        do 86 M=1,LEF
 86     H(M) = H(M) + 2*F0*F(M)         
C   F(M) is the H of the previous convolution
C This seems to be correct (18 July 1984)
        RETURN
        END
 
        subroutine NORMAL
        include 'cov.cmm'
 
        Y = CM1 * 2.0
        Z = CM2 * 2.0
                write (3,714) MIE,MIH,H0,xn,y,z
 714  format(' NORM',2i5,'  H0=',1pe12.5,'  xn,y,z=',0p2f12.3,1pe12.5)
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
                write (3,635) CM0,CM1,CMQ
 635    format  (7x,'area=',1pe12.5,3x,'straight mean=',
     1          1pe12.5,3X,'CM1/CM0=',0pf14.4)
        if (cm0 - H0 .ne. 0) T = (1. - H0) / (cm0 - H0)
        CM1 = CM1*T
        CM2 = 0.
        CM3 = 0.
        CM4 = 0.
        do 31  L=1,leh
 31     H(L) = H(L) * T
                print*,     ' N,leh=',N,leh
                write (3,*) ' N,leh=',N,leh
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
 34             write (3,332) XN,CM0,CM1,CM2,CM3,CM4,Y,Z
 332            format  (7X,'Precision control, ',1P6E13.5,
     1          /24x,'mean=',e12.5,'   variance=',e12.5)
        RETURN
        END
 
        subroutine SHRINK
        include 'cov.cmm'
 
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
