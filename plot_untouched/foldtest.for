C
C Test Bichsel's convolution method 
C [ intended for batch jobs ] 
C
C  Input 1   --- integer   particle type (1-6)  
C        2   --- real      partlcle momentum  (MeV) 
C        3   --- real      Silicon thickness  (micron) 
C
C Miscellaneous output goes to          stream   6
C Folded dE/dx spectrum table goes to   stream  21
C 
      subroutine foldtest 
C
      Include 'bichinc/cov.cmm'
      common/ BICHSP / OUTE(1252), OUTH(1252),                                
     &                 OUTASP(1252),OUTASS(1252)                              
C 
      Dimension PMASS(6)                                                      
      Data PMASS/938.256, 139.578, 3727.328 , 
     &          0.511004, 497.034, 105.865  / 
C
      Integer Lun
      data    Lun/21/ 
      character*64 fname  
      character*4  pname(6)
      data pname/'prtn','pion','alph','elec','kaon','muon'/
C
C Initialization 
C
C Particle types   NPM = 1    proton 
C                      = 2    pion   
C                      = 3    alpha  
C                      = 4    electron (poistron)
C                      = 5    kaon  
C                      = 6    muon
C        
      LEVDMP = 1 
      Write(*,'(a,$)') 'Particle (1=P 2=pi 3=alpha 4=e 5=K 6=mu) : '    
      read(5,*) NPM          ! Input particle type    
C
C Particle momentum     
C  bg = beta*gamme of incident particle 
C
      Write(*,'(a,$)') 'Momentum (MeV) : '                              
      read(5,*) iparmom      ! Input particle momentum in MeV 
      parmom = iparmom 
      pamass = Pmass(NPM)    ! MeV  
      bg     = parmom/pamass 
C
C Silicon thickness (cm) 
C read thickness in microns 
C
      Write(*,'(a,$)') 'Silicon thickness (Micron) : '          
      read(5,*) iexthmu      ! Input silicon thickness in microns 
      exth   = iexthmu/1.0E4 ! cm
C
C open output file
C
      write(fname,80) pname(NPM),iparmom,iexthmu 
 80   format('foldspec_',a4,'_',i4,'mev_',i3,'micron.dat')
      write(*,85) fname 
 85   format(/,'Output file: ',a,/) 
      open(Lun,file=fname,status='new')   
c 
      Call CONVLU 
C
      write(Lun,90)
90    format(/,1x,' L        E        ADC        H         ASP, ASS',/)
      Do 150 L=1,1252 
        Eprint = OUTE(L) 
        if( Eprint.eq.0.0) goto 150  
        ADC = OUTE(L)/(3.70*30.0) 
        write(Lun,100)L, OUTE(L), ADC, OUTH(L), OUTASP(L), OUTASS(L) 
100     format(1x,I4,E12.4,f8.3,3E12.4)
150   continue 
C
      close(Lun) 
c
      END
c
      include 'convolu.for'
      include 'cov3.for'
