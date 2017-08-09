        SUBROUTINE SIGCAL( IPTYPE, BETGAM ) 
C--------------------------------------------------------------------- 
C
C Interface subroutine to drive the Bichsel calculation of 
C single collision spectrum 
C INPUT: IPTYPE = 1   Proton
C               = 2   Pion 
C               = 3   Alpha 
C               = 4   Electron 
C        BETGAM = Incident Paerticle momentum factor  Beta*Gamma
C
C---------------------------------------------------------------------
C Author: Su Dong             Date: 28/June/89
C---------------------------------------------------------------------
C
        include 'bichinc/cov.cmm'
        logical first, newpar
        data first/.true./ 
        data npmold/0/, bgold/0.0/
C
C Interface particle type & momentum factor into internal common
C
        NPM = IPTYPE
        bg  = betgam    
C
C Check if the particle ID or momentum changed 
C
        newpar = .false.
        if( npm .ne. npmold ) newpar = .true. 
        if(  bg .ne. bgold  ) newpar = .true. 
C
C Set up basic parameters & kinematics parameters 
C
        If( first .or. newpar ) then                  
          call PREP
          call PREPE                        
          npmold = npm 
          bgold  = bg
        Endif 
C
c Read in 1/epsilon  & GOS data  (only do it once)  
C
C       Input files:
C       HEPS.TAB      in subroutine EPRED ( dielec.const. epsilon )
C       MACOM.TAB     in subroutine AERED ( GOS K,L shells  data  ) 
C       EMERC.TAB     in subroutine EMRED ( GOS  M  shell   data  )
C
c They are assumed to reside in directory with logical name BICHDAT: 
c which has to be assigned before running this program 
c
        If( First ) then 
          call EPRED 
          call AERED
          call EMRED
        Endif 
C
c Calculate single collision spectrum & corrected dE/dx  
c
        If( First .or. newpar ) call SPECT
C
        First = .false.
C
        RETURN  
        END
