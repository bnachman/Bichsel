C
C Original program from Hans Bichsel slightly modified I/O interface 
C to generate the energy loss distribution for a requested particle 
C and silicon thickness. 
C Output  COV.OPA     Control debugging output 
C         COV.SPE     Folded E-loss spectrum      
C
        program COV 
C
C Main drive 
C
C Convolution program with realistic energy loss spectrum for silicon
C Link COV.FOR & convolution subroutines: COV3.FOR
C
        include 'cov.cmm'
        character*22 TEE 
 
        open (unit=3, file='COV.OPA',status='new')
 
        call TIMDAT(TEE)
                write (3,600) 
 600            format (/9x,'COV: includes subroutines COV3.FOR')
                write (3,*) TEE
 
cc Don't know what ners does 
cc      print*, '$give ners (# of spectra printed):   '
cc              read*, ners
        ners = 1 
        call PREP
        call PREPE
        call EPRED
        call AERED
        call EMRED
        call SPECT
                if (ners .eq. 0) STOP
        open (unit=9, file='COV.SPE',status='new')
                write (9,*) ners,betasq,tdedx
                write (9,'(a)') TEE
                write (9,*) ' COV.SPE , t=',exth
        call CONV
C 
C       Input files:
C       14      heps.tab        in subr EPRED
C       15      MACOM.TAB       in subr AERED
C       16      EMERC.TAB       in subr EMRED
C
        END
