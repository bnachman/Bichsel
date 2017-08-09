C
C Modified to UNIX fortran
C
	SUBROUTINE TIMDAT(TEE)

C	returns time and date

 	CHARACTER*22 TEE
        INTEGER MDATE(3), MTIME(3) 
        character*3 chmon(12) 
        data chmon/'JAN','FEB','MAR','APR','MAY','JUN',
     &             'JUL','AUG','SEP','OCT','NOV','DEC'/  
        character*2 chday, chhour, chmin, chsec 

C       CHARACTER*10 timx,day  
C
c	CALL TIME( timx )
c	CALL DATE(  day )
C
        call itime( mtime ) 
        call idate( mdate )           
C
        TEE = ' ' 
        call chfzero(mdate(1),chday)
        call chfzero(mtime(1),chhour)
        call chfzero(mtime(2),chmin)
        call chfzero(mtime(3),chsec)
c
	WRITE (TEE,50) chmon(mdate(2)),chday,mdate(3), 
     &                 chhour,chmin,chsec 
50	FORMAT (a3,'-',a2,'-',i4,1x,a2,':',a2,':',a2) 

	RETURN
	END
C
C Print numbers into characters with filled leading zeros
C
        subroutine chfzero(ival,chout)
        integer ival
        character*2 chout  
        if(ival.lt.10) then 
          write(chout,10) ival
 10	  format('0',i1)
        else
          write(chout,20) ival
 20	  format(i2)
        endif 
        return
        end
