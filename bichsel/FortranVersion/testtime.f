c
c Test UNIX version of date/time
c
	program testtime

C	returns time and date

	CHARACTER*22 TEE
        CHARACTER*25 datetime
C
cc	CALL TIME( timx )
cc	CALL DATE(  day )
        call fdate( datetime )
        call timdat(TEE)        
C
	WRITE (*,10) datetime, TEE 
10	FORMAT (1x,'fDate : ',a25,
     &        /,1x,'TEE   : ',a22) 
c	
	END
