c**********************************************************************************************
	subroutine rootf(al,xf)	 !Procedimento para a funcao f!
      implicit real*8 (a-h,o-z)
	complex*16 f
c
c             	  ! Varredura de raizes !
	  hinc=1.d-1
	  n=0 ! Contador das raizes para a funcao f !
	  do ai=1.d-2,25.d0,hinc !*inicio do looping de procura de raizes de f *!
	    fa=dreal(f(al,ai))
	    fb=dreal(f(al,ai+hinc))
		  if(fa*fb.lt.0.d0)then
		   n=n+1 ! Contador das raizes !
c	       write(5,*)'para l=',al
c		   write(5,*)'raiz entre',' r1=',ai,' e r2=',ai+hinc
	         ! Refinamento da raiz !
		     jj=0  ! Contador das iteracoes !
	         a=ai
	         b=ai+hinc
	         c=b-a
	         !precisao!
	         d=1.d-12
		     do while(c.ge.d)
	          c=b-a
	          xf=(a+b)/2.d0
	          fga=dreal(f(al,a))
	          fgx=dreal(f(al,xf))
	          e=fga*fgx
	          if(e.lt.0.d0)then
	           b=xf
	          else
	           a=xf
	          endif
	          jj=jj+1	 ! Contador das iteracoes !
	         enddo
	         xf=(a+b)/2.d0
	         fgx=dreal(f(al,xf))
	                !Torna as raizes positivas !
                if(fgx.lt.0.d0)then
                do while(fgx.lt.0.d0)
                xf=xf+1.d-13
	          fgx=dreal(f(al,xf))
                enddo
                else
                xf=xf
                fgx=dreal(f(al,xf))
                endif
c
	    ! Condicao lógica para a contagem de raizes !
	      else
	      n=n	           
		  endif
c	      
	  enddo	 !***** fim do loooping de procura de raizes de f *****!
c
	           if(n.eq.0)then
	             xf=0.d0
	           endif
c ************** Fim do procedimento para a funcao F ***************
c
      return
      end	  !Fim do procedimento para a funcao F!
c******************************************************************************************
