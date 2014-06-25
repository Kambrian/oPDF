function func(r,vr,vt)
		implicit none
		real*8::func,r,vr,vt
		real(kind=8),external::func_model
		integer*4, parameter:: MaxNumPars=10
		real*8::pars(MaxNumPars)=(/10.**7.34d0,15.2d0,0.715d0,69.014d0,2.301d0,7.467d0,0.d0, 0.d0, 0.d0, 0.d0/)
		func=func_model(r,vr,vt,pars)
end function 
      
function func_model(r,vr,vt, pars)
        implicit none
        integer(kind=8)::i,tag 
        real(kind=8)::beta,rs,r,vr,vt,gamma1,gamma2,CC,rou0,func_model,pi,G,coeff,epsilon,term2,term2a,term2b,rm,rmtmp,rmtmp2,tmp
        real(kind=8),external::func2,midsql,midsqu,midinf
		integer*4, parameter:: MaxNumPars=10
		real*8,intent(in)::pars(MaxNumPars)
	integer istat
	
	if(r<1. .or. r>1000.)then
	  func_model=0.
	  return
	end if
	
		rou0=pars(1)
		rs=pars(2)
		beta=pars(3)
		CC=pars(4)
		gamma1=pars(5)
		gamma2=pars(6)
! 		write(*,*) rou0,rs,beta,CC,gamma1,gamma2
        pi=acos(-1.)
        G=4.3007*10.**(-6)

        coeff=-rs**(2.*beta-gamma1-gamma2)*(r*vt)**(-2.*beta)*2.**(beta-1.5)*(rs**2*4.*pi*G*rou0)**(beta-1.5)/pi**1.5/gamma(beta+0.5)/gamma(1.-beta)

        epsilon=log(1.+r/rs)/(r/rs)-vr**2/rs**2/8./pi/G/rou0-vt**2/rs**2/8./pi/G/rou0

        if(epsilon>0.d0)then
           tag=0
           rmtmp=1.
           rmtmp2=2000.
           do while(abs(rmtmp2-rmtmp)/rmtmp>0.0000001d0)
              tag=tag+1
              rmtmp=rmtmp2
              rmtmp2=log(1.+rmtmp)/epsilon
              if(tag>5000)then
!                  write(*,*)'warning: reaching maximum step of iteration2! for (r,vr,vt, epsilon, rmtmp)=', r,vr,vt, epsilon, rmtmp
!                  func_model=0.
!                  return
                 exit
              end if
           end do
           rm=rmtmp2

           tmp=10.d0*rm
           call qromo(func2,rm,tmp,term2a,midsql,beta,rs,gamma1,gamma2,CC,epsilon,istat)
!            if(istat==0)then
! 	      write(*,*) 'at ', r,vr,vt
! 	      func_model=0.
! 	      return
!            end if
           call qromo(func2,tmp,10.d0**60,term2b,midinf,beta,rs,gamma1,gamma2,CC,epsilon,istat)
!            if(istat==0)then
! 	      write(*,*) 'at2 ', r,vr,vt
! 	      func_model=0.
! 	      return
!            end if
           term2=term2a+term2b
        else
           term2=0.
        end if
        func_model=coeff*term2
	
! 	write(*,*) r,vr,vt,log(func_model)
      end

      function func2(x,beta,rs,gamma1,gamma2,CC,epsilon)
        implicit none
        real(kind=8)::x,rs,beta,gamma1,gamma2,CC,epsilon,func2,term1,term2,term3,term4

        term1=(2.*beta+1.)*x**(2.*beta)*(x/(1.+x)-log(1.+x))-(1./(1.+x)**2-1./(1.+x))*x**(2.*beta+1.)
        term1=term1/(x/(1.+x)-log(1.+x))**2
        term2=(2.*beta-gamma1)*(x/CC)**gamma1*rs**(-gamma2)+(2.*beta-gamma2)*(x/CC)**gamma2*rs**(-gamma1)
        term2=term2/((x/CC)**gamma1*rs**(-gamma2)+(x/CC)**gamma2*rs**(-gamma1))**2
        term3=x**(2.*beta+1.)/(x/(1.+x)-log(1.+x))
        term4=(2.*beta-gamma1)*rs**(-gamma1-gamma2)*(gamma1/CC-2.*gamma2/CC)*(x/CC)**(gamma1+gamma2-1.)+(2.*beta-gamma2)*rs**(-gamma1-gamma2)*(gamma2/CC-2.*gamma1/CC)*(x/CC)**(gamma1+gamma2-1.)
        term4=term4-(2.*beta-gamma1)*rs**(-2.*gamma2)*gamma1/CC*(x/CC)**(2.*gamma1-1.)-(2.*beta-gamma2)*rs**(-2.*gamma1)*gamma2/CC*(x/CC)**(2.*gamma2-1.)
        term4=term4/((x/CC)**gamma1*rs**(-gamma2)+(x/CC)**gamma2*rs**(-gamma1))**3
        func2=(term1*term2+term3*term4)
        func2=func2*(epsilon-log(1.+x)/x)**(beta-0.5)

      end


      !  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.                                                                                                                                                                                  
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END


      SUBROUTINE midsql(funk,aa,bb,s,n,para1,para2,para3,para4,para5,para6)
      INTEGER n
      REAL*8 aa,bb,s,funk,para1,para2,para3,para4,para5,para6
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b

      func(x,para1,para2,para3,para4,para5,para6)=2.*x*funk(aa+x**2,para1,para2,para3,para4,para5,para6)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b),para1,para2,para3,para4,para5,para6)
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,para1,para2,para3,para4,para5,para6)
          x=x+ddel
          sum=sum+func(x,para1,para2,para3,para4,para5,para6)
          x=x+del
          !write(*,*)func(x,para1,para2,para3,para4,para5,para6),x,aa,bb,para1,para2,para3,para4,para5,para6                                                                                                                                 
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END


      SUBROUTINE midinf(funk,aa,bb,s,n,para1,para2,para3,para4,para5,para6)
      INTEGER n
      REAL*8 aa,bb,s,funk,para1,para2,para3,para4,para5,para6
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x
      func(x,para1,para2,para3,para4,para5,para6)=funk(1./x,para1,para2,para3,para4,para5,para6)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b),para1,para2,para3,para4,para5,para6)
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,para1,para2,para3,para4,para5,para6)
          x=x+ddel
          sum=sum+func(x,para1,para2,para3,para4,para5,para6)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END


      SUBROUTINE qromo(func,a,b,ss,choose,para1,para2,para3,para4,para5,para6,istat)
      INTEGER JMAX,JMAXP,K,KM,istat
      REAL*8 a,b,func,ss,EPS,para1,para2,para3,para4,para5,para6
      EXTERNAL func,choose
      PARAMETER (EPS=1.e-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
!U    USES polint                                                                                                                                                                                                                            
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      istat=1
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j,para1,para2,para3,para4,para5,para6)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      istat=0
!       write(*,*) 'too many steps in qromo'
      END


        !  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.                                                                                                                                                                                
      SUBROUTINE qromb(func,a,b,ss,para1,para2,para3,para4,para5,para6)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS,para1,para2,para3
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!U    USES polint,trapzd                                                                                                                                                                                                                     
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j,para1,para2,para3,para4,para5,para6)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
!      pause 'too many steps in qromb'                                                                                                                                                                                                       
      END


      !  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.                                                                                                                                                                                  
      SUBROUTINE trapzd(func,a,b,s,n,para1,para2,para3,para4,para5,para6)
      INTEGER n
      REAL*8 a,b,s,func,para1,para2,para3,para4,para5,para6
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a,para1,para2,para3,para4,para5,para6)+func(b,para1,para2,para3,para4,para5,para6))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,para1,para2,para3,para4,para5,para6)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.                 



