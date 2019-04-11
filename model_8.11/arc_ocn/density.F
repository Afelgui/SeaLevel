       subroutine DENSity(dens,Temp,Sal)
      Include 'model.par'
      real *8 z,dz,dzk,dx,dy,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      Real*8 TR,SR,SSR, P, RHOR, Cr,
     *       temp(nodes,kgr),sal(nodes,kgr),dens(nodes,kgr)
      integer nu_diag
      common/iodiag/nu_diag
      kz1=2
      DO 8 K=KZ1,KKZ
           K1=K-1
        zr=z(k)
     
      DO lj=1,kxLine
       n=IxD(1,lj)
        i1=IxD(2,lj)
          i2=IxD(3,lj)
      DO  m=I1,I2      
       i=kp(m,n)
       if(kbt(m,n).ge.k) then
       TR=temp(i,k-1)
       SR=sal(i,k-1)
       
        if(SR.lt.100.and.(SR.ge.38.or.SR.lt.0)) then
c        Write(nu_diag,*)'salt!!!>38',m,n,k1,sal(i,k-1)
c      read *
       sal(i,k-1)=(
     *        	  sal(kp(m-1,n),k-1)*inbt(m-1,n,k)
     *             +sal(kp(m+1,n),k-1)*inbt(m+1,n,k) 
     *             +sal(kp(m,n-1),k-1)*inbt(m,n-1,k) 
     *             +sal(kp(m,n+1),k-1)*inbt(m,n+1,k)
     *            )/(inbt(m-1,n,k)+inbt(m+1,n,k)
     *                +inbt(m,n-1,k)+inbt(m,n+1,k)) 
       end if
	if (SR.lt.100)then
       SSR=SQRT(SR)
c      P=0.802*(SR-35.)
c      DENS=(28.14-0.0735*TR-0.00469*TR*TR+P)
      P=9.801D0*1.025D0*ZR*0.0001
c           p=0.
C
      RHOR = 999.842594 - 1025.0
      RHOR = RHOR + (6.793952E-2 +
     $        (- 9.095290E-3 + (1.001685E-4 +
     $        (- 1.120083E-6 + 6.536332E-9*TR)*TR)*TR)*TR)*TR
cC
      RHOR = RHOR + SR*(
     *0.824493+
     *(-4.0899E-3+(7.6438E-5+(-8.2467E-7+5.3875E-9*TR)*TR)*TR)*TR+
     *((-5.72466E-3+(1.0227E-4-1.6546E-6*TR)*TR)+4.8314E-4*SSR)*SSR)
cC
      CR=1449.1+.0821*P+(4.55-.045*TR)*TR
     $                              +1.34*(SR-35.)
      RHOR=RHOR + 1.E5*P/CR**2 
     $     *(1.-2.*P/CR**2)
	 else
	RHOR=0.
	endif
       else
      RHOR=0.
       end if
      DENS(i,k-1)=RHOR
      end do
      end do
   8  continue
      RETURN
      END

      Real*8 Function DENS(TR,SR,ZR)
      Real*8 TR, SR, P, RHOR, CR,SSR
c      P=0.802*(SR-35.)
c      DENS=(28.14-0.0735*TR-0.00469*TR*TR+P)
      P=9.801D0*1.025D0*ZR*0.0001
c           p=0.
C
      SSR=sqrt(SR)
      RHOR = 999.842594 - 1025.0
      RHOR = RHOR + (6.793952E-2 +
     $        (- 9.095290E-3 + (1.001685E-4 +
     $        (- 1.120083E-6 + 6.536332E-9*TR)*TR)*TR)*TR)*TR
cC
      RHOR = RHOR + SR*(
     *0.824493+
     *(-4.0899E-3+(7.6438E-5+(-8.2467E-7+5.3875E-9*TR)*TR)*TR)*TR+
     *((-5.72466E-3+(1.0227E-4-1.6546E-6*TR)*TR)+4.8314E-4*SSR)*SSR)
cC
      CR=1449.1+.0821*P+(4.55-.045*TR)*TR
     $                              +1.34*(SR-35.)
      RHOR=RHOR + 1.E5*P/CR**2 
     $     *(1.-2.*P/CR**2)
cC
c      DENS=(1025+RHOR)*1.e-3
      DENS=RHOR
      RETURN
      END
