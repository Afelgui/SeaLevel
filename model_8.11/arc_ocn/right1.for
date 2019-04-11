           SUBROUTINE RIGHT1
      Include 'model.par'
      COMMON/NH2/ NH2
       real *8 Z,Dz,Dzk,dx,dy,h,h1
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
	Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Cf/UMY(kh),UMX(kh),RU(kh),ayv(nh),axv(mh)
      real *8 TaxV,TayV
      common/Wind/TaxV(mh,nh),TayV(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      COMMON/V1/ROZ(mh,nh)
      COMMON/V4/Risl2
      common/v5/grad_x(mh,nh),grad_y(mh,nh)
      common/v6/Fui1(mh,nh),Fvi1(mh,nh)
      common/clame/hl1(mh,nh),hl2(mh,nh)
 
      real *8 Fui1,Fvi1,hl1,hl2
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      REAL *8 FTOK,PSI0,PSI1,DZET,GANT,Gant0,Gant1
     *       ,SL,SM,SMB,SLB,grad1_x(mh,nh),grad1_y(mh,nh)
     *       , ROZ,Risl2,grad_x,grad_y,grx,gry
     *,      roz1(mh,nh),roz2(mh,nh),roz3(mh,nh)
      REAL *8 TAX(mh,nh),TAY(mh,nh)
      REAL *8 Q1,Q2,Q3,Q4,Q5,Q6,FF,ff1,ff2
      REAL *8 DX0,DX1,DY0,DY1
      DO 101 K=1,kxLine
      N=IXD(1,K)
      I1=IXD(2,K)
      I2=IXD(3,K)
      DO10 M=I1,I2
       if(h1(m,n).gt.0) then
        
      TAX(M,N)    =  TAXV(M,N)/H1(M,N)*(hl1(m,n)+hl1(m,n+1))/2.
      grad1_x(m,n) =  grad_x(m,n)*(hl1(m,n)+hl1(m,n+1))/2.
      fui1(m,n)    = fui1(m,n)*(hl1(m,n)+hl1(m,n+1))/2.
      TAY(M,N)    =  TAYV(M,N)/H1(M,N)*(hl2(m,n)+hl2(m+1,n))/2.
      grad1_y(m,n) = grad_y(m,n)*(hl2(m,n)+hl2(m+1,n))/2.
      fvi1(m,n)    = fvi1(m,n)*(hl2(m,n)+hl2(m+1,n))/2.
       end if
   10 continue
  101  continue
       DO 88 K=1,kxLine
      N=IXD(1,K)
      I1=IXD(2,K)
      I2=IXD(3,K)
      DO88 M=I1,I2
      DX0=DX(m-1,n)
      DX1=DX(m,n)
      DY0=DY(m-1,n)
      DY1=DY(m,n)
           Z1=iPH(m,n)
             Z2=iPV(m,n)
                Z3=iPH(m-1,n)
                  Z4=iPV(m-1,n-1)
                    Z5=iPH(m-1,n-1)
                      Z6=iPV(m,n-1)
      Q1=Z1*(ROZ(m,n)+ROZ(m+1,n+1)+ROZ(m+1,n))
      Q2=Z2*(ROZ(m,n)+ROZ(m+1,n+1)+ROZ(m,n+1))
      Q3=Z3*(ROZ(m,n)+ROZ(m-1,n)+ROZ(m,n+1))
      Q4=Z4*(ROZ(m,n)+ROZ(m-1,n)+ROZ(m-1,n-1))
      Q5=Z5*(ROZ(m,n)+ROZ(m,n-1)+ROZ(m-1,n-1))
      Q6=Z6*(ROZ(m,n)+ROZ(m,n-1)+ROZ(m+1,n))
      FF=H(m+1,n)*(Q6-Q1)+H(m+1,n+1)*(Q1-Q2)+H(m,n+1)*(Q2-Q3)+
     *H(m-1,n)*(Q3-Q4)+H(m-1,n-1)*(Q4-Q5)+H(m,n-1)*(Q5-Q6)
      z16=z1*z2*z3*z4*z5*z6
       ff=0
      roz1(m,n)=Z1*DY1*3*grad1_Y(m,n)
     *       -z3*DY1*3*grad1_Y(m-1,n)
     *       +Z6*DY0*3*grad1_Y(m,n-1)
     *       -Z4*DY0*3*grad1_Y(m-1,n-1)
     *     -Z2*DX1*3*grad1_X(m,n)
     *     +Z6*DX1*3*grad1_X(m,n-1)
     *     -Z3*DX0*3*grad1_X(m-1,n)
     *     +Z5*DX0*3*grad1_X(m-1,n-1)
      ff=0
      
      roz2(m,n)=Z1*DY1*3*TAY(m,n)
     *       -z3*DY1*3*TAY(m-1,n)
     *       +Z6*DY0*3*TAY(m,n-1)
     *       -Z4*DY0*3*tay(m-1,n-1)
     *     -Z2*DX1*3*TAX(m,n)
     *     +Z6*DX1*3*TAX(m,n-1)
     *     -Z3*DX0*3*TAX(m-1,n)
     *     +Z5*DX0*3*TAX(m-1,n-1)

      roz3(m,n)=Z1*DY1*3*fvi1(m,n)
     *       -z3*DY1*3*fvi1(m-1,n)
     *       +Z6*DY0*3*fvi1(m,n-1)
     *       -Z4*DY0*3*fvi1(m-1,n-1)
     *     -Z2*DX1*3*fui1(m,n)
     *     +Z6*DX1*3*fui1(m,n-1)
     *     -Z3*DX0*3*fui1(m-1,n)
     *     +Z5*DX0*3*fui1(m-1,n-1)
      roz1(m,n)=roz1(m,n)*z16
      roz2(m,n)=roz2(m,n)*z16*0
      roz3(m,n)=roz3(m,n)*z16
      roz(m,n)=roz1(m,n)+roz2(m,n)+roz3(m,n)
	      roz(m,n)=roz1(m,n)+roz2(m,n)+roz3(m,n)
   88 CONTINUE
      DO 111 K=1,kxLine
      N=IXD(1,K)
      I1=IXD(2,K)
      I2=IXD(3,K)
      DO121 M=I1,I2
        
      TAX(M,N)    =  TAXV(M,N)*(hl1(m,n)+hl1(m,n+1))/2.
      TAY(M,N)    =  TAYV(M,N)*(hl2(m,n)+hl2(m+1,n))/2.
  121 continue
  111 continue

      do m=61,110

      roz(m,2)=
     *	(TAY(m,2)
     *       -TAY(m-1,2))/dlx(m,2)+
     *     (-TAX(m,2+1)
     *     +TAX(m,2))/dly(m,2)
c     roz(m,2)=0
	 end do

      Risl2=0
      ff1=0
      ff2=0
c      open(27,file='c:/users/elen/new_arctic/roz2.bin',form='binary')
c      write(27)roz1,fui,fvi
c      close(27)
c      stop
      RETURN
      END