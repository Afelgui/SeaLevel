

         SUBROUTINE advx_V1(uu,vv)
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      COMMON/RUNPARS/nmyear
       COMMON/NH2/ NH2
      real *8 z,dz,dzk,dx,dy,h,h1,sl,slb
      COMMON/botV/H1(mh,nh)
      COMMON/botTS/ H(mh,nh)
      real *8 hreal
      COMMON/botReal/hreal(mh,nh)
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
            Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt1
      common/Cfsmag/cfsmag(mh,nh,kgr)
      common/Cf/UMY(nh),UMX(mh),RU(kh),ayv(nh),axv(mh)        
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/v6/Fui(mh,nh),Fvi(mh,nh),Fui1(mh,nh),Fvi1(mh,nh)	
	common/brt/UBrT(nodes), VBrT(nodes)

     
      real *8 uu(nodes,kgr),vv(nodes,kgr)
      Real *8 A1(1000),B1(1000),C1(1000),Fu1(1000),U1(1000),
     *        Fv1(1000),v1(1000)         
      REAL *8 dcx(mh,nh),dcy(mh,nh),gam(mh,nh),curv_uv(mh,nh)
     *        ,Uhol(mh,nh),Vhol(mh,nh),FT_st(mh,nh) 
     *        ,u(mh,nh),v(mh,nh),u0(mh,nh),v0(mh,nh)
     *        ,rx(mh,nh),ry(mh,nh),advX(mh),advY(nh),rzpole(nh)
     *        ,Fui,Fvi,aa1,cc1,bb1
     *        ,adx4,adx5,ady1,ady3,dz1,vmod,ubot,vbot,hr1
      integer nu_diag
      logical arc_special
      common/iodiag/nu_diag

C*-------------------------------------------
      arc_special=.true.




      rtropic=5.0e+8   ! viscosity for latitudes from -20 to +20 
      rmiddle=2.0e+7   ! viscosity for latitudes from 30 to 60
      rpole=1.0e+6     ! viscosity for Arctic

      rtropic=5.0e+9  ! viscosity for latitudes from -20 to +20 
      rmiddle=5.0e+8  ! viscosity for latitudes from 30 to 60
      rpole=5.0e+8   ! viscosity for Arctic

      do k=1,2
      rzpole(k)=2.e+7
	end do
	do k=3,kh
	rzpole(k)=rpole
      end do
c      test 4

      do k=1,kh
      rzpole(k)=5.e+5
	end do
      rmiddle=5.0e+5   ! viscosity for latitudes from 30 to 60


      dt=dt1
      tt=0.0
c      if(myear.eq.1954.and.month.eq.1)then
c         tt=0.5*(1.0+cos(nday*3.14159265/31.0))
c         tt=tt*tt
c      elseif((myear.eq.1954.and.month.gt.1).or.(myear.gt.1954))then
c         tt=0.0
c      endif
c      fr=0.75*sqrt(tt)   !fraction for Holloway length scale
c      fr=1.0              !fraction for Holloway length scale
c      fr=3.5/4.5
 	fr=0             !fraction for Holloway length scale
      cnr=0.
                     do n=1,nh
                     do m=1,mh
                     fui(m,n)=0
                     fvi(m,n)=0
c           Ft_st(m,n)=sl(m,n)*(fr*4.e+5)**2*h(m,n)
           Ft_st(m,n)=-sl(m,n)*(fr*4.e+5)**2*hreal(m,n)
           Uhol(m,n)=0
           Vhol(m,n)=0
	dcx(m,n)=0
	dcy(m,n)=0
                     end do
                     end do
      DO 15 K=1,kxlu
       N=IXDu(1,K)
        M1=IXDu(2,K)
         M2=IXDu(3,K)
      DO 15 M=M1,M2
          K2=KP(M,N)
      IF(H1(M,N).EQ.0)GOTO 15
      hr1=0
c	4./(hreal(M,N)+hreal(m+1,n+1)+hreal(m,n+1)+hreal(m+1,n))
      Uhol(m,n)=-hr1*
     *           ((Ft_st(M,N+1)-Ft_st(M,N))+
     *            (Ft_st(M+1,N+1)-Ft_st(M+1,N)))
     *             *IPV(M,N)*IPH(M,N)/(dly(m,n)+dly(m+1,n))
  
      Vhol(m,n)=hr1*
     *          ((Ft_st(M+1,N)-Ft_st(M,N))+
     *            (Ft_st(M+1,N+1)-Ft_st(M,N+1)))
     *             *IPH(M,N)*IPV(M,N)/(dlx(m,n)+dlx(m,n+1))

      ertyuyu=1
   15 CONTINUE
      if(myear.le.nmyear)then
         tt=(myear-nmyear+1)*360+(month-1)*30+min(nday,30)-1
         tt=0.5*(1.0+cos(tt*3.14159265/720.0))
      else
         tt=0.0
      endif
       do m=1,mh
          if (arc_special) then
c             DO n=1,30
c                ry(m,n)=1.0e+8
c                rx(m,n)=1.0e+8
c             end do
c             DO n=31,100
c                ry(m,n)=0.7e+8
c                rx(m,n)=0.7e+8
c             end do
c             DO n=101,nh
c                ry(m,n)=1.e+9
c                rx(m,n)=1.e+9
c             end do




        tt=0

                Do n=1,nh
                   if(n.le.40)then
                      rx(m,n)=5.0e+8+rtropic*(1.-tt)
                      ry(m,n)=5.0e+8+rtropic*(1.-tt)
                   elseif(n.le.50)then
                      rinter=0.1*(rtropic*(50-n)+rmiddle*(n-40))
                      rx(m,n)=1.0e+8*tt+rinter*(1.-tt)
                      ry(m,n)=1.0e+8*tt+rinter*(1.-tt)
                   elseif(n.le.80)then
                      rx(m,n)=1.0e+8*tt+rmiddle*(1.-tt)
                      ry(m,n)=1.0e+8*tt+rmiddle*(1.-tt)
                   elseif(n.le.90)then
                      rinter=0.1*(rmiddle*(100-n)+rpole*(n-90))
                      rx(m,n)=1.0e+7*tt+rinter*(1.-tt)
                      ry(m,n)=1.0e+7*tt+rinter*(1.-tt)
                   else
                      rx(m,n)=5.0e+6*tt+rpole*(1.-tt)
                      ry(m,n)=5.0e+6*tt+rpole*(1.-tt)
                   endif
                enddo
c          else
c             DO n=1,nh
c               ry(m,n)=ayv(n)
c               rx(m,n)=axv(m)
c                ry(m,n)=(9*tt+1)*1.e+8
c                rx(m,n)=(9*tt+1)*1.e+8

c             end do
          end if
       end do 

c      kgr vertical levels
       do 999 k=1,kgr
             do m=1,mh
				   DO n=90,nh
					  ry(m,n)=rzpole(k)     
c			  cfsmag(m,n,k)
					  rx(m,n)=rzpole(k)
c			  cfsmag(m,n,k)
				   end do

    

c				   DO n=90,nh
c					  ry(m,n)=rzpole(k)     
c			  cfsmag(m,n,k)
c					  rx(m,n)=rzpole(k)
c			  cfsmag(m,n,k)
c				   end do

				   DO n=1,30
					  ry(m,n)=5.0e+9
					  rx(m,n)=5.0e+9
				   end do
				   DO n=31,100
					  ry(m,n)=5.0e+8
					  rx(m,n)=5.0e+8
				   end do
				   DO n=101,nh
					  ry(m,n)=5.0e+7
					  rx(m,n)=5.0e+7
				   end do



             end do




c      DO 40 K=1,kxLu
c      n=IxDu(1,K)
c      I1=IxDu(2,K)
c      I2=IxDu(3,K)
      do  n=2,nh-1
        do m=2,mh-1
c      DO m=I1,I2
      dy1=(dly(m,n)+dly(m+1,n)+dly(m+1,n+1)+dly(m,n+1))/4.
      dy0=(dly(m,n)+dly(m+1,n)+dly(m+1,n-1)+dly(m,n-1))/4.
       
      dx1=(dlx(m,n)+dlx(m+1,n)+dlx(m+1,n+1)+dlx(m,n+1))/4.
      dx0=(dlx(m,n)+dlx(m-1,n)+dlx(m-1,n+1)+dlx(m,n+1))/4.
       if(dx1.eq.0)dx1=dx(10,10)
       if(dy1.eq.0)dy1=dy(10,10)
      gam(m,n)=(dlx(m,n)+dlx(m,n+1))*(dly(m,n)+dly(m+1,n))
     *        *3./2.
      dcx(m,n)=-dy1/dx1*(rx(m,n)+rx(m+1,n)+rx(m,n+1))
     *         -dy0/dx1*(rx(m,n)+rx(m+1,n)+rx(m,n-1))
      dcy(m,n)=-dx1/dy1*(ry(m,n)+ry(m,n+1)+ry(m+1,n))
     *         -dx0/dy1*(ry(m,n)+ry(m,n+1)+ry(m-1,n))
      end do
	dcx(2,n)=dcx(1,n)
      end do	

c       do 999 k=1,kgr
   40 CONTINUE
      P2=0.
      P3=0.
          do n=1,nh-1
           do m=1,mh-1
           k2=kp(m,n)
              u(m,n)=uu(k2,k)
              v(m,n)=vv(k2,k)
           if(k.gt.kbtv(m,n)-1) then
          u(m,n)=0.
          v(m,n)=0.
          Uhol(m,n)=0.
          Vhol(m,n)=0.
            end if
      DX1=(dlx(m,n)+dlx(m,n+1))/2.
      DY1=(dly(m,n)+dly(m+1,n))/2.
          curv_uv(m,n)=(v(m,n)*(dly(m+1,n)-dly(m,n))/dx1/dy1-
     *               u(m,n)*(dlx(m,n+1)-dlx(m,n))/dx1/dy1)
           u0(m,n)=u(m,n)
           v0(m,n)=v(m,n)
           end do
          end do
c      U- x Advection a diffusion
c      goto 777
      DO1 Kj=1,kxLu
      N=IXDu(1,Kj)
      I1=IXDu(2,Kj)
      I2=IXDu(3,Kj)
       dy0=dy(i1-1,n-1)
       dy1=dy(i1-1,n)

      m=i1-1 
      advX(m)=3*
     *    	(u0(m,n)*(dlx(m+1,n)+dlx(m+1,n+1))
     *    	+u0(m+1,n)*(dlx(m,n)+dlx(m,n+1)))
     *         /(dlx(m+1,n)+dlx(m+1,n+1)+dlx(m,n)+dlx(m,n+1))
     *          *dly(m+1,n)
      DO2 M=I1,I2
c******************************
c*   8 *----- 3* ------*7
c*   4 *----- 2* ------*5
c*   6 *----- 1* ------*8
c******************************
       dy0=dy(m,n-1)
       dy1=dy(m,n)

      adx4=-advX(m-1)
      advX(m)=3*
     *    	(u0(m,n)*(dlx(m+1,n)+dlx(m+1,n+1))
     *    	+u0(m+1,n)*(dlx(m,n)+dlx(m,n+1)))
     *         /(dlx(m+1,n)+dlx(m+1,n+1)+dlx(m,n)+dlx(m,n+1))
     *          *dly(m+1,n)
      adx5=advX(m)

	  a1(m)=adx4-cnr*DABS(adx4)+dcx(m-1,n)
        c1(m)=adx5-cnr*DABS(adx5)+dcx(m,n)
       a1(m)=0*adx4-0*cnr*DABS(adx4)+dcx(m-1,n)
        c1(m)=0*adx5-0*cnr*DABS(adx5)+dcx(m,n)

         b1(m)=gam(m,n)/dt-(A1(M)+C1(M))
          fu1(m)=gam(m,n)/dt*u(m,n)
     *		+(dcx(m-1,n)*Uhol(m-1,n)+dcx(m,n)*Uhol(m+1,n)
     *          -(dcx(m-1,n)+dcx(m,n))*Uhol(m,n))
     *           +curv_uv(m,n)*gam(m,n)*v0(m,n)
      if(ipv(m,n)*iph(m,n).eq.0.or.k.gt.kbtv(m,n)-1) then
      a1(m)=0
      c1(m)=0
      b1(m)=gam(m,n)/dt
      fu1(m)=0
      end if
    2 CONTINUE
      P2=0.
      P3=0.
      ps2=0.
         Call Fact1(A1,B1,C1,Fu1,U1,p2,p3,i1,i2)
                   do m=I1,I2       
                      u(m,n)=u1(m)
                   end do
      do m=i1,i2
         dz1=dzk(k+1)
c               if(k.eq.kbtv(m,n)-1) dz1=dz(k)/2.
         if(k.le.kbtv(m,n)-1) then
           aa1=a1(m)
      	  cc1=c1(m)
             bb1=-aa1-cc1
         fui(m,n)=fui(m,n)+(-(aa1*u(m-1,n)
     *	   +bb1*u(m,n)+cc1*u(m+1,n))
     *     	+(dcx(m-1,n)*Uhol(m-1,n)+dcx(m,n)*Uhol(m+1,n)
     *          -(dcx(m-1,n)+dcx(m,n))*Uhol(m,n))
     *           +curv_uv(m,n)*gam(m,n)*v0(m,n)

     *            )*dz1
     *               /h1(m,n)/gam(m,n)

          end if
       end do

c*    V-x diffusion-advection
      DO22 M=I1,I2
c******************************
c*   8 *----- 3* ------*7
c*   4 *----- 2* ------*5
c*   6 *----- 1* ------*8
c******************************

      adx4=-advX(m-1)
      adx5=advX(m)

      a1(m)=(adx4-cnr*DABS(adx4)+dcx(m-1,n))
        c1(m)=(adx5-cnr*DABS(adx5)+dcx(m,n))
      a1(m)=(0*adx4-0*cnr*DABS(adx4)+dcx(m-1,n))
        c1(m)=(0*adx5-0*cnr*DABS(adx5)+dcx(m,n))


         b1(m)=gam(m,n)/dt-(A1(M)+C1(M))
           fv1(m)=gam(m,n)/dt*v(m,n)
     *		+(dcx(m-1,n)*Vhol(m-1,n)+dcx(m,n)*Vhol(m+1,n)
     *          -(dcx(m-1,n)+dcx(m,n))*Vhol(m,n))
     *           -curv_uv(m,n)*gam(m,n)*u0(m,n)
      if(ipv(m,n)*iph(m,n).eq.0.or.k.gt.kbtv(m,n)-1) then
      a1(m)=0
      c1(m)=0
      b1(m)=gam(m,n)/dt
      fv1(m)=0
      end if
   22 CONTINUE
      P2=0.
      P3=0.
      ps2=0.
         Call Fact1(A1,B1,C1,Fv1,V1,p2,p3,i1,i2)
                    do m=I1,I2       
                       v(m,n)=v1(m)
                    end do




c        if(n.ge.nc1.and.n.le.nc2) then
c             u(mc1-1,n)=u(mc2,n)
c              u(mc2+1,n)=u(mc1,n)
c               v(mc1-1,n)=v(mc2,n)
c                v(mc2+1,n)=v(mc1,n)
c         end if
      do m=i1,i2
         dz1=dzk(k+1)
c               if(k.eq.kbtv(m,n)-1) dz1=dz(k)/2.
         if(k.le.kbtv(m,n)-1) then
            aa1=a1(m)
	      cc1=c1(m)
            bb1=-aa1-cc1
            fvi(m,n)=fvi(m,n)+(-(aa1*v(m-1,n)+bb1*v(m,n)+cc1*v(m+1,n))
     *     	+(dcx(m-1,n)*Vhol(m-1,n)+dcx(m,n)*Vhol(m+1,n)
     *          -(dcx(m-1,n)+dcx(m,n))*Vhol(m,n))
     *           -curv_uv(m,n)*gam(m,n)*u0(m,n)
     *            )*dz1
     *               /h1(m,n)/gam(m,n)
		  
          end if
      end do

    1 CONTINUE
c     U - y advection-diffusion
      DO 5 Kj=1,kyLu
      M=IYDu(1,Kj)
      I1=IYDu(2,Kj)
      I2=IYDu(3,Kj)
      dx1=dx(m,i1)
      dx0=dx(m-1,i1)
      ady3=0.
      n=i1-1
        advY(n)=3*
     *	  (v0(m,n)*(dly(m,n+1)+dly(m+1,n+1))
     *      +v0(m,n+1)*(dly(m,n)+dly(m+1,n)))
     *     /(dly(m,n+1)+dly(m+1,n+1)+dly(m,n)+dly(m+1,n))
     *     *dlx(m,n+1)
      DO7 N=I1,I2
      z2=ipv(m,n)
      z3=iph(m-1,n)
      dx1=dx(m,n)
      dx0=dx(m-1,n)
        ady1=-advY(n-1)
        advY(n)=3*
     *	  (v0(m,n)*(dly(m,n+1)+dly(m+1,n+1))
     *      +v0(m,n+1)*(dly(m,n)+dly(m+1,n)))
     *     /(dly(m,n+1)+dly(m+1,n+1)+dly(m,n)+dly(m+1,n))
     *     *dlx(m,n+1)
        ady3=advY(n)

        A1(n)=(ady1-cnr*DABS(ady1)+dcy(m,n-1))
          C1(n)=(ady3-cnr*DABS(ady3)+dcy(m,n))
	 A1(n)=(0*ady1-0*cnr*DABS(ady1)+dcy(m,n-1))
         C1(n)=(0*ady3-0*cnr*DABS(ady3)+dcy(m,n))


            B1(n)=gam(m,n)/dt-(A1(n)+C1(n))
             Fu1(n)=gam(m,n)/dt*u(m,n)
     *		+(dcy(m,n-1)*Uhol(m,n-1)+dcy(m,n)*Uhol(m,n+1)
     *          -(dcy(m,n-1)+dcy(m,n))*Uhol(m,n))
      if(ipv(m,n)*iph(m,n).eq.0.or.k.gt.kbtv(m,n)-1) then
      a1(n)=0
      c1(n)=0
      b1(n)=gam(m,n)/dt
       fu1(n)=0
      end if
    7 CONTINUE
         Call Fact1(A1,B1,C1,Fu1,U1,p2,p3,i1,i2)
         do n=I1,I2
            u(m,n)=u1(n)
         end do
          do n=i1,i2
              dz1=dzk(k+1)
c               if(k.eq.kbtv(m,n)-1) dz1=dz(k)/2.
               if(k.le.kbtv(m,n)-1) then
              	aa1=a1(n)
	            cc1=c1(n)
                  bb1=-aa1-cc1
        fui(m,n)=fui(m,n)+(-(aa1*u(m,n-1)+bb1*u(m,n)+cc1*u(m,n+1))
     *		+(dcy(m,n-1)*Uhol(m,n-1)+dcy(m,n)*Uhol(m,n+1)
     *          -(dcy(m,n-1)+dcy(m,n))*Uhol(m,n))

     *	  )*dz1
     *           /h1(m,n)/gam(m,n)
	end if

        end do

c****** V-y diffusion-advection
      DO77 N=I1,I2
        ady1=-advY(n-1)
        ady3=advY(n)

         A1(n)=ady1-cnr*DABS(ady1)+dcy(m,n-1)
          C1(n)=ady3-cnr*DABS(ady3)+dcy(m,n)
	A1(n)=0*ady1-0*cnr*DABS(ady1)+dcy(m,n-1)
          C1(n)=0*ady3-0*cnr*DABS(ady3)+dcy(m,n)

            B1(n)=gam(m,n)/dt-(A1(n)+C1(n))
              Fv1(n)=gam(m,n)/dt*v(m,n)
     *		+(dcy(m,n-1)*Vhol(m,n-1)+dcy(m,n)*Vhol(m,n+1)
     *          -(dcy(m,n-1)+dcy(m,n))*Vhol(m,n))
      if(ipv(m,n)*iph(m,n).eq.0.or.k.gt.kbtv(m,n)-1) then
      a1(n)=0
      c1(n)=0
      b1(n)=gam(m,n)/dt
      fv1(n)=0
      end if
   77 CONTINUE
         Call Fact1(A1,B1,C1,Fv1,V1,p2,p3,i1,i2)
              do n=I1,I2
                 v(m,n)=v1(n)
              end do
          do n=i1,i2
              dz1=dzk(k+1)
c               if(k.eq.kbtv(m,n)-1) dz1=dz(k)/2.
               if(k.le.kbtv(m,n)-1) then 
               	aa1=a1(n)
	            cc1=c1(n)
                  bb1=-aa1-cc1
	  fvi(m,n)=fvi(m,n)+(-(aa1*v(m,n-1)+bb1*v(m,n)+cc1*v(m,n+1))
     *		+(dcy(m,n-1)*Vhol(m,n-1)+dcy(m,n)*Vhol(m,n+1)
     *          -(dcy(m,n-1)+dcy(m,n))*Vhol(m,n))
     *	  )*dz1
     *              /h1(m,n)/gam(m,n)
               end if
        end do


    5 CONTINUE
 777  continue
c            do n=nc1,nc2
c             u(mc1-1,n)=u(mc2,n)
c              u(mc2+1,n)=u(mc1,n)
c               v(mc1-1,n)=v(mc2,n)
c                v(mc2+1,n)=v(mc1,n)
c            end do
          do n=1,nh
           do m=1,mh
              uu(kp(m,n),k)=u(m,n)
              vv(kp(m,n),k)=v(m,n)
              u(m,n)=0
              v(m,n)=0
           end do
           end do
        
  999 CONTINUE
      do k=1,kgr
      uu(1,k)=0.
      vv(1,k)=0.
      end do
  666  format(1x,7f5.1,2x,7f5.1)
   80 FORMAT(2(7(E14.7)/)/)
c       read *
      RETURN
      END
