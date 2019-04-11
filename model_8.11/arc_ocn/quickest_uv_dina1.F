

	SUBROUTINE qst_UV(uu,vv,uu_ts,vv_ts,wwvel,ww,lvl,lvl0 )
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      COMMON/RUNPARS/nmyear
      COMMON/NH2/ NH2
      real *8 z,dz,dzk,dx,dy,h,h1,sl,slb,dx4,dy4
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
      integer *2 ikbt,inbt, kz
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      common/Cfsmag/cfsmag(mh,nh,kgr)
      common/Cf/UMY(nh),UMX(mh),RU(kh),ayv(nh),axv(mh)        
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/v6/Fui(mh,nh),Fvi(mh,nh)
	common/brt/UBrT(nodes), VBrT(nodes)
	common/brkl/Uu_Brkl(nodes,kgr), Vv_brkl(nodes,kgr)
	common/int/Uu_int(nodes), Vv_int(nodes)
	common/cwz/   cwz00(mh,nh),	cwzvel00(mh,nh)


      real *8 uu(nodes,kgr),vv(nodes,kgr),  WWvel(mh,nh,kgr)
      Real *8 A1(1000),B1(1000),C1(1000),Fu1(1000),U1(1000),
     *        Fv1(1000),v1(1000) 
     *    , uxflux(mh,kh), vxflux(mh,kh)
     *     , uyflux(mh,nh,kh), vyflux(mh,nh,kh)
     *     , uzflux(mh,kh), vzflux(mh,kh), hgm(mh,kh)
     *    ,Temp(nodes,kgr),Temp1(nodes,kgr),ww(nodes,kgr)
     *      ,uu_ts(mh,nh,kgr) ,vv_ts(mh,nh,kgr)
        
      REAL *8 dcx(mh,nh),dcy(mh,nh),gam(mh,nh),curv_uv(mh,nh)
     *        ,Uhol(mh,nh),Vhol(mh,nh),FT_st(mh,nh) 
     *        ,u(nh,kh),v(nh,kh), w(mh,kh), u0(mh,kh),v0(mh,kh)
     *  ,grad,gradV, curvV, curv curz, curm, curvz 
     *       , sfl, sfll, sflV,sfllV,uz1,vz1
     *    ,u01(mh,kh),v01(mh,kh),u2(nh,kh),v2(nh,kh),uz,vz
     *   ,wer1,cwz00
     *        ,lvl(mh,nh),lvl0(mh,nh),lvl_uv
            
      integer nu_diag,  kgr1

      common/iodiag/nu_diag

C*-------------------------------------------
c     	dt=3*3600

	ll=0
c     UV - y advection
      DO M=M0,MM


	   DO N=1,Nh

            K2=KP(M,N)
		  DO K=kz0,KKZ			 
		     u2(N,K)=uu(K2,K-1)
			 v2(N,K)=vv(K2,K-1)
			 uyflux(m,n,k)=0.
			 vyflux(m,n,k)=0.
			 u(N,K)=uu(K2,K-1)
			 v(N,K)=vv(K2,K-1)
			 if(k.gt.kbtv(m,n)) then
				u(n,k)=0.
				v(n,k)=0.  
				u2(n,k)=0.
				v2(n,k)=0.        
               end if
            end do
	   end do
    
         do ll=1,kylu

	      if (iydu(1,ll).eq.m) then


               I1=IYDu(2,LL)
               I2=IYDu(3,LL)
               n=i1
			 dy4=(dly(m,n)+dly(m+1,n)+dly(m,n-1)+dly(m+1,n-1))/4
			 kzv=max(kbtv(m,n),kbtv(m,n-1))
	         do k=kz0,kzv

c				 curm=(vv(kp(m,n-1),k-1)*(dly(m,n)+dly(m+1,n))
c     *                   +vv(kp(m,n),k-1)*(dly(m,n-1)+dly(m+1,n-1)))
c     *              /(dly(m,n)+dly(m+1,n)+dly(m,n-1)+dly(m+1,n-1))

	curm=vv_ts(m,n,k-1)

                    curm=curm/dy4*dt !!!!courant	


c		    	curm=(v2(n,k)+v2(n-1,k))/dy4/2.D+0 *dt  !!!!courant

			    sfl=0.5*(u(n,k)+u(n-1,k))

	if (m.eq.56.and.n.eq.20.and.k.eq.2) then
	continue
	end if


			    uyflux(m,n-1,k)=sfl*curm* dy4*dlx(m,n) *dzk(k)
                  sflV=0.5*(v(n,k)+v(n-1,k))
			    vyflux(m,n-1,k)=sflV*curm
     *                * dy4*dlx(m,n) *dzk(k)
	         end do

          	n=i2
		if (m.eq.66.and.n.eq.125) then
	continue
	end if


          	dy4=(dly(m,n+1)+dly(m+1,n+1)+dly(m,n)+dly(m+1,n))/4
          	kzv=max(kbtv(m,n),kbtv(m,n+1))

	     	do k=kz0,kzv
		if (m.eq.54.and.n.eq.124.and.k.eq.5) then
	continue
	end if
		       
			    curm=vv_ts(m,n+1,k-1)
                    curm=curm/dy4*dt !!!!courant	
				  
c				  curm=(v2(n+1,k)+v2(n,k))/dy4/2.D+0  *dt !!!!courant

                 sfl=0.5*(u(n+1,k)+u(n,k))
			   uyflux(m,n,k)=sfl*curm     * dy4*dlx(m,n+1) *dzk(k) 
			   sflV=0.5*(v(n+1,k)+v(n,k))
			   vyflux(m,n,k)=sflV*curm
     *                * dy4*dlx(m,n+1) *dzk(k)   
  	        end do

              DO N=I1+1,I2
          	   kzv=max(kbtv(m,n),kbtv(m,n-1))
                 DO K=kz0,kzv    !!!!!kzv+1  ?
				  dy4=(dly(m,n)+dly(m+1,n)+dly(m,n-1)+dly(m+1,n-1))/4
				  dy3=(dly(m,n-1)+dly(m+1,n-1)
     *			  +dly(m,n-2)+dly(m+1,n-2))/4
				  dy5=(dly(m,n+1)+dly(m+1,n+1)+dly(m,n)+dly(m+1,n))/4
                    grad=(u(n,k)-u(n-1,k))/DY4 !!!gradient
		          gradV=(v(n,k)-v(n-1,k))/DY4 !!!gradient

                    curm=vv_ts(m,n,k-1)

                    curm=curm/dy4*dt !!!!courant	  
				  curv=0.D+0 
	              curvV=0.D+0           !!!!curvature
		          if (curm.gt.0) then
	                 curv=2.D+0*
     *		         ((u(n,k)-u(n-1,k))/dy4
     *                     -(u(n-1,k)-u(n-2,k))/dy3)
     *			         /(dy3+dy4)
					 curvV=2.D+0*
     *		           ((v(n,k)-v(n-1,k))/dy4
     *                     -(v(n-1,k)-v(n-2,k))/dy3)
     *			            /(dy3+dy4)
			      end if

			     if(curm.lt.0) then                    
			          curv=2.D+0*((u(n+1,k)-u(n,k))/ dy5
     *		           -(u(n,k)-u(n-1,k))/ dy4)/(dy4+dy5)
			          curvV=2.D+0*((v(n+1,k)-v(n,k))/ dy5
     *		            -(v(n,k)-v(n-1,k))/ dy4)/(dy4+dy5)
			     end if


	              sfll=0.5D+0*(u(n-1,k)+u(n,k))
	              sfl=sfll 
     *                -curm*grad* dy4/2.D+0
     *               -(1-curm**2)*curv*(dy4**2)/6.D+0

	              uyflux(m,n-1,k)=curm*sfl* dy4*dlx(m,n)*dzk(k) 

 
	  

          	      sfllV=0.5D+0*(v(n-1,k)+v(n,k))
	              sflV=sfllV 
     *                -curm*gradV* dy4/2.D+0
     *                -(1-curm**2)*curvV*(dy4**2)/6.D+0

	              vyflux(m,n-1,k)=curm*sflV* dy4*dlx(m,n) *dzk(k) 

		if (m.eq.56.and.n.eq.20.and. k.eq.2) then
	continue
	end if
		     


                 end do
    
              end do

		
	      end if
	   end do  



     
      end do
		

	ll=0
c     UV - xz advection
      DO N=N0,NN
	   DO m=1,mh-1
            K2=KP(M,N)
			k3=KP(M+1,N)
			k4=KP(M,N+1)
			k5=KP(M+1,N+1)
	      DO K=2,KKZ 
		     U01(M,K)=uu(K2,K-1)
			 V01(M,K)=vv(K2,K-1)
	         W(M,K)=WWvel(m,n,k-1)

			 U0(M,K)=uu(K2,K-1)
			 V0(M,K)=vv(K2,K-1)
	         
			 uxflux(m,k)=0.
			 vxflux(m,k)=0.
			 uzflux(m,k)=0.
			 vzflux(m,k)=0.
			 if(k.gt.kbtv(m,n)) then
			 	u0(m,k)=0.
				v0(m,k)=0.
					U01(M,K)=0.
					 V01(M,K)=0.

			 end if
            end do
	   end do
         do ll=1,kxlu
	continue
		  if (ixdu(1,ll).eq.n) then	
              I1=IxDu(2,LL)
              I2=IxDu(3,LL)  	 
          	do m=i1,i2
	           k=kz0	
	           hgm(m,k)=dlx(m,n)*dly(m,n)*dzk(k)
	           kzv=kbtv(m,n)
c			kzv=max(kbt(m,n),kbt(m,n+1),kbt(m+1,n),kbt(m+1,n+1))
	           do k=kz0,kzv
	              hgm(m,k)=dlx(m,n)*dly(m,n)*dzk(k)
			   end do 		
	        end do 
c     UV - x advection
              m=i1
              dx4=(dlx(m-1,n)+dlx(m,n)+dlx(m-1,n+1)+dlx(m,n+1))/4
	        kzv=max(kbtv(m,n),kbtv(m-1,n))
	        do k=kz0,kzv

				curm=uu_ts(m,n,k-1)

		       curm= curm/dx4*dt


 		       sfl=0.5*(u0(m,k)+u0(m-1,k))
      	       sflV=0.5*(v0(m,k)+v0(m-1,k))

	if (m.eq.53.and.n.eq.124.and. k.eq.9) then
	continue
	end if

                 uxflux(m-1,k)=sfl*curm
     *               *(dly(m,n))
     *			  *dx4*dzk(k)

                 vxflux(m-1,k)=sflV*curm
     *               *(dly(m,n))
     *			  *dx4*dzk(k)
              end do
	        m=i2
	        dx4=(dlx(m+1,n)+dlx(m,n)+dlx(m+1,n+1)+dlx(m,n+1))/4
 	        kzv=max(kbtv(m,n),kbtv(m+1,n))
	        do k=kz0,kzv

					curm=uu_ts(m+1,n,k-1)
  

		       curm= curm/dx4*dt


			   sfl=0.5*(u0(m+1,k)+u0(m,k))
		       sflV=0.5*(v0(m,k)+v0(m+1,k))

                 uxflux(m,k)=sfl*curm
     *               *(dly(m+1,n))*dzk(k)
     *			  *dx4
	           vxflux(m,k)=sflV*curm
     *              *(dly(m+1,n))*dzk(k)
     *			  *dx4
              end do

              DO M=I1+1,I2
                 dx4=(dlx(m-1,n)+dlx(m,n)+dlx(m-1,n+1)+dlx(m,n+1))/4
	           dx3=(dlx(m-2,n)+dlx(m-1,n)+dlx(m-2,n+1)+dlx(m-1,n+1))/4
	           dx5=(dlx(m,n)+dlx(m+1,n)+dlx(m,n+1)+dlx(m+1,n+1))/4
			   kzv=max(kbtv(m,n),kbtv(m-1,n))

                 DO K=kz0,kzv

	              grad=(u0(m,k)-u0(m-1,k))/dx4 !!!gradient
		          gradV=(v0(m,k)-v0(m-1,k))/dx4 !!!gradient

	              curm=uu_ts(m,n,k-1)

   

		          curm= curm/dx4*dt

                    curv=0.D+0 !curvature
	              curvV=0.D+0 !curvature	

c-----------fix curvature	--------------------
				 if  (curm.gt.0) then
					 curv=2.D+0*((u0(m,k)-u0(m-1,k))/dx4
     *					  -(u0(m-1,k)-u0(m-2,k))/dx3)
     *					   /(dx4+dx3)
					 curvV=2.D+0*((v0(m,k)-v0(m-1,k))/dx4
     *					  -(v0(m-1,k)-v0(m-2,k))/dx3)
     *					   /(dx4+dx3)
				 end if


				 if (curm.lt.0) then
					 curv=2.D+0*((u0(m+1,k)-u0(m,k))/dx5
     *			                  -(u0(m,k)-u0(m-1,k))/dx4)
     *				         /(dx5+dx4)
					 curvV=2.D+0*((v0(m+1,k)-v0(m,k))/dx5
     *			                  -(v0(m,k)-v0(m-1,k))/dx4)
     *				         /(dx5+dx4)

				 end if

c-------------------------fix fluxes---------------
	              sfll=0.5D+0*(u0(m,k)+u0(m-1,k))
                    sfl=0.5D+0*(u0(m,k)+u0(m-1,k))
     *        	   -curm*grad*dx4/2.D+0
     *        	   -(1-curm**2)*curv*(dx4**2)/6.D+0


                    uxflux(m-1,k)=sfl*curm
     *			    *dly(m,n)*dx4*dzk(k)

	              sfllV=0.5D+0*(v0(m,k)+v0(m-1,k))
                    sflV=0.5D+0*(v0(m,k)+v0(m-1,k))
     *	             -curm*gradV*dx4/2.D+0
     *	              -(1-curm**2)*curvV*(dx4**2)/6.D+0

                    vxflux(m-1,k)=sflV*curm
     *			    *dly(m,n)*dx4*dzk(k)

			
	if (m.eq.126.and.n.eq.132.and. k.eq.11) then
	continue
	end if
		        
		        
                 end do

              end do

                
c     UV - z advection
              DO M=I1,I2

c	           kzv=kbtv(m,n)
		kzv=max(kbt(m,n),kbt(m+1,n),kbt(m,n+1),kbt(m+1,n+1))
                 DO K=kz0+2,kzv-1
                    grad=(u0(m,k)-u0(m,k-1))/dz(k-1)
	              gradV=(v0(m,k)-v0(m,k-1))/dz(k-1)		
				  curm= (w(m,k-1))/dz(k-1)*dt
             			
               	  curv=0.D+0 !curvature	

				  curvV=0.D+0 !curvature
c-----------fix curvature	--------------------
				  if(curm.gt.0) then
		             curv=2.D+0*((u0(m,k)-u0(m,k-1))/dz(k-1)
     *					  -(u0(m,k-1)-u0(m,k-2))/dz(k-2))
     *					   /(dz(k-1)+dz(k-2))
	                 curvV=2.D+0*((v0(m,k)-v0(m,k-1))/dz(k-1)
     *					  -(v0(m,k-1)-v0(m,k-2))/dz(k-2))
     *					   /(dz(k-1)+dz(k-2))
	
				  end if


				  if(curm.lt.0) then	
		
					  curv=2.D+0*((u0(m,k+1)-u0(m,k))/dz(k)
     *			      -(u0(m,k)-u0(m,k-1))/dz(k-1))
     *				  /(dz(k-1)+dz(k))

					  curvV=2.D+0*((v0(m,k+1)-v0(m,k))/dz(k)
     *			      -(v0(m,k)-v0(m,k-1))/dz(k-1))
     *				  /(dz(k-1)+dz(k))
				   end if

c-------------------------fix fluxes---------------
				   sfll=(u0(m,k)*dzk(k)+u0(m,k-1) *dzk(k-1))
     *				  /(dzk(k)+dzk(k-1))	           
				   sfl=sfll
     *			    -curm*grad*dz(k-1)/2.D+0
     *				 -(1-curm**2)*curv*(dz(k-1)**2)/6.D+0

				   sfll=(v0(m,k)*dzk(k)+v0(m,k-1)*dzk(k-1))
     *				  /(dzk(k)+dzk(k-1))
		     	   sflV=sfll
     *			    -curm*gradV*dz(k-1)/2.D+0
     *				 -(1-curm**2)*curvV*(dz(k-1)**2)/6.D+0

		    	   uzflux(m,k-1)=sfl*curm
     *			    *(dly(m,n)+dly(m+1,n))/2.D+0
     *			     *dz(k-1)*(dlx(m,n)+dlx(m,n+1))/2.
				   vzflux(m,k-1)=sflV*curm
     *			      *(dly(m,n)+dly(m+1,n))/2.D+0
     *			      *dz(k-1)*(dlx(m,n)+dlx(m,n+1))/2.

		

		        
                 end do

  	           uzflux(m,kz0-1)=0
	           uzflux(m,kzv)=0
	           vzflux(m,kz0-1)=0
	           vzflux(m,kzv)=0

	           k=kz0+1

	    	   curm= (w(m,k-1))/dz(k-1)*dt
			   sfl=(u0(m,k)*dzk(k)+u0(m,k-1)*dzk(k-1))
     *				  /(dzk(k)+dzk(k-1))
     	           uzflux(m,k-1)=sfl*curm
     *			    *(dly(m,n)+dly(m+1,n))/2.D+0
     *			  *dz(k-1) *(dlx(m,n)+dlx(m,n+1))/2.
			   sflV=(v0(m,k)*dzk(k)+v0(m,k-1)*dzk(k-1))
     *				  /(dzk(k)+dzk(k-1))
	           vzflux(m,k-1)=sflV*curm
     *			    *(dly(m,n)+dly(m+1,n))/2.D+0
     *			  *dz(k-1) *(dlx(m,n)+dlx(m,n+1))/2.

	           k=kzv

		       curm= (w(m,k-1))/dz(k-1)*dt
		       sfl=(u0(m,k)*dzk(k)+u0(m,k-1)*dzk(k-1))
     *				  /(dzk(k)+dzk(k-1))
	           uzflux(m,k-1)=sfl*curm
     *			    *(dly(m,n)+dly(m+1,n))/2.D+0
     *			  *dz(k-1)*(dlx(m,n)+dlx(m,n+1))/2. 
		       sflV=(v0(m,k)*dzk(k)+v0(m,k-1)*dzk(k-1))
     *				  /(dzk(k)+dzk(k-1))
	           vzflux(m,k-1)=sflv*curm
     *			    *(dly(m,n)+dly(m+1,n))/2.D+0
     *			  *dz(k-1) *(dlx(m,n)+dlx(m,n+1))/2.
              end do

			DO  m=i1,i2
				kzv=kbtv(m,n) 
			    DO  K=kz0,kzv	
	
	 lvl_uv=(lvl0(m,n)+lvl0(m+1,n)+lvl0(m,n+1)+lvl0(m+1,n+1))/4.

	if (k.eq.2) then
       wer=
     *    -(uxflux(m,k)-uxflux(m-1,k))
     *    -(uyflux(m,n,k)-uyflux(m,n-1,k))
     *    -uzflux(m,k) *dzk(k)/(dzk(k)-lvl_uv)
	else
       wer=
     *    -(uxflux(m,k)-uxflux(m-1,k))
     *    -(uzflux(m,k)-uzflux(m,k-1))
     *    -(uyflux(m,n,k)-uyflux(m,n-1,k))
      end if
			if (k.eq.2) then
         
				wer1=
     *		dt*u0(m,k)*cwz00(m,n)*dzk(k)/(dzk(k)-lvl_uv)/4.D+0
			else
				wer1=0
			end if
    		
		            u0(m,k)=u0(m,k)*hgm(m,k)+wer+wer1
				    u0(m,k)=u0(m,k)/hgm(m,k)
		if (isNan(u0(m,k))) then
	print *,'qst',m,n,k,u0(m,k),cwz00(m,n)
      read *
	end if

	if (k.eq.2) then
       wer=
     *    -(vxflux(m,k)-vxflux(m-1,k))
     *    -(vyflux(m,n,k)-vyflux(m,n-1,k))
     *    -vzflux(m,k)  *dzk(k)/(dzk(k)-lvl_uv)
	else
       wer=
     *    -(vxflux(m,k)-vxflux(m-1,k))
     *    -(vzflux(m,k)-vzflux(m,k-1))
     *    -(vyflux(m,n,k)-vyflux(m,n-1,k))
      end if

	if (m.eq.107.and.n.eq.2.and.k.eq.2) then
	continue
	end if


			if (k.eq.2) then
				wer1=
     *		dt*v0(m,k)*cwz00(m,n)*dzk(k)/(dzk(k)-lvl_uv)/4.D+0
			else     
				wer1=0
			end if

	                v0(m,k)=v0(m,k)*hgm(m,k)+wer+wer1
	          	    v0(m,k)=v0(m,k)/hgm(m,k)	   
				    u0(m,k)=u0(m,k)
c     *	   	*inbt(m,n,k)*inbt(m+1,n,k)*inbt(m,n+1,k)*inbt(m+1,n+1,k)
				    v0(m,k)=v0(m,k)
c    *		*inbt(m,n,k)*inbt(m+1,n,k)*inbt(m,n+1,k)*inbt(m+1,n+1,k) 
			    end do
			end do

		 end if	

	continue
	  end do

	  DO  K=1,kgr-1
        DO  m=m0,mm
	      k2=kp(m,n) 	
	      uu(K2,K)=U0(m,K+1)
	      vv(K2,K)=V0(m,K+1)
	  end do
	  end do


      end do


      RETURN

	END
