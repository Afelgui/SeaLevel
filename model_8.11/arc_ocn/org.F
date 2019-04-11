      SUBROUTINE Org(restart,renew)
      Include 'model.par'
      logical restart, renew
      integer *2 Kp,Ipv,Iph,mask,Ixd,Iyd,Ixyd,Ixdu,Iydu,Isl,Is_pnt
      Common/Gr1/Kp(mh,nh),Ipv(mh,nh),Iph(mh,nh),mask(mh,nh)
     *          ,iXd(3,kxLine),iYd(3,kyLine),iXYd(3,kxyLin)
     *          ,iXdu(3,kxLu),iYdu(3,kyLu)
      COMMON/Gr2/ M0,MM,N0,NN,KKP,KZ0,KKZ
      real *8 Z,Dz,Dzk,dx,dy,h,h1
      COMMON/Gr3/Z(kh),Dz(kh),Dzk(kh),dx(mh,nh),dy(mh,nh)
      common/Gr4/kbt(mh,nh),kbtv(mh,nh)
	
      integer *2 ikbt,inbt
      common/Gr5/ikbt(mh,nh,kh),inbt(mh,nh,kh)
      real *8 dlx,dly
	real *8 z1(kh)

      common/d3/dlx(mh,nh),dly(mh,nh)
      COMMON/botV/h1(mh,nh)
      COMMON/botTS/h(mh,nh)
      real *8 hreal
      COMMON/botReal/hreal(mh,nh)
      COMMON/d1/ Rad,Gro,R,dt
      real *8 hl1,hl2,sl,slb
      COMMON/d2/SL(mh,nh),SLB(mh,nh)
      common/clame/hl1(mh,nh),hl2(mh,nh)
	integer*4 idate, isec
	common/date/idate,isec
      dimension ih(mh,nh)
      logical home, oldmode
      real *8 xx(mh,nh),yy(mh,nh),zz(kh),xxu(mh,nh),yyu(mh,nh)
      Common/d33/ xx,yy,zz,xxu,yyu
      integer nu_diag
      character*80 rname
      integer maxmon(12), myr, mmo, mdy
      data maxmon/31,28,31,30,31,30,31,31,30,31,30,31/
      common/iodiag/nu_diag
      DATA TEO,TEOV,omega,PI
	
     */-2.25,0.,0.729E-4,3.14159265/
       rad=6.37e+8
      R=1.e-6
      dt=3600.*2.
      oldmode=.false.
	idate=19480101



      KKP=nodes
         M0=2
           MM=mh-1
             N0=2
               NN=nh-1
                 KZ0=2
                    KKZ=kh-2

      Open(10,file='c:\dina\input data\domain\grid.dat'
     *        ,form='unformatted'
     *        ,access='sequential'
     *        ,status='unknown')
      read(10)xx,yy,zz,dlx,dly,xxu,yyu
      close(10,STATUS='KEEP')
	Write(nu_diag,*)mh,nh


      open (12,file='c:\dina\input data\domain\kbt_new4.dat'
     *        ,form='unformatted'
     *        ,access='sequential'
     *        ,status='unknown')

      read (12)kbt,kbtv
      close(12)
      
	

      Open(10,file='c:\dina\input data\domain\bound.dat'
     *        ,form='unformatted'
     *        ,access='sequential'
     *        ,status='unknown')

      read(10)((KP(i,j),i=1,mh),j=1,nh),
     *         ((IPV(i,j),i=1,mh),j=1,nh),
     *         ((IPH(i,j),i=1,mh),j=1,nh),
     *         ((mask(i,j),i=1,mh),j=1,nh)
     *        ,((IXD(i,j),i=1,3),j=1,kxLine),
     *         ((IYD(i,j),i=1,3),j=1,kyLine),
     *         ((IXdu(i,j),i=1,3),j=1,kxLu),
     *         ((IYdu(i,j),i=1,3),j=1,kyLu)
      
      close(10,STATUS='KEEP')


            do m=1,mh
            do n=1,nh
               dx(m,n)=1.
               dy(m,n)=1.
c========in new version kbt and kbtv  less by 1============					
					
					if (kbt(m,n).ne.1) then 
					kbt(m,n)=kbt(m,n)-1
c					kbt(m,n)=5
					end if

					if (kbtv(m,n).ne.1) then 
					kbtv(m,n)=kbtv(m,n)-1
c					kbtv(m,n)=5
					end if

            end do
            end do
      do m=1,mh
      do n=1,nh
           h(m,n)=h(m,n)*100
           h1(m,n)=h1(m,n)*100
           hl1(m,n)=dlx(m,n)*rad/dx(m,n)
           hl2(m,n)=dly(m,n)*rad/dy(m,n)
              dlx(m,n)=dlx(m,n)*rad*PI/180.
              dly(m,n)=dly(m,n)*rad*PI/180.

      end do
      end do
           do m=1,mh
           do n=1,nh
              dx(m,n)=dx(m,n)*PI/180.
              dy(m,n)=dy(m,n)*PI/180.
           end do
           end do
c=============old dz============
c       do  k=1,kh
c        z(k)=zz(k)*100.
c      end do 


c=============new dz============
      do  k=1,kh
        z1(k)=zz(k)*100.
      end do 

	do  k=1,kh-1
        z(k)=(z1(k)+z1(k+1))/2.
      end do 
	z(kh)=z(kh-1)

	      do  k=1,kh-1
               dz(k)=z(k+1)-z(k)
            end do
            dz(kh)=0.

c=============old dz============

c      do  k=2,KKZ
c        dzk(k)=(dz(k)+dz(k-1))*0.5
c      end do

c=============new dz============
      dzk(1)=0.
      do  k=2,KKZ		
        dzk(k)=z1(k+1)-z1(k)
      end do

      

      print *,z
	print *,dz
	print *,dzk
      read *
       do m=60,111
	kbt(m,2)=kbt(m,3)
	end do
	do m=60,111
	kbtv(m,2)=min(kbt(m,2),kbt(m+1,2))
	end do
	do m=60,111
	kbt(m,2)=max(kbtv(m,2),kbtv(m-1,2))
	end do


	


      do m=1,mh
      do n=1,nh
           kt=kbt(m,n)
	if(kt.gt.1)then
	     hn=0.
	        do k=2,kt
	           hn=hn+dzk(k)
	        end do
	      h(m,n)=hn
	end if
           kv=kbtv(m,n)
	if(kv.gt.1)then
	     hn=0.
	        do k=2,kv
	           hn=hn+dzk(k)
	        end do
	      h1(m,n)=hn
		  end if
	end do
	end do 


C*-------------------------------------------
        DO  n=1,nd
         do m=1,md
           tet=yy(m,n)*pi/180.
           tetv=yyu(m,n)*pi/180.
           SL(m,n)=2.*omega*sin(tet)
           SLB(m,n)=2.*omega*sin(tetv)
         end do
       end do   





         do k=2,kkz
           do n=1,nh-1 
              do m=1,mh-1
               SQ=(KBT(M,N)/K)*(KBT(M+1,N)/K)*
     *            (KBT(M,N+1)/K)*(KBT(M+1,N+1)/K)
c               sq=kbtv(m,n)/k
                    IF(SQ.GT.0.) then
                      ikbt(m,n,k)=1
                    else
                      ikbt(m,n,k)=0
                    end if
                 end do
              end do
          end do
         do k=1,kh
           do n=1,nh 
              do m=1,mh
                    IF(k.le.kbt(m,n)) then
                      inbt(m,n,k)=1
                    else
                      inbt(m,n,k)=0
                    end if
                 end do
              end do
          end do

	do k=1,kh
           do n=1,nh 
              do m=1,mh
                    IF(k.le.kbtv(m,n)) then
                      ikbt(m,n,k)=1
                    else
                      ikbt(m,n,k)=0
                    end if
                 end do
              end do
          end do


      return
      end
