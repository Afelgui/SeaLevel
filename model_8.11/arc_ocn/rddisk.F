       SUBROUTINE rddisk(Temp,Sal,UU,VV,WW,Ftok,restart)
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      integer*4 idate, isec
      Common/date/idate,isec
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
     *       ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr)
     * ,Ftok(mh,nh)
      integer nu_diag
      common/iodiag/nu_diag
      integer *2 Kp
      Common/Gr1/Kp(mh,nh)
      Logical restart, oldmode
      Character*80 rname
      Character*80 rpath
      oldmode=.true.
      rpath='c:/users/elen/new_arctic/'
      Open(15,file='z:\dina\input data\start_TSUV\nT_3y.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((Temp(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file='z:\dina\input data\start_TSUV\nS_3y.dat',
     *   form='binary',status='unknown')
            read(15)((Sal(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file='z:\dina\input data\start_TSUV\nU_3y.dat',
     *   form='binary',status='unknown')
            read(15)((UU(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file='z:\dina\input data\start_TSUV\nV_3y.dat',
     *   form='binary',status='unknown')
            read(15)((VV(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file='z:\dina\input data\start_TSUV\nW_3y.dat',
     *   form='binary',status='unknown')
            read(15)((WW(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file='z:\dina\input data\start_TSUV\nStream_3y.dat'
     * 	,form='binary'
     *,status='unknown')
c            len=mh*nh*8
         read(15)((Ftok(m,n),m=1,mh),n=1,nh)
      CLOSE(unit=15,status='keep')
       myear=int(idate/10000)
       month=int((idate-10000*myear)/100)
       nday=idate-10000*myear-100*month
       Write(nu_diag,*) 'Starting date: year=',myear,' month=',month,
     * ' day=',nday
C*
      RETURN
      END
       SUBROUTINE rddisk1(Temp,Sal,UU,VV,WW,Ftok,restart)
C*-------------------------------------------
      Include 'model.par'
      COMMON/DAY/myear,month,nday
      integer*4 idate, isec
      Common/date/idate,isec
      real *8 Temp(nodes,kgr),Sal(nodes,kgr)
     *       ,UU(nodes,kgr),VV(nodes,kgr),WW(nodes,kgr)
     * ,Ftok(mh,nh)
      integer nu_diag
      common/iodiag/nu_diag
      integer *2 Kp
      Common/Gr1/Kp(mh,nh)
      Logical restart, oldmode
      Character*80 rname
      Character*80 rpath
      oldmode=.false.
      rpath='c:/users/elen/new_arctic/'
      Open(15,file='c:\dina\input data\start_TSUV\nTTT_17.dat',
     *   form='binary',status='unknown')
            len=nodes*kgr*8
            read(15)((Temp(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
      Open(15,file='c:\dina\input data\start_TSUV\nSSS_17.dat',
     *   form='binary',status='unknown')
            read(15)((Sal(i,k),i=1,nodes),k=1,kgr)
      CLOSE(unit=15,status='keep')
       myear=int(idate/10000)
       month=int((idate-10000*myear)/100)
       nday=idate-10000*myear-100*month
       Write(nu_diag,*) 'Starting date: year=',myear,' month=',month,
     * ' day=',nday

      RETURN
      END


