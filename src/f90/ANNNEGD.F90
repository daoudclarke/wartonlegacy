       program annnegd
       dimension idisarr(120),irecarr(1000),ireclen(150)
       open (unit=1,file='discrim',access='direct',form=&
       'formatted',recl=965,status='old')
       read (1,1001,rec=1)icont,(idisarr(jf),jf=1,120)
       print *,'discrims',(idisarr(jf),jf=1,icont)
       close (unit=1)
       stop
       open (unit=3,file='modjlen',access='direct',form=&
       'formatted',recl=600,status='old')
       read (3,1005,rec=1)(ireclen(jf),jf=1,79)
       print *,'reclen',(ireclen(jf),jf=1,120)
       stop


1005   format (1000i4)       
       open (unit=1,file='discrim',access='direct',form=&
       'formatted',recl=965,status='old')
       open (unit=2,file='modjpol',access='sequential')
       iccc=0
       iddd=0
200    iccc=iccc+1   
       limt=ireclen(iccc)
       read (2,*,end=400)(irecarr(jf),jf=1,limt)
       if (irecarr(3).ne.4)goto 201
       iddd=iddd+1
       print *,'irecarr',(irecarr(jf),jf=1,limt)
       print *,'recno',irecarr(1),'disc',irecarr(2),'deg',irecarr(3),iccc
201    goto 200
400    close (unit=2)       
       print *,'iccc=',iccc,'iddd',iddd
       stop
       read (1,1001,rec=1)icont,(idisarr(jf),jf=1,120)
       print *,'discrims',(idisarr(jf),jf=1,icont)
       close (unit=1)
       stop
1001   format(i5,120i8)
       do i=1,120
       idisarr(i)=0
       end do
       idis=51
       icont=0
55     ddis=idis
       rtd3=sqrt(ddis/3)
       irtd3=rtd3
       iib=mod(idis,2)
       iit=(iib*iib+idis)/4
       iia=iib
       idegm=1
2      if (iia.le.1)goto 1
3      if (mod(iit,iia).ne.0)goto 20
       if (iia.eq.iib)goto 6
       if (iia*iia.eq.iit)goto 6
       if (iib.eq.0)goto 6
       goto 10
6      idegn=1
       idegm=idegm+idegn
       goto 20
1      iia=1
       goto 20
63     iit=(iib*iib+idis)/4
       iia=iib
       goto 2
10     idegn=2
       idegm=idegm+idegn
20     iia=iia+1
       if (iia*iia.le.iit)goto 3
       iib=iib+2
       if (iib.le.irtd3)goto 63
       print *,'class no.=',idegm ,'for d=-',idis
       if (idegm.gt.4)goto 49
       icont=icont+1
       idisarr(icont)=idis
49     if (idis.eq.50000)goto 100
       if (mod(idis,4).eq.0)goto 50
       idis=idis+1
       goto 55
50     idis=idis+3
       goto 55
100    print *,'total no of discriminants whose class is <5 =',icont       
       write (1,1001,rec=1)icont,(idisarr(jf),jf=1,120)
       
       end
