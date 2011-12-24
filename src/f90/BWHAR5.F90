       program bwhar5
       dimension littr(20),normar(20),iabp(20),iabpn(20),ity(20)
       dimension n(60)
       open (unit=1,file='gretpar',access='direct',form=&
       'formatted',recl=303,status='old')
       read (1,10111,rec=1)iconr,ispecp,limprm,limprm2,mimprm,lenb2,&
       (n(jf),jf=1,60),izer1,izer2,lenb1,nizz
10111  format (2i6,i8,i8,i3,i6,60i4,3i6,i8)
       print *,'limprm',limprm,'iconr',iconr
       close (unit=1)
       stop
       open(unit=11,file='kernel',access='direct',form=&
       'formatted',recl=1000,status='old')
       open(unit=12,file='matches',access='direct',form=&
       'formatted',recl=2000,status='old')
       open(unit=13,file='newhits',access='sequential')
       open(unit=14,file='hits',access='sequential')
       open(unit=15,file='compdata',access='direct',form=&
       'formatted',recl=212,status='old')
       
       
       
       
       
       
       open (unit=1,file ='llrnel2',access='direct',form=&
       'formatted',recl=1000,status='old')
       open (unit=2,file='match4',access='direct',form=&
       'formatted',recl=2000,status='old')
       
       open(unit=4,file='tnx',access='sequential',position='rewind')
       icon=0
       open(unit=5,file='ternel2',access='direct',form=&
       'formatted',recl=15820,status='old')
1      read(4,*,end=9)irecnn,narc,kkb,ihitn,icur,(littr(i),i=1,20)
       
       print *,'irecnn=',irecnn,'kia=',kkb,'ihitn=',ihitn,'icur=',icur
       
       print *,'narc',narc,'littr',(littr(i),i=1,20)
       
       read(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
       do i=1,icur
!       print *,'i',i,'iabp',iabp(i)
       print *,'i',i,'iabp',iabp(i),'iabpn',iabpn(i),'ity(i)',ity(i)
       
       end do
       icon=icon+1 
       if (iabpn(1).eq.1)goto 1

!       if (icon.ne.28)goto 1
       
       
       stop
       goto 1
9      close(unit=4)
       end
