       program bwar5
!      to print out file for quadratic sieving       
       dimension littr(20),normar(20),iabp(20),iabpn(20),ity(20)
       open(unit=4,file='rnx',access='sequential',position='rewind')
       icon=0
       


1      read(4,*,end=9)irecnn,kkb,ihitn,icur,(littr(i),i=1,20)
       
       
       print *,'irecnn=',irecnn,'kkb',kkb,'ihitn=',ihitn,'icur=',icur
       
       print *,'littr',(littr(i),i=1,20)
       
       read(4,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
       
       do i=1,icur
       print *,'i',i,'iabp',iabp(i),'iabpn',iabpn(i),'ity',ity(i)
       
       end do
111    icon=icon+1 
       if (icon.ne.21)goto 1
       goto 9
       
       stop
       goto 1
9      close(unit=4)
       end
