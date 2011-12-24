      program bobmerg
!     for merging quadratic sieve files      
      dimension littr(20),iabp(50),iabpn(50),ity(50)
      
      
      
      open (unit=1,file='faxm',access='sequential')
      open (unit=2,file='fax2',access='sequential')
      open (unit=3,file='faxm2',access='sequential')
1     read (1,*,end=100)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
      read (1,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      write (3,*)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
      write (3,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      goto 1
100   rewind(unit=1)
2     read (2,*,end=200)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
      read (2,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      irecnn=irecnn+2701
      write(3,*)irecnn,kkb,ihitn,icur,(littr(j1),j1=1,20)
      write(3,*)irecnn,(iabp(i),iabpn(i),ity(i),i=1,icur)
      goto 2
200   close (unit=1)
      close (unit=2)
      close (unit=3)
      end
