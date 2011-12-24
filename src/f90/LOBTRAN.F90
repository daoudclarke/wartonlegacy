      program lobtran 
!     program to create long records for quadratic sieve null space rtn.
      dimension ibigrec(1000000),littr(20),nar(5000)
      open(unit=1,file ='compdatc',access='direct',&
      form='formatted',recl=3025,status='old')
      open(unit=3,file='bigrec',access='direct',&
      form='formatted',recl=880500,status='old')
      read(1,1,rec=1)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,2934)
      print *,kkb,ihitn
      do jk=1,1000000
      ibigrec(jk)=0
      end do
      mm=3164
      nn=2935
      open (unit=2,file='compdate',access='direct',form=&
      'formatted',recl=3025,status='old')
      do i=1,mm
      read(1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,2934)
      write(2,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,2934)
      end do
      close(unit=1)
      close(unit=2)
      stop


      
      
      
      nk1=0
      nk2=1
      do i=1,mm
      read(1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      ind =nk2
      ind2=(ind-1)*nn
      ibigrec(ind2+1)=kkb
      do jk=1,nn-1
      ibigrec(ind2+1+jk)=nar(jk)
      end do
      if (nk2.eq.300)goto 10
      nk2=nk2+1
      goto 12
10    nk1=nk1+1
      write(3,3,rec=nk1)(ibigrec(jk),jk=1,100000)
      nk2=1
12    end do
1     format(i1,i10,20i4,2934i1)
3     format(880500i1)      
      end
      
      
      
      
      
      
      
      
      
      
