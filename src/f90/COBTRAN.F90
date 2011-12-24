      program cobtran 
!     program to transpose matrix from quadratic sieve :long records
      dimension itran(1000000),littr(20),nar(5000)
      open(unit=1,file ='compdatc',access='direct',&
      form='formatted',recl=3025,status='old')
      open(unit=2,file='frandat1',access='direct',&
      form='formatted',recl=949200,status='OLD')
      read(1,1,rec=1)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,2934)
      print *,kkb,ihitn
      
      
      mm=3164
      nn=2935
      do i=1,mm
      read(1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      
      itran(i)=kkb
      do j=2,300
      mark=(j-1)*mm +i
      itran(mark)=nar(j-1)
      end do
      
      end do
      
      write(2,2,rec=1)(itran(jk),jk=1,mm*300)
      
      
      
      
      
      
      
      
      
      
      
      
      do k=1,8
      do i=1,mm
      read (1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      
      js=k*300
      
      do j=1,300
      mark=(j-1)*mm+i
      itran(mark)=nar(js+j-1)
      end do
      end do
      
      
      write(2,2,rec=k+1)(itran(jk),jk=1,mm*300)
      print *,'k=',k
      
      end do
      do jf=1,mm*300
      itran(jf)=0
      end do
      
      do i=1,mm
      read (1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      do j=2701,2935
      mark=(j-2701)*mm+i
      itran(mark)=nar(j-1)
      end do
      end do
      
      write(2,2,rec=10)(itran(jk),jk=1,mm*300)
      print *,'k=10'
      
      
1     format(i1,i10,20i4,2934i1)
2     format(949200i1)
      end
