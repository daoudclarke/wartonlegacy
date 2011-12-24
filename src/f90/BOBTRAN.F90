      program bobtran 
!     program to transpose matrix from quadratic sieve greta
! or gretan suite files
      dimension littr(20),nar(30000),n(60),iabn(30000)
      integer,allocatable,dimension (:,:)::itran 
      open (unit=7,file='gretpar',access='direct',form=&
      'formatted',recl=303,status='old')
      read (7,6,rec=1)iconr,ispecp,limprm,limprm2,mimprm,lenb2,&
      (n(jf),jf=1,60),irecnn,kkll,lenb1,nizz
      narc9=kkll+97
      icurmax=0
      open(unit=1,file ='gretf4',access='direct',&
      form='formatted',recl=narc9,status='old')
      nunit=1500000/irecnn
      ntem=kkll/nunit
      allocate (itran(nunit,irecnn))
6     format(2i6,i8,i8,i3,i6,60i4,3i6,i8)      
      open(unit=2,file='trandat1',access='sequential')
      open (unit=3,file='trandir',access='direct',form=&
      'formatted',recl=irecnn,status='old')
      goto 201
      
      
      
      
199   read(2,*,end=200)j,icur,(iabn(jk),jk=1,icur+1)
      print *,'j',j,'icur',icur,'iabn',iabn(1)
      goto 199
200   close (unit=2)
      close (unit=1)
      close (unit=7)
      stop
201   do i=1,10000
      iabn(i)=0
      end do
      
!      icur=2
!      iabn(1)=1
!      iabn(2)=2
!      j=1
!      write (2,*)j,icur,(iabn(jf),jf=1,icur)
!      rewind (unit=2)
!      read (2,*)i,idur,(iabn(jf),jf=1,idur)
!      print *,'iabn',(iabn(jf),jf=1,idur)
!      stop
      mm=irecnn
      nn=kkll
      
      
      
      
      read(1,1,rec=1)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      print *,kkb,ihitn
      
      
      
      do i=1,mm
      read(1,1,rec=i)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      itran(1,i)=kkb
      do j=2,nunit
      itran(j,i)=nar(j-1)
      end do
      
      end do
      do i=1,100
      write (3,2,rec=i)(itran(i,jk),jk=1,mm)
      end do
      do j=101,nunit
      icur=0
      do jf=1,mm
      if (itran(j,jf).eq.0)goto 101
      icur=icur+1
      iabn(icur)=jf
101   end do      
      write(2,*)j,icur,(iabn(jk),jk=1,icur+1)
      print *,'j',j,'icur',icur
      if (icur.le.icurmax)goto 102
      icurmax=icur
      jmax=j
102   end do
      print *,'icurmax',icurmax,'jmax',jmax
      
      do k=1,ntem-1
      print *,'k=',k
      js=k*nunit
      do i=1,mm
      read (1,1,rec=i)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      
      do j=1,nunit
      itran(j,i)=nar(js+j-1)
      end do
      end do
      
      do j=1,nunit
      icur=0
      do jf=1,mm
      if (itran(j,jf).eq.0)goto 111
      icur=icur+1
      iabn(icur)=jf
111   end do
      j2=js+j
      print *,'j2',j2,'icur',icur
      write(2,*)j2,icur,(iabn(jk),jk=1,icur+1)
      end do
      end do
      js=nunit*ntem
      if (js.eq.kkll)goto 11
      do i=1,mm
      read (1,1,rec=i)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
      do j=1,kkll-js
      itran(j,i)=nar(js+j-1)
      end do
      end do
      do j=js+1,kkll
      icur=0
      do jf=1,mm
      if (itran(j-js,jf).eq.0)goto 121
      icur=icur+1
      iabn(icur)=jf
121   end do
      
      
      print *,'j',j,'icur',icur
      write(2,*)j,icur,(iabn(jf),jf=1,icur+1)
      
      end do
      close (unit=2)
      close (unit=7)
      close (unit=1)
1     format(i1,i6,i10,20i4,30000i1)
2     format(30000i1)
11    end
