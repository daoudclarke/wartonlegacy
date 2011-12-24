      program bbpcont
      dimension ibbp(30000),nbbp(30000),iconar(600),jconar(600),n(60)
      dimension kconar(600)
      do jf=1,30
      iconar(jf)=0
      jconar(jf)=0
      kconar(jf)=0
      end do
      open (unit=10,file='gretpar',access='direct',form=&
      'formatted',recl=303,status='old')
      read (10,10111,rec=1)iconr,ispecp,limprm,limprm2,mimprm,lenb2,&
      (n(jf),jf=1,60),izer1,izer2,lenb1,nizz
      print *,'lenb2',lenb2,'lenb1',lenb1
10111 format (2i6,i8,i8,i3,i6,60i4,3i6,i8)
      nkkl=3*lenb2
      ikkl=9*lenb2
      open (unit=8,file='nbbp',access='direct',form=&
      'formatted',recl=nkkl,status='old')
81111 format (20000i3)
      read (8,81111,rec=1)(nbbp(jf),jf=1,lenb2)
      open (unit=9,file='gretibbp',access='direct',form=&
      'formatted',recl=ikkl,status='old')
      read (9,91111,rec=1)(ibbp(jf),jf=1,lenb2)
91111 format (20000i9)
      do i=1,lenb2
      indic=ibbp(i)/100000+1
!      print *,'indic',indic
      iconar(indic)=iconar(indic)+1
      end do
      print *,'singles',(iconar(jf),jf=1,12)
      do i=1,lenb2
      if (nbbp(i).eq.1)goto 1
      indic=ibbp(i)/100000+1
      jconar(indic)=jconar(indic)+1
1     end do
      print *,'plurals',(jconar(jf),jf=1,12)
!      open (unit=1,file='gretabcd',access='direct',form=&
!      'formatted',recl=72,status='old')
!      write (1,1111,rec=1)(iconar(jf),jf=1,12)
1111  format (12i6)
!      close (unit=1,status='delete')
      do i=1,lenb2
      if (nbbp(i).lt.5)goto 2
      indic=ibbp(i)/100000+1
      kconar(indic)=kconar(indic)+1
2     end do      
      print *,'threes',(kconar(jf),jf=1,12)
      end
