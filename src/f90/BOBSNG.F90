      program bobsng
      dimension nar(2000),nar2(8000)
      open(unit=1,file='compdata',access='direct',form=&
      'formatted',recl =182,status='old')
      open(unit=2,file='singles',access='direct',form=&
      'formatted',recl=1,status='old')
      ihit =125
      
      
      
      do ii=1,ihit
      read(1,1,rec=ii)kia,kib,isgn,(nar(i),i=1,54),&
      (nar2(j),j=1,115)
1     format(i6,i6,i1,54i1,115i1)
      print *,'row no.=',ii
      imark=ii
      write(2,2,rec=imark)isgn
      
      do jj=2,55
      imark=(jj-1) *ihit +ii
      write(2,2,rec=imark)nar(jj-1)
      
      end do
      
      do jj=56,170
      imark=(jj-1)*ihit +ii
      
      write(2,2,rec=imark)nar2(jj-55)
2     format(i1)
      
      end do
      
      end do
      end
