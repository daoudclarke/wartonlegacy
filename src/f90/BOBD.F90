      program bobd
      dimension nar(2000),nar2(4000)
      open(unit=1,file='compdata',access='direct',form=&
      'formatted',recl=196,status='old')
      open(unit=2,file='singles',access='direct',form=&
      'formatted',recl=1,status='old')
      
      
      
      
      read(1,1,rec=151)kia,kib,isgn,(nar(i),i=1,64),(nar2(j),j=1,119)
1     format(i6,i6,i1,64i1,119i1)
      if (isgn.eq.0)goto 2
      print *,'isgn',isgn
2     do i=1,64
      if (nar(i).eq.0)goto 10
      print *,'i',i,'nar',nar(i)
10    end do
      do i=1,119
      if (nar2(i).eq.0)goto 12
      print *,'i',i,'nar2',nar2(i)
12    end do
      stop
      jpos=0
      do i=1,2187
      jpos=jpos+500
      read(2,20,rec=jpos)jar
      if (jar.eq.0)goto 14
      print *,'i',i,'jar',jar
14    end do      
20    format(i1)      
      close(unit=1)
      close(unit=2)
      end
