      program bobnewd
      dimension iarray(200)
      kmax=20
      open(unit=1,file='bobnewb',access='direct',form=&
      'formatted',recl=kmax,status='old')
      do jf=1,200
      iarray(jf)=0
      end do
      do jf=1,20
      iarray(jf)=1
      end do
      kmax=20
1     format (200i1)
      write(1,1,rec=1)(iarray(jf),jf=1,kmax)
      end
