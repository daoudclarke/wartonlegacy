      program xxxx     
      dimension kr(10000,10)
      do iz=1,10000
      do iz2=1,10
      kr(iz,iz2)=iz *133
      end do
      end do
      
      open (unit=2,file='nargkr',access='sequential',form=&
      'unformatted')
      write (2)kr
      rewind(unit=2)
      read (2)kr
      print *,((kr(iz,iz2),iz=1,10000),iz2=1,10)
      close(unit=2)
      stop
1     open (unit=1,file='margkr',access='sequential')
      read (1,*)((kr(iz,iz2),iz=1,4500),iz2=1,10)
      print *,((kr(iz,iz2),iz=1,4500),iz2=1,10)
      close(unit=1)
      
      end
