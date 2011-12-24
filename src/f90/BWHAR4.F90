       program a
       dimension ipr(100000)
       open(unit=3,file='volts.dat')
       read(3, *)ia,ib
       print *,ia,ib
       close(unit=3)
       open(unit=4,file ='recl.dat',access='direct',form=&
       'formatted',recl=390000,status='old')
       n=1
       read(4,5,rec=n) (ipr(i),i=1,65000)
5      format(65000i6)
       print *,(ipr(i),i=10832,10842)
       print *,'no of primes=',ipr(1)
       close(4)
       end
