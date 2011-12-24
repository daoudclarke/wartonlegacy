       program a
       dimension ipr(200000)
       open(unit=3,file = 'volts.dat')
       do i=1,65002
       ipr(i)=0
       end do
       ia = 10
       ib = 23
       write(3, *)ia,ib
       print *,ia,ib
       close(unit=3)
       open(unit=3,file ='recl.dat',access='direct',form =&
       'formatted',recl=390000,status ='old')
       ipr(1)=2
       ipr(2)=3
       ipr(3)=5
       ipr(4)=7
       jend=4
       do i=11,800001,2
       do j=2,jend
       itst=mod(i,ipr(j))
       if(itst.eq.0)goto 10
       end do
       jend=jend+1
       ipr(jend)=i
10     end do
       print *,'primes',(ipr(j),j=1,jend)
       print *,'no of primes <800001=',jend
       
       do i=65000,1,-1
       ipr(i+1)=ipr(i)
       end do
       ipr(1)=jend
       print *,ipr(1),ipr(2),ipr(jend +1)

       
       ic = 9 *2743
       id=10*2743
       n=1
       write(3,5,rec=n) (ipr(i),i=1,65000)
5      format(65000i6)
       close(3)
       end
