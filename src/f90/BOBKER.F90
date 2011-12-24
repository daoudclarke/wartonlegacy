       program bobker
!      disk version of kernel routine        
       
       
       dimension indic(3500),iqi(3500),iv(3500,20)
       dimension ic(3500),jar(3500),jjar(3500)
       dimension iquin(3500),ipt(3500)
       open(unit=2,file='trandat1',access='direct',form=&
       'formatted',recl=3164,status='old')
       open(unit=3,file='qernel',access='direct',form=&
       'formatted',recl=15820,status='old')
       mm=3164
       nn=2935
       ipd=2
       do i=1,nn
       ipt(i)=0
       iquin(i)=0
       end do
       
       do i =1,mm
       indic(i) =0
       end do
              
       ir = 0
       do i=1,nn
       
       ic(i)=0
       end do
       do k = 1,mm
       print *,'k=',k,'ir=',ir
       
       do j2=1,nn
       if (k.eq.1)goto 1000
       if (k.eq.mm)goto 1000
       if (ipt(j2).ne.k)goto 1000
       iqi(j2)=iquin(j2)
       goto 1002

1000   read(2,2,rec=j2)(jjar(jk),jk=1,mm)
       iqi(j2)=jjar(k)
1002   end do
       
       
       
2      format(3164i1)       
3      format(3164i5)       
       do j =1,nn
       if((iqi(j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
34     read (2,2,rec=j)(jar(jk),jk=1,mm)
       do jj=1,nn          
       if (jj.eq.j)goto 38
       
       if (iqi(jj).eq.0)goto 38
       read(2,2,rec=jj)(jjar(jk),jk=1,mm)
       
       do i =k,mm
       itt=jjar(i)+jar(i)  
       jjar(i)=mod(itt,ipd)  
        
       
36     end do       
       write(2,2,rec=jj)(jjar(jk),jk=1,mm)
       if (k.eq.mm)goto 38
       iquin(jj)=jjar(k+1)
       ipt(jj)=k+1
38     end do
  


       ic(j) =k   
       indic(k) =j
       goto 50

40     ir =ir +1
       do is =1,mm
       iv(is,ir) = 0
       end do
       iv(k,ir) =1
       do is =1,mm

       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is,ir) =iqi(jv)
42     end do
       print *,'vector found',(iv(nq,ir),nq=1,mm)
       
       write(3,3,rec=ir)(iv(nq,ir),nq=1,mm)
       if (ir.lt.20)goto 50
       close(unit=2)
       close(unit=3)
       stop

50     end do
       close(unit=2)
       close(unit=3)
       end
