       program cobker
!      disk version of kernel routine using very long records       
!      obtained from cobtran       
       
       dimension indic(3500),iqi(3500),iv(3500,20)
       dimension ic(3500),jar(1000000),jjar(1000000)
       
       open(unit=2,file='frandat1',access='direct',form=&
       'formatted',recl=949200,status='old')
       open(unit=3,file='rernel',access='direct',form=&
       'formatted',recl=15820,status='old')
       mm=3164
       nn=2935
       ipd=2
       
       
       do i =1,mm
       indic(i) =0
       end do
              
       ir = 0
       do i=1,nn
       
       ic(i)=0
       end do
       lengr=mm*300
       do k = 1,100
       print *,'k=',k,'ir=',ir
       
       do kc=1,10
       read(2,2,rec=kc)(jjar(jk),jk=1,lengr)
       print *,'kc',kc
       do j2=1,300
       mark=(kc-1)*300+j2
       mark2=(j2-1)*mm+k
       iqi(mark)=jjar(mark2)
       end do
       end do

       
2      format(949200i1)       
3      format(3164i5)       
       do j =1,nn
       if((iqi(j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
34     jint=int(j/300)
       jrem=mod(j,300)
       
       read (2,2,rec=jint+1)(jar(jk),jk=1,lengr)
       do kr=1,10
       read(2,2,rec=kr)(jjar(jk),jk=1,lengr)
       do jr=1,300          
       jj=(kr-1)*300+jr
       if (jj.eq.j)goto 38
       
       if (iqi(jj).eq.0)goto 38
       iind=(jrem-1)*mm+k-1
       iind2=(jr-1)*mm+k-1
       
       do i =k,mm
       iind=iind+1
       iind2=iind+1
       itt=jjar(iind2)+jar(iind)  
       jjar(iind2)=mod(itt,ipd)  
        
       
       end do       
38     end do
       write(2,2,rec=kr)(jjar(jk),jk=1,lengr)
       print *,'kr=',kr
       end do
  


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
