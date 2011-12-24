       program hobxt
!      computation of null space for quadratic sieve method        
       
       
       dimension indic(1000),iqi(700,2263),iv(1000,10)
       dimension ic(2600)
       dimension nar(2600)
       dimension littr(20)
       open(unit=1,file='compdatc',access='direct',form=&
       'formatted',recl=2353,status='old')
       open(unit=2,file='kernel3',access='direct',form=&
       'formatted',recl=4810,status='old')
       open(unit=3,file='matches3',access='direct',form=&
       'formatted',recl=17460,status='old')
       
       mm=700
       nn=2263
       ipd=2
       do i=1,mm 
       do j=1,nn 
       
       
       iqi(i,j)=0
       end do
       end do
       
       read(1,1,rec=151)kkb,ihitn,(littr(jf),jf=1,20),(nar(ii),ii=1,nn-1)
       
       print *,'nar',(nar(ii),ii=1,nn-1)
       
       icon=1
       do i=1,mm
       
       read(1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(ii),ii=1,nn-1)
       
       iqi(i,1)=kkb
       do k1=1,nn-1
       iqi(i,k1+1)=nar(k1)
       end do
       
       
       
       
       
       end do
1      format(i1,i10,20i4,2262i1)       
              
       do i =1,mm
       indic(i) =0
       end do
              
       ir = 0
       do i=1,nn
       
       ic(i)=0
       end do
       do k = 1,mm
       print *,'k=',k,'ir=',ir
       
       
              
       do j =1,nn
       if((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40

34     do jj=1,nn
       if(jj.eq.j)goto 38
       
       
       mul=iqi(k,jj)
       do i=k,mm
       itt = iqi(i,jj) + iqi(i,j) *mul
       iqi(i,jj) =mod(itt,ipd)
       if(iqi(i,jj).ge.0)goto 36
       iqi(i,jj)=iqi(i,jj) +ipd

36     end do
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
       iv(is,ir) =iqi(k,jv)
42     end do
       print *,'vector found','ir',ir,'ivs',(iv(nq,ir),nq=1,mm)
       write(2,2,rec=ir)(iv(nq,ir),nq=1,mm)
2      format(962i5)     
       
       if (ir.lt.10)goto 50
       close(unit=1)
       close(unit=2)
       close(unit=3)
       stop
50     end do
       end

