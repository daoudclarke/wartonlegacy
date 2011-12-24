       program cobxt
        
       
       common  ibarray(100),isarray(100),igarray(100),inv(20)
       common  mult1(100),mult2(100),mult3(200),mult4(200),mult5(100)
       dimension indic(600),iq(20,20),iqi(600,600),iv(600,10)
       dimension irowp(20),irown(20),ip(20),ic(1000),jar(600),jjar(600)
       dimension jpre(600),nar(300),nar2(300),match(2000),izzy(100)
       open(unit=1,file='compdatb',access='direct',form=&
       'formatted',recl=528,status='old')
       open(unit=2,file='kernel2',access='direct',form=&
       'formatted',recl=2610,status='old')
       open(unit=3,file='matches2',access='direct',form=&
       'formatted',recl=5220,status='old')
       
       mm=522
       nn=516
       ipd=2
       do i=1,mm +10
       do j=1,nn +10
       
       
       iqi(i,j)=0
       end do
       end do
       
       read(1,1,rec=151)kia,kib,isgn,(nar(ii),ii=1,191),&
       (nar2(jj),jj=1,284),(izzy(jl),jl=1,39),lloc
       print *,'nar',(nar(ii),ii=1,191),'nar2',(nar2(jj),jj=1,284)
       
       icon=1
       do i=1,mm
       
       read(1,1,rec=i)kia,kib,isgn,(nar(ii),ii=1,191),(nar2(jj),jj=1,284)&
       ,(izzy(jl),jl=1,39),lloc
       iqi(i,1)=isgn
       do k1=1,191
       iqi(i,k1+1)=nar(k1)
       end do
       do k2=1,284
       iqi(i,k2+192)=nar2(k2)
       end do
       do k3=1,39
       iqi(i,k3+476) =izzy(k3)
       end do
       iqi(i,516)=1
       
       match(icon)=kia
       icon=icon+1
       match(icon)=kib
       icon=icon+1
       end do
1      format(i6,i6,i1,191i1,284i1,39i1,i1)       
       write(3,3,rec=1)(match(i2),i2=1,mm*2)
       print *,(match(jf),jf=1,mm*2)
       
3      format(1044i5)       
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
       
       do n1=k,mm
       
       jpre(n1)=iqi(n1,jj)
       end do
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
2      format(522i5)     
       
       if (ir.lt.10)goto 50
       close(unit=1)
       close(unit=2)
       close(unit=3)
       stop
50     end do
       end

