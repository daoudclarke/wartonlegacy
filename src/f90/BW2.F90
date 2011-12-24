       program bbw2
        
       
       common  ibarray(100),isarray(100),igarray(100),inv(20)
       common  mult1(100),mult2(100),mult3(200),mult4(200),mult5(100)
       dimension indic(20),iq(20,20),iqi(20,20),iv(20,20)
       dimension irowp(20),irown(20),ip(20),ic(20)
       do i =1,20
       indic(i) =0
       end do
       ipd =13
       inv(1) =1
       inv(2) =7
       inv(3)=9
       inv(4) =10
       inv(5) =8
       inv(6) =11
       inv(7) =2
       inv(8) = 5
       inv(9) =3
       inv(10)=4
       inv(11) =6
       inv(12)=12
       ip(1) =0
       ip(2)=1
       ip(3) =0
       ip(4)=10
       ip(5)=10
       ip(6) =8
       ip(7) =2
       ip(8)=8
       
       do i=1,20
       do j=1,20
       iq(i,j) = 0
       end do
       end do
       do i =1,20
       irowp(i) =0
       end do
       irowp(8) =1
       do k = 2,8
       do j =1,ipd
       do i =1,8
       irown(i) =irowp(i+1) - ip(i) *irowp(1)
       end do
       do i=1,8
       irowp(i) =irown(i)
       end do
       end do
       do i =1,8
       irowp(9-i) = mod(irowp(9-i),ipd)
       if(irowp(9-i).ge.0)goto 10
       irowp(9-i)=irowp(9-i)+ipd

10     iq(k,i)=irowp(9-i)
       end do
       end do
       iq(1,1) =1
       do k = 1,8
       print *,iq(k,1),iq(k,2),iq(k,3),iq(k,4),iq(k,5),iq(k,6),iq(k,7),iq(k,8)
       end do
       do i=1,8
       do j =1,8
       iqi(i,j)=iq(i,j)
       end do
       end do
       do i=1,8
       iqi(i,i)=iqi(i,i)-1
       if(iqi(i,i).ge.0)goto 12
       iqi(i,i) =iqi(i,i) +ipd 
12     end do       
       ir = 0
       do i=1,20
       
       ic(i)=0
       end do
       do k = 1,8
       do j =1,8
       if((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
34     mul =inv(iqi(k,j))       
       do i =k,8
       iqi(i,j) =(-1) * iqi(i,j)  * mul
       iqi(i,j)=mod(iqi(i,j),ipd) 
       if(iqi(i,j).ge.0)goto 33
       iqi(i,j) =iqi(i,j) +ipd
33     end do
       do jj=1,8
       if(jj.eq.j)goto 38
        
       mul =iqi(k,jj)
       do i=k,8
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
       do is =1,8
       iv(is,ir) = 0
       end do
       iv(k,ir) =1
       do is =1,8
       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is,ir) =iqi(k,jv)
42     end do
50     end do
       do i=1,8
       print *,'iqis',iqi(i,1),iqi(i,2),iqi(i,3),iqi(i,4),iqi(i,5),&
       iqi(i,6),iqi(i,7),iqi(i,8)
       end do
       do jr =1,ir
       print *,'ivs',iv(1,jr),iv(2,jr),iv(3,jr),iv(4,jr),iv(5,jr),&
       iv(6,jr),iv(7,jr),iv(8,jr)
       end do
       stop
       
       ibarray(1) =1
       ibarray(2) =0
       ibarray(3) =1
       ibarray(4)=0
       ibarray(5) =10
       ibarray(6)=10
       ibarray(7)=8
       ibarray(8)=2
       ibarray(9) =8
       idegb = 8
       isarray(1) =1
       isarray(2)=5
       isarray(3)=9
       isarray(4)=0
       isarray(5)=5
       isarray(6)=5
       isarray(7)=0
       isarray(8)=7
       
       idegs=6
       
       
       
       
       do i = 1,idegs+1
       mult1(i)=isarray(i)
       mult5(i)=mult1(i)
       end do
       ipd =13
       idegm=idegs
       nn =2
       ie =4
230    if(ie.eq.1)goto 250       
       ie =ie/2
       do i=1,idegm +1
       mult2(i) =mult1(i)
       end do
       idegn = idegm 
       call multy(idegm,idegn,ipd)
       idegm =idegm +idegn
       do i = 1,idegm +1
       mult1(i) =mult3(i)
       end do


       
       if(nn.ge.ie)goto 240
       goto 230
240    nn = nn-ie   
       do i = 1,idegs
       mult2(i) = mult5(i)
       end do
       
       idegn =idegs
       call multy(idegm,idegn,ipd)
       idegm =idegm +idegn
       do i = 1,idegm +1
       mult1(i)=mult3(i)
       end do
       goto 230
       
250    do i =1,idegb+1   
       isarray(i) =ibarray(i)
       end do
       idegs = idegb
       do i = 1,idegm+1
       ibarray(i) =mult1(i)
       end do
       idegb =idegm
       ibarray(idegb+1) =ibarray(idegb +1) + 1
       if(ibarray(idegb+1).ge.0)goto 255
       ibarray(idegb +1) =ibarray(idegb+1) -ipd
255    call subgcd(idegs,idegb,idegg,ipd)
       print *,'idegg',idegg,idegs,idegb,isarray(1),ibarray(1)
       print *,igarray(1),igarray(2),igarray(3),igarray(4),&
       igarray(5),igarray(6),igarray(7),igarray(8),igarray(9)
       stop
       
       
       
       
       
       
       
       
       print *,mult3(1),mult3(2),mult3(3),mult3(4)
       print *,mult3(11),mult3(12),mult3(13),mult3(14),mult3(15),mult3(16)
       do i = 1,18
       print *,i,mult3(i)
       end do
       
       stop
       
       call subgcd(6,8,idegg,13)
       print *,igarray(1),igarray(2),igarray(3),igarray(4),igarray(5),&
       igarray(6),igarray(7),igarray(8),igarray(9)
       print *,'idegg',idegg
       end
       subroutine subgcd(idegs,idegb,idegg,ipd)
       
       common ibarray(100),isarray(100),igarray(100),inv(20)
       common mult1(100),mult2(100),mult3(200),mult4(200),mult5(100)
       dimension iws(100),iwb(100),itempb(100)
       icon=0
       do i= 1,idegs +1
       iws(i) = isarray(i)
       end do
       do i=1,idegb+1
       iwb(i) = ibarray(i)
       end do
       iwsd =idegs
       iwbd =idegb
10     loopl = iwbd-iwsd +1
       icon = icon+1
       mul = inv(iws(1))
       do i =1,loopl
       iqt= iwb(i) * mul
       iwb(i) =0
       do j =2,iwsd +1
       iwb(i+j-1)=iwb(i+j-1)-iws(j) * iqt
       iwb(i+j-1) =mod(iwb(i+j-1),ipd)
       print *,'imodd',iwb(i+j-1)
       
       if(iwb(i+j-1).ge.0)goto 12
       iwb(i+j-1)=iwb(i+j-1) +ipd
12     end do
       print *,iwb(1),iwb(2),iwb(3),iwb(4),iwb(5),iwb(6),iwb(7),iwb(8),iwb(9)
       
       end do
       print *,iwb(1),iwb(2),iwb(3),iwb(4),iwb(5),iwb(6),iwb(7),iwb(8),iwb(9)
       
       

       do i= 1,iwbd +1

       if(iwb(i).ne.0)goto 20
       end do
       
       goto 100
20     itempbd=iwbd+2 -i
       do jj =1,itempbd
       itempb(jj) = iwb(i+jj-1)
       
       
       end do
       print *,'itempb',itempbd,itempb(1),itempb(2),itempb(3),itempb(4),itempb(5)&
       ,itempb(6),itempb(7),itempb(8)
       
       do i= 1,iwsd + 1
       iwb(i) =iws(i)
       end do
       iwbd = iwsd
       do jj=1,itempbd
       iws(jj) = itempb(jj)
       end do
       iwsd =itempbd-1
       print *,'iw',iwsd,iwbd
       print *,'iws',iws(1),iws(2),iws(3),iws(4),iws(5),iws(6),iws(7)
       print *,iwb(1),iwb(2),iwb(3),iwb(4),iwb(5),iwb(6),iwb(7),iwb(8)
       if (icon .ne.4)goto 10
       
       goto 10
100    idegg=iwsd
       do i=1,10
       igarray(i)=0
       end do
       do i = 1,iwsd +1
       igarray(i) = iws(i)
       end do
       return
       end
       subroutine multy(idegm,idegn,ipd)
       common ibarray(100),isarray(100),igarray(100),inv(20)
       common mult1(100),mult2(100),mult3(200),mult4(200),mult5(100)
       do i=1,200
       mult3(i) =0
       mult4(i)=0
       end do
        
       do i =1,idegm +1
       do j =1,idegn+1
       mult3(i+j-1)=mult3(i+j-1) +mult1(i) *mult2(j)
       mult3(i+j-1) =mod(mult3(i+j-1),ipd)
       if(mult3(i+j-1).ge.0)goto 10
       mult3(i+j-1)= mult3(i+j-1) +ipd
10     end do 
       end do
       return
       end
