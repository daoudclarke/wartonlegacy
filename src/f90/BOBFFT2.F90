       program fftmp2
       double precision uni(40000,2),fta(40000,2),ftb(40000,2)
       double precision kcarr(40000,2),ftc(40000,2) 
       double precision root,bnum,pi
       dimension karr(40000),kbarr(40000),numb(40000)
       dimension inn1(40000),inn2(40000)
       root =2.0**0.1
       bnum=root**10
       print *,'root',root
       print *,'bnum=',bnum
       
       
       pi = 3.14159265389324d0
       print *,pi
       

       do i=1,380
       karr(i) =0
       kbarr(i)=0
       kcarr(i,1)=0
       kcarr(i,2)=0
       end do
       inn1(1)=0
       inn1(2)=17
       inn1(3)=5
       inn1(4)=1
       inn1(5)=2
       inn1(6)=3
       inn1(7)=4
       inn1(8)=5
       inn1(9)=6
       inn1(10)=7
       inn1(11)=8
       inn1(12)=9
       inn1(13)=1
       inn1(14)=2
       inn1(15)=3
       inn1(16)=4
       inn1(17)=6
       inn1(18)=3
       inn1(19)=1

       



       inn2(1)=0
       inn2(2)=29
       inn2(3)=1
       inn2(4)=4
       inn2(5)=7
       inn2(6)=3
       inn2(7)=2
       inn2(8)=2
       inn2(9)=6
       inn2(10)=5
       inn2(11)=3
       inn2(12)=2
       inn2(13)=1
       inn2(14)=1
       inn2(15)=4
       inn2(16)=5
       inn2(17)=3
       inn2(18)=1
       inn2(19)=7
       inn2(20)=3
       inn2(21)=3
       inn2(22)=1
       inn2(23)=3
       inn2(24)=5
       inn2(25)=3
       inn2(26)=2
       inn2(27)=8
       inn2(28)=2
       inn2(29)=3
       inn2(30)=8
       inn2(31)=3



       goto 2
       inn1(1)=0
       inn1(2)=15000
       inn2(1)=0
       inn2(2)=15000
       do jf=3,15002
       inn1(jf)=1
       inn2(jf)=2
       end do
2      ilens=inn2(2)
       do jf=3,ilens+2
       karr(ilens+3-jf)=inn2(jf)
       end do
       ilens=inn2(2)
       do jf=3,ilens+2
       kbarr(ilens+3-jf)=inn2(jf)
       end do
       goto 1

       
       
       
       
       karr(1)=3
       karr(2)=2
       karr(3)=1
       karr(4)=7
       kbarr(1)=4
       kbarr(2)=3
       kbarr(3)=2
       goto 1
       do i = 1,10
       karr(i) =1
       kbarr(i) =2
       end do
       
       
1      ilen =70
       uni(ilen,1) =1.0d0
       uni(ilen,2)=0.0d0
       do i =1,ilen-1
       uni(i,1) =cos(2*pi *i/ilen)
       uni(i,2)=sin(2*pi*i/ilen)
       
       end do
       goto 3
       uni(1,1)=0.7071067811865475241d0
       uni(3,1)=uni(1,1)*(-1.0)
       uni(1,2)=uni(1,1)
       uni(3,2)=uni(1,1)
       uni(5,1)=uni(3,1)
       uni(5,2)=uni(3,1)
       uni(7,1)=uni(1,1)
       uni(7,2)=uni(3,1)
       uni(4,1)=-1.0d0
       uni(6,2)=-1.0d0
       uni(4,2)=0.0d0
       uni(2,1)=0.0d0
       uni(2,2)=1.0d0
       uni(6,1)=0.0d0
       do i=1,ilen
       print *,uni(i,1),uni(i,2)
       end do
       stop










       
3      do i =1,ilen
       fta(i,1)= 0
       ftb(i,1)= 0
       ftc(i,1) =0
       fta(i,2)=0
       ftb(i,2)=0
       ftc(i,2)=0
       
       
       end do
       
       do i =1,ilen
       do j =1,ilen
       ipow = mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 10
       ipow =ipow +ilen
10     fta(i,1) =fta(i,1) +uni(ipow,1)*karr(j)
       fta(i,2)=fta(i,2)+uni(ipow,2)*karr(j)
       ftb(i,1) =ftb(i,1) +uni(ipow,1) *kbarr(j)
       ftb(i,2)=ftb(i,2)+uni(ipow,2)*kbarr(j)
       end do
       end do
       print *,'fta',(fta(jf,1),jf=1,ilen )
       print *,'ftb',(ftb(jf,1),jf=1,ilen)
       
       
       do i=1,ilen
       ftc(i,1) =fta(i,1) *ftb(i,1)- fta(i,2)*ftb(i,2)
       ftc(i,2) =fta(i,1)*ftb(i,2)+fta(i,2)*ftb(i,1)
       end do
       print *,'ftc1',(ftc(jf,1),jf=1,ilen)
       print *,'ftc2',(ftc(jf,2),jf=1,ilen)
       
       
       
       
       do i =1,ilen
       do j=1,ilen
       ipow =mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 20
       ipow =ilen
20     kcarr(i,1) =kcarr(i,1) +uni(ipow,1) *ftc(j,1) -uni(ipow,2)&
       *ftc(j,2)





       end do
       end do
       
       
       
30     do i=1,ilen-1
       if (i.ne.1)goto 32
       compt =kcarr(1,1)
32     kcarr(i,1)=kcarr(i+1,1)/ilen
       end do
       kcarr(ilen,1) =compt/ilen
       print *,'ans',(kcarr(i,1),i=1,ilen)
       
       rum=kcarr(16,1) +0.5
       num=int(rum)
       print *,'num',num
       do jf=1,ilen
       numb(jf)=0
       
       end do
       do jf=ilen,2,-1
       rum=kcarr(jf,1)+0.5
       numb(jf)=int(rum)
       end do
       
       
       do jf=ilen,2,-1
       itemp=numb(jf)
       numb(jf)=mod(itemp,10)
       
       numb(jf-1)=numb(jf-1)+int(itemp/10)
       end do
       
       do jf=ilen,1,-1
       numb(jf+2)=numb(jf)
       end do
       do jf=3,ilen+2
       if (numb(jf).ne.0)goto 40
       end do
       print *,'number zero'
       stop
40     ilen2=jf-2        
       do jk=ilen2+2,ilen+2
       numb(jk-ilen2+1)=numb(jk)
       end do
       numb(2)=ilen+1-ilen2
       print *,'numb',(numb(jk),jk=1,numb(2)+2)
       
       
       end
