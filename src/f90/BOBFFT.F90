       program fftmp
       complex  uni(10000),fta(10000),ftb(10000),ftc(10000),kcarr(10000)
       
       dimension karr(10000),kbarr(10000),numb(10000)
       dimension inn1(10000),inn2(10000)
       complex ci,compt
       ci =(0,1)
       pi = 3.14159265
       do i=1,2000
       karr(i) =0
       kbarr(i)=0
       kcarr(i)=(0,0)
       
       end do
       inn1(1)=0
       inn1(2)=9
       inn1(3)=1
       inn1(4)=2
       inn1(5)=3
       inn1(6)=4
       inn1(7)=5
       inn1(8)=6
       inn1(9)=7
       inn1(10)=9
       inn1(11)=1
       



       inn2(1)=0
       inn2(2)=9
       inn2(3)=1
       inn2(4)=2
       inn2(5)=3
       inn2(6)=4
       inn2(7)=5
       inn2(8)=6
       inn2(9)=8
       inn2(10)=0
       inn2(11)=3
       inn1(1)=0
       inn1(2)=200
       inn2(1)=0
       inn2(2)=200
       do jf=3,202
       inn1(jf)=1
       inn2(jf)=2
       end do
       ilens=inn1(2)
       do jf=3,ilens+2
       karr(ilens+3-jf)=inn1(jf)
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
       
       
1      ilen =402
       uni(ilen) =(1,0)
       
       do i =1,ilen-1
       uni(i) =exp(2*pi *ci*i/ilen)
       end do
       print *,'unis',(uni(i),i=1,ilen)
       
       do i =1,ilen
       fta(i)= (0,0)
       ftb(i)= (0,0)
       ftc(i) =(0,0)
       
       end do
       
       do i =1,ilen
       do j =1,ilen
       ipow = mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 10
       ipow =ipow +ilen
10     fta(i) =fta(i) +uni(ipow )*karr(j)
       ftb(i) =ftb(i) +uni(ipow) *kbarr(j)
       end do
       end do
       print *,(fta(i),i=1,ilen)
       print *,(ftb(i),i=1,ilen)
       
       do i=1,ilen
       ftc(i) =fta(i) *ftb(i)
       end do
       print *,'ftcs',(ftc(i),i=1,ilen)
       do i =1,ilen
       do j=1,ilen
       ipow =mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 20
       ipow =ilen
20     kcarr(i) =kcarr(i) +uni(ipow) *ftc(j)
       

       end do
       end do
       
       
       
30     do i=1,ilen-1
       if (i.ne.1)goto 32
       compt =kcarr(1)
32     kcarr(i)=kcarr(i+1)/ilen
       end do
       kcarr(ilen) =compt/ilen
       print *,'ans',(kcarr(i),i=1,ilen)
       rum=real(kcarr(16)) +0.5
       num=int(rum)
       print *,'num',num
       do jf=1,ilen
       numb(jf)=0
       
       end do
       do jf=ilen,2,-1
       rum=real(kcarr(jf))+0.5
       numb(jf)=int(rum)
       end do
       print *,'numb',(numb(jf),jf=1,ilen)
       stop
       do jf=ilen,2,-1
       itemp=numb(jf)
       numb(jf)=mod(itemp,10)
       
       numb(jf-1)=numb(jf-1)+int(itemp/10)
       end do
       print *,'numb',(numb(jf),jf=1,ilen)
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
