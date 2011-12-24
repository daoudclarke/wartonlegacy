       program fftmp3
       double precision uni(40000,2),fta(40000,2),ftb(40000,2)
       double precision kcarr(40000,2),ftc(40000,2) 
       double precision triar(2,40000,2)
       double precision argr,pi,x1
       dimension karr(40000),kbarr(40000),numb(40000)
       dimension inn1(40000),inn2(40000),ibarr(60),ibarr2(60),ibarr3(60)
       dimension ibarr4(60)
       
       pi = 3.14159265389324d0
       print *,pi
       
       do i=1,20000
       karr(i) =0
       kbarr(i)=0
       kcarr(i,1)=0
       kcarr(i,2)=0
       end do
       inn1(1)=0
       inn1(2)=4
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
       inn2(2)=4
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
       inn1(2)=1000
       inn2(1)=0
       inn2(2)=1000
       do jf=3,1002
       inn1(jf)=1
       inn2(jf)=2
       end do
2      ilens=inn1(2)
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
       
       
1      ilen =8
       klen=3
       uni(ilen,1) =1.0d0
       uni(ilen,2)=0.0d0
       
       do i =1,ilen-1
       
       argr=(2*pi*i)/ilen
       uni(i,1) =dcos(argr)
       uni(i,2)=dsin(argr)
       print *,uni(i,1),uni(i,2)
       end do
       goto 3
7      uni(1,1)=0.7071067811329052087d0
       uni(1,2)=uni(1,1)
       uni(2,1)=0.0d0
       uni(2,2)=1.0d0
       uni(3,1)=-1.0*uni(1,1)
       uni(3,2)=uni(1,1)
       uni(4,1)=-1.0d0
       uni(4,2)=0.0d0
       uni(5,1)=uni(3,1)
       uni(5,2)=uni(3,1)
       uni(6,1)=0.0d0
       uni(6,2)=-1.0d0
       uni(7,1)=uni(1,1)
       uni(7,2)=0.0d0
       do i=1,ilen
       print *,i,uni(i,1),uni(i,2)
       end do
       


3      do i =1,ilen
       fta(i,1)= 0
       ftb(i,1)= 0
       ftc(i,1) =0
       fta(i,2)=0
       ftb(i,2)=0
       ftc(i,2)=0
       
       
       end do
       
       isw1=1
       do jf=1,ilen
       triar(1,jf,1)=kbarr(jf)
       triar(1,jf,2)=0.0d0
       end do
       do i=1,ilen
       do j=1,ilen
       ipow=mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 10
       ipow=ipow+ilen
10     fta(i,1)=fta(i,1)+uni(ipow,1)*karr(j)
       fta(i,2)=fta(i,2)+uni(ipow,2)*karr(j)
       ftb(i,1)=ftb(i,1)+uni(ipow,1)*kbarr(j)
       ftb(i,2)=ftb(i,2)+uni(ipow,2)*kbarr(j)
       
       end do
!       print *,'ftas',fta(i,1),fta(i,2)
       print *,'ftbs',ftb(i,1),ftb(i,2)
       end do
       
       isw1=2
       goto 250




99     do k=1,klen
       do jf=1,ilen
       jz=jf-1
       
       
       
       do jm=1,klen
       ibarr(jm)=0
       ibarr2(jm)=0
       ibarr3(jm)=0
       ibarr4(jm)=0
       end do
       ind1=1
       jq=jz
136    jp=jq/2
       itemp=jq-jp*2
       ibarr(ind1)=itemp
       jq=jp
       if (jp.eq.0)goto 138 
       ind1=ind1+1
       goto 136
138    do jv=1,ind1       
       ibarr2(jv)=ibarr(jv)
       ibarr3(jv)=ibarr(jv)
       ibarr4(jv)=0
       end do
       inj=klen-k+1
       ibarr2(inj)=0
       ibarr3(inj)=1
       do jv=1,k
       ibarr4(klen-jv+1)=ibarr(klen-k+jv)
       end do
       ind31=0
       ind32=0
       indn=0
       ipow=1
       do jv=1,klen
       
       ind31=ind31+ibarr2(jv)*ipow
       ind32=ind32+ibarr3(jv)*ipow
       indn=indn+ibarr4(jv)*ipow
       ipow=ipow*2
       end do
       
       
       
       if (mod(indn,ilen).eq.0)goto 2000
       indp=mod(indn,ilen)
       goto 2002
2000   indp=ilen

2002   triar(2,jf,1)=triar(1,ind31+1,1)+uni(indp,1)*triar(1,ind32+1,1)&
       -uni(indp,2)*triar(1,ind32+1,2)
       triar(2,jf,2)=triar(1,ind31+1,2)+uni(indp,1)*triar(1,ind32+1,2)&
       +uni(indp,2)*triar(1,ind32+1,1)
       i31=ind31+1
       i32=ind32+1
!       print *,'jf',jf,'indp',indp,'i31',i31,'i32',i32
       
!       print *,'fft3 fir triar',triar(2,jf,1),triar(2,jf,2)
       
       if (isw1.ne.3)goto 2005
       if (jf.ne.8)goto 2005
       if (k.ne.klen)goto 2005
       
       
       
2005   if (k.ne.klen)goto 991
       kcarr(indn+1,1)=triar(2,jf,1)
       kcarr(indn+1,2)=triar(2,jf,2)
       print *,'fft3 triar',triar(2,jf,1),triar(2,jf,2)
       



991    end do
       
       
       do jf=1,ilen
       do kf=1,2
       triar(1,jf,kf)=triar(2,jf,kf)
       end do
       end do
       
       
       end do
       print *,'fft3 kcarrs',(kcarr(jf,1),jf=1,ilen)
       stop
       
       
       
       if (isw1.eq.2)goto 250
       if (isw1.eq.3)goto 300
       
       isw1=2
       
       
       do jf=1,ilen
       do kf=1,2
       fta(jf,kf)=kcarr(jf,kf)
       end do
       end do
       do jf=1,ilen
       
       triar(1,jf,1)=kbarr(jf)
       triar(1,jf,2)=0.0d0
       end do
       goto 99

200    isw1=3       
       print *,'kcars',(kcarr(jf,1),jf=1,ilen)
       print *,'kcar2',(kcarr(jf,2),jf=1,ilen)
       
       do jf=1,ilen
       do kf=1,2
       ftb(jf,kf)=kcarr(jf,kf)
       end do
       end do
250    isw1=3
       do i=1,ilen
       ftc(i,1) =fta(i,1) *ftb(i,1)- fta(i,2)*ftb(i,2)
       ftc(i,2) =fta(i,1)*ftb(i,2)+fta(i,2)*ftb(i,1)
       end do
       print *,'ftc1',(ftc(jf,1),jf=1,ilen)
       print *,'ftc2',(ftc(jf,2),jf=1,ilen)
       
       do jf=1,ilen
       triar(1,jf,1)=ftc(jf,1)
       triar(1,jf,2)=ftc(jf,2)
       end do
       goto 99
       do i=1,ilen
       kcarr(i,1)=0.0
       end do
       
       do i=1,ilen
       do j=1,ilen
       ipow=mod((i-1)*(j-1),ilen)
       if (ipow.ne.0)goto 20
       ipow=ilen
20     kcarr(i,1)=kcarr(i,1)+uni(ipow,1) *ftc(j,1)-uni(ipow,2)&
       *ftc(j,2)
       end do
       end do

300    do i=1,ilen-1
       if (i.ne.1)goto 32
       compt =kcarr(1,1)
32     kcarr(i,1)=kcarr(i+1,1)/ilen
       end do
       kcarr(ilen,1) =compt/ilen
       print *,'kcarrs',(kcarr(jf,1),jf=1,ilen)
       
       
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
       numb(1)=0
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
       print *,'karr',(karr(jf),jf=1,ilen)
       print *,'kbarr',(kbarr(jf),jf=1,ilen)
       end
