       
       program bbrcur1
!      variation of program using linear rcurrences -uses truncated
!      multiplication of polynomials
!      Recurrence algorithm to find kernel of matrix
       common karr(1600),kbarr(1600),kcarr(3200)
       
       common mat2(201,201),isol(200),ipol(500),mpol(500)
       common ibarray(1000),isarray(1000),igarray(1000)
       dimension matt(2,200,200),nmatt(200,200),ibarr(200),iyarr(200)
       dimension jbarr(200),iseqarr(401,200)
       dimension jsol(200),kchek(200),npol(500)
       dimension ibpol(500)
       goto 100
       do i=1,800
       karr(i)=1
       end do
       do i=1,800
       kbarr(i)=1
       end do
       karr(2)=0
       karr(3)=0
       karr(4)=0
       kbarr(2)=0
       kbarr(3)=0
       kbarr(4)=0
       ilen=750
       ilen2=750
       do jj=1,1000
       call mpmul0(ilen,ilen2,ilen3)
       end do
       print *,'prod1',(kcarr(jf),jf=1,ilen3)
       
       stop



       ibarray(1)=1
       ibarray(2)=0
       ibarray(3)=0
       ibarray(4)=0
       ibarray(5)=1
       ibarray(6)=1
       ibarray(7)=0
       ibarray(8)=0
       ibarray(9)=0
       ibarray(10)=1
       
       idegb=9
       isarray(1)=1
       isarray(2)=0
       isarray(3)=0
       isarray(4)=0
       isarray(5)=1
       idegs=4
       call polgcd(idegb,idegs,idegg)
!       print *,'idegg',idegg
!       print *,'gcd',(igarray(jf),jf=1,idegg+1)
       


       call bbw4(idegb,idegs,idegg,idegr)
100    maxpol=0
       ndeg=0
       ncpol=0
       npol(1)=1
       mpol(1)=0
       isw1=0
       isw2=0
       mmr=0
       nn=100

       npow=2*nn
       npowt=npow
       ibarr(1)=1
       ibarr(2)=1
       ibarr(3)=0
       ibarr(4)=0
       ibarr(5)=0
       ibarr(6)=1
       ibarr(7)=1
       do i=1,nn/2
       ibarr(i)=1
       end do
       do i=nn/2+1,nn
       ibarr(i)=0
       end do
       ibarr(3)=0
       ibarr(97)=1
       do i=1,nn
       iyarr(i)=0
       jbarr(i)=ibarr(i)
       end do
       imark=0
       
       do i=1,nn
       do j=1,nn
       nmatt(i,j)=0
       matt(1,i,j)=0
       
       end do
       end do

       do i=1,nn
       matt(1,i,i)=1
       end do

       do i=1,nn
       do j=1,i
       nmatt(i,j)=1
       end do
       end do
       
       do i=1,nn
       do j=1,nn
       matt(2,i,j)=nmatt(i,j)
       end do
       end do

       ibbig=1
5      do jf=1,nn       
       if (jbarr(jf).eq.1)goto 6
       end do
       goto 20
6      if (isw1.eq.1)goto 7
       isw1=1
       do ibig=1,2
       do i=1,nn
       isum=0
       do k=1,nn
       isum=isum+matt(ibig,i,k)*jbarr(k)
       end do
       iseqarr(ibig,i)=mod(isum,2)
       end do
       end do
       do ibig=2,npow
       
       do i=1,nn
       isum=0
       do k=1,nn
       isum=isum+matt(2,i,k)*iseqarr(ibig,k)
       end do
       iseqarr(ibig+1,i)=mod(isum,2)
       end do
       end do
7      if (isw2.eq.0)goto 107
       
       npow=npowt


       do i=1,npow+1
       if (iseqarr(i,ibbig).eq.0) goto 101
       ided1=i-1
       goto 102
101    end do
       print *,'gross sequence error'
       stop
107    do i=1,npow+1
       ipol(i)=iseqarr(i,ibbig)
       end do
       isw2=1
       goto 711

102    do i=ided1+1,npow+1
       karr(i-ided1)=iseqarr(i,ibbig)
       end do
       ilen=npow+1-ided1
       do i=1,ndeg+1
       if (npol(i).eq.0)goto 103
       ided2=i-1
       goto 104
103    end do
       print *,'nsus',(npol(jf),jf=1,ndeg+1)
       print *,'setup error'
       stop
104    ilen2=ndeg+1-ided2 
       do i=ided2+1,ndeg+1
       kbarr(i-ided2)=npol(i)
       end do
       call mpmul0(ilen,ilen2,ilen3)
       do i=1,npow+nn
       ibpol(i)=0
       end do
       do i=1,ilen3
       ibpol(i+ided1+ided2)=kcarr(i)
       end do
       iteml=ilen3+ided1+ided2
       print *,'product',(ibpol(jf),jf=1,iteml)

!       print *,'iseq',(ipol(jf),jf=1,npow),'ibbig',ibbig
!       if (ncpol.ne.5)goto 711
!       if (isw2.eq.0)goto 711
!       isw2=1
!       npow=npowt
!       do i=1,npow +nn
!       ibpol(i)=0
!       end do
!       do i=1,npow+1
!       do j=1,ndeg+1
!       if ((ipol(i).eq.0).or.(npol(j).eq.0))goto 714
!       ibpol(i+j-1)=ibpol(i+j-1)+1
! 714    end do
!       end do
       do i=ndeg+1,iteml
!       do i=ndeg+1,npow+ndeg+1
!       ipol(i-ndeg)=ibpol(i)
       ipol(i-ndeg)= ibpol(i)
       end do
       npow=npow-ndeg

711    ndegpr=ndeg
       iend=npow/2+1
       do i=1,iend
       do j=1,iend
       ind=i+j-1
       mat2(i,j)=ipol(ind)
       end do
!       print *,'i',i,'mat',(mat2(i,jk),jk=1,iend)
       end do
       if (ncpol.ne.5)goto 712
       
712    npow=npow
       call minpol(npow,mmr)
       
!       print *,'ipol',(ipol(jf),jf=1,npow+1)
       npow=npowt
      ncpol=ncpol+1 
       mpol(mmr+1)=1
       mdeg=mmr
       print *,'mpol',(mpol(jf),jf=1,mdeg+1),'ncpol',ncpol
!       print *,'npol',(npol(jf),jf=1,ndeg+1)
       if (mdeg.le.maxpol)goto 43
       maxpol=mdeg
       
!      instructions for gcds of polynomials       
43     a=a       
       goto 42
       if (mdeg.gt.ndeg)goto 40
       do i=1,ndeg+1
       ibarray(ndeg+2-i)=npol(i)
       end do
       idegb=ndeg
       do i=1,mdeg+1
       isarray(mdeg+2-i)=mpol(i)
       end do
       idegs=mdeg
       goto 41
40     do i=1,ndeg+1       
       isarray(ndeg+2-i)=npol(i)
       end do
       idegs=ndeg
       do i=1,mdeg+1
       ibarray(mdeg+2-i)=mpol(i)
       end do
       idegb=mdeg
41     call polgcd(idegb,idegs,idegg)
       print *,'idegg',idegg,'igarray',(igarray(jf),jf=1,idegg+1)
       do i=1,idegb+1
       npol(i)=ibarray(idegb+2-i)
       end do
       ndeg=idegb
       do i=1,idegs+1
       ibarray(i)=isarray(i)
       end do
       idegb=idegs
       do i=1,idegg+1
       isarray(i)=igarray(i)
       end do
       idegs=idegg
       call bbw4(idegb,idegs,idegg,idegr)
       print *,'degquo',idegg,'quo',(igarray(jf),jf=1,idegg+1)
       do i=1,idegg+1
       mpol(i)=igarray(idegg+2-i)
       end do
       mdeg=idegg
       mmr=mdeg





       
42     do i=1,2*nn+1
       ibpol(i)=0
       end do
       do i=1,mdeg+1
       if (mpol(i).eq.0)goto 121
       ided1=i-1
       goto 122
121    end do       
       print *,'gross error type 2'
       stop
122    do i=ided1+1,mdeg+1       
       karr(i-ided1)=mpol(i)
       end do
       ilen=mdeg+1-ided1
       do i=1,ndeg+1
       if (npol(i).eq.0)goto 131
       ided2=i-1
       goto 132
131    end do
       print *,'setup error 2'
       stop
132    do i=ided2+1,ndeg+1
       kbarr(i-ided2)=npol(i)
       end do
       ilen2=ndeg+1-ided2
       call mpmul0(ilen,ilen2,ilen3)
       do i=1,ilen3
       npol(i+ided1+ided2)=kcarr(i)
       end do
       print *,'secnpol',(npol(jf),jf=1,ilen3+ided1+ided2)
!       do i=1,mdeg+1
!       do j=1,ndeg+1
!       if ((mpol(i).eq.0).or.(npol(j).eq.0))goto 38
       
!       ibpol(i+j-1)=ibpol(i+j-1)+1
!       mind=i+j-1
!       print *,'ij',i,j,'ibpol',ibpol(i+j-1),'mulg',mulg,'mind',mind,&
!       'ibpol4',ibpol(4)
! 38     end do
!       end do
       ndeg=mdeg+ndeg
!       do i=1,ndeg+1
!       npol(i)=mod(ibpol(i),2)
!       end do
       
!       print *,'ndeg',ndeg
       print *,'npol',(npol(jf),jf=1,ndeg+1),'ibbig',ibbig
       ibbig=ibbig+1
       if ((ndeg.lt.nn).and.(ncpol.lt.nn))goto 5



!       if (ncpol.ne.5)goto 391
!       stop
391    a=a
       do i=1,nn
       iyarr(i)=0
       end do
       do k=2,ndeg+1
       if (npol(k).eq.0)goto 39
       do i=1,nn
       
       iyarr(i)=iyarr(i)+iseqarr(k-1,i)
       iyarr(i)=mod(iyarr(i),2)
       end do
39     end do
       print *,'iyarr',(iyarr(jf),jf=1,nn)
       

2      do i=1,nn
       isum=0
       do k=1,nn
       isum=isum+matt(2,i,k)*iyarr(k)
       end do
       isum=mod(isum,2)
       jbarr(i)=isum
       end do
       if (npol(1).eq.0)goto 13
       do jf=1,nn
       jbarr(jf)=jbarr(jf)+ibarr(jf)
       jbarr(jf)=mod(jbarr(jf),2)
       end do
13     print *,'jbarr',(jbarr(jk),jk=1,nn),'imark',imark
       do jf=1,nn

       end do
       do jf=1,nn
       if (jbarr(jf).eq.1)goto 9
       end do
       goto 20
! 9      end do
9      print *,'no solution'
       stop
10     print *,'matrix singular' 
       ntpol=ndeg+1
       print *,'ntpol',(npol(jf),jf=1,ndeg+1)
       print *,'ysols',(iyarr(jf),jf=1,nn)
       stop
       do i=1,nn
       isum=0
       do k=1,ndeg+1
       isum=isum+npol(k)*iseqarr(k,i)
       end do
       isum=mod(isum,2)
       jsol(i)=isum
       end do
       
       print *,'solutions',(jsol(jf),jf=1,nn)
       print *,'matrix singular'
       goto 280
20     do jf=1,nn
       jsol(jf)=iyarr(jf)
       end do
       if (npol(1).eq.0)goto 10
       do i=1,nn
       isum=0
       do k=1,nn
       isum=isum+nmatt(i,k)*jsol(k)
       end do
       kchek(i)=mod(isum,2)
       end do
       do i=1,nn
       if (kchek(i).eq.ibarr(i))goto 50
       goto 10
50     end do


       print *,'matrix nonsingular'
       ntpol=ndeg+1
       print *,'ntpol',ntpol,(npol(jf),jf=1,ndeg+1)
       print *,'solutions',(jsol(jf),jf=1,nn)
280    print *,'ncpol',ncpol,'ndeg',ndeg,'maxpol',maxpol,'ndegpr',ndegpr
       end
       subroutine minv(ising,mmr)
       common karr(1600),kbarr(1600),kcarr(3200)
       common mat2(201,201),isol(200),ipol(500),mpol(500)
       common ibarray(1000),isarray(1000),igarray(1000)
       dimension mat1(200,200),irhs1(200),jmark(200)
       dimension jx(200),markr(200),ichek(200)
       n=mmr
       
       do i =1,n
       
       markr(i) =0
       jmark(i) = 0
       jx(i) =0
       isol(i) =0
       end do
       
       do i =1,n
       irhs1(i) = mat2(i,mmr+1)
       end do
!       print *,'mmr',mmr,'irhs1',(irhs1(jf),jf=1,mmr)
       do i=1,n
       do j = 1,n
       mat1(i,j) = mat2(i,j)
       end do
       end do
       j =0
110    j = j+1
       do i =1,n
       if(i.eq.n)goto 140
       if(mat1(i,j).eq.0)goto 137
       if(markr(i).eq.1)goto 137
       jmark(j)=i
       markr(i) =1
       do ik =i+1,n
       if (markr(ik).eq.1)goto 134
       if(mat1(ik,j).eq.0)goto 134
       do jj =j,n
       if(mat1(i,jj).eq.0)goto 132
       
       mat1(ik,jj) =mat1(ik,jj) +mat1(i,jj)
132    end do
!       print *,'irhses', 'ik',ik,'i',i,'irhsik',irhs1(ik),'irhsi',irhs1(i)
       irhs1(ik) = irhs1(ik) +irhs1(i)
134    end do
       do li =1,n
       do lj =j,n
       mat1(li,lj)=mod(mat1(li,lj),2)
       if(mat1(li,lj).ge.0)goto 2000
       mat1(li,lj)=mat1(li,lj) +2
       
2000   end do 
       end do
       goto 150
137    end do
       goto 150
140    if(markr(i).eq.1)goto 149
!       print *,'matlast','ij',i,j,mat1(i,j)
       if(mat1(i,j).eq.1)goto 146
       jx(j) =1
       goto 150
146    jmark(j) =i       
       markr(i)=1
       goto 150
149    jx(j) =1       
150    a=a
!       print *,'j',j,'irhs',(irhs1(jf),jf=1,n)

       if(j.lt.n)goto 110
!       print *,'rhs',(irhs1(jf),jf=1,n)
       do iix =1,n
!       print *,'matrix',mat1(iix,1),mat1(iix,2),mat1(iix,3),mat1(iix,4)
       end do
!       print *,'jmarks',(jmark(jk),jk=1,n)
       if(jx(n).eq.1)goto 156
       im=jmark(n)
       isol(n) = irhs1(im)
156    do j =1,n-1
       ind =n-j+1
       ind2 =ind-1
       if(jx(ind2).eq.1)goto 180
       im=jmark(ind2)
       isol(ind2)=irhs1(im)
       do jj =ind,n
       isol(ind2) =isol(ind2)+isol(jj)*mat1(im,jj)
       
       end do
180    end do
!       print *,'sols',(isol(jf),jf=1,mmr),'mmr',mmr
       do j =1,n
       if(jx(j).eq.1)goto 190
       isol(j) =mod(isol(j),2)
       if(isol(j).ge.0)goto 190
       isol(j) =isol(j) +2
190    end do    
!       print *,'solutions',(isol(jf),jf=1,n)
!       print *,'jx',(jx(jf),jf=1,n)
       do i=1,n
       isum=0
       do k=1,n
       isum=isum+mat2(i,k)*isol(k)
       end do
       ichek(i)=mod(isum,2)
       end do
       do i=1,n
!       print *,'ichek',ichek(i),'irhs1',mat2(i,mmr+1)
       if (ichek(i).eq.mat2(i,mmr+1))goto 300
!       print *,'matrix singular'
       ising=1
       goto 302
300    end do
!       print *,'solution ok'
       ising=0
302    if (mmr.ne.3)goto 303
       
303    return



       end
       subroutine minpol(npow,mmr)
       common karr(1600),kbarr(1600),kcarr(3200)
       common mat2(201,201),isol(200),ipol(500),mpol(500)
       common ibarray(1000),isarray(1000),igarray(1000)
       mgap=1
       if (npow.lt.40)goto 1111
       mgap=npow/20
1111   injj=2
       if (ipol(1).eq.ipol(2))goto 10
! 11     kend=npow/2
!       do kk1=2,kend
!       mmr=kk1
      jincr=mgap
11    mlow=2       
      mmr=mlow
111    mmr=mmr+jincr
       call minv(ising,mmr)
       print *,'ising',ising,'mmr',mmr
       if (ising.eq.1)goto 99
       if (2*mmr.gt.npow)goto 87 
       do kk2=2*mmr+1,npow
       
       irh=ipol(kk2)
       isum=0
       ind=1
       do  kk3=kk2-mmr,kk2-1
       isum=isum+ipol(kk3)*isol(ind)
       ind=ind+1
       end do
       itst=mod(isum,2)
       if (mmr.ne.3)goto 791
!       print *,'itst',itst,'irh',irh,'kk2',kk2
       
791    if (itst.eq.irh)goto 79
!       goto 99
       goto 991
79     end do
       do i=1,mmr
       mpol(i)=isol(i)
       end do
       goto 100
       do i=1,mmr
       if (isol(i).eq.1)goto 671
       end do
       print *,'error polynomial'
       stop
671    injj2=i

       do i=injj2,mmr
       mpol(i+1-injj2)=isol(i)
       end do
       mmr=mmr+1-injj2
91     a=a       
!       print *,'minpol found',' mmr',mmr,(mpol(jf),jf=1,mmr)
       if (mmr.ne.3)goto 100
       


       goto 100
! 99     end do
99     if (injj.eq.2)goto 661       
       
       
       goto 111
661    injj=1
       mlow=mmr-mgap
       mmr=mlow
       jincr=1
       goto 111

991    if (injj.eq.2)goto 111       
       mlow=mmr
       jincr=mgap
       injj=2

       goto 111
          
       
       print *,'looks like limit reached'
       stop
       
       goto 100
10     do kk2=3,npow       
       if (ipol(kk2).ne.ipol(1))goto 11
       end do
       mmr=1
       mpol(1)=1
!       print *,'minpol found',' mmr',mmr,(mpol(jf),jf=1,mmr)
       goto 100
87     do i=1,mmr
       mpol(i)=isol(i)
       end do
       print *,'limit reached',' mmr',mmr,(mpol(jf),jf=1,mmr)
       

100    return
       end
!       subroutine polsq
!       dimension ipol(110),ipolt(110),ipoltp(110),ipoltd(100)
!       dimension ipolv(110),,ipolw(110),ipola(10,110)
!       ipcon=0
!       iee=1
!       do i=1,ipol(2)+2
!       ipoltp(i)=ipol(i)
!       end do
!       if (ipoltp(2).eq.1)goto 100
!       
!       do i=1,ipoltp(2)+2
!       ipolpd(i)=0
!       end do
!       ipdc=ipoltp(2)/2
!       ist=mod(ipoltp(2),2)
!       if (ist.eq.0)goto 6
!!      ist=1 
!       goto 7
!6      ist=2
!7      isw1=0
!       markd=0
!       do i=1,ipdc
!       if (ipoltp(ist).eq.0)goto 10
!       ipoltd(ist+1)=1
!       if (isw1.eq.1)goto 10
!       isw1=1
!       markd=ist+1
!10     ist=ist+2
!       end do
!       if (markd.eq.0)goto 14 
!       ipoltd(2)=ipoltp(2)+1-markd
!       idegs=ipoltd(2)-1
!       do i=1,idegs+1
!       isarray(i)=ipoltd(i+2)
!       end do
!       idegb=ipoltp(2)
!       do i=1,idegb+1
!       ibarray(i)=ipoltp(i+2)
!       end do
!       call polgcd(idegb,idegs,idegg)
!       ipolt(2)=idegg+1
!       do i=1,idegg+1
!       ipolt(i+2)=igarray(i)
!       end do
!       do i=1,idegg+1
!       isarray(i)=igarray(i)
!       end do
!       idegs=idegg
!       call bbw4(idegb,idegs,idegg,idegr)
!       do i=1,idegg+1
!       ipolv(i+2)=igarray(i)
!       end do
!       ipolv(2)=idegg+1
       
       
!       goto 16
! 14     do i=1,ipoltp(2)+2
!       ipolt(i)=ipoltp(i)
!       end do
!       ipolv(1)=0
!       ipolv(2)=1
!       ipolv(3)=1
! 16     k=0
!       if (ipolv(2).gt.1)goto 20
!       do i=1,ipolt(2)+2
!       ipoltd(i)=0
!       end do
!       do i=3,ipolt(2)+2
!       if (mod(ipolt(2)-i,2).eq.1)goto 19
!       idq=(ipolt(2)-i)/2
!       ipoltd(ipolt(2)-idq)=ipolt(i)
!       end do







      subroutine bbw4(idegb,idegs,idegg,idegr)
      common karr(1600),kbarr(1600),kcarr(3200)
      common mat2(201,201),isol(200),ipol(500),mpol(500)
      common ibarray(1000),isarray(1000),igarray(1000)
      dimension irarray(1000),iwb(1000),iws(1000)
      


      do i =1,idegs+1
      iws(i) =isarray(i)
      end do
      do i=1,idegb+1
      iwb(i) =ibarray(i)
      end do
      iwsd = idegs
      iwbd = idegb
      loopl = iwbd -iwsd +1
      
      do i =1,loopl
      igarray(i) =iwb(i) 
      iwb(i) =0
      if (igarray(i).eq.0)goto 19
      do j =2,iwsd + 1
      iwb(i+j-1) =iwb(i+j-1) +iws(j) 

      iwb(i+j-1) =mod(iwb(i+j-1),2)
      end do
19    end do
      do i =1,iwbd +1
      if(iwb(i).ne.0)goto 20
      end do
      goto 100
20    itempbd =iwbd +2 -i
      do jj =1,itempbd
      irarray(jj) =iwb(i+jj-1)
      end do
      idegr = itempbd -1
      goto 110
100   do i=1,100
      irarray(1) = 0
      end do
      idegr =-1
110   print *,'quos',(igarray(jf),jf=1,loopl)
      idegg=loopl-1
      if (idegr.eq.-1)goto 111
      print *,'rem',(irarray(jf),jf=1,itempbd)
111   return
      end
      subroutine polgcd(idegb,idegs,idegg)
      common karr(1600),kbarr(1600),kcarr(3200)
      common mat2(201,201),isol(200),ipol(500),mpol(500)
      common ibarray(1000),isarray(1000),igarray(1000)
      
      dimension iwb(1000),iws(1000),itempb(1000)
      
      
      do i=1,idegs+1
      iws(i)=isarray(i)
      end do
      do i=1,idegb+1
      iwb(i)=ibarray(i)
      end do
      iwsd=idegs
      iwbd=idegb
10    loopl=iwbd-iwsd+1 
      do i=1,loopl
      iqt=iwb(i)
      iwb(i)=0
      if (iqt.eq.0)goto 19
      do j=2,iwsd+1
      iwb(i+j-1)=iwb(i+j-1)+iws(j)
      iwb(i+j-1)=mod(iwb(i+j-1),2)
      end do
19    end do
      do i=1,iwbd+1
      if (iwb(i).ne.0)goto 20
      end do
      goto 100
20    itempbd=iwbd+2-i
      do jj=1,itempbd
      itempb(jj)=iwb(i+jj-1)
      end do
      do i=1,iwsd+1
      iwb(i)=iws(i)
      end do
      iwbd=iwsd
      do jj=1,itempbd
      iws(jj)=itempb(jj)
      end do
      iwsd=itempbd-1
      goto 10
100   idegg=iwsd
      do i=1,iwsd+1
      igarray(i)=iws(i)
      end do
      return
      end
      
      
      
      
      
      subroutine mpmul(ilen,ilen2,ilen3)
      common karr(1600),kbarr(1600),kcarr(3200)
      common mat2(201,201),isol(200),ipol(500),mpol(500)
      common ibarray(1000),isarray(1000),igarray(1000)
      
      dimension karrp(800),kbarrp(800),karra(2,800),kbarra(2,800)
      dimension ilena(2),ilena2(2),ibdf(1600)
      dimension iarr(1600),iadf(1600)
      ilenp=ilen
      ilenp2=ilen2
      
      do i=1,ilen
      karrp(i)=karr(i)
      end do
      do i=1,ilen2
      kbarrp(i)=kbarr(i)
      end do
      if (ilenp.gt.ilenp2)goto 1
      ilenas=ilenp2/2
!      print *,'p1 ilenas',ilenas
      if (mod(ilenp2,2).eq.0)goto 2
      ilenas=ilenas+1
      goto 2
1     ilenas=ilenp/2
      if (mod(ilenp,2).eq.0)goto 2
      ilenas=ilenas+1

2    a=a      
     if (ilenp.lt.ilenas+1)goto 210
      ilena(1)=ilenp-ilenas
      
      do jf=1,ilena(1)
      karra(1,jf)=karr(jf)
      end do
      
      
      
211   ired=ilena(1)+1 
      ijk=1
      do jf=ired,ired+ilenas-1
      if (karr(jf).eq.0)goto 204
      ilena(2)=ilenas+1-ijk
      do jk=1,ilena(2)
      karra(2,jk)=karr(jf-1+jk)
      end do
      goto 212
204   ijk=ijk+1
      end do
      ilena(2)=0

      goto 212
210   ilena(1)=0
      
      ilena(2)=ilenp
      do jf=1,ilenp
      karra(2,jf)=karr(jf)
      end do
      
      
212   a=a
!      print *,ilena(1),ilena(2),'ilenas',ilenas
      do i=1,2
      if (ilena(i).eq.0)goto 213
!      print *,'i',i,'karra',(karra(i,jf),jf=1,ilena(i))
213   end do
      
      
      if (ilenp2.lt.ilenas+1)goto 310
      ilena2(1)=ilenp2-ilenas
      do jf=1,ilena2(1)
      kbarra(1,jf)=kbarrp(jf)
      end do
      
311   ired=ilena2(1)+1
      ijk=1
      do jf=ired,ired+ilenas-1
      if (kbarr(jf).eq.0)goto 304
      ilena2(2)=ilenas+1-ijk
      do jk=1,ilena2(2)
      kbarra(2,jk)=kbarr(jf-1+jk)
      end do
      goto 312
304   ijk=ijk+1 
      end do
      ilena2(2)=0
      goto 312
310   ilena2(1)=0
      
      ilena2(2)=ilenp2
      do jf=1,ilenp2
      kbarra(2,jf)=kbarr(jf)
      end do





      
      
312   a=a
      do i=1,2
      if (ilena2(i).eq.0)goto 313
!      print *,'ilen2a',ilena2(1),ilena2(2)
!      print *,'i',i,'kbarra',(kbarra(i,jf),jf=1,ilena2(i))
313   end do
      
      if ((ilena(1).eq.0).or.(ilena2(1).eq.0))goto 535
      do jf=1,ilena(1)
      karr(jf)=karra(1,jf)
      end do
      ilen=ilena(1)
      do jf=1,ilena2(1)
      kbarr(jf)=kbarra(1,jf)
      end do
      ilen2=ilena2(1)
      isub=2
      call mpmul2(ilen,ilen2,ilen3)
!      print *,'ilen3',ilen3,'mul1',(kcarr(jf),jf=1,ilen3)
      
      goto 540
535   iarr(1)=0
      iarr(2)=0
      goto 545
540   do jf=1,ilen3
      iarr(jf+2)=kcarr(jf)
      end do
      

!      print *,'point  1'
      iarr(1)=0
      iarr(2)=ilen3
      karr(1)=iarr(1)
      karr(2)=iarr(2)+ilenas
      do jf=1,ilen3
      karr(jf+2)=iarr(jf+2)
      end do
      do jf=ilen3+3,ilen3+2 +ilenas 
      karr(jf)=0
      end do
      kbarr(1)=0 
      kbarr(2)=karr(2)+ilenas
      do jf=3,karr(2)+2
      kbarr(jf)=karr(jf)
      end do
      do jf=karr(2)+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'iarr1',iarr(2),'shift1',(iarr(jf),jf=1,iarr(2)+2)
      

545   if (ilena(1).eq.0)goto 546
      do jf=1,ilena(1)
      karr(jf+2)=karra(1,jf)
      end do
      karr(1)=0
      karr(2)=ilena(1)
      goto 547
546   karr(1)=0
      karr(2)=0
      
547   if (ilena(2).eq.0)goto 548
      do jf=1,ilena(2)
      kbarr(jf+2)=karra(2,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena(2)
      goto 549
548   kbarr(1)=0
      kbarr(2)=0
549   call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      iadf(jf)=kcarr(jf)
      end do
      if (ilena2(2).eq.0)goto 5495
      do jf=1,ilena2(2)
      karr(jf+2)=kbarra(2,jf)
      end do
      karr(1)=0
      karr(2)=ilena2(2)
      goto 5497
5495  karr(1)=0       
      karr(2)=0
5497  if (ilena2(1).eq.0)goto 5498
      
      do jf=1,ilena2(1)
      kbarr(jf+2)=kbarra(1,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena2(1)
      goto 5499
5498  kbarr(1)=0      
      kbarr(2)=0
5499  call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int2',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=3,kcarr(2)+2
      kbarr(jf-2)=kcarr(jf)
      end do 
      ilen2=kcarr(2)
      isgn=kcarr(1)
      do jf=3,iadf(2)+2
      karr(jf-2)=iadf(jf)
      end do
      ilen=iadf(2)
      isub=3
      call mpmul2(ilen,ilen2,ilen3)
!      print *,'2ilen3',ilen3,'kcarr',(kcarr(jf),jf=1,ilen3)
      
550   do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      
      
      
      kbarr(1)=mod(iadf(1)+isgn,2)
!      print *,'kbarr1',kbarr(1)
      kbarr(2)=ilen3+ilenas
      
      
      do jf=ilen3+3,ilen3+2+ilenas
      kbarr(jf)=0
      end do
!      print *,'kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
      
!      print *,'iarr2',iarr(2)
      
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'second iarr',iarr(2),'arr',(iarr(jf),jf=1,iarr(2)+2)
      
552   if (ilena(2).eq.0)goto 553      
      do jf=1,ilena(2)
      karr(jf)=karra(2,jf)
      end do
      ilen=ilena(2)
      if (ilena2(2).eq.0)goto 553


      do jf=1,ilena2(2)
      kbarr(jf)=kbarra(2,jf)
      end do
      ilen2=ilena2(2)
      isub=4 
      call mpmul2(ilen,ilen2,ilen3)
!      print *,'3ilen3',ilen3,'kcarr',(kcarr(jf),jf=1,ilen3)
      
560   do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      kbarr(jf+2)=kcarr(jf)
      end do
      
      karr(1)=0
      kbarr(1)=0
      karr(2)=ilen3
      
      kbarr(2)=ilen3+ilenas
      do jf=ilen3+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kbarr(jf)=kcarr(jf)
      end do
!      print *,'last',(kcarr(jf),jf=1,kcarr(2)+2)
      
!      print *,'iarrlen',iarr(2)
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      ilen3=kcarr(2)
!      print *,'fir ilen3',ilen3,'arr',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=3,ilen3 +2
      kcarr(jf-2)=kcarr(jf)
      
      end do
!      print *,'kcarrml',(kcarr(jf),jf=1,ilen3)
      
      goto 110
553   ilen3=iarr(2)
      do jf=1,ilen3
      kcarr(jf)=iarr(jf+2)
      end do
      goto 110
      stop
      
      jdf=ilenp-(i-1)*ilenas
      

!      if ((ilen.gt.80).or.(ilen2.gt.80))goto 9
      goto 10
9     a=a
10    do i=1,ilen+ilen2
      kcarr(i) =0
      end do
!      print *,'two',ilen,ilen2
      do i =1,ilen
      do j=1,ilen2
      itemp =karr(i) *kbarr(j)
      
      itemp2=int(itemp/10000)
      irem1 =itemp-itemp2 *10000
      
      
      kcarr(i+j)=kcarr(i+j) +irem1
      kcarr(i+j-1)=kcarr(i+j-1) +itemp2
      
      if (kcarr(i+j).lt.10000)goto 20
      kcarr(i+j)=kcarr(i+j)-10000
      kcarr(i+j-1) =kcarr(i+j-1)+1
20    do k =1,i+j-2
      if (kcarr(i+j-k).lt.10000)goto 22
      kcarr(i+j-k)=kcarr(i+j-k) -10000
      kcarr(i+j-k-1)=kcarr(i+j-k-1)+1
      end do
22    end do
      end do
      ilen3 =ilen +ilen2
!      print *,'four ilen3',ilen3
      if(kcarr(1).ne.0)goto 100
      ilen3 =ilen3-1
!      print *,'three ilen3',ilen3
      do i=1,ilen3
      kcarr(i)=kcarr(i+1)
      end do

100   if (isub.eq.2)goto 540
      if (isub.eq.3)goto 550
      if (isub.eq.4)goto 560
110   ilen=ilenp
      ilen2=ilenp2
      do jf=1,ilen
      karr(jf)=karrp(jf)
      end do
      do jf=1,ilen2
      kbarr(jf)=kbarrp(jf)
      end do
      return


      end

      
      subroutine mpmul2(ilen,ilen2,ilen3)
      common karr(1600),kbarr(1600),kcarr(3200)
      common mat2(201,201),isol(200),ipol(500),mpol(500)
      common ibarray(1000),isarray(1000),igarray(1000)
      dimension karrp(800),kbarrp(800),karra(2,800),kbarra(2,800)
      dimension ilena(2),ilena2(2),ibdf(1600)
      dimension iarr(1600),iadf(1600)
      ilenp=ilen
      ilenp2=ilen2
      
      
      do i=1,ilen
      karrp(i)=karr(i)
      end do
      do i=1,ilen2
      kbarrp(i)=kbarr(i)
      end do
      if (ilenp.gt.ilenp2)goto 1
      ilenas=ilenp2/2
      if (mod(ilenp2,2).eq.0)goto 2
      ilenas=ilenas+1
      goto 2
1     ilenas=ilenp/2
      if (mod(ilenp,2).eq.0)goto 2
      ilenas=ilenas+1

2    a=a      
     if (ilenp.lt.ilenas+1)goto 210
      ilena(1)=ilenp-ilenas
      
      do jf=1,ilena(1)
      karra(1,jf)=karr(jf)
      end do
      
      
      
211   ired=ilena(1)+1 
      ijk=1
      do jf=ired,ired+ilenas-1
      if (karr(jf).eq.0)goto 204
      ilena(2)=ilenas+1-ijk
      do jk=1,ilena(2)
      karra(2,jk)=karr(jf-1+jk)
      end do
      goto 212
204   ijk=ijk+1
      end do
      ilena(2)=0

      goto 212
210   ilena(1)=0
      
      ilena(2)=ilenp 
      do jf=1,ilenp
      karra(2,jf)=karr(jf)
      end do
      
      
      
      
212   a=a
!      print *,ilena(1),ilena(2),'ilenasmul2',ilenas
      do i=1,2
      if (ilena(i).eq.0)goto 213
!      print *,'i',i,'karra',(karra(i,jf),jf=1,ilena(i))
213   end do
      
!      print *,'mul2',' ilenp',ilenp,'ilenp2',ilenp2
      if (ilenp2.lt.ilenas+1)goto 310
      ilena2(1)=ilenp2-ilenas
      do jf=1,ilena2(1)
      kbarra(1,jf)=kbarrp(jf)
      end do
      
      
311   ired=ilena2(1)+1
      ijk=1
      do jf=ired,ired+ilenas-1
      if (kbarr(jf).eq.0)goto 304
      ilena2(2)=ilenas+1-ijk
      do jk=1,ilena2(2)
      kbarra(2,jk)=kbarr(jf-1+jk)
      end do
      goto 312
304   ijk=ijk+1 
      end do
      ilena2(2)=0
      goto 312
310   ilena2(1)=0
      
      ilena2(2)=ilenp2
      do jf=1,ilenp2
      kbarra(2,jf)=kbarr(jf)
      end do
      
      
      
312   a=a
      do i=1,2
      if (ilena2(i).eq.0)goto 313
!      print *,'ilen2a',ilena2(1),ilena2(2)
!      print *,'i',i,'kbarra',(kbarra(i,jf),jf=1,ilena2(i))
313   end do
      
      if ((ilena(1).eq.0).or.(ilena2(1).eq.0))goto 535
      do jf=1,ilena(1)
      karr(jf)=karra(1,jf)
      end do
      ilen=ilena(1)
      do jf=1,ilena2(1)
      kbarr(jf)=kbarra(1,jf)
      end do
      ilen2=ilena2(1)
      isub=2
      goto 10
535   iarr(1)=0
      iarr(2)=0
      goto 545
540   do jf=1,ilen3
      iarr(jf+2)=kcarr(jf)
      end do
      

!      print *,'point  1'
      iarr(1)=0
      iarr(2)=ilen3
      karr(1)=iarr(1)
      karr(2)=iarr(2)+ilenas
      do jf=1,ilen3
      karr(jf+2)=iarr(jf+2)
      end do
      do jf=ilen3+3,ilen3+2 +ilenas 
      karr(jf)=0
      end do
      kbarr(1)=0 
      kbarr(2)=karr(2)+ilenas
      do jf=3,karr(2)+2
      kbarr(jf)=karr(jf)
      end do
      do jf=karr(2)+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'iarr1',iarr(2)
      

545   if (ilena(1).eq.0)goto 546
      do jf=1,ilena(1)
      karr(jf+2)=karra(1,jf)
      end do
      karr(1)=0
      karr(2)=ilena(1)
      goto 547
546   karr(1)=0
      karr(2)=0
      
547   if (ilena(2).eq.0)goto 548
      do jf=1,ilena(2)
      kbarr(jf+2)=karra(2,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena(2)
      goto 549
548   kbarr(1)=0
      kbarr(2)=0
549   call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      iadf(jf)=kcarr(jf)
      end do
      if (ilena2(2).eq.0)goto 5495
      do jf=1,ilena2(2)
      karr(jf+2)=kbarra(2,jf)
      end do
      karr(1)=0
      karr(2)=ilena2(2)
      goto 5497
5495  karr(1)=0       
      karr(2)=0
5497  if (ilena2(1).eq.0)goto 5498
      
      do jf=1,ilena2(1)
      kbarr(jf+2)=kbarra(1,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena2(1)
      goto 5499
5498  kbarr(1)=0      
      kbarr(2)=0
5499  call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int2',(kcarr(jf),jf=1,kcarr(2)+2)
      do jf=3,kcarr(2)+2
      kbarr(jf-2)=kcarr(jf)
      end do 
      ilen2=kcarr(2)
      isgn=kcarr(1)
      do jf=3,iadf(2)+2
      karr(jf-2)=iadf(jf)
      end do
      ilen=iadf(2)
      isub=3
      goto 10
550   do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      
      
      
      kbarr(1)=mod(iadf(1)+isgn,2)
!      print *,'kbarr1',kbarr(1)
      kbarr(2)=ilen3+ilenas
      
      
      do jf=ilen3+3,ilen3+2+ilenas
      kbarr(jf)=0
      end do
!      print *,'kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
      
!      print *,'iarr2',iarr(2)
      
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'second iarr',iarr(2)
      
552   if (ilena(2).eq.0)goto 553      
      do jf=1,ilena(2)
      karr(jf)=karra(2,jf)
      end do
      ilen=ilena(2)
      if (ilena2(2).eq.0)goto 553


      do jf=1,ilena2(2)
      kbarr(jf)=kbarra(2,jf)
      end do
      ilen2=ilena2(2)
      isub=4 
      goto 10
560   do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      kbarr(jf+2)=kcarr(jf)
      end do
      
      karr(1)=0
      kbarr(1)=0
      karr(2)=ilen3
      
      kbarr(2)=ilen3+ilenas
      do jf=ilen3+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kbarr(jf)=kcarr(jf)
      end do
!      print *,'last',(kcarr(jf),jf=1,kcarr(2)+2)
!      print *,'iarrlen',iarr(2)
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      ilen3=kcarr(2)
!      print *,'fir ilen3',ilen3
      
      do jf=3,ilen3 +2
      kcarr(jf-2)=kcarr(jf)
      
      end do
!      print *,'kcarrml',(kcarr(jf),jf=1,ilen3)
      
      goto 110
553   ilen3=iarr(2)
      do jf=1,ilen3
      kcarr(jf)=iarr(jf+2)
      end do
      goto 110
      stop





      
      
      jdf=ilenp-(i-1)*ilenas




!      if ((ilen.gt.80).or.(ilen2.gt.80))goto 9
      goto 10
9     a=a
10    do i=1,ilen+ilen2
      kcarr(i) =0
      end do
!      print *,'two',ilen,ilen2
      do i =1,ilen
      do j=1,ilen2
      if ((karr(i).eq.0).or.(kbarr(j).eq.0))goto 600
      kcarr(i+j-1)=kcarr(i+j-1)+1
600   end do      
      end do
      ilen3=ilen+ilen2-1
      do i=1,ilen3-1
      kcarr(i)=mod(kcarr(i),2)
      end do
      
!      print *,'four ilen3',ilen3
      

100   if (isub.eq.2)goto 540
      if (isub.eq.3)goto 550
      if (isub.eq.4)goto 560
110   ilen=ilenp
      ilen2=ilenp2
      do jf=1,ilen
      karr(jf)=karrp(jf)
      end do
      do jf=1,ilen2
      kbarr(jf)=kbarrp(jf)
      end do
      return


      end

      
       



      subroutine mpadd(isora)
      common karr(1600),kbarr(1600),kcarr(3200)
      common mat2(201,201),isol(200),ipol(500),mpol(500)
      common ibarray(1000),isarray(1000),igarray(1000)
      if(karr(2).eq.0)goto 50
      if(kbarr(2).eq.0)goto 55
      if(kbarr(2).gt.karr(2))goto 2
      do i=1,karr(2)+2
      kcarr(i)=karr(i)
      end do
      indic=karr(2)-kbarr(2)
      do i=kbarr(2)+2,3,-1
      kcarr(i+indic)=kcarr(i+indic)+kbarr(i)
      kcarr(i+indic)=mod(kcarr(i+indic),2)
      end do
      kcarr(1)=0
9     do jf=3,kcarr(2)+2
      if (kcarr(jf).eq.0)goto 10
      goto 11
10    end do      
      kcarr(2)=0
      goto 100
11    idiff=jf-3      
      if (idiff.eq.0)goto 100
      do i=3,kcarr(2)+2-idiff
      kcarr(i)=kcarr(i+idiff)
      end do
      kcarr(2)=kcarr(2)-idiff
      goto 100

      
      goto 100
2     do i=1,kbarr(2)+2       
      kcarr(i)=kbarr(i)
      end do
      indic=kbarr(2)-karr(2)
      do i=karr(2)+2,3,-1
      kcarr(i+indic)=kcarr(i+indic)+ karr(i)
      kcarr(i+indic)=mod(kcarr(i+indic),2)
      end do
      kcarr(1)=0
      goto 9
      
      
50    do i=2,kbarr(2)+2
      kcarr(i)=kbarr(i)
      end do
      
      kcarr(1)=0
      goto 100
55    do i=1,karr(2)+2
      kcarr(i)=karr(i)
      end do
      kcarr(1)=0
100   return
      end
      subroutine mpmul0(ilen,ilen2,ilen3)
      common karr(1600),kbarr(1600),kcarr(3200)
      common mat2(201,201),isol(200),ipol(500),mpol(500)
      common ibarray(1000),isarray(1000),igarray(1000)
      
      dimension karrp(800),kbarrp(800),karra(2,800),kbarra(2,800)
      dimension ilena(2),ilena2(2),ibdf(1600)
      dimension iarr(1600),iadf(1600)
      ilenp=ilen
      ilenp2=ilen2
      
      do i=1,ilen
      karrp(i)=karr(i)
      end do
      do i=1,ilen2
      kbarrp(i)=kbarr(i)
      end do
      if (ilenp.gt.ilenp2)goto 1
      ilenas=ilenp2/2
!      print *,'p1 ilenas',ilenas
      if (mod(ilenp2,2).eq.0)goto 2
      ilenas=ilenas+1
      goto 2
1     ilenas=ilenp/2
      if (mod(ilenp,2).eq.0)goto 2
      ilenas=ilenas+1

2    a=a      
     if (ilenp.lt.ilenas+1)goto 210
      ilena(1)=ilenp-ilenas
      
      do jf=1,ilena(1)
      karra(1,jf)=karr(jf)
      end do
      
      
      
211   ired=ilena(1)+1 
      ijk=1
      do jf=ired,ired+ilenas-1
      if (karr(jf).eq.0)goto 204
      ilena(2)=ilenas+1-ijk
      do jk=1,ilena(2)
      karra(2,jk)=karr(jf-1+jk)
      end do
      goto 212
204   ijk=ijk+1
      end do
      ilena(2)=0

      goto 212
210   ilena(1)=0
      
      ilena(2)=ilenp
      do jf=1,ilenp
      karra(2,jf)=karr(jf)
      end do
      
      
212   a=a
!      print *,ilena(1),ilena(2),'ilenas',ilenas
      do i=1,2
      if (ilena(i).eq.0)goto 213
!      print *,'i',i,'karra',(karra(i,jf),jf=1,ilena(i))
213   end do
      
      
      if (ilenp2.lt.ilenas+1)goto 310
      ilena2(1)=ilenp2-ilenas
      do jf=1,ilena2(1)
      kbarra(1,jf)=kbarrp(jf)
      end do
      
311   ired=ilena2(1)+1
      ijk=1
      do jf=ired,ired+ilenas-1
      if (kbarr(jf).eq.0)goto 304
      ilena2(2)=ilenas+1-ijk
      do jk=1,ilena2(2)
      kbarra(2,jk)=kbarr(jf-1+jk)
      end do
      goto 312
304   ijk=ijk+1 
      end do
      ilena2(2)=0
      goto 312
310   ilena2(1)=0
      
      ilena2(2)=ilenp2
      do jf=1,ilenp2
      kbarra(2,jf)=kbarr(jf)
      end do





      
      
312   a=a
      do i=1,2
      if (ilena2(i).eq.0)goto 313
!      print *,'ilen2a',ilena2(1),ilena2(2)
!      print *,'i',i,'kbarra',(kbarra(i,jf),jf=1,ilena2(i))
313   end do
      
      if ((ilena(1).eq.0).or.(ilena2(1).eq.0))goto 535
      do jf=1,ilena(1)
      karr(jf)=karra(1,jf)
      end do
      ilen=ilena(1)
      do jf=1,ilena2(1)
      kbarr(jf)=kbarra(1,jf)
      end do
      ilen2=ilena2(1)
      isub=2
      call mpmul(ilen,ilen2,ilen3)
!      print *,'ilen3',ilen3,'mul1',(kcarr(jf),jf=1,ilen3)
      
      goto 540
535   iarr(1)=0
      iarr(2)=0
      goto 545
540   do jf=1,ilen3
      iarr(jf+2)=kcarr(jf)
      end do
      

!      print *,'point  1'
      iarr(1)=0
      iarr(2)=ilen3
      karr(1)=iarr(1)
      karr(2)=iarr(2)+ilenas
      do jf=1,ilen3
      karr(jf+2)=iarr(jf+2)
      end do
      do jf=ilen3+3,ilen3+2 +ilenas 
      karr(jf)=0
      end do
      kbarr(1)=0 
      kbarr(2)=karr(2)+ilenas
      do jf=3,karr(2)+2
      kbarr(jf)=karr(jf)
      end do
      do jf=karr(2)+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'iarr1',iarr(2),'shift1',(iarr(jf),jf=1,iarr(2)+2)
      

545   if (ilena(1).eq.0)goto 546
      do jf=1,ilena(1)
      karr(jf+2)=karra(1,jf)
      end do
      karr(1)=0
      karr(2)=ilena(1)
      goto 547
546   karr(1)=0
      karr(2)=0
      
547   if (ilena(2).eq.0)goto 548
      do jf=1,ilena(2)
      kbarr(jf+2)=karra(2,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena(2)
      goto 549
548   kbarr(1)=0
      kbarr(2)=0
549   call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=1,kcarr(2)+2
      iadf(jf)=kcarr(jf)
      end do
      if (ilena2(2).eq.0)goto 5495
      do jf=1,ilena2(2)
      karr(jf+2)=kbarra(2,jf)
      end do
      karr(1)=0
      karr(2)=ilena2(2)
      goto 5497
5495  karr(1)=0       
      karr(2)=0
5497  if (ilena2(1).eq.0)goto 5498
      
      do jf=1,ilena2(1)
      kbarr(jf+2)=kbarra(1,jf)
      end do
      kbarr(1)=0
      kbarr(2)=ilena2(1)
      goto 5499
5498  kbarr(1)=0      
      kbarr(2)=0
5499  call mpadd(1)
      if (kcarr(2).eq.0)goto 552
!      print *,'int2',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=3,kcarr(2)+2
      kbarr(jf-2)=kcarr(jf)
      end do 
      ilen2=kcarr(2)
      isgn=kcarr(1)
      do jf=3,iadf(2)+2
      karr(jf-2)=iadf(jf)
      end do
      ilen=iadf(2)
      isub=3
      call mpmul(ilen,ilen2,ilen3)
!      print *,'2ilen3',ilen3,'kcarr',(kcarr(jf),jf=1,ilen3)
      
550   do jf=1,ilen3
      kbarr(jf+2)=kcarr(jf)
      end do
      
      
      
      kbarr(1)=mod(iadf(1)+isgn,2)
!      print *,'kbarr1',kbarr(1)
      kbarr(2)=ilen3+ilenas
      
      
      do jf=ilen3+3,ilen3+2+ilenas
      kbarr(jf)=0
      end do
!      print *,'kbarr',(kbarr(jf),jf=1,kbarr(2)+2)
      
!      print *,'iarr2',iarr(2)
      
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      iarr(jf)=kcarr(jf)
      end do
!      print *,'second iarr',iarr(2),'arr',(iarr(jf),jf=1,iarr(2)+2)
      
552   if (ilena(2).eq.0)goto 553      
      do jf=1,ilena(2)
      karr(jf)=karra(2,jf)
      end do
      ilen=ilena(2)
      if (ilena2(2).eq.0)goto 553


      do jf=1,ilena2(2)
      kbarr(jf)=kbarra(2,jf)
      end do
      ilen2=ilena2(2)
      isub=4 
      call mpmul(ilen,ilen2,ilen3)
!      print *,'3ilen3',ilen3,'kcarr',(kcarr(jf),jf=1,ilen3)
      
560   do jf=1,ilen3
      karr(jf+2)=kcarr(jf)
      kbarr(jf+2)=kcarr(jf)
      end do
      
      karr(1)=0
      kbarr(1)=0
      karr(2)=ilen3
      
      kbarr(2)=ilen3+ilenas
      do jf=ilen3+3,kbarr(2)+2
      kbarr(jf)=0
      end do
      call mpadd(0)
      do jf=1,kcarr(2)+2
      kbarr(jf)=kcarr(jf)
      end do
!      print *,'last',(kcarr(jf),jf=1,kcarr(2)+2)
      
!      print *,'iarrlen',iarr(2)
      do jf=1,iarr(2)+2
      karr(jf)=iarr(jf)
      end do
      call mpadd(0)
      ilen3=kcarr(2)
!      print *,'fir ilen3',ilen3,'arr',(kcarr(jf),jf=1,kcarr(2)+2)
      
      do jf=3,ilen3 +2
      kcarr(jf-2)=kcarr(jf)
      
      end do
!      print *,'kcarrml',(kcarr(jf),jf=1,ilen3)
      
      goto 110
553   ilen3=iarr(2)
      do jf=1,ilen3
      kcarr(jf)=iarr(jf+2)
      end do
      goto 110
      stop
      
      jdf=ilenp-(i-1)*ilenas
      

!      if ((ilen.gt.80).or.(ilen2.gt.80))goto 9
      goto 10
9     a=a
10    do i=1,ilen+ilen2
      kcarr(i) =0
      end do
!      print *,'two',ilen,ilen2
      do i =1,ilen
      do j=1,ilen2
      itemp =karr(i) *kbarr(j)
      
      itemp2=int(itemp/10000)
      irem1 =itemp-itemp2 *10000
      
      
      kcarr(i+j)=kcarr(i+j) +irem1
      kcarr(i+j-1)=kcarr(i+j-1) +itemp2
      
      if (kcarr(i+j).lt.10000)goto 20
      kcarr(i+j)=kcarr(i+j)-10000
      kcarr(i+j-1) =kcarr(i+j-1)+1
20    do k =1,i+j-2
      if (kcarr(i+j-k).lt.10000)goto 22
      kcarr(i+j-k)=kcarr(i+j-k) -10000
      kcarr(i+j-k-1)=kcarr(i+j-k-1)+1
      end do
22    end do
      end do
      ilen3 =ilen +ilen2
!      print *,'four ilen3',ilen3
      if(kcarr(1).ne.0)goto 100
      ilen3 =ilen3-1
!      print *,'three ilen3',ilen3
      do i=1,ilen3
      kcarr(i)=kcarr(i+1)
      end do

100   if (isub.eq.2)goto 540
      if (isub.eq.3)goto 550
      if (isub.eq.4)goto 560
110   ilen=ilenp
      ilen2=ilenp2
      do jf=1,ilen
      karr(jf)=karrp(jf)
      end do
      do jf=1,ilen2
      kbarr(jf)=kbarrp(jf)
      end do
      return
      end


      
