       program bbrecur
!      first program using linear rcurrences 
!      Recurrence algorithm to find kernel of matrix
       common mat2(201,201),isol(200),ipol(500),mpol(500)
       common ibarray(1000),isarray(1000),igarray(1000)
       dimension matt(2,200,200),nmatt(200,200),ibarr(200),iyarr(200)
       dimension jbarr(200),iseqarr(401,200)
       dimension jsol(200),kchek(200),npol(500)
       dimension ibpol(500)
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
       maxpol=0
       ndeg=0
       ncpol=0
       npol(1)=1
       mpol(1)=0
       isw1=0
       mmr=0
!       nn=50
!  nnwas always 40       
       nn=15
!       nn=40
       if (nn.ne.3)goto 80
       npow=10
       goto 81
80     npow=2*nn
81     a=a       
       ibarr(1)=1
       ibarr(2)=1
       ibarr(3)=1
       ibarr(4)=1
       ibarr(5)=0
       ibarr(6)=1
       ibarr(7)=1
       do i=1,nn/2
       ibarr(i)=1
       end do
       do i=nn/2+1,nn
       ibarr(i)=1
       end do
       ibarr(2)=0
       ibarr(3)=0
!       ibarr(12)=0
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
       
       nmatt(1,3)=1
       nmatt(2,3)=1
       nmatt(3,3)=1
       nmatt(4,4)=1
       nmatt(1,6)=1
       nmatt(2,6)=1
       do i=1,nn
       do j=1,nn
       matt(2,i,j)=nmatt(i,j)
       end do
       end do

       do ibbig=1,nn
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
       print *,'ibig',ibig,'iseq',(iseqarr(ibig,jf),jf=1,nn)
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



7      do i=1,npow+1
       ipol(i)=iseqarr(i,ibbig)
       end do
       print *,'iseq',(ipol(jf),jf=1,npow),'ibbig',ibbig
       if (ncpol.ne.5)goto 711
       
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
       do j=1,ndeg+1
       if ((mpol(i).eq.0).or.(npol(j).eq.0))goto 38
       
       ibpol(i+j-1)=ibpol(i+j-1)+1
!       mind=i+j-1
!       print *,'ij',i,j,'ibpol',ibpol(i+j-1),'mulg',mulg,'mind',mind,&
!       'ibpol4',ibpol(4)
38     end do
       end do
       ndeg=mdeg+ndeg
       do i=1,ndeg+1
       npol(i)=mod(ibpol(i),2)
       if (npol(i).ge.0)goto 381
       npol(i)=npol(i)+2
381    end do
       
!       print *,'ndeg',ndeg
       print *,'npol',(npol(jf),jf=1,ndeg+1),'ibbig',ibbig
       
       if (ncpol.ne.5)goto 391
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
!       print *,'iyarr',(iyarr(jf),jf=1,nn)
       

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
13     a=a
!       print *,'jbarr',(jbarr(jk),jk=1,nn),'imark',imark
       do jf=1,nn

       end do
       do jf=1,nn
       if (jbarr(jf).eq.1)goto 9
       end do
       goto 20
9      end do
       goto 10
       print *,'no solution'
       do i=1,nn
       iyarr(i)=0
       end do
       do k=3,ndeg+1
       if (npol(k).eq.0)goto 59
       do i=1,nn
       
       iyarr(i)=iyarr(i)+iseqarr(k-2,i)
       iyarr(i)=mod(iyarr(i),2)
       end do
59     end do
       print *,'iyarr',(iyarr(jf),jf=1,nn),'ndeg',ndeg
       do i=1,nn
       isum=0
       do k=1,nn
       isum=isum+matt(2,i,k)*iyarr(k)
       end do
       isum=mod(isum,2)
       jbarr(i)=isum
       end do
!       if (npol(1).eq.0)goto 13
       print *,'jbarr',(jbarr(jk),jk=1,nn),'imark',imark
       stop
       do jf=1,nn
       jbarr(jf)=jbarr(jf)+ibarr(jf)
       jbarr(jf)=mod(jbarr(jf),2)
       end do
       a=a
       print *,'jbarr',(jbarr(jk),jk=1,nn),'imark',imark
       stop
       do jf=1,nn
       jsol(jf)=iyarr(jf)
       end do
       if (npol(1).eq.0)goto 10
       do i=2,nn
       isum=0
       do k=1,nn
       isum=isum+nmatt(i,k)*jsol(k)
       end do
       kchek(i)=mod(isum,2)
       end do
       print *,'kchek',(kchek(jf),jf=2,nn)
       stop
10     print *,'matrix singular','ibbig',ibbig,'jbarr',(jbarr(jf),& 
       jf=1,nn),'ibarr',(ibarr(jf),jf=1,nn)
       do i=1,nn
       print *,'in mat',(nmatt(i,jf),jf=1,nn)
       end do

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
       print *,'ibarr',(ibarr(jf),jf=1,nn)
       do i=1,nn
       print *,'in mat',(nmatt(i,jf),jf=1,nn)
       end do
       ntpol=ndeg+1
       print *,'ntpol',ntpol,(npol(jf),jf=1,ndeg+1)
       print *,'solutions',(jsol(jf),jf=1,nn)
280    print *,'ncpol',ncpol,'ndeg',ndeg,'maxpol',maxpol,'ndegpr',ndegpr
       end
       subroutine minv(ising,mmr)
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
       print *,'mmr',mmr,'irhs1',(irhs1(jf),jf=1,mmr)
       do i=1,n
       do j = 1,n
       mat1(i,j) = mat2(i,j)
       end do
       print *,'mat1',(mat1(i,jf),jf=1,mmr)
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
       print *,'solutions',(isol(jf),jf=1,n)
!       print *,'jx',(jx(jf),jf=1,n)
       do i=1,n
       isum=0
       do k=1,n
       isum=isum+mat2(i,k)*isol(k)
       end do
       ichek(i)=mod(isum,2)
       end do
       do i=1,n
       print *,'ichek',ichek(i),'irhs1',mat2(i,mmr+1)
       if (ichek(i).eq.mat2(i,mmr+1))goto 300
       print *,'matrix singular'
       ising=1
       goto 302
300    end do
       print *,'solution ok'
       ising=0
302    if (mmr.ne.3)goto 303
       
303    return



       end
       subroutine minpol(npow,mmr)
       common mat2(201,201),isol(200),ipol(500),mpol(500)
       common ibarray(1000),isarray(1000),igarray(1000)
       mgap=1
       if (npow.lt.40)goto 1111
       mgap=npow/20
1111   jincr=mgap
       injj=2
       if (ipol(1).eq.ipol(2))goto 10
! 11     kend=npow/2
!       do kk1=2,kend
!       mmr=kk1
11      mlow=2       
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
       goto 991
79     end do
!       print *,'firspol',(isol(jk),jk=1,mmr)
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
       print *,'mmr',mmr,'injj2',injj2
       
       do i=injj2,mmr
       mpol(i+1-injj2)=mod(isol(i),2)
       end do
       mmr=mmr+1-injj2
91     a=a       
       print *,'minpol found',' mmr',mmr,(mpol(jf),jf=1,mmr)
       
       if (mmr.ne.3)goto 100
       


       goto 100
! 99     end do
99    if (injj.eq.2)goto 661       
      goto 111 
661   injj=1       
      mlow=mmr-mgap 
      mmr=mlow
      jincr=1
      goto 111
991   if (injj.eq.2)goto 111
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



      subroutine bbw4(idegb,idegs,idegg,idegr)
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
