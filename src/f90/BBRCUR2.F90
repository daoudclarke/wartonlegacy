       program bbrcur2
!      Variation of bbrecur
!      Recurrence algorithm to find kernel of matrix
       common mat2(601,601),isol(1000),ipol(10000),mpol(2000)
       common ibarray(1000),isarray(1000),igarray(1000)
       dimension ibarr(10000),iyarr(10000)
       dimension jbarr(10000)
       dimension jsol(10000),kchek(10000),npol(21000),ivec(30000),iabn(200)
       dimension ibpol(21000),n(60),nar(30000),littr(20),nar2(30000)
       dimension memarr(100000),memind(10000)
       open (unit=7,file='gretpar',access='direct',form=&
       'formatted',recl=303,status='old')
       read (7,107,rec=1)iconr,ispecp,limprm,limprm2,mimprm,lenb2,&
       (n(jf),jf=1,60),irecnn,kkll,lenb1,nizz
       
       narc9=kkll+97
       open (unit=1,file='gretf4',access='direct',form=&
       'formatted',recl=narc9,status='old')
       open (unit=2,file='trandat1',access='sequential')
       open (unit=3,file='trandir',access='direct',form=&
       'formatted',recl=irecnn,status='old')
       open (unit=4,file='gretfs',access='direct',form=&
       'formatted',recl=kkll,status='old')
101    format (i1,i6,i10,20i4,30000i1)       
103    format (30000i1)       
104    format (30000i1)       
107    format (2i6,i8,i8,i3,i6,60i4,3i6,i8)       
!       read (4,104,rec=kkll+1)(ivec(jk),jk=1,200)
!       goto 520


       goto 210
       mmcon=1
       memind(1)=1
       do ii=1,kkll-100
       
       read(2,*)ireco,icur,(iabn(jk),jk=1,icur+1)
       
       memind(mmcon+1)=memind(mmcon)+icur
       if (icur.eq.0)goto 501
       membeg=memind(mmcon)
       do i=1,icur
       memarr(membeg-1+i)=iabn(i)
       end do
       
501    mmcon=mmcon+1
       print *,'row no.',ii,'size so far',memind(mmcon)
       
       end do
       mmcon=mmcon-1
       print *,'size of array=',memind(mmcon+1)
       print *,'indexmem',(memind(jk),jk=1,20)
       print *,'memarray',(memarr(jk),jk=1,200)
       
!       goto 520
       read(1,101,rec=kkll+1)nar(1),narc,ihitn,(littr(jf),jf=1,20),&
       (nar(jk),jk=2,kkll)
       
       write (4,104,rec=1)(nar(jf),jf=1,kkll)
       do k=2,kkll*2
       do i=1,100
       read (3,103,rec=i)(nar2(jf),jf=1,kkll)
       isum=0
       do j=1,kkll
       if ((nar(j).eq.0).or.(nar2(j).eq.0))goto 200
       isum=isum+1
200    end do
       ivec(i)=mod(isum,2)
       end do
       do ii=1,mmcon
       isum=0
       if (memind(ii).eq.memind(ii+1))goto 510
       do jf=memind(ii),memind(ii+1)-1
       indic=memarr(jf)
       isum=isum+nar(indic)
       end do
       ivec(ii+100)=mod(isum,2)
       goto 512
510    ivec(ii+100)=0       
512    end do
!      remove disk access instructions 
!       do i=101,kkll
!       read (2,*)ireco,icur,(iabn(jk),jk=1,icur+1)
!       if (icur.eq.0)goto 201
!       isum=0
!       do jf=1,icur
!       indic=iabn(jf)
!       isum=isum+nar(indic)
!       end do
!       ivec(i)=mod(isum,2)
!       goto 202
!201    ivec(i)=0
!202    end do
       write (4,104,rec=k)(ivec(jf),jf=1,kkll)
!      restoring vector
       do jk=1,kkll
       nar(jk)=ivec(jk)
       end do
!       rewind (unit=2)
       print *,'k',k
       end do
520    close (unit=2)
       close (unit=1)
       close (unit=3)
       close (unit=4)
       close (unit=7)
       stop
210    read (1,101,rec=kkll+1)ibarr(1),narc,ihitn,(littr(jf),jf=1,20),&
       (ibarr(jk),jk=1,kkll)

       print *,'rhs first 100 elements',(ibarr(jk),jk=1,100)

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
       print *,'idegg',idegg
       print *,'gcd',(igarray(jf),jf=1,idegg+1)
       


       call bbw4(idegb,idegs,idegg,idgr)
       maxpol=0
       ndeg=0
       ncpol=0
       npol(1)=1
       mpol(1)=0
       isw1=0
       mmr=0
       nn=kkll

       npow=2*kkll
       
       
       do i=1,nn
       iyarr(i)=0
       jbarr(i)=ibarr(i)
       end do
       
       
       

       do ibbig=1,nn
5      do jf=1,nn       
       if (jbarr(jf).eq.1)goto 7
       end do
       goto 20

       



7      do i=1,npow
       read (4,104,rec=i)(ivec(jf),jf=1,kkll)



       
       ipol(i)=ivec(ibbig)
       end do
!       print *,'ipol',(ipol(jf),jf=1,npow),'ibbig',ibbig
       print *,'ibbig',ibbig
       iend=601
       do i=1,iend
       do j=1,iend
       ind=i+j-1
       mat2(i,j)=ipol(ind)
       end do
       end do
       npow=npow
       call minpol(npow,mmr)
!       print *,'ipol',(ipol(jf),jf=1,npow+1)
      ncpol=ncpol+1 
       mpol(mmr+1)=1
       mdeg=mmr
       print *,'mpol',(mpol(jf),jf=1,mdeg+1),'ncpol',ncpol
!       print *,'npol',(npol(jf),jf=1,ndeg+1)
       if (mdeg.le.maxpol)goto 42
       maxpol=mdeg
       goto 42
!      instructions for gcds of polynomials       
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
       end do
       
!       print *,'ndeg',ndeg
       print *,'npol',(npol(jf),jf=1,100)
       
!       if (ncpol.ne.5)goto 391
!       stop
391    a=a
       do i=1,nn
       iyarr(i)=0
       end do
       do k=2,ndeg+1
       if (npol(k).eq.0)goto 39
       isum=0
       read (4,104,rec=k-1)(ivec(jf),jf=1,kkll)
       do i=1,nn
       isum=isum+ivec(i)
       
       end do
       iyarr=mod(isum,2)
39     end do
       print *,'iyarr',(iyarr(jf),jf=1,nn)
       

2      do i=1,100
       isum=0
       read (3,103,rec=i)(nar2(jf),jf=1,kkll)
       do k=1,kkll
       if ((nar2(k).eq.0).or.(iyarr(k).eq.0))goto 250
       isum=isum+1
       
       
250    end do
       isum=mod(isum,2)
       jbarr(i)=isum
       end do
       rewind (unit=2)
       do i=101,kkll
       read (2,*)ireco,icur,(iabn(jk),jk=1,icur+1)
       if (icur.eq.0)goto 251
       isum=0
       do jf=1,icur
       indic=iabn(jf)
       isum=isum+iyarr(indic)
       end do
       jbarr(i)=mod(isum,2)
       goto 252
251    jbarr(i)=0       
252    end do
       if (npol(1).eq.0)goto 13
       do jf=1,nn
       jbarr(jf)=jbarr(jf)+ibarr(jf)
       jbarr(jf)=mod(jbarr(jf),2)
       end do
13     print *,'jbarr',(jbarr(jk),jk=1,nn),'ncpol',ncpol,'ndeg',ndeg
       
       do jf=1,nn
       if (jbarr(jf).eq.1)goto 9
       
       end do
       goto 20
       
       
9      end do
       print *,'no solution'
       stop
10     print *,'matrix singular' 
       print *,'ysols',(iyarr(jf),jf=1,nn)
       ntpol=ndeg+1
       print *,'ntpol',(npol(jf),jf=1,ndeg+1)
       
       print *,'matrix singular'
       goto 280
       
       
       
20     do jf=1,nn
       jsol(jf)=iyarr(jf)
       end do
       rewind (unit=2)
       if (npol(1).eq.0)goto 10
       do i=1,100
       
       read (3,103,rec=i)(nar2(jf),jf=1,kkll)
       isum=0
       do k=1,nn
       if ((nar2(k).eq.0).or.(jsol(k).eq.0))goto 260
       isum=isum+1
260    end do
       kchek(i)=mod(isum,2)
       end do
       do i=101,nn
       read (2,*)ireco,icur,(iabn(jf),jf=1,icur+1)
       if (icur.eq.0)goto 261
       isum=0
       do jk=1,icur
       indic=iabn(jk)
       isum=isum+jsol(indic)
       end do
       kchek(i)=mod(isum,2)
       goto 262
261    kchek(i)=0
262    end do




       
       do i=1,nn
       if (kchek(i).eq.ibarr(i))goto 50
       goto 10
50     end do


       print *,'matrix nonsingular'
       ntpol=ndeg+1
       print *,'ntpol',ntpol,(npol(jf),jf=1,ndeg+1)
       print *,'solutions',(jsol(jf),jf=1,nn)
280    print *,'ncpol',ncpol,'ndeg',ndeg,'maxpol',maxpol
       
       
       write (4,104,rec=2*kkll+1)(jsol(jf),jf=1,nn)
       close (unit=1)
       close (unit=2)
       close (unit=3)
       close (unit=4)
       close (unit=7)
       
       
       end
       subroutine minv(ising,mmr)
       common mat2(601,601),isol(1000),ipol(10000),mpol(2000)
       common ibarray(1000),isarray(1000),igarray(1000)
       dimension mat1(600,600),irhs1(1000),jmark(1000)
       dimension jx(1000),markr(1000),ichek(1000)
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
       if(mat1(ik,j).eq.0)goto 134
       if (markr(ik).eq.1)goto 134
       do jj =j,n
       if(mat1(i,jj).eq.0)goto 132
       
       mat1(ik,jj) =mat1(ik,jj) +mat1(i,jj)
132    end do
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
150    if(j.lt.n)goto 110
!       print *,'rhs',irhs1(1),irhs1(2),irhs1(3),irhs1(4)
       do iix =1,n
!       print *,'matrix',mat1(iix,1),mat1(iix,2),mat1(iix,3),mat1(iix,4)
       end do
!       print *,'jmarks',jmark(1),jmark(2),jmark(3),jmark(4)
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
!       print *,'sols',isol(1),isol(2),isol(3),isol(4)
       do j =1,n
       if(jx(j).eq.1)goto 190
       isol(j) =mod(isol(j),2)
       if(isol(j).ge.0)goto 190
       isol(j) =isol(j) +2
190    end do    
!       print *,'solutions',isol(1),isol(2),isol(3),isol(4)
!       print *,'jx',jx(1),jx(2),jx(3),jx(4)
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
302    return
       end
       subroutine minpol(npow,mmr)
       common mat2(601,601),isol(1000),ipol(10000),mpol(2000)
       common ibarray(1000),isarray(1000),igarray(1000)
       if (ipol(1).eq.ipol(2))goto 10
11     kend=600
       do kk1=2,kend
       mmr=kk1
       if (mmr.gt.599)goto 87
       call minv(ising,mmr)
       if (ising.eq.1)goto 99
        
       print *,'mmr',mmr
       do kk2=2*mmr+1,npow
       
       irh=ipol(kk2)
       isum=0
       ind=1
       do  kk3=kk2-mmr,kk2-1
       isum=isum+ipol(kk3)*isol(ind)
       ind=ind+1
       end do
       itst=mod(isum,2)
       if (itst.eq.irh)goto 79
       goto 99
79     end do

       do i=1,mmr
       mpol(i)=isol(i)
       end do
       
91     a=a       
!       print *,'minpol found',' mmr',mmr,(mpol(jf),jf=1,mmr)
       goto 100
99     end do
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
       stop

100    return
       end



      subroutine bbw4(idegb,idegs,idegg,idegr)
      common mat2(601,601),isol(1000),ipol(10000),mpol(2000)
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
      common mat2(601,601),isol(1000),ipol(10000),mpol(2000)
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
