       program bwbig
       common ibarray(100),isarray(100),igarray(100),inv(20)
       common mult1(100),mult2(100),mult3(200),mult4(200),mult5(100)
       common mat1(20,20),irhs1(20)
       common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
       common isol(20)
       dimension mpol(6),ipol(6)
       dimension ip(20)
       idd =3
       ipold =3
       idegm =2
       idegn =2
       ipd =13
       ipol(1) =0
       ipol(2) =0
       ipol(3) =-2
       ii =1
       mult1(1) =3
       mult1(2) =5
       mult1(3) =8
       do i =1,idd
       mat2(i,idd-1) = mult1(i)
       end do
       do i =1,idd
       mult2(i) =mult1(i)
       end do
10     call multy(idegm,idegn,ipd)
       idegg=idegm +idegn
       ii = ii +1
       if (ii.lt.idd )goto 12
       goto 14

12     do j =1,idegg-idegm
       mul=mult3(j)
       do i =1,idd-1
       
       mult3(j) =0
       mult3(i+j)=mult3(i+j) +mul *ipol(i+1) *(-1)
       end do
       end do
       do i = 1,idd
       
       mat2(i,idd -ii)=mult3(i+idegm)
       
       
       
       mult1(i)=mult3(i+idegm)
       
       end do
       goto 10
14     do i =1,idd -1
       mat2(i,idd)=0
       end do
       print *,'idegg',idegg,mult3(1),mult3(2),mult3(3),mult3(4),&
       mult3(5)
       
       mat2(i,idd) =1
       do i =1,ipold
       mpol(i) = -1 * ipol(i)
       end do
       do j = 1,idegg-idegm
       mul =mult3(j)
       print *,'mul',mul,mpol(1),mpol(2)
       do i =1,idd-1
       mult3(j) =0
       mult3(i+j) =mult3(i+j) +mul *mpol(i+1)
       end do
       end do
       do i =1,idd
       irhs2(i) =mult3(i+idegg-idegm)
       end do
       do i =1,idd
       do j =1,idd
       mat2(i,j) =mod(mat2(i,j),ipd)
       if (mat2(i,j).ge.0)goto 20
       mat2(i,j) =mat2(i,j) +ipd
20     end do       
       irhs2(i)=mod(irhs2(i),ipd)
       if (irhs2(i).ge.0)goto 22
       irhs2(i) =irhs2(i) +ipd
22     end do
       print *,mat2(1,1),mat2(1,2),mat2(1,3),irhs2(1)
       print *,mat2(2,1),mat2(2,2),mat2(2,3),irhs2(2)
       print *,mat2(3,1),mat2(3,2),mat2(3,3),irhs2(3)
       
       call subbw3(idd,ipd)
       do i =1,20
       ip(I) =0
       end do
       ip(1)=isol(1)
       do i =2,idd+1
       ip(2*i -1) =isol(i)
       end do
       stop
       
       end
       
       
       
       
       subroutine subbw2
        
       
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
       return
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
       subroutine subbw6
       print *,'modulus?'
       read *,ia
       print *,'number?'
       read *,ib
       
       ib =mod(ib,ia)
       
       print *,ib
       
       if(ib.ge.0)goto 10
       ib =ia + ib
10     iu =1
       print *,'iaib',ia,ib
       id = ia
       if(ib.eq.0)goto 888
       iv1=0
       iv3 =ib
815    if(iv3.eq.0)goto 830
       iqq = int(id/iv3)
       it3 =id -iqq*iv3
       it1 =iu -iqq*iv1
       iu =iv1
       id = iv3
       iv1 = it1
       iv3 = it3
       goto 815
830    iv =(id -ia *iu)/ib
       if(iu.le.0)goto 870
       iv =iv *(1-ia)
       iv =mod(iv,ia)
       if(iv.ge.0)goto 870
       iv=iv+ia
870    print *,'gcd=',id
       if (id.gt.1)goto 890
       print *,'inverse=',iv
       goto 890
888    iv =ib
890    return
       end
       subroutine subbw4
       dimension ibarray(100),isarray(100),irarray(100),inv(20)
      dimension iqt(100),iwb(100),iws(100)
      ipd =13
      inv(1) =1
      inv(2) =7
      inv(3) =9
      inv(4) =10
      inv(5) =8
      inv(6) =11
      inv(7) =2
      inv(8) =5
      inv(9) =3
      inv(10) =4
      inv(11) =6
      inv(12) =12
      ibarray(1) =1
      ibarray(2) =3
      ibarray(3)=3
      ibarray(4)=7
      idegb =3
      isarray(1) =1
      isarray(2) =2
      isarray(3) =4
      idegs=2


      do i =1,idegs+1
      iws(i) =isarray(i)
      end do
      do i=1,idegb+1
      iwb(i) =ibarray(i)
      end do
      iwsd = idegs
      iwbd = idegb
      loopl = iwbd -iwsd +1
      mul = inv(iws(1))
      do i =1,loopl
      iqt(i) =iwb(i) *mul
      iwb(i) =0
      do j =2,iwsd + 1
      iwb(i+j-1) =iwb(i+j-1) -iws(j) *iqt(i)
      iwb(i+j-1) =mod(iwb(i+j-1),ipd)
      if(iwb(i+j-1).ge.0)goto 12
      iwb(i+j-1) =iwb(i+j-1) +ipd
12    end do
      end do
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
      idegr =0
110   print *,'quos',iqt(1),iqt(2),iqt(3)
      print *,'rem',irarray(1),irarray(2)
      return
      end
      subroutine subbw3(idd,ipd)
      common ibarray(100),isarray(100),igarray(100),inv(20)
      common mult1(100),mult2(100),mult3(200),mult4(200),mult5(100)
      common mat1(20,20),irhs1(20)
      common mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      common isol(20)
      
      idd2 =idd 
      inv(1) =1
      inv(2) =7
      inv(3) =9
      inv(4) =10
      inv(5) =8
      inv(6) =11
      inv(7) =2
      inv(8)=5
      inv(9) =3
      inv(10) =4
      inv(11) =6
      inv(12) =12
      do i = 1,idd2
      irhs1(i) =irhs2(i)
      end do
      do i =1,idd2
      do j = 1,idd2
      mat1(i,j) =mat2(i,j)
      end do
      end do
      do j = 1,idd2
      jmark(j) =0
      markr(j) =0
      jx(j) =0
      end do
      do j =1,idd2
      do i =1,idd2
      
      if(mat1(i,j) .eq.0)goto 137 
      if(markr(i).eq.1)goto 137
      jmark(j) =i
      markr(i) =1
      print *,'firstmarksj',jmark(1),jmark(2),jmark(3),jmark(4)
      print *,'firstmarksr',markr(1),markr(2),markr(3),markr(4)
      mul = inv(mat1(i,j))
      do j1 = j,idd2
      itt = mat1(i,j1) *mul
      mat1(i,j1) = mod(itt,ipd)
      if (mat1(i,j1).ge.0)goto 10
      mat1(i,j1)=mat1(i,j1) +ipd
10    print *,'mully',mul,i,j,j1
      end do
      irhs1(i) =irhs1(i) *mul
      itt =irhs1(i)
      irhs1(i) = mod(itt,ipd)
      if(irhs1(i).ge.0)goto 12
      irhs1(i) =irhs1(i) +ipd
      if(i.eq.idd2)goto 150
12    do ik = i+1,idd2
      if(markr(ik).eq.1)goto 134
      if(mat1(ik,j).eq.0)goto 134
      mul = mat1(ik,j)
      do jj =j,idd2
      
      
      itt = mat1(ik,jj) -mat1(i,jj) *mul
      mat1(ik,jj) = mod(itt,ipd)
      if(mat1(ik,jj).ge.0)goto 132
      mat1(ik,jj) =mat1(ik,jj) +ipd
      print *,'mulijikjj',mul,i,j,ik,jj

132   end do  
      irhs1(ik) =irhs1(ik) - mul * irhs1(i)
      itt = irhs1(ik)
      irhs1(ik) = mod(itt,ipd)
      if(irhs1(ik).ge.0)goto 134
      irhs1(ik) = irhs1(ik) + ipd
134   end do
      goto 150
137   end do

      
      jx(j) =1 
      
      
      
      
      
      



150   end do
      do i =1,idd2
      print *,'matrix',mat1(i,1),mat1(i,2),mat1(i,3),mat1(i,4),&
      irhs1(i)
      print *,'marks',jmark(i),markr(i)
      end do
      
      do i = 1,idd2
      print *,jx(i)
      end do
      
      do i =1,idd2
      if (jx(i).eq.0)goto 152 
      print *,'matrix nonsingular'
      goto 200
152   end do
      
      isol(idd2) =irhs1(idd2)
      do j =1,idd
      ind =idd2 -j+1
      ind2 =ind-1
      im=jmark(ind2)
      isol(ind2)=irhs1(im)
      do jj =ind,idd2
      isol(ind2)=isol(ind2)-isol(jj) *mat1(im,jj)
      end do
      end do
      do i =1,idd2
      itt = isol(i)
      isol(i) =mod(itt,ipd)
      if(isol(i).ge.0)goto 160
      isol(i) =isol(i) +ipd
160   end do
      print *,'solutions',isol(1),isol(2),isol(3)
200   return
      end

       subroutine subbw7
       dimension mat1(20,20),mat2(20,20),irhs1(20),irhs2(20),jmark(20)
       dimension jx(20),markr(20),isol(20)
       n=4
       do i =1,n
       irhs2(i) = 1
       markr(i) =0
       jmark(i) = 0
       jx(i) =0
       isol(i) =0
       do j = 1,n
       mat2(i,j) =0
       end do
       end do
       mat2(2,1) =1
       mat2(2,2) =1
       mat2(3,2)= 1
       mat2(3,3) =1
       mat2(4,3)=1
       mat2(4,1)=1
       mat2(1,1)=1
       mat2(3,1) =1
       mat2(4,1)=1
       irhs2(1)=0
       irhs2(3) =0
       do i =1,n
       irhs1(i) = irhs2(i)
       do j = 1,n
       mat1(i,j) = mat2(i,n+1-j)
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
       do jj =j,n
       if(mat1(ik,jj).eq.0)goto 132
       if(markr(ik).eq.1)goto 132
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
       if(mat1(i,j).eq.1)goto 146
       jx(j) =1
       goto 150
146    jmark(j) =i       
       goto 150
149    jx(j) =1       
150    if(j.lt.n)goto 110
       print *,'rhs',irhs1(1),irhs1(2),irhs1(3),irhs1(4)
       do iix =1,n
       print *,'matrix',mat1(iix,1),mat1(iix,2),mat1(iix,3),mat1(iix,4)
       end do
       print *,'jmarks',jmark(1),jmark(2),jmark(3),jmark(4)
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
       print *,'sols',isol(1),isol(2),isol(3),isol(4)
       do j =1,n
       if(jx(j).eq.1)goto 190
       isol(j) =mod(isol(j),2)
       if(isol(j).ge.0)goto 190
       isol(j) =isol(j) +2
190    end do    
       print *,'solutions',isol(1),isol(2),isol(3),isol(4)
       print *,'jx',jx(1),jx(2),jx(3),jx(4)
       return
       end
       
       subroutine subbw5
       common iprod(200)
       print *,'number?'
       read *,iaa
       print *,'mod'
       read *,ip
       iaa =mod(iaa,ip)
       if(iaa.ge.0)goto 10
       iaa = iaa + ip
       
10     iprecod =ip -1
       i=0
26     itemp = int(iprecod/2)
       irem1 =iprecod -itemp*2
       if(itemp.eq.0)goto 46
       if(irem1.gt.0)goto 40
       iprecod = itemp
       i =i+1
       if(i.lt.200)goto 26
40     iq =iprecod
       ie = i
       goto 48
46     iq =1
       ie = i
48     i =1
       print *,'qs',iq
       n = 1
52     n = n*607
54     itemp=int(n/1000)
56     irem1 =n-1000 *itemp
       n = irem1
       id = n
       call sub400(id,ip,k)
       if(k.eq.-1)goto 68
       i =i +1
       if(i.lt.1000)goto 52
68     ipn =iq
       iaas = n
       call sub516(ibprod,iaas,ipn,ip)
       iz = ibprod
       print *,'zz',iz
       iy = iz
       ir = ie
       ipn = (iq-1)/2
       iaas = iaa
       call sub516(ibprod,iaas,ipn,ip)
       ix = ibprod
       print *,ix
       itemp =iaa *ix *ix
       itemp2=int(itemp/ip)
       irem1 =itemp-itemp2*ip
       ib =irem1
       print *,'bbb',ib
       itemp = iaa *ix
       itemp2 =int(itemp/ip)
       ix = itemp - itemp2 *ip
100    itemp =int(ib/ip)
       irem1 = ib - ip *itemp
       if(irem1.eq.1)goto 200
       i=1
       m=1
110    ipow =2**m
       goto 700
112    irem2 =ibprod
       print *,'irem2',irem2
       if(irem2.eq.1)goto 130
       m =m+1
       goto 110
130    if(m.eq.ir)goto 180
       ipow =2**(ir-m-1)
       goto 800
134    it = ibprod
       iy =it*it
       itemp =int(iy/ip)
       irem1 =iy -itemp *ip
       iy =irem1
       itemp =int(m/ip)
       irem1 =m -ip *itemp
       m =irem1
       ir = irem1
       itemp =ix * it
       itemp2 =int(itemp/ip)
       irem1 =itemp -itemp2 *ip
       ix = irem1
       itemp =ib *iy
       itemp2 =int(itemp/ip)
       irem1 =itemp -itemp2 *ip
       ib = irem1
       goto 100
180    print *,'no square root exists'
       
       goto 900
200    print *,'square root=',ix
       
       goto 900
700    ipn =ipow
       iaas =ib
       call sub516(ibprod,iaas,ipn,ip)
       goto 112
800    ipn =ipow
       iaas=iy
       call sub516(ibprod,iaas,ipn,ip)
       goto 134
900    return       
       end

       subroutine sub400(id,ip,k)
       ide =int(id/2)
       ipe= int(ip/2)
       ide2 =id -ide *2
       ipe2 =ip -ipe*2
       if(ide2.eq.0)goto 414
       goto 416
414    if(ipe2.eq.0)goto 512
416    iv = 0
       ipe = ip
       ii = 0
419    print *,'pefirst',ipe
       ipe2 =int(ipe/2)
       ipe3 =ipe -ipe2 *2
       if(ipe3.eq.1)goto 432
       ipe =ipe2
       iv = iv+1
       ii =ii +1
       if(ii.lt.50)goto 419
432    ive = int(iv/2)
       ive2 =iv -ive*2
       if(ive2.eq.0)goto 450
       iae =(id **2 -1)/8
       iae1 =int(iae/2)
       iae2 =iae -iae1 *2
       iae3 =iae2 +2
       k=(-1)**iae3
       goto 451
450    k =1
451    ide = id
452    if(ide.eq.0)goto 510 
       ive = 0
       ii = 0
456    ide1 =int(ide/2)
       ide2 =ide -ide1 *2
       if(ide2.eq.1)goto 470
       iv = iv+1
       ide =ide1
       ii =ii+1
       if(ii.lt.50)goto 456
470    ive =int(iv/2)
       ive2=iv -ive *2
       if(ive2.eq.0)goto 486
       iae=(ipe **2 -1)/8
       iae2 =int(iae/2)
       iae3 =iae -iae2 *2
       iae3 = iae3 +2
       
       
       
       
       k = (-1)**iae3*k
486    iae2 =((ide -1)*(ipe-1))/4       
       iae3 =int(iae2/2)
       iae4 =iae2-iae3 *2
       iae4 =iae4 +2
       k =(-1) **iae4 *k
       ir=abs(ide)
       print *,'r',ir,ipe
       itemp=int(ipe/ir)
       ide =ipe -itemp*ir
       ipe = ir
       print *,'pe',ipe,ide
       goto 452
510    if(ipe.eq.1)goto 513
512    k =0
513    print *,k,ipe
       return
       end
       
       subroutine sub516(ibprod,iaas,ipn,ip)
       common iprod(200)
       icoun =1
       iq1 = ipn
       i =1
522    iprod(i) =iaas
       i =i +1
       if(i.lt.201)goto 522
       i =1
527    ibigi =1
       j =1
530    if(ibigi.gt.ipn)goto 560
       if(ibigi.eq.ipn)goto 590
       iprep =iprod(i)
       iprod(i) =iprod(i) *iprod(i)
       itemp =int(iprod(i)/ip)
       irem1 =iprod(i) -ip *itemp
       iprod(i)=irem1
       ibigi =ibigi*2
       j = j+1
       if(j.lt.201)goto 530
       print *,'jaze',j,i,iprep
       goto 604
560    ipn = ipn-int(ibigi/2)
       iprod(i) =iprep
       icoun =icoun +1
       i =i +1
       goto 527
590    ibprod =iprod(1)
       if(icoun.eq.1)goto 604
       i =2
594    ibprod =ibprod *iprod(i)
       itemp= int(ibprod/ip)
       irem1 =ibprod -itemp *ip
       ibprod = irem1
       i =i +1
       if(i.lt.icoun+1)goto 594
604    return
       end
