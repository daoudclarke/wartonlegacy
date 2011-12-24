      program bbw4
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
      end
