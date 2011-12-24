      dimension mat1(20,20),irhs1(20),inv(20)
      dimension mat2(20,20),irhs2(20),jmark(20),markr(20),jx(20)
      dimension isol(20)
      
      ipd =13
      idd =2
      idd2 = idd+1
      do i = 1,idd2
      irhs2(i) = 3
      end do
      irhs2(3) = 0
      mat2(1,1) =2
      mat2(1,2) =1
      mat2(1,3) =0
      
      mat2(2,1) =3
      mat2(2,2) =2
      mat2(2,3) =0
      
      mat2(3,1) =4
      mat2(3,2) =4
      mat2(3,3) =1
      
      
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
200   end
