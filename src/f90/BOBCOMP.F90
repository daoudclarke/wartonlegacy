      program bobcomp
      double precision complex i2,i3,i4
      double precision x1,x2,x3
      dimension iloc(10)
      i2 =(1234567893,4987654)
      i3 =(56784,567)
      i4 = i2 *i3
      
      print *,'i4',i4
      stop
      x1 =9234567890123456
      print *,'number?'
      mul =12
      read *,i5
      
      do i =1,1
      i6 =i5*mul
      print *,'i6',i6
      end do
      do i =1,1

      iapd =int((1 *100000000 +2 *1000000 + 3 *10000 +4 *100 +5)/&
      (6 *1000000 +7 *10000 + 8 *100 +9))
      print *,'iapd',iapd
      end do
      do i =1,1
      n1 = 9 * 100000000
      print *,n1
      end do
      do i =1,10
      iloc(i) =i
      end do
      do j =1,1500
      do i=1,9
      iloc(i) =iloc(i+1)
      
      end do
      print *,j

      end do
      do i=1,1500
      do j =1,10
      iloc(j) =iloc(j) -1
      end do
      print *,'i',i
      end do
      end
