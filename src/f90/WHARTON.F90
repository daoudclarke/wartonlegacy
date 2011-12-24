      dimension varb(450000)
      varb(2000) = 1
      sum=0
      do 10 i=1,10
      sum=sum+i
10    continue
      print *, sum
      a = 53
      call mpabrt
      


      end



