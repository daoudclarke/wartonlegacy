       program bb
       double precision b
       integer::icount,icount2
       dimension icount(10),icount2(10)
       icount(1) = 123456789
       icount2(1) = icount(1) * 2
       print *,icount2(1)
       j = dprod(icount1(1),icount2(1))
       
       print *,j
       end
