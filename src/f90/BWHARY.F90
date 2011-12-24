       program bb
       double precision b
       integer(kind=selected_int_kind(9))::icount,icount2
       dimension icount(10),icount2(10)
       icount(1) = 1234567890987
       icount2(1) = icount(1) * 2
       print *,icount2(1)
       
       
       end
