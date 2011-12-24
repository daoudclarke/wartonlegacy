       program bbw6
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
890    end
