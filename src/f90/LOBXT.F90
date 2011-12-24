       program lobxt
        
       
       
       
       dimension indic(3500),iqi(301,2950),iqi2(301,2950),iv(3500,10)
       dimension ic(3000),jpiv(300),kct(300)
       dimension nar(3000),littr(20)
       dimension iharr(300),larr(300,20)
       open(unit=1,file='compdate',access='direct',form=&
       'formatted',recl=3025,status='old')
       open(unit=2,file='ternel2',access='direct',form=&
       'formatted',recl=15820,status='old')
       read (1,1,rec=301)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,2934)
       print *,kkb,(nar(jf),jf=1,20)
       
       
       mm=3164
       nn=2935
       do i=1,mm
       indic(i)=0
       end do
       do i=1,nn
       ic(i)=0
       end do
       ir=0
       
       do nnn=1,10
       ibeg=(nnn-1)*300+1
       iend =ibeg+299
       nk1=1
       do k=ibeg,iend
       print *,'k=',k
       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       iqi(nk1,1)=kkb
       do jk=1,nn-1
       iqi(nk1,jk+1)=nar(jk)
       end do
       nk1=nk1+1
       end do
       do jf=1,300
       kct(jf)=0
       end do
       do k=1,300
       do j=1,nn
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
       
       
       
       
       
1      format(i1,i10,20i4,2934i1)       
       
       
34     jpiv(k)=j    
       
              
       

       do jj=1,nn
       if(jj.eq.j)goto 38
       if (iqi(k,jj).eq.0)goto 38
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       do i=k,300
       itt = iqi(i,jj) + iqi(i,j) 
       iqi(i,jj) =mod(itt,2)
       if(iqi(i,jj).ge.0)goto 36
       iqi(i,jj)=iqi(i,jj) +2

36     end do
38     end do  
       knd=(nnn-1)*300+k
       ic(j) =knd   
       indic(knd) =j
       goto 50

40     ir =ir +1
       
       
       
       do is =1,mm
       iv(is,ir) = 0
       end do
       knd=(nnn-1)*300+k
       iv(knd,ir) =1
       do is =1,mm

       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is,ir) =iqi(k,jv)
42     end do
       print *,'vector found','ir',ir,'ivs',(iv(nq,ir),nq=1,mm)
       write(2,2,rec=ir)(iv(nq,ir),nq=1,mm)
2      format(3164i5)     
       
       if (ir.lt.10)goto 50
       close(unit=1)
       close(unit=2)
       
       stop
50     end do
       if (nnn.ne.1)goto 501
       
       
       
501    if (nnn.eq.10)goto 100
       
       do nnn2=nnn+1,10
       print *,'nnn2=',nnn2
       ibeg=(nnn2-1)*300+1
       iend=ibeg+299
       nk2=1
       do k=ibeg,iend
       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       
       
       iharr(nk2)=ihitn
       
       
       do jf=1,20
       larr(nk2,jf)=littr(jf)
       end do
       
       
       iqi(nk2,1)=kkb
       do jf=1,nn-1
       iqi(nk2,jf+1)=nar(jf)
       
       end do
       
       nk2=nk2+1 
       end do
       
       do k=1,300
       
       jnd=jpiv(k)
       do jf=1,kct(k)
       jnd2=iqi2(k,jf)
       
       do k2=1,300
       itt=iqi(k2,jnd2)+iqi(k2,jnd)
       iqi(k2,jnd2)=mod(itt,2)
       end do
       end do
       if (k.lt.280)goto 601
       

601    end do
       
       
       
       nnd=(nnn2-1)*300
       do k=1,300
       kkb=iqi(k,1)
       ihitn=iharr(k)
       do jf=1,20
       littr(jf)=larr(k,jf)
       end do
       
       
       do jf=1,nn-1
       nar(jf)=iqi(k,jf+1)
       end do
       ivan=nnd+k
       write(1,1,rec=ivan)kkb,ihitn,(littr(jf),jf=1,20),&
       (nar(jk),jk=1,nn-1)
       end do
       end do
       end do
100    print *,'not enough rows'
       close(unit=1)
       close(unit=2)



       
       
       end

