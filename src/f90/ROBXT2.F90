       program robxt2
!      computation of kernel of matrix quad sieve 60 digit nos.        
!      for use by quadratic sieve program and MPQS       
       
       
       dimension indic(6000),iqi(201,5300),iqi2(201,5300),iv(5600)
       dimension ic(6000),jpiv(6000),kct(201)
       dimension nar(6000),littr(20)
       dimension iharr(200),larr(200,20)
       open(unit=1,file='rompdatc',access='direct',form=&
       'formatted',recl=5391,status='old')
       open(unit=2,file='vernel22',access='direct',form=&
       'formatted',recl=28000,status='old')
       read (1,1,rec=20)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,5300)
       print *,kkb,(nar(jf),jf=1,300)
       do i=1,5300
       nar(i)=0
       end do
       kkb=0
       ihitn=0
       do i=1,20
       littr(i)=0
       end do












       do i=5232,5600
       write (1,1,rec=i)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,5300)
       end do
       
       
       
       
       mm=5600
       nn=5300
       do i=1,mm
       indic(i)=0
       end do
       do i=1,nn
       ic(i)=0
       end do
       ir=0
       
       do nnn=1,27
       ibeg=(nnn-1)*200+1
       iend =ibeg+199
       nk1=1
       do k=ibeg,iend
       print *,'k=',k,'ir=',ir
       read(1,1,rec=k)kkb,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       iqi(nk1,1)=kkb
       do jk=1,nn-1
       iqi(nk1,jk+1)=nar(jk)
       end do
       nk1=nk1+1
       end do
       do jf=1,200
       kct(jf)=0
       end do
       do k=1,200
       do j=1,nn
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40
       
       
       
       
       
1      format(i1,i10,20i4,5300i1)       
       
       
34     jpiv(k)=j    
       
              
       

       do jj=1,nn
       if(jj.eq.j)goto 38
       if (iqi(k,jj).eq.0)goto 38
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       do i=k,200
       itt = iqi(i,jj) + iqi(i,j) 
       iqi(i,jj) =mod(itt,2)
       if(iqi(i,jj).ge.0)goto 36
       iqi(i,jj)=iqi(i,jj) +2

36     end do
38     end do  
       knd=(nnn-1)*200+k
       ic(j) =knd   
       indic(knd) =j
       goto 50

40     ir =ir +1
       print *,'k',k,(iqi(k,jzz),jzz=1,500)
       print *,'ic',(ic(jzz),jzz=1,100)
       
       
       do is =1,mm
       iv(is) = 0
       end do
       knd=(nnn-1)*200+k
       iv(knd) =1
       do is =1,mm

       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is) =iqi(k,jv)
42     end do
       print *,'vector found','ir',ir,'ivs',(iv(nq),nq=1,mm)
       write(2,2,rec=ir)(iv(nq),nq=1,mm)
2      format(5600i5)     
       
       if (ir.lt.50)goto 50
       close(unit=1)
       close(unit=2)
       print *, 'ir=',ir
       stop
50     end do
       if (nnn.ne.1)goto 501
       
       
       
501    if (nnn.eq.28)goto 100
       
       do nnn2=nnn+1,28
       print *,'nnn2=',nnn2
       
       ibeg=(nnn2-1)*200+1
       iend=ibeg+199
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
       
       do k=1,200
       
       jnd=jpiv(k)
       do jf=1,kct(k)
       jnd2=iqi2(k,jf)
       
       do k2=1,200
       itt=iqi(k2,jnd2)+iqi(k2,jnd)
       iqi(k2,jnd2)=mod(itt,2)
       end do
       end do
       
       

601    end do
       
       
       
       nnd=(nnn2-1)*200
       do k=1,200
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

