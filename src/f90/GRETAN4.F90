       program gretan4
!      this program is a development of greta4 exploiting the existence
!      of zero submatrices where possible and reordering rows
!      double large prime variant
!      third generaslized MPQS program
!      computation of kernel of matrix quad sieve 66 digit nos.        
!      for use by quadratic sieve program and MPQS       
       
       
       dimension indic(40000),iv(40000),iarfin(40000),iv2(40000)
       dimension ic(40000),jpiv(40000),kct(1000)
       dimension nar(40000),littr(20)
       dimension iharr(1000),larr(1000,20),n(60)
       integer,allocatable,dimension (:,:)::iqi,iqi2
       open(unit=7,file='gretpar',access='direct',form=&
       'formatted',recl=321,status='old')
       read(7,6,rec=1)iconr,ispecp,limprm,limprm2,mimprm,lenb2,&
       (n(jf),jf=1,60),irecnn,kkll,lenb1,nizz,lmat1,nmat1,nmat2
       print *,'irecnn',irecnn,'ispecp',ispecp,'lmat1',lmat1,'nmat1',nmat1,&
       'nmat2',nmat2
!       stop
       narc9=kkll+97
6      format(2i6,i8,i8,i3,i6,60i4,3i6,i8,3i6)       
       open(unit=1,file='gretf4',access='direct',form=&
       'formatted',recl=narc9,status='old')
       
       nunit=950000/kkll
       ntem=kkll/nunit
       nzz=(ntem+2)*nunit
       nardiv=ntem+1
       allocate (iqi(nunit,kkll))
       allocate (iqi2(nunit,kkll))
       if (irecnn.gt.nzz)goto 742
       mm=nzz
       goto 743
742    mm=irecnn       
       
743    write(7,6,rec=1)iconr,ispecp,limprm,limprm2,mimprm,lenb2,&
       (n(jf),jf=1,60),irecnn,kkll,lenb1,mm,lmat1,nmat1,nmat2
       print *,'mm',mm,'kkll',kkll
       
       ireccl=7*irecnn
       open (unit=12,file='grettran',access='direct',form=&
       'formatted',recl=ireccl,status='old')
       read (12,120,rec=1)(iarfin(jk),jk=1,irecnn)
!       print *,'iarfin',(iarfin(jk),jk=786,1000)
       
120    format (30000i7)       
       open(unit=2,file='gretf5',access='direct',form=&
       'formatted',recl=nzz,status='old')
       read (1,1,rec=20)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
       print *,kkb,(nar(jf),jf=1,300)
       do i=1,kkll
       nar(i)=0
       end do
       kkb=0
       ihitn=0
       do i=1,20
       littr(i)=0
       end do
       narc=0






       if (irecnn.ge.nzz)goto 74




       do i=irecnn+1,nzz
       write (1,1,rec=i)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,kkll)
       end do
       
       
       

74     nn=kkll
       do i=1,mm
       indic(i)=0
       end do
       do i=1,nn
       ic(i)=0
       end do
       ir=0
       
       do nnn=1,nardiv
       ibeg=(nnn-1)*nunit+1
       iend =ibeg+nunit-1
       nmult=ibeg-1
       nk1=1
       do k=ibeg,iend
       print *,'k=',k,'ir=',ir
       read(1,1,rec=k)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       iqi(nk1,1)=kkb
       do jk=1,nn-1
       iqi(nk1,jk+1)=nar(jk)
       end do
       nk1=nk1+1
       end do
       do jf=1,nunit
       kct(jf)=0
       end do
       do k=1,nunit
       do j=1,nmat1+1
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       if (nmult+k.lt.lmat1)goto 210
       do j=nmat1+2,nmat2

       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
210    do j=nmat2+1,nn       
       if ((iqi(k,j).ne.0).and.(ic(j).eq.0))goto 34
       end do
       goto 40







       
1      format(i1,i6,i10,20i4,50000i1)       
       
       
34     jpiv(k)=j    
       
              
       

       do jj=1,nmat1+1
       if(jj.eq.j)goto 38
       if (iqi(k,jj).eq.0)goto 38
       
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       do i=k,nunit
       itt = iqi(i,jj) + iqi(i,j) 
       iqi(i,jj) =mod(itt,2)
       

36     end do
38     end do  
       if (nmult+k.lt.lmat1)goto 220
       do jj=nmat1+2,nmat2
       if (jj.eq.j)goto 218
       if (iqi(k,jj).eq.0)goto 218
       
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       do i=k,nunit
       itt=iqi(i,jj)+iqi(i,j)
       iqi(i,jj)=mod(itt,2)
       end do
218    end do
220    do jj=nmat2+1,nn
       if (jj.eq.j)goto 228
       if (iqi(k,jj).eq.0)goto 228
       
       kct(k)=kct(k)+1
       ktem=kct(k)
       iqi2(k,ktem)=jj
       do i=k,nunit
       itt=iqi(i,jj)+iqi(i,j)
       iqi(i,jj)=mod(itt,2)
       end do
228    end do




       
       
       knd=(nnn-1)*nunit+k
       ic(j) =knd   
       indic(knd) =j
       goto 50

40     ir =ir +1
       
       
       
       do is =1,mm
       iv(is) = 0
       end do
       knd=(nnn-1)*nunit+k
       iv(knd) =1
       do is =1,mm

       jv = indic(is)
       if (jv.eq.0)goto 42
       iv(is) =iqi(k,jv)
42     end do
       print *,'vector found','ir',ir,'ivs',(iv(nq),nq=1,mm)
       do jdd=1,mm
       if (jdd.gt.irecnn)goto 250
       injdd=iarfin(jdd)
       iv2(injdd)=iv(jdd)
       goto 251
250    iv2(jdd)=iv(jdd)
251    end do       
       write(2,2,rec=ir)(iv2(nq),nq=1,mm)
2      format(50000i1)     
       
       if (ir.lt.80)goto 50
       close(unit=1)
       close(unit=2)
       close(unit=12)
       print *, 'ir=',ir
       stop  

50     end do
       if (nnn.ne.1)goto 501
       
       
       
501    if (nnn.eq.nardiv+1)goto 100
       
       do nnn2=nnn+1,nardiv+1
       print *,'nnn2=',nnn2
       
       ibeg=(nnn2-1)*nunit+1
       iend=ibeg+nunit-1
       nk2=1
       do k=ibeg,iend
       read(1,1,rec=k)kkb,narc,ihitn,(littr(jf),jf=1,20),(nar(jk),jk=1,nn-1)
       
       
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
       
       do k=1,nunit
       if (kct(k).eq.0)goto 601
       jnd=jpiv(k)
       do jf=1,kct(k)
       jnd2=iqi2(k,jf)
       
       do k2=1,nunit
       itt=iqi(k2,jnd2)+iqi(k2,jnd)
       iqi(k2,jnd2)=mod(itt,2)
       end do
       end do
       
       

601    end do
       
       
       
       nnd=(nnn2-1)*nunit
       do k=1,nunit
       kkb=iqi(k,1)
       ihitn=iharr(k)
       do jf=1,20
       littr(jf)=larr(k,jf)
       end do
       
       
       do jf=1,nn-1
       nar(jf)=iqi(k,jf+1)
       end do
       ivan=nnd+k
       write(1,1,rec=ivan)kkb,narc,ihitn,(littr(jf),jf=1,20),&
       (nar(jk),jk=1,nn-1)
       end do
       end do
       end do
100    print *,'ir=',ir
       close(unit=1)
       close(unit=2)
       close(unit=12)


       
       
       end

