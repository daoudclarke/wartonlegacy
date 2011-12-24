      program a
      common lgpr(10)
      dimension karr(1600),kbarr(1600),kcarr(3200)
      dimension kz(60000),loca(1,100),locb(1,100)
      real lgpr
      double precision x,x2,x3,x4 ,x5,x6
      umb=12345678.0
      rel=log10(umb)
      numb=int(rel)
      umb2=numb
      drem=rel-umb2
      drem2=drem*10000.0
      numb2=int(drem2)
      print *,'rel',rel,'numb',numb,'umb2',umb2,'derem',drem
      print *,'drem2',drem2,'numb2',numb2
      stop

      loca(1,1)=7
      locb(1,1)=13
      loca(1,1)=loca(1,1)*locb(1,1)
      print *,'loca',loca(1,1)
      stop
      numb=1234567
      do i=1,60000
      kz(i)=i
      end do
      open (unit=1,file='hazy',access='direct',form=&
      'formatted',recl=540000,status='old')
!      write (1,1111,rec=1)(kz(jk),jk=1,60000)
      read (1,1111,rec=1)(kz(jk),jk=1,60000)
      print *,'kzs',(kz(jk),jk=2001,2100)
      close (unit=1,status='delete')
      stop
1111  format(60000i9)
      
      do i=1,750
      karr(i)=1
      kbarr(i)=1
      end do
      do i=2,4
      karr(i)=0
      kbarr(i)=0
      end do
      ilen=750
      ilen2=750
      do jj=1,1000
      do i=1,ilen+ilen2
      kcarr(i)=0
      end do
      do i=1,ilen
      do j=1,ilen2
      if ((karr(i).eq.0).or.(kbarr(j).eq.0))goto 6
      kcarr(i+j-1)=kcarr(i+j-1)+1
6     end do
      end do
      ilen3=ilen+ilen2-1
      do i=1,ilen3
      if (kcarr(i).eq.0)goto 7
      kcarr(i)=mod(kcarr(i),2)
7     end do
      end do
      print *,'prodb',(kcarr(jf),jf=1,ilen3)
      print *,'ilen3',ilen3,'jj',jj
      stop




      
      do i=1,100000000
      ians=numb*2
      end do
      print *,'ians',ians,'i',i
      stop
      
      
      ix=-5
      iy=-6
      jx=-7
      jy=-8
      k1=mod(ix+iy,3)
      k2=mod(ix+jx,7)
      k3=mod(jx+jy,5)
      print *,'k1',k1,'k2',k2,'k3',k3
      stop
      x=834333222111d-304*2
      x2=10d-560*2
      print *,'x2=',x2
      print *,'x=',x
      x3=cos(6.2831858/10000)
      x4=sin(6.2831858/10000)
      x5=0.0000000000000000000000000001*2
      x6=10.0**12
      print *,'x3=',x3,'x4=',x4,'x5=',x5,'x6=',x6


      i1=31/7
      i2=int(31/7)
      print *,'i1',i1,'i2',i2
      item=int(-31/7)
      print *,'item',item
      
      ib = 100
      print *,ib
      lgpr(1)=8.67*3.2
      lgpr(2)=7.67*3.0
      lgpr(3)=17.3
      
      call sub1
      ig=int(1/2)
      print *,'ig=',ig
      
      end
      subroutine sub1
      common lgpr(10)
      real lgpr
      print *,'lgpr',lgpr(1),lgpr(2),lgpr(3)
      if (lgpr(3).eq.17.3)goto 2
      goto 1
2     print *,'equal'
      if (lgpr(1)+lgpr(2).gt.(lgpr(3)+33.444))goto 3
      goto 1
3     print *,'greater'
1     return
      end
      
