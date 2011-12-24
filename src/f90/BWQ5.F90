       program bwq5
!      faster program for square root mod p       
       common iprod(200)
       dimension ipr(65000)
       open (unit=3,file='recl.data',access='direct',form=&
       'formatted',recl=390000,status='old')
       read (3,5,rec=1)(ipr(i),i=1,65000)
5      format(65000i6)       
       
       
       
       
       print *,'number?'
       read *,iaa
       print *,'mod'
       read *,ip
       icong=0
       
       
       
       
       
       
       ix = 0
       
       iaa =mod(iaa,ip)
       
       if(iaa.gt.0)goto 10
       if(iaa.eq.0)goto 200
       iaa = iaa + ip
       
10     iprecod =ip -1
       i=0
26     itemp = int(iprecod/2)
       irem1 =iprecod -itemp*2
       if(itemp.eq.0)goto 46
       if(irem1.gt.0)goto 40
       iprecod = itemp
       i =i+1
       if(i.lt.200)goto 26
40     iq =iprecod
       ie = i
       goto 48
46     iq =1
       ie = i
48     i =1
       print *,'qs',iq,'ie',ie
       
       n = 1
52     n = n*607
54     itemp=int(n/1000)
56     irem1 =n-1000 *itemp
       n = irem1
       id = n
       call sub400(id,ip,k)
       if(k.eq.-1)goto 68
       i =i +1
       if(i.lt.1000)goto 52
68     ipn =iq
       print *,'n',n,'i',i
       stop
       iaas = n
       call sub516(ibprod,iaas,ipn,ip)
       iz = ibprod
       print *,'zz',iz
       
       iy = iz
       ir = ie
       ipn = (iq-1)/2
       iaas = iaa
       call sub516(ibprod,iaas,ipn,ip)
       ix = ibprod
       print *,ix
       itemp=mod((iaa *ix),ip)
       
       irem1=mod((itemp*ix),ip)
       
       ib =irem1
       print *,'bbb',ib
       
       itemp = iaa *ix
       ix =mod(itemp,ip)
       
100    irem1 =mod(ib,ip)
       
       if(irem1.eq.1)goto 200
       i=1
       m=1
110    ipow =2**m
       goto 700
112    irem2 =ibprod
       print *,'irem2',irem2
       if(irem2.eq.1)goto 130
       m =m+1
       goto 110
130    if(m.eq.ir)goto 180
       ipow =2**(ir-m-1)
       goto 800
134    it = ibprod
       iy =it*it
       itemp =int(iy/ip)
       irem1 =iy -itemp *ip
       iy =irem1
       itemp =int(m/ip)
       irem1 =m -ip *itemp
       m =irem1
       ir = irem1
       itemp =ix * it
       itemp2 =int(itemp/ip)
       irem1 =itemp -itemp2 *ip
       ix = irem1
       itemp =ib *iy
       itemp2 =int(itemp/ip)
       irem1 =itemp -itemp2 *ip
       ib = irem1
       goto 100
180    print *,'no square root exists'
       
       goto 900
200    print *,'square root=',ix
       icong=icong+1
       if (icong.eq.6500)goto 1000
       print *,'icong=',icong
       
       goto 900
700    ipn =ipow
       iaas =ib
       call sub516(ibprod,iaas,ipn,ip)
       goto 112
800    ipn =ipow
       iaas=iy
       call sub516(ibprod,iaas,ipn,ip)
       goto 134
900    goto 1000       
       print *,'loop too short'
       stop
1000   end       
       subroutine sub400(id,ip,k)
       ide =int(id/2)
       ipe= int(ip/2)
       ide2 =id -ide *2
       ipe2 =ip -ipe*2
       if(ide2.eq.0)goto 414
       goto 416
414    if(ipe2.eq.0)goto 512
416    iv = 0
       ipe = ip
       ii = 0
419    print *,'pefirst',ipe
       ipe2 =int(ipe/2)
       ipe3 =ipe -ipe2 *2
       if(ipe3.eq.1)goto 432
       ipe =ipe2
       iv = iv+1
       ii =ii +1
       if(ii.lt.50)goto 419
432    ive = int(iv/2)
       ive2 =iv -ive*2
       if(ive2.eq.0)goto 450
       iae =(id **2 -1)/8
       iae1 =int(iae/2)
       iae2 =iae -iae1 *2
       iae3 =iae2 +2
       k=(-1)**iae3
       goto 451
450    k =1
451    ide = id
452    if(ide.eq.0)goto 510 
       ive = 0
       ii = 0
456    ide1 =int(ide/2)
       ide2 =ide -ide1 *2
       if(ide2.eq.1)goto 470
       iv = iv+1
       ide =ide1
       ii =ii+1
       if(ii.lt.50)goto 456
470    ive =int(iv/2)
       ive2=iv -ive *2
       if(ive2.eq.0)goto 486
       iae=(ipe **2 -1)/8
       print *,'newipe',ipe,'iae',iae
       iae2 =int(iae/2)
       iae3 =iae -iae2 *2
       iae3 = iae3 +2
       
       
       
       
       k = (-1)**iae3*k
486    iae2 =((ide -1)*(ipe-1))/4       
       print *,'newiae2',iae2,'ide',ide,'ipe',ipe
       iae3 =int(iae2/2)
       iae4 =iae2-iae3 *2
       iae4 =iae4 +2
       k =(-1) **iae4 *k
       ir=abs(ide)
       print *,'r',ir,ipe
       itemp=int(ipe/ir)
       ide =ipe -itemp*ir
       ipe = ir
       print *,'pe',ipe,ide
       goto 452
510    if(ipe.eq.1)goto 513
512    k =0
513    print *,k,ipe
       return
       end
       
       subroutine sub516(ibprod,iaas,ipn,ip)
       common iprod(200)
       iy=iaas
       n=ipn
       ipow=2
14     if (ipow.gt.n)goto 20
       ipow=ipow*2
       goto 14
20     ie=int(ipow/2)
       n1=n
       n1=n1-ie
26     if (ie.eq.1)goto 100
       ie =int(ie/2)
       iy=iy*iy
       iy=mod(iy,ip)
       if (n1.lt.ie)goto 26
       n1=n1-ie
       iy=iy*iaas
       iy=mod(iy,ip)
       goto 26
100    ibprod=iy




620    return
       end
