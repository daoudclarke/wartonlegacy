       program bbw7
       dimension mat1(20,20),mat2(20,20),irhs1(20),irhs2(20),jmark(20)
       dimension jx(20),markr(20),isol(20),ichek(20)
       n=4
       do i =1,n
       irhs2(i) = 1
       markr(i) =0
       jmark(i) = 0
       jx(i) =0
       isol(i) =0
       do j = 1,n
       mat2(i,j) =0
       end do
       end do
       mat2(4,2)=1
       mat2(2,1) =1
       mat2(2,2) =1
       mat2(3,2)= 1
       mat2(3,3) =1
       mat2(4,3)=1
       mat2(4,1)=1
       mat2(1,1)=1
       mat2(3,1) =1
       mat2(4,1)=1
       mat2(4,4)=1
       irhs2(1)=1
       irhs2(2)=1
       irhs2(3)=1
       irhs2(4)=1
       do i =1,n
       irhs1(i) = irhs2(i)
       do j = 1,n
       mat1(i,j) = mat2(i,j)
       end do
       end do
       j =0
110    j = j+1
       do i =1,n
       if(i.eq.n)goto 140
       if(mat1(i,j).eq.0)goto 137
       if(markr(i).eq.1)goto 137
       jmark(j)=i
       markr(i) =1
       do ik =i+1,n
       if(mat1(ik,j).eq.0)goto 134
       if (markr(ik).eq.1)goto 134
       do jj =j,n
       if(mat1(i,jj).eq.0)goto 132
       
       mat1(ik,jj) =mat1(ik,jj) +mat1(i,jj)
132    end do
       irhs1(ik) = irhs1(ik) +irhs1(i)
134    end do
       do li =1,n
       do lj =j,n
       mat1(li,lj)=mod(mat1(li,lj),2)
       if(mat1(li,lj).ge.0)goto 2000
       mat1(li,lj)=mat1(li,lj) +2
       
2000   end do 
       end do
       goto 150
137    end do
       goto 150
140    if(markr(i).eq.1)goto 149
       print *,'matlast','ij',i,j,mat1(i,j)
       if(mat1(i,j).eq.1)goto 146
       jx(j) =1
       goto 150
146    jmark(j) =i       
       markr(i)=1
       goto 150





149    jx(j) =1       
150    if(j.lt.n)goto 110
       print *,'rhs',irhs1(1),irhs1(2),irhs1(3),irhs1(4)
       do iix =1,n
       print *,'matrix',mat1(iix,1),mat1(iix,2),mat1(iix,3),mat1(iix,4)
       end do
       print *,'jmarks',jmark(1),jmark(2),jmark(3),jmark(4)
       if(jx(n).eq.1)goto 156
       im=jmark(n)
       isol(n) = irhs1(im)
156    do j =1,n-1
       ind =n-j+1
       ind2 =ind-1
       if(jx(ind2).eq.1)goto 180
       im=jmark(ind2)
       isol(ind2)=irhs1(im)
       do jj =ind,n
       isol(ind2) =isol(ind2)+isol(jj)*mat1(im,jj)
       
       end do
180    end do
       print *,'sols',isol(1),isol(2),isol(3),isol(4)
       do j =1,n
       if(jx(j).eq.1)goto 190
       isol(j) =mod(isol(j),2)
       if(isol(j).ge.0)goto 190
       isol(j) =isol(j) +2
190    end do    
       print *,'solutions',isol(1),isol(2),isol(3),isol(4)
       print *,'jx',jx(1),jx(2),jx(3),jx(4)
       do i=1,n
       isum=0 
       do k=1,n
       isum=isum+mat2(i,k)*isol(k)
       end do
       ichek(i)=mod(isum,2)
       end do
       do i=1,n
       print *,'ichek',ichek(i),'irhs2',irhs2(i)
       if (ichek(i).eq.irhs2(i))goto 300
       print *,'matrix singular'
       goto 302
300    end do
       print *,'solution ok' 
       
302    end
