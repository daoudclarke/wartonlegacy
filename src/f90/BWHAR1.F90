       
      
      
 !  This is the test progrm for the F-90 based MP translation modules.
 
 !  David H. Bailey   96-03-15
 
 use mpmodule
 implicit type (mp_real) (a-h, o-z)
 type (mp_integer) ia, ib, ic
 type (mp_complex) c, d, e
 parameter (n = 25)
 dimension a(n), b(n)
 
 call mpinit
 
 !   Character-to-MP assignment, generic MPREAL function, pi and e.
 
 x = '1.234567890 1234567890 1234567890 D-100'
 ee = exp (mpreal (1.d0))
 call mpwrite (6, x, mppic, ee)
 s = 0.d0
 
 !   Loop with subscripted MP variables.
 do i = 1, n
   a(i) = 2 * i + 1
   b(i) = 2.d0 * a(i) * (a(i) + 1.d0)
   s = s + b(i) ** 2
 enddo
 call mpwrite (6, s)
 
 !   Expression with mixed MPI and MPR entities.
 
 ia = s
 ib = 262144
 s = (s + 327.25d0) * mod (ia, 4 * ib)
call mpwrite (6, s)
 
 !   A complex square root reference.
 
 e = sqrt (mpcmpl (2.d0 * s, s))
 call mpwrite (6, e)
 
 !   External and intrinsic MP function references in expressions.
 
 s = dot (n, a, b)
 t = 2.d0 * sqrt (s) ** 2
call mpwrite (6, s, t)
 
 s = s / 1048576.d0
 t = s + 2.d0 * log (s)
 x = 3 + nint (t) * 5
 call mpwrite (6, s, t, x)
 
 !   A deeply nested expression with function references.
 
 x = (s + (2 * (s - 5) + 3 * (t - 5))) * exp (cos (log (s)))
 call mpwrite (6, x)
 
 !   A special MP subroutine call (computes both cos and sin of S).
 
 call mpcssnf (s, x, y)
 t = 1.d0 - (x ** 2 + y ** 2)
 
 !   IF-THEN-ELSE construct involving MP variables.
 
 if (abs (t) .lt. mpeps) then
   call mpwrite (6, t)
 else
   call mpwrite (6, mpeps)
 endif
 
 stop
 end
 
 function dot (n, a, b)
 
 !   MP function subprogram.
 
 use mpmodule
 type (mp_real) a(n), b(n), dot, s
 
 s = 0.d0
 
 do i = 1, n
   s = s + a(i) * b(i)
 enddo
 
 dot = s
 return
 end
//E*O*F tmpmod90.f//

echo x - tmpmod.out
sed 's=^X ==' > "tmpmod.out" << '//E*O*F tmpmod.out//'
X 10 ^      -100 x  1.2345678901234567890123456789,
 10 ^         0 x  3.14159265358979323846264338327950288419716939937510582097,
       0 x  2.71828182845904523536028747135266249775724709369995957496,
 10 ^         8 x  1.5929408,
 10 ^        14 x  1.52779903171104,
10 ^         7 x  1.79886916621058384791816469746598239276010631272205670324,
 10 ^         6 x  4.24655405854065559427127871323097782191862421687286286356,
 10 ^         6 x  1.8734,
 10 ^         6 x  3.7468,
 10 ^         0 x  1.78661346435546875,
 10 ^         0 x  2.94725728147199086310081336757253060017559122444062761706,
 10 ^         1 x  1.8,
 10 ^         1 x -2.49203075346978107784578072984228366112142524386642299644,
 10 ^         0 x  0.,
//E*O*F tmpmod.out//

exit 0

