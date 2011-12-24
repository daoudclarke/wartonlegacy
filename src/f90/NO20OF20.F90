X     return
X   end function
X 
X   function mp_qqtoz (qa, qb)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_complex):: mp_qqtoz
X     intent (in):: qa, qb
X     call mpmmpc (qa%mpr, qb%mpr, mp4, mp_qqtoz%mpc)
X     return
X   end function
X 
X   function mp_iitoz (ia, ib)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_complex):: mp_iitoz
X     intent (in):: ia, ib
X     xa = cmplx (ia, ib, kdb)
X     call mpxzc (xa, mp_iitoz%mpc)
X     return
X   end function
X 
X   function mp_rrtoz (ra, rb)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_complex):: mp_rrtoz
X     intent (in):: ra, rb
X     xa = cmplx (ra, rb, kdb)
X     call mpxzc (xa, mp_rrtoz%mpc)
X     return
X   end function
X 
X   function mp_ddtoz (da, db)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_complex):: mp_ddtoz
X     intent (in):: da, db
X     xa = cmplx (da, db, kdb)
X     call mpxzc (xa, mp_ddtoz%mpc)
X     return
X   end function
X 
X   function mp_aatoz (aa, ab)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     character*(*), intent (in):: aa, ab
X     type (mp_complex):: mp_aatoz
X     l = len (aa)
X     do i = 1, l
X       az(i) = aa(i:i)
X     enddo
X     call mpinpc (az, l, mp_aatoz%mpc)
X     l = len (ab)
X     do i = 1, l
X       az(i) = ab(i:i)
X     enddo
X     call mpinpc (az, l, mp_aatoz%mpc(mp41))
X     return
X   end function
X 
X   subroutine mp_cssh (qa, qb, qc)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: qa
X     intent (out):: qb, qc
X     call mpcssh (qa%mpr, mpl02%mpr, qb%mpr, qc%mpr)
X     return
X   end subroutine
X 
X   subroutine mp_cssn (qa, qb, qc)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: qa
X     intent (out):: qb, qc
X     call mpcssn (qa%mpr, mppic%mpr, qb%mpr, qc%mpr)
X     return
X   end subroutine
X 
X   function mp_qtoj (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_qtoj
X     intent (in):: qa
X     call mpeq (qa%mpr, mpt1)
X     call mpinfr (mpt1, mp_qtoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_ztoj (za)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_ztoj
X     intent (in):: za
X     call mpeq (za%mpc, mpt1)
X     call mpinfr (mpt1, mp_ztoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_itoj (ia)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_itoj
X     intent (in):: ia
X     da = ia
X     call mpdmc (da, 0, mp_itoj%mpi)
X     return
X   end function
X 
X   function mp_rtoj (ra)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_rtoj
X     intent (in):: ra
X     da = ra
X     call mpdmc (da, 0, mpt1)
X     call mpinfr (mpt1, mp_rtoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_ctoj (ca)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_ctoj
X     intent (in):: ca
X     da = ca
X     call mpdmc (da, 0, mpt1)
X     call mpinfr (mpt1, mp_ctoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_dtoj (da)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_dtoj
X     intent (in):: da
X     call mpdmc (da, 0, mpt1)
X     call mpinfr (mpt1, mp_dtoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_xtoj (xa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_xtoj
X     intent (in):: xa
X     da = xa
X     call mpdmc (da, 0, mpt1)
X     call mpinfr (mpt1, mp_xtoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_atoj (aa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     character*(*), intent (in):: aa
X     type (mp_integer):: mp_atoj
X     l = len (aa)
X     do i = 1, l
X       az(i) = aa(i:i)
X     enddo
X     call mpinpc (az, l, mpt1)
X     call mpinfr (mpt1, mp_atoj%mpi, mpt2)
X     return
X   end function
X 
X   function mp_nrt (qa, ib)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_nrt
X     intent (in):: qa, ib
X     call mpnrt (qa%mpr, ib, mp_nrt%mpr)
X     return
X   end function
X 
X   function mp_rand ()
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_rand
X     call mprand (mp_rand%mpr)
X     return
X   end function
X 
X   subroutine mp_inpj (iu, j1, j2, j3, j4, j5, j6, j7, j8, j9)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: j1, j2, j3, j4, j5, j6, j7, j8, j9
X     optional:: j2, j3, j4, j5, j6, j7, j8, j9
X     call mpinp (iu, j1%mpi, az)
X     if (present (j2)) call mpinp (iu, j2%mpi, az)
X     if (present (j3)) call mpinp (iu, j3%mpi, az)
X     if (present (j4)) call mpinp (iu, j4%mpi, az)
X     if (present (j5)) call mpinp (iu, j5%mpi, az)
X     if (present (j6)) call mpinp (iu, j6%mpi, az)
X     if (present (j7)) call mpinp (iu, j7%mpi, az)
X     if (present (j8)) call mpinp (iu, j8%mpi, az)
X     if (present (j9)) call mpinp (iu, j9%mpi, az)
X     return
X   end subroutine
X 
X   subroutine mp_inpq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: q1, q2, q3, q4, q5, q6, q7, q8, q9
X     optional:: q2, q3, q4, q5, q6, q7, q8, q9
X     call mpinp (iu, q1%mpr, az)
X     if (present (q2)) call mpinp (iu, q2%mpr, az)
X     if (present (q3)) call mpinp (iu, q3%mpr, az)
X     if (present (q4)) call mpinp (iu, q4%mpr, az)
X     if (present (q5)) call mpinp (iu, q5%mpr, az)
X     if (present (q6)) call mpinp (iu, q6%mpr, az)
X     if (present (q7)) call mpinp (iu, q7%mpr, az)
X     if (present (q8)) call mpinp (iu, q8%mpr, az)
X     if (present (q9)) call mpinp (iu, q9%mpr, az)
X     return
X   end subroutine
X 
X   subroutine mp_inpz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: z1, z2, z3, z4, z5, z6, z7, z8, z9
X     optional:: z2, z3, z4, z5, z6, z7, z8, z9
X     call mpinp (iu, z1%mpc, az)
X     call mpinp (iu, z1%mpc(mp41), az)
X     if (present (z2)) call mpinp (iu, z2%mpc, az)
X     if (present (z2)) call mpinp (iu, z2%mpc(mp41), az)
X     if (present (z3)) call mpinp (iu, z3%mpc, az)
X     if (present (z3)) call mpinp (iu, z3%mpc(mp41), az)
X     if (present (z4)) call mpinp (iu, z4%mpc, az)
X     if (present (z4)) call mpinp (iu, z4%mpc(mp41), az)
X     if (present (z5)) call mpinp (iu, z5%mpc, az)
X     if (present (z5)) call mpinp (iu, z5%mpc(mp41), az)
X     if (present (z6)) call mpinp (iu, z6%mpc, az)
X     if (present (z6)) call mpinp (iu, z6%mpc(mp41), az)
X     if (present (z7)) call mpinp (iu, z7%mpc, az)
X     if (present (z7)) call mpinp (iu, z7%mpc(mp41), az)
X     if (present (z8)) call mpinp (iu, z8%mpc, az)
X     if (present (z8)) call mpinp (iu, z8%mpc(mp41), az)
X     if (present (z9)) call mpinp (iu, z9%mpc, az)
X     if (present (z9)) call mpinp (iu, z9%mpc(mp41), az)
X     return
X   end subroutine
X 
X   function mp_jtoq (ja)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_jtoq
X     intent (in):: ja
X     call mpeq (ja%mpi, mp_jtoq%mpr)
X     return
X   end function
X 
X   function mp_ztoq (za)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_ztoq
X     intent (in):: za
X     call mpeq (za%mpc, mp_ztoq%mpr)
X     return
X   end function
X 
X   function mp_itoq (ia)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_itoq
X     intent (in):: ia
X     da = ia
X     call mpdmc (da, 0, mp_itoq%mpr)
X     return
X   end function
X 
X   function mp_rtoq (ra)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_rtoq
X     intent (in):: ra
X     da = ra
X     call mpdmc (da, 0, mp_rtoq%mpr)
X     return
X   end function
X 
X   function mp_ctoq (ca)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_ctoq
X     intent (in):: ca
X     da = ca
X     call mpdmc (da, 0, mp_ctoq%mpr)
X     return
X   end function
X 
X   function mp_dtoq (da)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_dtoq
X     intent (in):: da
X     call mpdmc (da, 0, mp_dtoq%mpr)
X     return
X   end function
X 
X   function mp_xtoq (xa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_xtoq
X     intent (in):: xa
X     da = xa
X     call mpdmc (da, 0, mp_xtoq%mpr)
X     return
X   end function
X 
X   function mp_atoq (aa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     character*(*), intent (in):: aa
X     type (mp_real):: mp_atoq
X     l = len (aa)
X     do i = 1, l
X       az(i) = aa(i:i)
X     enddo
X     call mpdexc (az, l, mp_atoq%mpr)
X     return
X   end function
X 
X   subroutine mp_outj (iu, j1, j2, j3, j4, j5, j6, j7, j8, j9)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: j1, j2, j3, j4, j5, j6, j7, j8, j9
X     optional:: j2, j3, j4, j5, j6, j7, j8, j9
X     call mpout (iu, j1%mpi, mpoud, az)
X     if (present (j2)) call mpout (iu, j2%mpi, mpoud, az)
X     if (present (j3)) call mpout (iu, j3%mpi, mpoud, az)
X     if (present (j4)) call mpout (iu, j4%mpi, mpoud, az)
X     if (present (j5)) call mpout (iu, j5%mpi, mpoud, az)
X     if (present (j6)) call mpout (iu, j6%mpi, mpoud, az)
X     if (present (j7)) call mpout (iu, j7%mpi, mpoud, az)
X     if (present (j8)) call mpout (iu, j8%mpi, mpoud, az)
X     if (present (j9)) call mpout (iu, j9%mpi, mpoud, az)
X      return
X   end subroutine
X 
X   subroutine mp_outq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: q1, q2, q3, q4, q5, q6, q7, q8, q9
X     optional:: q2, q3, q4, q5, q6, q7, q8, q9
X     call mpout (iu, q1%mpr, mpoud, az)
X     if (present (q2)) call mpout (iu, q2%mpr, mpoud, az)
X     if (present (q3)) call mpout (iu, q3%mpr, mpoud, az)
X     if (present (q4)) call mpout (iu, q4%mpr, mpoud, az)
X     if (present (q5)) call mpout (iu, q5%mpr, mpoud, az)
X     if (present (q6)) call mpout (iu, q6%mpr, mpoud, az)
X     if (present (q7)) call mpout (iu, q7%mpr, mpoud, az)
X     if (present (q8)) call mpout (iu, q8%mpr, mpoud, az)
X     if (present (q9)) call mpout (iu, q9%mpr, mpoud, az)
X      return
X   end subroutine
X 
X   subroutine mp_outz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: z1, z2, z3, z4, z5, z6, z7, z8, z9
X     optional:: z2, z3, z4, z5, z6, z7, z8, z9
X     call mpout (iu, z1%mpc, mpoud, az)
X     call mpout (iu, z1%mpc(mp41), mpoud, az)
X     if (present (z2)) call mpout (iu, z2%mpc, mpoud, az)
X     if (present (z2)) call mpout (iu, z2%mpc(mp41), mpoud, az)
X     if (present (z3)) call mpout (iu, z3%mpc, mpoud, az)
X     if (present (z3)) call mpout (iu, z3%mpc(mp41), mpoud, az)
X     if (present (z4)) call mpout (iu, z4%mpc, mpoud, az)
X     if (present (z4)) call mpout (iu, z4%mpc(mp41), mpoud, az)
X     if (present (z5)) call mpout (iu, z5%mpc, mpoud, az)
X     if (present (z5)) call mpout (iu, z5%mpc(mp41), mpoud, az)
X     if (present (z6)) call mpout (iu, z6%mpc, mpoud, az)
X     if (present (z6)) call mpout (iu, z6%mpc(mp41), mpoud, az)
X     if (present (z7)) call mpout (iu, z7%mpc, mpoud, az)
X     if (present (z7)) call mpout (iu, z7%mpc(mp41), mpoud, az)
X     if (present (z8)) call mpout (iu, z8%mpc, mpoud, az)
X     if (present (z8)) call mpout (iu, z8%mpc(mp41), mpoud, az)
X     if (present (z9)) call mpout (iu, z9%mpc, mpoud, az)
X     if (present (z9)) call mpout (iu, z9%mpc(mp41), mpoud, az)
X      return
X   end subroutine
X 
X   function mp_nint (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_nint
X     intent (in):: qa
X     call mpnint (qa%mpr, mp_nint%mpi)
X     return
X   end function
X 
X   function mp_jtor (ja)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: ja
X     real:: mp_jtor
X     call mpmdc (ja%mpi, da, ia)
X     mp_jtor = da * 2.d0 ** ia
X     return
X   end function
X 
X   function mp_qtor (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: qa
X     real:: mp_qtor
X     call mpmdc (qa%mpr, da, ia)
X     mp_qtor = da * 2.d0 ** ia
X     return
X   end function
X 
X   function mp_ztor (za)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     intent (in):: za
X     real:: mp_ztor
X     call mpmdc (za%mpc, da, ia)
X     mp_ztor = da * 2.d0 ** ia
X     return
X   end function
X 
X   function mp_signj (ja, jb)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_integer):: mp_signj
X     intent (in):: ja, jb
X     call mpeq (ja%mpi, mp_signj%mpi)
X     mp_signj%mpi(1) = sign (mp_signj%mpi(1), jb%mpi(1))
X     return
X   end function
X 
X   function mp_signq (qa, qb)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_signq
X     intent (in):: qa, qb
X     call mpeq (qa%mpr, mp_signq%mpr)
X     mp_signq%mpr(1) = sign (mp_signq%mpr(1), qb%mpr(1))
X     return
X   end function
X 
X   function mp_sin (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_sin
X     intent (in):: qa
X     call mpcssn (qa%mpr, mppic%mpr, mpt1, mp_sin%mpr)
X     return
X   end function
X 
X   function mp_sinz (za)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_complex):: mp_sinz
X     intent (in):: za
X     call mpeq (za%mpc(mp41), mpt2)
X     mpt2(1) = - mpt2(1)
X     call mpexp (mpt2, mpl02%mpr, mpt1)
X     call mpdmc (1.d0, 0, mpt3)
X     call mpdiv (mpt3, mpt1, mpt2)
X     call mpcssn (za%mpc, mppic%mpr, mpt3, mpt4)
X     call mpadd (mpt1, mpt2, mpt5)
X     call mpmuld (mpt5, 0.5d0, 0, mpt6)
X     call mpmul (mpt6, mpt4, mp_sinz%mpc)
X     call mpsub (mpt1, mpt2, mpt5)
X     call mpmuld (mpt5, -0.5d0, 0, mpt6)
X     call mpmul (mpt6, mpt3, mp_sinz%mpc(mp41))
X     return
X   end function
X 
X   function mp_sinh (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_sinh
X     intent (in):: qa
X     call mpcssh (qa%mpr, mpl02%mpr, mpt1, mp_sinh%mpr)
X     return
X   end function
X 
X   function mp_sqrtq (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_sqrtq
X     intent (in):: qa
X     call mpsqrt (qa%mpr, mp_sqrtq%mpr)
X     return
X   end function
X 
X   function mp_sqrtz (za)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_complex):: mp_sqrtz
X     intent (in):: za
X     call mpcsqr (mp4, za%mpc, mp_sqrtz%mpc)
X     return
X   end function
X 
X   function mp_tan (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_tan
X     intent (in):: qa
X     call mpcssn (qa%mpr, mppic%mpr, mpt1, mpt2)
X     call mpdiv (mpt1, mpt2, mp_tan%mpr)
X     return
X   end function
X 
X   function mp_tanh (qa)
X     implicit complex (c), double precision (d), type (mp_integer) (j), &
X       type (mp_real) (q), complex (kdb) (x), type (mp_complex) (z)
X     type (mp_real):: mp_tanh
X     intent (in):: qa
X     call mpcssh (qa%mpr, mpl02%mpr, mpt1, mpt2)
X     call mpdiv (mpt1, mpt2, mp_tanh%mpr)
X     return
X   end function
X 
X end module
X 
X module mpmodule
X use mpfunmod
X use mpintmod
X use mprealmod
X use mpcmpmod
X use mpgenmod
X end module
//E*O*F mpmod90.f//

echo x - tmpmod90.f
sed 's=^X ==' > "tmpmod90.f" << '//E*O*F tmpmod90.f//'
X program tmod
X 
X !  This is the test progrm for the F-90 based MP translation modules.
X 
X !  David H. Bailey   96-03-15
X 
X use mpmodule
X implicit type (mp_real) (a-h, o-z)
X type (mp_integer) ia, ib, ic
X type (mp_complex) c, d, e
X parameter (n = 25)
X dimension a(n), b(n)
X 
X call mpinit
X 
X !   Character-to-MP assignment, generic MPREAL function, pi and e.
X 
X x = '1.234567890 1234567890 1234567890 D-100'
X ee = exp (mpreal (1.d0))
X call mpwrite (6, x, mppic, ee)
X s = 0.d0
X 
X !   Loop with subscripted MP variables.
X 
X do i = 1, n
X   a(i) = 2 * i + 1
X   b(i) = 2.d0 * a(i) * (a(i) + 1.d0)
X   s = s + b(i) ** 2
X enddo
X call mpwrite (6, s)
X 
X !   Expression with mixed MPI and MPR entities.
X 
X ia = s
X ib = 262144
X s = (s + 327.25d0) * mod (ia, 4 * ib)
X call mpwrite (6, s)
X 
X !   A complex square root reference.
X 
X e = sqrt (mpcmpl (2.d0 * s, s))
X call mpwrite (6, e)
X 
X !   External and intrinsic MP function references in expressions.
X 
X s = dot (n, a, b)
X t = 2.d0 * sqrt (s) ** 2
X call mpwrite (6, s, t)
X 
X s = s / 1048576.d0
X t = s + 2.d0 * log (s)
X x = 3 + nint (t) * 5
X call mpwrite (6, s, t, x)
X 
X !   A deeply nested expression with function references.
X 
X x = (s + (2 * (s - 5) + 3 * (t - 5))) * exp (cos (log (s)))
X call mpwrite (6, x)
X 
X !   A special MP subroutine call (computes both cos and sin of S).
X 
X call mpcssnf (s, x, y)
X t = 1.d0 - (x ** 2 + y ** 2)
X 
X !   IF-THEN-ELSE construct involving MP variables.
X 
X if (abs (t) .lt. mpeps) then
X   call mpwrite (6, t)
X else
X   call mpwrite (6, mpeps)
X endif
X 
X stop
X end
X 
X function dot (n, a, b)
X 
X !   MP function subprogram.
X 
X use mpmodule
X type (mp_real) a(n), b(n), dot, s
X 
X s = 0.d0
X 
X do i = 1, n
X   s = s + a(i) * b(i)
X enddo
X 
X dot = s
X return
X end
//E*O*F tmpmod90.f//

echo x - tmpmod.out
sed 's=^X ==' > "tmpmod.out" << '//E*O*F tmpmod.out//'
X 10 ^      -100 x  1.2345678901234567890123456789,
X 10 ^         0 x  3.14159265358979323846264338327950288419716939937510582097,
X 10 ^         0 x  2.71828182845904523536028747135266249775724709369995957496,
X 10 ^         8 x  1.5929408,
X 10 ^        14 x  1.52779903171104,
X 10 ^         7 x  1.79886916621058384791816469746598239276010631272205670324,
X 10 ^         6 x  4.24655405854065559427127871323097782191862421687286286356,
X 10 ^         6 x  1.8734,
X 10 ^         6 x  3.7468,
X 10 ^         0 x  1.78661346435546875,
X 10 ^         0 x  2.94725728147199086310081336757253060017559122444062761706,
X 10 ^         1 x  1.8,
X 10 ^         1 x -2.49203075346978107784578072984228366112142524386642299644,
X 10 ^         0 x  0.,
//E*O*F tmpmod.out//

exit 0

