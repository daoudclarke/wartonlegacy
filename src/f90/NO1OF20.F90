#! /bin/sh
# This is a shell archive.  Remove anything before this line, then
# unpack it by saving it in a file and typing "sh file".  (Files
# unpacked will be owned by you and have default permissions.)
#
# This archive contains:
# readme-f90-ieee mpfun90.f testmp90.f mpmod90.f tmpmod90.f tmpmod.out

echo x - readme-f90-ieee
sed 's=^X ==' > "readme-f90-ieee" << '//E*O*F readme-f90-ieee//'
X README File for the Fortran-90 Multiprecision System
X 
X IEEE version
X Update as of 96-03-20
X 
X The following files are included in the "shar" file accompanying this file:
X 
X Name            Description
X 
X mpfun90.f       MPFUN library file
X testmp90.f      Test program for MPFUN
X mpmod90.f       Translator modules
X tmpmod90.f      Test program for translator modules
X tmpmod.out      Reference output of tmpmod90.f
X 
X These programs should work correctly on any computer system based on
X the IEEE floating-point standard.  For Cray vector systems or others,
X request separate versions of these programs from the author.
X 
X Because Fortran-90 processors are still scarce, I will list only the
X instructions for compiling and testing on an IBM RS6000 workstation.
X For other systems, you will need to modify these instructions
X accordingly.  I would appreciate learning the command sequences that
X work on other systems.
X 
X Command                           Notes
X 
X xlf90 -O3 -c mpfun90.f            There should be no fatal diagnostics.
X 
X xlf90 -O3 testmp90.f mpfun90.o    There should be no fatal diagnostics.
X 
X a.out > testmp.out                Check testmp.out to make sure all tests
X                                     passed.
X 
X xlf90 -O3 -c mpmod90.f            There should be no fatal diagnostics.
X 
X xlf90 tmpmod90.f mpmod90.o mpfun90.o   There should be no fatal diagnostics.
X 
X a.out > test.out                  This should complete normally.
X 
X diff test.out tmpmod90.out        There should be no differences here.
X 
X On an IBM RS6000/590 workstation, and on any other IBM system that
X employs the Power-2 processor, it may be necessary to add the flag
X -qarch=pwr2 to the xlf90 commands.  This flag is not needed if the
X pwr2 option has been "hardwired" in the xlf90 command on your system.
X 
X These codes have been tested quite thoroughly, but a few bugs may
X remain.  If you encounter any, please let me know and I will fix them
X as soon as possible.
X 
X David H. Bailey
X NASA Ames Research Center
X Mail Stop T27A-1
X Moffett Field, CA 94035-1000
X Tel.: 415-604-4410
X Fax:  415-604-3957
X E-mail: dbailey@nas.nasa.gov
X 
//E*O*F readme-f90-ieee//

echo x - mpfun90.f
sed 's=^X ==' > "mpfun90.f" << '//E*O*F mpfun90.f//'
X !*****************************************************************************
X 
X !   MPFUN: A MULTIPLE PRECISION FLOATING POINT COMPUTATION PACKAGE
X 
X !   IEEE Fortran-90 version
X !   Version Date:  96-03-15
X 
X !   Author:
X 
X !      David H. Bailey                 Telephone:   415-604-4410
X !      NASA Ames Research Center       Facsimile:   415-604-3957
X !      Mail Stop TA27-1                Internet:    dbailey@nas.nasa.gov
X !      Moffett Field, CA 94035
X !      USA
X 
X !   Restrictions:
X 
X !   This software has now been approved by NASA for unrestricted distribution.
X !   However, usage of this software is subject to the following:
X 
X !   1. This software is offered without warranty of any kind, either expressed
X !      or implied.  The author would appreciate, however, any reports of bugs
X !      or other difficulties that may be encountered.
X !   2. If modifications or enhancements to this software are made to this
X !      software by others, NASA Ames reserves the right to obtain this enhanced
X !      software at no cost and with no restrictions on its usage.
X !   3. The author and NASA Ames are to be acknowledged in any published paper
X !      based on computations using this software.  Accounts of practical
X !      applications or other benefits resulting from this software are of
X !      particular interest.  Please send a copy of such papers to the author.
X 
X !   Description:
X 
X !   The following information is a brief description of this program.  For
X !   full details and instructions for usage, see the paper "A Portable High
X !   Performance Multiprecision Package", available from the author.
X 
X !   This package of Fortran subroutines performs multiprecision floating point
X !   arithmetic.  If sufficient main memory is available, the maximum precision
X !   level is at least 16 million digits.  The maximum dynamic range is at
X !   least 10^(+-14,000,000).  It employs advanced algorithms, including an
X !   FFT-based multiplication routine and some recently discovered
X !   quadratically convergent algorithms for pi, exp and log.  The package also
X !   features extensive debug and self-checking facilities, so that it can be
X !   used as a rigorous system integrity test.  All of the routines in this
X !   package have been written to facilitate vector, parallel processing and/or
X !   RISC processing.
X 
X !   Most users will not wish to manually write code that directly calls these
X !   routines, since this is a tedious and error-prone process.  To assist
X !   such users, the author has prepared Fortran-90 modules that permit one
X !   to write ordinary Fortran program, yet have the required MPFUN routines
X !   automatically called.  Contact the author for details.
X 
X !   Machine-specific tuning notes may be located by searching for the text
X !   string !> in this program file.  It is highly recommended that these notes
X !   be read before running this package on a specific system.  Certain
X !   vectorizable DO loops that are often not recognized as such by vectorizing
X !   compilers are prefaced with Cray !dir$ ivdep directives.  On other vector
X !   systems these directives should be replaced by the appropriate equivalents.
X !   Double precision should be disabled when compiling this code on Cray
X !   systems (use the -dp flag).
X 
X !   Instructions for compiling and testing this program on various specific
X !   systems are included in the readme file that accompanies this file.
X 
X !*****************************************************************************
X 
X module mpfunmod
X 
X !   This section initializes global parameters and arrays default values.
X !   Note that this Fortran-90 version completely dispenses with common blocks.
X !   The global parameters previously in common block MPCOM1 are now global
X !   parameters in this module and can be accessed in any user program with a
X !   USE MPFUNMOD statement.  Note however the names now start with 'MP'.  The
X !   user also no longer needs to allocate scratch space in common MPCOM3,
X !   MPCOM4 or MPCOM5, since the required space is now allocated dynamically
X !   by each routine as required.
X 
X double precision, private:: bbx, bdx, bx2, rbx, rdx, rx2, rxx
X private kdp, nbt, npr, mcr, nrow, nsp1, nsp2
X double precision, allocatable, private:: dd1, dd2
X parameter (kdp = kind (0.d0))
X complex (kdp), allocatable, private:: dc1, dc2, dc3, uu1, uu2
X !>
X !   The parameters BBX, NBT, NPR and MCR depend on system word size.
X !   On IEEE systems and most other 32 bit systems, set BBX = 4096.D0,
X !   NBT = 24, NPR = 32, and MCR = 7.  On Cray systems, set BBX = 2048.D0,
X !   NBT = 22, NPR = 16, and MCR = 8.  The parameters NROW, NSP1 and NSP2 are
X !   spacing parameters to avoid bank and cache conflicts in the FFT routines.
X !   NROW = 16, NSP1 = 2 and NSP2 = 9 appear to work well on most systems.
X !   Set NROW = 64 on Crays.
X 
X parameter (bbx = 4096.d0, nbt = 24, npr = 32, mcr = 7, &
X   bdx = bbx ** 2, bx2 = bdx ** 2, rbx = 1.d0 / bbx, &
X   rdx = rbx ** 2, rx2 = rdx ** 2, rxx = 16.d0 * rx2, &
X   nrow = 16, nsp1 = 2, nsp2 = 9)
X dimension mpker(72), dc1(:), dc2(:), dc3(:), dd1(:), dd2(:), uu1(:), uu2(:)
X 
X data mpnw, mpidb, mpldb, mpndb, mpier, mpmcr, mpird / 16, 0, 6, 22, 99, mcr, 1/
X data mpker /72 * 2/
X 
X contains
X 
X subroutine dpadd (a, na, b, nb, c, nc)
X 
X !   This adds the DPE numbers (A, NA) and (B, NB) to yield the sum (C, NC).
X 
X implicit double precision (a-h, o-z)
X dimension pt(64)
X save pt
X data pt/ 64 * 0.d0/
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   c = 0.d0
X   nc = 0
X   return
X endif
X 
X !   If this is the first call to DPADD, initialize the PT table.
X 
X if (pt(1) .eq. 0.d0) then
X   pt(1) = 0.5d0
X 
X   do i = 2, 64
X     pt(i) = 0.5d0 * pt(i-1)
X   enddo
X endif
X 
X !   This operation reduces to five cases.
X 
X if (b .eq. 0.d0) then
X   c = a
X   nc = na
X else if (a .eq. 0.d0) then
X   c = b
X   nc = nb
X else if (na .eq. nb) then
X   c = a + b
X   nc = na
X else if (na .gt. nb) then
X   k = na - nb
X   nc = na
X   if (k .gt. 64) then
X     c = a
X   else
X     c = a + b * pt(k)
X   endif
X else
X   k = nb - na
X   nc = nb
X   if (k .gt. 64) then
X     c = b
X   else
X     c = b + a * pt(k)
X   endif
X endif
X if (c .eq. 0.d0) then
X   nc = 0
X   goto 130
X endif
X 
X !   Normalize the result to a decent range if it is not.
X 
X 110  if (abs (c) .ge. bdx) then
X   c = rdx * c
X   nc = nc + nbt
X   goto 110
X endif
X 
X 120  if (abs (c) .lt. 1.d0) then
X   c = bdx * c
X   nc = nc - nbt
X   goto 120
X endif
X 
X 130  return
X end subroutine
X 
X subroutine dpdec (a, na, b, nb)
X 
X !   This converts the DPE number (A, NA) to decimal form, i.e. B * 10^NB,
X !   where |B| is between 1 and 10.
X 
X implicit double precision (a-h, o-z)
X parameter (xlt = 0.3010299956639812d0)
X 
X if (a .ne. 0.d0) then
X   t1 = xlt * na + log10 (abs (a))
X   nb = t1
X   if (t1 .lt. 0.d0) nb = nb - 1
X   b = sign (10.d0 ** (t1 - nb), a)
X else
X   b = 0.d0
X   nb = 0
X endif
X 
X return
X end subroutine
X 
X subroutine dpdiv (a, na, b, nb, c, nc)
X 
X !   This divides the DPE number (A, NA) by (B, NB) to yield the quotient
X !   (C, NC).
X 
X implicit double precision (a-h, o-z)
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   c = 0.d0
X   nc = 0
X   return
X endif
X if (b .eq. 0.d0) then
X   if (mpker(1) .ne. 0) then
X     write (mpldb, 1)
X 1   format ('*** DPDIV: Divisor is zero.')
X     mpier = 1
X     if (mpker(mpier) .eq. 2) call mpabrt
X   endif
X   return
X endif
X 
X !   Divide A by B and subtract exponents, unless A is zero.
X 
X if (a .eq. 0.d0) then
X   c = 0.d0
X   nc = 0
X   goto 120
X else
X   c = a / b
X   nc = na - nb
X endif
X 
X !   Normalize the result to a decent range if it is not.
X 
X 100  if (abs (c) .ge. bdx) then
X   c = rdx * c
X   nc = nc + nbt
X   goto 100
X endif
X 
X 110  if (abs (c) .lt. 1.d0) then
X   c = bdx * c
X   nc = nc - nbt
X   goto 110
X endif
X 
X 120  return
X end subroutine
X 
X subroutine dpmul (a, na, b, nb, c, nc)
X 
X !   This multiplies the DPE number (A, NA) by (B, NB) to yield the product
X !   (C, NC).
X 
X implicit double precision (a-h, o-z)
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   c = 0.d0
X   nc = 0
X   return
X endif
X 
X !   Multiply A by B and add exponents, unless either is zero.
X 
X if (a .eq. 0.d0 .or. b .eq. 0.d0) then
X   c = 0.d0
X   nc = 0
X   goto 120
X else
X   c = a * b
X   nc = na + nb
X endif
X 
X !   Normalize the result to a decent range if it is not.
X 
X 100  if (abs (c) .ge. bdx) then
X   c = rdx * c
X   nc = nc + nbt
X   goto 100
X endif
X 
X 110  if (abs (c) .lt. 1.d0) then
X   c = bdx * c
X   nc = nc - nbt
X   goto 110
X endif
X 
X 120  return
X end subroutine
X 
X subroutine dppwr (a, na, b, nb, c, nc)
X 
X !   This raises the DPE number (A, NA) to the (B, NB) power and places the
X !   result in (C, NC).
X 
X implicit double precision (a-h, o-z)
X parameter (cl2 = 1.4426950408889633d0)
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   c = 0.d0
X   nc = 0
X   return
X endif
X if (a .le. 0.d0) then
X   if (mpker(2) .ne. 0) then
X     write (mpldb, 1)
X 1   format ('*** DPPWR: Argument is less than or equal to zero.')
X     mpier = 2
X     if (mpker(mpier) .eq. 2) call mpabrt
X   endif
X   return
X endif
X 
X if (b .eq. 0.d0) then
X   c = 1.d0
X   nc = 0
X   goto 120
X endif
X 
X if (b .eq. 1.d0 .and. nb .eq. 0) then
X   c = a
X   nc = na
X   goto 120
X endif
X 
X !   Compute the base 2 logarithm of A and multiply by B.
X 
X al = cl2 * log (a) + na
X call dpmul (al, 0, b, nb, t1, n1)
X 
X !   Check for possible overflow or underflow.
X 
X if (n1 .gt. 6) then
X   if (t1 .gt. 0.d0) then
X     if (mpker(3) .ne. 0) then
X       write (mpldb, 2)
X 2     format ('*** DPPWR: Overflow')
X       mpier = 3
X       if (mpker(mpier) .eq. 2) call mpabrt
X     endif
X     return
X   else
X     c = 0.d0
X     nc = 0
X     goto 120
X   endif
X endif
X 
X !   Compute 2 raised to the power B * Log_2 (A).
X 
X t1 = t1 * 2.d0 ** n1
X nc = int (t1)
X c = 2.d0 ** (t1 - nc)
X 
X !   Normalize the result to a decent range if it is not.
X 
X 100  if (abs (c) .ge. bdx) then
X   c = rdx * c
X   nc = nc + nbt
X   goto 100
X endif
X 
X 110  if (abs (c) .lt. 1.d0) then
X   c = bdx * c
X   nc = nc - nbt
X   goto 110
X endif
X 
X 120  return
X end subroutine
X 
X subroutine dpsqrt (a, na, b, nb)
X 
X !   This computes the square root of the DPE number (A, NA) and places the
X !   result in (B, NB).
X 
X implicit double precision (a-h, o-z)
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   b = 0.d0
X   nb = 0
X   return
X endif
X if (a .lt. 0.d0) then
X   if (mpker(4) .ne. 0) then
X     write (mpldb, 1)
X 1   format ('*** DPSQRT: Argument is negative.')
X     mpier = 4
X     if (mpker(mpier) .eq. 2) call mpabrt
X   endif
X   return
X endif
X 
X if (a .eq. 0.d0) then
X   b = 0.d0
X   nb = 0
X   goto 120
X endif
X 
X !   Divide the exponent of A by two and then take the square root of A.  If
X !   NA is not an even number, then we have to multiply A by 10 before taking
X !   the square root.
X 
X nb = na / 2
X if (na .eq. 2 * nb) then
X   b = sqrt (a)
X else
X   b = sqrt (2.d0 * a)
X   if (na .lt. 0) nb = nb - 1
X endif
X 
X !   Normalize the result to a decent range if it is not.
X 
X 100  if (abs (b) .ge. bdx) then
X   b = rdx * b
X   nb = nb + nbt
X   goto 100
X endif
X 
X 110  if (abs (b) .lt. 1.d0) then
X   b = bdx * b
X   nb = nb - nbt
X   goto 110
X endif
X 
X 120  return
X end subroutine
X 
X subroutine dpsub (a, na, b, nb, c, nc)
X 
X !   This subtracts the DPE number (B, NB) from (A, NA) to yield the difference
X !   (C, NC).
X 
X implicit double precision (a-h, o-z)
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   c = 0.d0
X   nc = 0
X   return
X endif
X 
X bb = -b
X call dpadd (a, na, bb, nb, c, nc)
X 
X return
X end subroutine
X 
X subroutine mpabrt
X !>
X !   This routine terminates execution.  Many users will want to replace the
X !   default STOP with a call to a system routine that provides a traceback.
X !   Examples of code that produce traceback are included here (commented out)
X !   for some systems.
X 
X if (mpier .eq. 99) then
X   write (mpldb, 1)
X 1 format ('*** The MPFUN library has not been initialized.  If you are using'/&
X   'the Fortran-90 translation modules, you must insert the following line'/ &
X   'at the start of execution in your main program:'/ ' '/ 'CALL MPINIT')
X else
X   write (mpldb, 2) mpier
X 2 format ('*** MPABRT: Execution terminated, error code =',i4)
X endif
X 
X !   Use this line on Cray systems.
X 
X ! call abort
X 
X !   Use this line plus the C routine TRACBK (available from author) on
X !   Silicon Graphics IRIS systems.
X 
X ! call tracbk
X 
X !   On other systems, merely terminate execution.
X 
X stop
X end subroutine
X 
X subroutine mpadd (a, b, c)
X 
X !   This routine adds MP numbers A and B to yield the MP sum C.  It attempts
X !   to include all significance of A and B in the result, up to the maximum
X !   mantissa length MPNW.  Debug output starts with MPIDB = 9.
X 
X !   Max SP space for C: MPNW + 4 cells.
X 
X double precision d
X dimension a(mpnw+2), b(mpnw+2), c(mpnw+4), d(mpnw+4)
X 
X if (mpier .ne. 0) then
X   if (mpier .eq. 99) call mpabrt
X   c(1) = 0.
X   c(2) = 0.
X   return
X endif
X if (mpidb .ge. 9) then
X   no = min (int (abs (a(1))), mpndb) + 2
X   write (mpldb, 1) (a(i), i = 1, no)
X 1 format ('MPADD I'/(6f12.0))
X   no = min (int (abs (b(1))), mpndb) + 2
X   write (mpldb, 1) (b(i), i = 1, no)
X endif
X 
X ia = sign (1., a(1))
X ib = sign (1., b(1))
X na = min (int (abs (a(1))), mpnw)
X nb = min (int (abs (b(1))), mpnw)
X 
X !   This first IF block checks for zero inputs.
X 
X if (na .eq. 0) then
X 
X !   A is zero -- the result is B.
X 
X   c(1) = sign (nb, ib)
X 
X   do i = 2, nb + 2
X     c(i) = b(i)
X   enddo
X 
X   goto 420
X elseif (nb .eq. 0) then
X 
X !   B is zero -- the result is A.
X 
X   c(1) = sign (na, ia)
X 
X   do i = 2, na + 2
X     c(i) = a(i)
X   enddo
X 
X   goto 420
X endif
X ma = a(2)
X mb = b(2)
X 
X !   This IF block breaks the problem into different branches depending on
X !   the relative sizes of the exponents of A and B.
X 
X if (ma .eq. mb) then
X 
X !   A and B have the same exponent.
X 
X   nm = min (na, nb)
X   nx = max (na, nb)
X   if (ia .eq. ib) then
X 
X !   A and B have the same exponent and sign.
X 
X     d(1) = sign (nx, ia)
X     d(2) = ma
X     d(nx+3) = 0.d0
X     d(nx+4) = 0.d0
X 
X     do i = 3, nm + 2
X       d(i) = dble (a(i)) + dble (b(i))
X     enddo
X 
X     if (na .gt. nb) then
X 
X !   A is longer than B -- include extra words of A in C.
X 
X       do i = nm + 3, na + 2
X         d(i) = a(i)
X       enddo
X 
X     elseif (nb .gt. na) then
X 
X !   B is longer than A -- include extra words of B in C.
X 
X       do i = nm + 3, nb + 2
X         d(i) = b(i)
X       enddo
X     endif
X   else
X 
X !   A and B have the same exponent but the opposite sign.  It is thus
X !   necessary to scan through each vector until we find an unequal word.
X 
X     do i = 3, nm + 2
X       if (a(i) .ne. b(i)) goto 180
X     enddo
X 
X !   All words up to the common length are equal.
X 
X     if (na .eq. nb) then
X 
X !   The length of A is the same as B -- result is zero.
X 
X       c(1) = 0.d0
X       c(2) = 0.d0
X       goto 420
X     elseif (na .gt. nb) then
X 
X !   A is longer -- thus trailing words of A are shifted to start of C.
X 
X       nn = na - nb
X       d(1) = sign (nn, ia)
X       d(2) = a(2) - nb
X       d(nn+3) = 0.d0
X       d(nn+4) = 0.d0
X 
X       do i = 3, nn + 2
X         d(i) = a(i+nb)
X       enddo
X     elseif (nb .gt. na) then
X 
X !   B is longer -- thus trailing words of B are shifted to start of C.
X 
X       nn = nb - na
X       d(1) = sign (nn, ib)
X       d(2) = b(2) - na
X       d(nn+3) = 0.d0
X       d(nn+4) = 0.d0
X 
X       do i = 3, nn + 2
X         d(i) = b(i+na)
X       enddo
X     endif
X     goto 410
X 
X !   An unequal word was found.
X 
X 180 k = i - 3
X     if (a(k+3) .gt. b(k+3)) then
X 
X !   A is larger -- subtract B (shifted) from A.
X 
X       d(1) = sign (nx - k, ia)
X       d(2) = a(2) - k
X       d(nx-k+3) = 0.d0
X       d(nx-k+4) = 0.d0
X 
X       do i = 3, nm - k + 2
X         d(i) = dble (a(i+k)) - dble (b(i+k))
X       enddo
X 
X       do i = nb - k + 3, na - k + 2
X         d(i) = a(i+k)
X       enddo
X 
X       do i = na - k + 3, nb - k + 2
X         d(i) = - b(i+k)
X       enddo
X     else
X 
X !   B is larger -- subtract A (shifted) from B.
X 
X       d(1) = sign (nx - k, ib)
X       d(2) = b(2) - k
X       d(nx-k+3) = 0.d0
X       d(nx-k+4) = 0.d0
X 
X       do i = 3, nm - k + 2
X         d(i) = dble (b(i+k)) - dble (a(i+k))
X       enddo
X 
X       do i = nb - k + 3, na - k + 2
X         d(i) = - a(i+k)
X        enddo
X 
X       do i = na - k + 3, nb - k + 2
X         d(i) = b(i+k)
X        enddo
X     endif
X   endif
X elseif (ma .gt. mb) then
X 
X !   Exponent of A is greater.  In other words, A has a larger magnitude.
X 
X   mc = ma - mb
X   la = min (mc, na)
X   lb = min (mc + nb, mpnw + 2)
X   lm = min (na, lb)
X   lx = min (max (na, lb), mpnw)
X   d(1) = sign (lx, ia)
X   d(2) = a(2)
X   d(lx+3) = 0.d0
X   d(lx+4) = 0.d0
X 
X   do i = 3, la + 2
X     d(i) = a(i)
X   enddo
X 
X !   If B is shifted MPNW + 2 or more words to the right of A then C = A.
X 
X   if (mc .ge. mpnw + 2) then
X     d(1) = sign (na, ia)
X     d(la+3) = 0.d0
X     d(la+4) = 0.d0
X     goto 410
X   endif
X   if (mc .gt. na) then
X 
X !   There is a gap between A and the shifted B.  Fill it with zeroes.
X 
X     do i = na + 3, mc + 2
X       d(i) = 0.d0
X     enddo
X 
X     lm = mc
X   endif
X   if (ia .eq. ib) then
X 
X !   A and B have the same sign -- add common words with B shifted right.
X 
X     do i = mc + 3, lm + 2
X       d(i) = dble (a(i)) + dble (b(i-mc))
X     enddo
X 
X !   Include tail of A or B, whichever is longer after shift.
X 
X     if (na .gt. lb) then
X       do i = lm + 3, na + 2
X         d(i) = a(i)
X       enddo
X     else
X 
X       do i = lm + 3, lb + 2
X         d(i) = b(i-mc)
X       enddo
X     endif
X   else
X 
X !   A and B have different signs -- subtract common words with B shifted right.
X 
X     do i = mc + 3, lm + 2
X       d(i) = dble (a(i)) - dble (b(i-mc))
X     enddo
X 
X !   Include tail of A or B, whichever is longer after shift.
X 
X     do i = lm + 3, na + 2
X       d(i) = a(i)
X     enddo
X 
X     do i = lm + 3, lb + 2
X       d(i) = - b(i-mc)
X     enddo
X   endif
X else
X 
X !   Exponent of B is greater.  In other words, B has a larger magnitude.
X 
X   mc = mb - ma
X   lb = min (mc, nb)
X   la = min (mc + na, mpnw + 2)
X   lm = min (nb, la)
X   lx = min (max (nb, la), mpnw)
X   d(1) = sign (lx, ib)
X   d(2) = b(2)
X   d(lx+3) = 0.d0
X   d(lx+4) = 0.d0
X 
X   do i = 3, lb + 2
X     d(i) = b(i)
X   enddo
X 
X !   If A is shifted MPNW + 2 or more words to the right of B then C = B.
X 
X   if (mc .ge. mpnw + 2) then
X     d(1) = sign (nb, ib)
X     d(lb+3) = 0.d0
X     d(lb+4) = 0.d0
X     goto 410
X   endif
X   if (mc .gt. nb) then
X 
X !   There is a gap between B and the shifted A.  Fill it with zeroes.
X 
X     do i = nb + 3, mc + 2
X       d(i) = 0.d0
X      enddo
X 
X     lm = mc
X   endif
X   if (ib .eq. ia) then
X 
X !   B and A have the same sign -- add common words with A shifted right.
X 
