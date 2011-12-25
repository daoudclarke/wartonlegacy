*In the name of God, the Merciful, the Compassionate.*

#The Warton Legacy

This repository contains code written by my late uncle, Robert Warton. He devoted the latter part of his life to the study of prime numbers and the problem of factorisation of the products of large primes. This is the result of his many years of effort.

Unfortunately, he left behind very little documentation, and the code is poorly commented (if at all). Hopefully, with time, we can document what he has made and make it a useful resource for researchers and students.

What follows is some of the documentation he gave me.

# Report on Cryptography System

## Preliminaries

The programs are created from two compilers: UBASIC, a multi precision extension of BASIC and FTN90, a Fortran 90 developed by Numerical Algorithms Group.

## Overview

The system consists of more than one hundred programs in all performing the following functions: RSA code scrambling and unscrambling, code generating via robust prime construction, code breaking and secret sharing. Routines in number theory and algebraic numbers exist both within programs and separately.

## RSA Code Scrambling and Unscrambling on UBASIC

BOBTOT2 determines second public key minimising number of fixed points. On typing RUN the computer will ask for first public key (product of two integers) and then separately the two constituents (which will be known by sendee). It also asks for a trial key, which must be an odd number. The output will be the second public key of lowest number of fixed points. BOBRSA1 scrambles message which is a number less than public key 1. It requests from the user the message, first public key without its factors which are unknown by the sender and the second public key.

BOBRSA3 unscrambles transmission from BOBRSA1. This program can only be operated by the sendee who is the only one who knows the breakdown of public key 1. Requested information: scrambled message, public key 1 and its constituents and public key 2 output is the scrambled message.

## The Code Generators

BOBPRM2 determines whether number input is composite or a strong pseudo primte to twenty random bases. If number is found to be composite the program will find the nearest strong pseudo prime greater than number input. The probability of a number found in this manner not being prime is 1 in 10 to the power 60 approximately.

ANNPRM determines with certainty whether any number of the form h*2^s + 1 with h < 2^s is prime (^ denotes "to the power of").

ANNPRM2 determines with certainty whether any number of the form h*10^s + 1 with h < 10^s is prime.
One number found to be prime in this way is 7630*(10^1000) + 1.

The above three programs are on UBASIC. Still under construction on the Fortran90 system is a primality testing procedure for all numbers. It is a rigorous procedure which exploits the properties of elliptic curves and its running time is considerably longer than the above three programs.

## The Factorizers

These can be divided into two categories --- special number factorizers and general number factorizers. Because the analyst is unlikely to know anything about the number he has to decompose, it is recommended that he take advantage of the very much faster special factorizers and run them in the first instance for a suitable length of time. Then should the numbers happen to be be of the appropriate structure the message will succumb to this method without necessitating the use of the general factorizers.

### The Special Factorizers

* **P-1 Procedure.** Should any prime factor - 1 of a number be smooth (i.e. all factors of this number less than a prescribed number, e.g. 1,000,000, the number will be factorized using this approach. On UBASIC the appropriate program is BOBPQ2 and on FTN90, BOBPQ1.F90. Note that on all Fortran applications the input expected will be in lots of 4 digits where leading zeroes can be replaced by spaces. The first information requested will always be the length of the number radix 10000 i.e. number of batches of 4 digits.

* **P+1 Procedure.** Should any prime factor + 1 of a given number be smooth as defined above this procedure will factorize it. Relevant programs on UBASIC are BOBPLUS2 and ANNPQ2.

* **P^2 + 1 Procedure.** As above with any prime whose P^2 + 1 is smooth. Relevant program is ANN4 on UBASIC.

* **P^2 + P + 1.** As above with any prime whose P^2 + P + 1 is smooth. BOBFER enables numbers of any length whose factors have the same length and leading halves identical to be factorized instantaneously, e.g. 12345679876543 * 12345671234567.

### General Factorizers

In this category there are three main approaches - Elliptic Curves, (Multiple) Polynomial Quadratic Sieves wherein same prime base is used for each polynomial, and (Multiple) General Number Field Sieves. Much Factorization time will be saved if it is discovered early that numbers have the structure r*s^t+u where r and s are small integers. These numbers are most efficiently handled as Aurifeullian and by the Special Number Field Sieve Procedures.

Mention should be made of two Binary Quadratic Forms techniques - BOBCLAS and BOBGR1 (and extension BOBGR2), and Pollard Rho Factorizer (BOBRHO) all on UBASIC. Although these techniques have been fully developed, they are largely obsolete and have been replaced by the more powerful Elliptic Curve routines. They are however used **within** other techniques.

In practice after having tried the special factorizers above, it is a good idea to try Elliptical Curve Methods next. This is because their running times are determined not by the size of the composite but mainly by the size of the smallest factor. Thus much time is saved by snipping off small factors (under say 10^25) where they exist.

For both the UBASIC and FTN90 Elliptic Curve Methods the information that will be requested will be length of minimum factor in powers of 10. A good practical rule here is to input length in denary digits **plus** 1 for safety. Thus if you suspected that the smallest factor was between 10^20 and 10^21 it is a good idea to respond here with the number 22. For the stage 2 multiplier requested a value of 40 is usually near optimal, but in FTN90's BOBECM8, which has the record-breaking special K factor of 250, values as high as 80 should be employed. Due to the fact that UBASIC has not been extended to take advantage of the large memory available it was not possible to incorporate these innovations there.

Of the some dozen Elliptic Curve Techniques on UBASIC the fastest is BOBONEC7 for all but the very big runs which should use BOBONEC4. The information requested in order: "code", "length etc.", "stage 2 multiplier".

On UBASIC's BOBONEC7 provision is made for restarting incomplete jobs. The usual response to system query "elliptic curve number" is the number 3. However, should a run have executed upto curve x without terminating, the run can be resumed by replying with number x+1.

On FTN90 the most commonly used Elliptic Curve Techniques are BOBECM8.F90 and BOBECM4 which potentiates this approach by employing special polynomials to compute powers of points. The information requested in both cases is the same as in the UBASIC counterpoint except that information is expected in batches of 4 digits. 

###Multiple Polynomial Quadratic Sieve Techniques

Unlike the ECM methods above these techniques consist of suites or batches of programs:

**Ann Suite**. On FTN90 load and go with ANN1.F90, when completed load and go on ANN2.F90, when completed load and go on ANN3.F90, finally when this is completed load and go on ANN4.F90. This suite can be prepared as a batch so that it will only be necessary to load and go once. In either case apart from the run instruction the only input will be length of code radiz 10000 and the code itself in batches of 4 digits. This program will be used to decompose "hard" numbers of between 40 and 68 digits.

**Greta Suite**. As above except that 5 programs are involved - GRETA1.F90, GRETA2.F90, GRETA3.F90, GRETA4.F90, and GRETA5.F90. This suite is more powerful than the ANN suite as it employs the Double Large Prime Variant. Recommended maximum length of input code is 72 decimal digits. For longer codes (see below) a special linear dependence determiner should be used.