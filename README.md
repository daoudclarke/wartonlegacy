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

ANNPRM2 determines with certainty whether any number of the form h*10^s + 1 with h < 10^s is prime. One number found to be prime in this way is 7630*(10^1000) + 1.

The above three programs are on UBASIC. Still under construction on the Fortran90 system is a primality testing procedure for all numbers. It is a rigorous procedure which exploits the properties of elliptic curves and its running time is considerably longer than the above three programs.


