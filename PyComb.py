# **********************************************************************
#   Module: Combinatorial Functions
#
#   Language: Python 3.4
#   Author: Saied H. Khayat
#   Date:   March 2016
#   URL: https://github.com/saiedhk/PyComb
#
#   Copyright Notice: Free use of this library is permitted under
#   the guidelines and in accordance with the MIT License (MIT).
#   http://opensource.org/licenses/MIT
#
# **********************************************************************

"""
Module Description:

This Python package provides a number of special functions that are used in
Combinatorics (Combinatorial Analysis).
Python is well-suited to combinatorial calculations as they normally produce
outputs that are huge in terms of the number of digits representing them.
The Python's built-in unlimited precision integers make such calculations
effortlessly easy.

In particular, the following functions are implemented in this package:
    Binomial, Multinomial, Catalan, Stirling of first kind, Stirling of second kind, Bell, Lah, Narayana,
    Delannoy, Motzkin, Schroder, Eulerian.

Users should refer to the following reference for definitions, and applications
of the above functions.

Reference:
    NIST's Digital Library of Mathematical Functions
    Chapter 26: Combinatorial Analysis
    Written by: D. M. Bressoud
    url: http://dlmf.nist.gov/26

These functions have been tested and verified against known results from reference tables.
There are no known bugs. However the users are encouraged to test before use.
The author would appreciate any bug reports or comments.

"""


import math

RP = 16  # reduction period


# --------------------------------------------------------------------
def gcd(a, b):
    """Compute greatest common divisor of a and b.
       This function is used in some of the functions in PyComb module.
    """
    r = a % b
    while r != 0:
        a = b
        b = r
        r = a % b
    return b


# --------------------------------------------------------------------
def prod(m, n):
    """Compute the product of consecutive integers in range [m,n] inclusive."""
    if n < m:
        return 0
    elif n == m:
        return n

    result = 1
    for i in range(m, n + 1):
        result *= i
    return result


# --------------------------------------------------------------------
def binomial(n, k):
    """Compute binomial coefficient of (n,k), a.k.a. C(n,k)"""
    if (n < k) or (n < 0) or (k < 0):
        return 0
    elif n == k:
        return 1
    elif k == 1:
        return n
    elif k == 0:
        return 1
    elif n - k < k:
        k = n - k

    num = n - k + 1  # numerator
    den = 1  # denominator

    for i in range(2, k + 1):
        num *= n - (k - i)
        den *= i
        if i % RP == 0:  # reduce numerator and denominator at every reduction period
            g = gcd(num, den)
            num = num // g
            den = den // g

    return num // den


# --------------------------------------------------------------------
def multinomial(nlist):
    """Compute multinomial coefficient given an list of integers (n1,n2,n3,...,nk).
       a.k.a. C(n;n1,n2,n3,...,nk), where n = sum(n1,n2,n3,...,nk).
    """
    n = nlist.copy()

    if len(n) == 0:
        return 0
    for i in range(0, len(n)):
        if n[i] < 0:
            return 0

    nsum = sum(n)
    if nsum == 0:
        return 0

    nmax = max(n)
    n.remove(nmax)

    num = 1  # numerator
    den = 1  # denominator
    m = n[0]
    n.remove(n[0])

    for i in range(nmax + 1, nsum + 1):
        num *= i
        if m > 1:
            den *= m
            m -= 1
        else:
            if len(n) > 0:
                m = n[0]
                n.remove(n[0])
        if i % RP == 0:  # reduce numerator and denominator
            g = gcd(num, den)
            num = num // g
            den = den // g

    return num // den


# --------------------------------------------------------------------
def catalan(n):
    """Compute Catalan number of n."""
    if n < 0:
        return 0
    return binomial(2 * n, n) // (n + 1)


# --------------------------------------------------------------------
def stirling2(n, m):
    """Compute Stirling number of second kind of (n,m)"""
    if (n < m) or (n < 0) or (m < 0):
        return 0

    result = 0
    for i in range(0, m + 1):
        if (m - i) % 2 == 0:
            sign = 1
        else:
            sign = -1
        result += (sign * binomial(m, i) * (i ** n))

    return result // math.factorial(m)


# --------------------------------------------------------------------
def stirling1(n, m):
    """Compute Stirling number of first kind of (n,m)"""
    if (n < m) or (n < 0) or (m < 0):
        return 0

    result = 0
    sign = 1
    for k in range(0, n - m + 1):
        product = binomial(n - 1 + k, n - m + k) * binomial(2 * n - m, n - m - k)
        product *= stirling2(n - m + k, k)
        result += (sign * product)
        sign = -sign

    return result


# --------------------------------------------------------------------
def bell(n):
    "Compute Bell number of n"
    if n < 0:
        return 0

    result = 0
    for k in range(0, n + 1):
        result += stirling2(n, k)

    return result


# --------------------------------------------------------------------
def lah(n, k):
    """Compute Lah number of (n,k)"""
    if (n < k) or (n < 1) or (k < 1):
        return 0

    result = 1
    for i in range(k + 1, n + 1):
        result *= i
    result *= binomial(n - 1, k - 1)
    return result


# --------------------------------------------------------------------
def delannoy(m, n):
    """Compute Delannoy number of (m,n)."""
    if (n < 0) or (m < 0):
        return 0

    result = 0
    for k in range(0, n + 1):
        result += (2**k) * binomial(m, k) * binomial(n,k)

    return result


# --------------------------------------------------------------------
def motzkin(n):
    """Compute Motzkin number of n."""
    if n < 0:
        return 0

    result = 0
    sign = 1
    for k in range(0, n + 1):
        temp = binomial(n,k) * binomial(2*n + 2 - 2*k, n+1-k)
        result += sign * (temp // (n+2-k))
        sign = -sign

    return result

# --------------------------------------------------------------------
def narayana(n, k):
    """Compute Narayana number of (n,k)."""
    if (n < 0) or (k < 0) or (n < k):
        return 0
    elif n == 0:
        return 1
    elif (k == 0) and (n > k):
        return 0
    elif (n == k) or (k == 1):
        return 1

    result = binomial(n,k) * binomial(n,k-1)
    return result // n


# --------------------------------------------------------------------
def schroder(n):
    """Compute Schroder number of n."""
    if n < 0:
        return 0

    return delannoy(n,n) - delannoy(n+1,n-1)


# --------------------------------------------------------------------
def eulerian(n,k):
    """Compute Eulerian number of n."""
    if (n < 0) or (k < 0):
        return 0
    elif (n == 0) and (k == 0):
        return 1
    elif n <= k:
        return 0

    result = 0
    sign = 1
    for j in range(0, k+1):
        result += sign * binomial(n+1,j) * ((k+1-j)**n)
        sign = -sign

    return result

# --------------------------------------------------------------------
