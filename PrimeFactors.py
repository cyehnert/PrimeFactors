
#!/usr/bin/env python

import random
#  gcd only - sqrt is deliberately not used anywhere in this file, see is_prime
from math import gcd

#  the first 12 primes, used both to peel off easy factors and as Miller-Rabin
#  witnesses, for which the test is proven exact below 3317044064679887385961981
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
MR_PROVEN_BELOW = 3317044064679887385961981

def is_prime(n):
    #  Miller-Rabin. Exact for every n below MR_PROVEN_BELOW using the fixed
    #  witnesses above; past that we add random witnesses, which leaves a
    #  chance of error under 4**-20 per call
    if n < 2:
        return 0
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d = d // 2
        r = r + 1
    witnesses = list(SMALL_PRIMES)
    if n >= MR_PROVEN_BELOW:
        i = 0
        while i < 20:
            witnesses.append(random.randrange(2, n - 1))
            i = i + 1
    for a in witnesses:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        composite = 1
        i = 0
        while i < r - 1:
            x = x * x % n
            if x == n - 1:
                composite = 0
                break
            i = i + 1
        if composite:
            return 0
    return 1

def find_factor(n):
    #  Pollard's rho with Brent's cycle detection. Returns some non-trivial
    #  divisor of composite n - not necessarily a prime one. Expected work is
    #  about n**0.25, versus n**0.5 for trial division
    #
    #  Rho is randomized, so the SAME input takes wildly different times from
    #  run to run while the answer stays exact. Measured on this machine:
    #    2^101-1  (30 runs)  min 0.6s   median 2.5s   max 5.3s      8.6x spread
    #    2^128+1  ( 5 runs)  min 90s    median 211s   max 577s      6.4x spread
    #  and one earlier run of 2^128+1 finished in about 60s, so the real range
    #  is wider still - the distribution has a long right tail, because each
    #  retry draws a fresh (y, c) and we are waiting on the first lucky draw.
    #  A call that looks hung is usually just an unlucky draw, not a bug
    if n % 2 == 0:
        return 2
    while 1:
        y = random.randrange(1, n)
        c = random.randrange(1, n)
        m = 128
        g = 1
        r = 1
        q = 1
        ys = y
        x = y
        while g == 1:
            x = y
            i = 0
            while i < r:
                y = (y * y + c) % n
                i = i + 1
            k = 0
            while k < r and g == 1:
                ys = y
                bound = min(m, r - k)
                i = 0
                while i < bound:
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                    i = i + 1
                g = gcd(q, n)
                k = k + m
            r = r * 2
        if g == n:
            #  the batched gcd overshot; walk back one step at a time
            g = 1
            y = ys
            while g == 1:
                y = (y * y + c) % n
                g = gcd(abs(x - y), n)
        if g != n:
            return g
        #  g == n means this (y, c) pair failed; pick new ones and retry

def collect_factors(n, factors):
    #  split n until every piece is prime, appending the primes to factors
    if n == 1:
        return
    if is_prime(n):
        factors.append(n)
        return
    d = find_factor(n)
    collect_factors(d, factors)
    collect_factors(n // d, factors)

def prime_factors(n):
    if n < 2:
        return ""
    factors = []
    #  peel off the small primes by trial division first: it is cheaper than rho
    #  on the common cases and leaves rho only the genuinely large factors
    d = 2
    while d < 1000 and d * d <= n:
        while n % d == 0:
            factors.append(d)
            n = n // d
        if d == 2:
            d = 3
        else:
            d = d + 2
    if n > 1:
        collect_factors(n, factors)
    factors.sort()
    #  collapse the sorted primes into (p) and (p**e) groups
    result = ""
    i = 0
    while i < len(factors):
        p = factors[i]
        e = 0
        while i < len(factors) and factors[i] == p:
            e = e + 1
            i = i + 1
        if e == 1:
            result = result + "(" + str(p) + ")"
        else:
            result = result + "(" + str(p) + "**" + str(e) + ")"
    return result


def exp_loc(s):
    len_str = len(s)
    eloc = 0
    eflag = 0
    if len_str >= 3:
        for e in s:
            eloc += 1
            if e == '^':
                eloc -= 1
                eflag = 1
                break
        eloc = eloc * eflag
    return eloc

def sum_loc(s):
    len_str = len(s)
    sloc = 0
    sflag = 0
    if len_str >= 3:
        for e in s:
            sloc += 1
            if e == '+':
                sloc -= 1
                sflag = 1
                break
        sloc = sloc * sflag
    return sloc

def dif_loc(s):
    len_str = len(s)
    sloc = 0
    sflag = 0
    if len_str >= 3:
        for e in s:
            sloc += 1
            if e == '-':
                sloc -= 1
                sflag = 1
                break
        sloc = sloc * sflag
    return sloc

def tokenize(s):
    tokens = []
    number = ""
    for ch in s:
        if ch == '^' or ch == '+' or ch == '-':
            if number != "":
                tokens.append(number)
                number = ""
            tokens.append(ch)
        elif ch == ' ' or ch == '\t' or ch == '\r' or ord(ch) == 65279:
            if number != "":
                tokens.append(number)
                number = ""
        else:
            number = number + ch
    if number != "":
        tokens.append(number)
    return tokens

def is_valid(s1):
    for t in s1:
        if t != '^' and t != '+' and t != '-' and not t.isdigit():
            return 0
    return 1

if __name__ == '__main__':
    stop = 1
    while (int(stop) > 0):
        s = input("Input string: ")
        s1 = tokenize(s)
        #  spaces around the operators are optional: "a^n-b" and "a ^ n - b" both work
        #  all '^' are resolved first, then '+' and '-' in the order in which they occur,
        #  until one token remains, which is converted to integer and factored
        #  a leading '-' is not supported, and no ( or ) are involved at this point
        if len(s1) == 0:
            print("No input.")
        elif is_valid(s1) == 0:
            print("Cannot read: " + s)
        else:
            expression = len(s1) > 1
            eloc = exp_loc(s1)
            while eloc > 0:
                j1 = int(s1[eloc-1])
                j2 = int(s1[eloc+1])
                if (j2 <= 0):
                    c = 1
                if (j2 > 0):
                    c = j1
                n = 1
                while (n < j2):
                    c = c * j1
                    n += 1
                s1[eloc-1] = str(c)
                del s1[eloc]
                del s1[eloc]
                eloc = exp_loc(s1)
            while len(s1) > 1:
                sloc = sum_loc(s1)
                dloc = dif_loc(s1)
                if sloc > 0 and (dloc == 0 or sloc < dloc):
                    loc = sloc
                elif dloc > 0:
                    loc = dloc
                else:
                    break
                j1 = int(s1[loc-1])
                j2 = int(s1[loc+1])
                if s1[loc] == '+':
                    c = j1 + j2
                else:
                    c = j1 - j2
                s1[loc-1] = str(c)
                del s1[loc]
                del s1[loc]
            if len(s1) != 1:
                print("Cannot read: " + s)
            else:
                k = int(s1[0])
                if expression:
                    print(k)
                if k < 2:
                    print(str(k) + " has no prime factorization")
                else:
                    print(prime_factors(k))
        s = input("Continue (q or Q = quit): ")
        if (s == 'q' or s == 'Q'):
            stop = 0
