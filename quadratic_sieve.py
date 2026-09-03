
#!/usr/bin/env python

#  Multiple-polynomial quadratic sieve (MPQS).
#
#  Pollard's rho in PrimeFactors.py costs about sqrt(p), where p is the SMALLEST
#  prime factor of n. That makes it excellent whenever n has any modest factor,
#  and hopeless on a balanced semiprime, where the smallest factor is near
#  sqrt(n). The sieve's cost depends only on the size of n, not on the size of
#  its factors, so the two are complements rather than replacements. See
#  find_factor in PrimeFactors.py for the dispatch between them.
#
#  The idea: collect many u with u*u - n smooth over a fixed set of small primes,
#  then combine them so the product of the smooth parts is a perfect square. That
#  gives u*u == y*y (mod n) with luck u != +-y, and then gcd(u - y, n) splits n.
#  Sieving is what makes the smooth ones cheap to find - rather than test every
#  candidate, add log(p) at every position p divides, then look only at the
#  positions whose running total is close to the expected size of the value.
#
#  isqrt is used freely here: it is exact integer arithmetic at any size. The
#  float sqrt stays banned repo-wide (see CLAUDE.md); log on a big int is safe
#  because CPython reads the exponent instead of converting to a double.

from math import gcd, isqrt, log

#  witnesses for the local primality test, exact below 3317044064679887385961981
MR_BASES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

#  (max digits, factor base bound B, sieve half-width M).
#  B follows the textbook exp(0.5 * sqrt(ln n * ln ln n)) scaled down, which is
#  what pure Python wants: a smaller base means fewer relations to collect and a
#  much cheaper matrix, at the cost of more polynomials. Tuned by measurement -
#  the figures are in the comment on qs_factor
QS_PARAMS = [
    (24,    1200,    15000),
    (30,    3000,    25000),
    (36,    8000,    40000),
    (42,   20000,    65000),
    (48,   50000,   100000),
    (54,  110000,   160000),
    (60,  250000,   250000),
    (66,  500000,   400000),
]

#  the sieve adds log(p) once per prime, so prime powers are under-counted, and
#  the logs are rounded to whole bits. This is how much slack the threshold
#  allows for that. Too small loses relations, too large wastes trial division
SIEVE_ERROR_BITS = 12

#  collect this many relations beyond the width of the matrix, so elimination
#  yields several independent dependencies to try
RELATION_MARGIN = 12

#  large prime variation: a value that trial division reduces to a single prime
#  below LARGE_PRIME_MULT * B is kept as a "partial". Two partials sharing that
#  prime multiply into a usable relation, which is far cheaper than sieving up
#  another full one. The bound stays below B*B, and a residue under B*B with no
#  factor in the base is necessarily prime - that is what makes this sound
LARGE_PRIME_MULT = 128

#  how far the sieve threshold is lowered to let partials through. Lowering it by
#  the large prime's full width catches every partial, but the extra candidates
#  cost more trial division than the extra relations are worth - measured 2.2x
#  SLOWER at 43 digits that way. Tuned by measurement; see the comment on qs_factor
LP_TOLERANCE_BITS = 8

#  give up rather than sieve forever; find_factor falls back to rho
MAX_POLYNOMIALS = 100000

def is_small_prime(x):
    #  deterministic Miller-Rabin. Only ever called on the polynomial primes q,
    #  which are far below the proven bound. Kept local so this module imports
    #  nothing from PrimeFactors.py, which imports this one
    if x < 2:
        return 0
    for p in MR_BASES:
        if x % p == 0:
            return x == p
    d = x - 1
    r = 0
    while d % 2 == 0:
        d = d // 2
        r = r + 1
    for a in MR_BASES:
        t = pow(a, d, x)
        if t == 1 or t == x - 1:
            continue
        composite = 1
        i = 0
        while i < r - 1:
            t = t * t % x
            if t == x - 1:
                composite = 0
                break
            i = i + 1
        if composite:
            return 0
    return 1

def int_root(n, k):
    #  largest r with r**k <= n, by Newton iteration on integers. The seed is at
    #  or above the true root, so the iteration descends onto it
    if n < 2:
        return n
    r = 1 << ((n.bit_length() + k - 1) // k)
    while 1:
        t = ((k - 1) * r + n // pow(r, k - 1)) // k
        if t >= r:
            break
        r = t
    while r > 1 and pow(r, k) > n:
        r = r - 1
    while pow(r + 1, k) <= n:
        r = r + 1
    return r

def perfect_power(n):
    #  if n == b**k for some k >= 2 return b, else 0.
    #  This is a correctness requirement, not an optimization: on a perfect power
    #  the sieve produces only trivial gcds and never terminates
    k = 2
    while k <= n.bit_length():
        r = int_root(n, k)
        if r < 2:
            break
        if pow(r, k) == n:
            return r
        k = k + 1
    return 0

def sieve_primes(limit):
    #  every prime up to limit, by Eratosthenes
    flags = bytearray(b'\x01') * (limit + 1)
    flags[0] = 0
    if limit >= 1:
        flags[1] = 0
    p = 2
    while p * p <= limit:
        if flags[p]:
            j = p * p
            while j <= limit:
                flags[j] = 0
                j = j + p
        p = p + 1
    primes = []
    p = 2
    while p <= limit:
        if flags[p]:
            primes.append(p)
        p = p + 1
    return primes

def legendre(a, p):
    #  1 if a is a non-zero square mod p, p-1 if it is not, 0 if p divides a
    return pow(a % p, (p - 1) // 2, p)

def sqrt_mod(a, p):
    #  Tonelli-Shanks: r with r*r % p == a % p, assuming a is a square mod p
    a = a % p
    if a == 0:
        return 0
    if p == 2:
        return a
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q = q // 2
        s = s + 1
    z = 2
    while legendre(z, p) != p - 1:
        z = z + 1
    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)
    while t != 1:
        i = 0
        t2 = t
        while t2 != 1:
            t2 = t2 * t2 % p
            i = i + 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = b * b % p
        t = t * c % p
        r = r * b % p
    return r

def qs_parameters(n):
    #  bit_length rather than len(str(n)): str() on a very large int is quadratic
    #  and Python caps it at 4300 digits by default
    digits = int(n.bit_length() * 0.30102999566398120) + 1
    for row in QS_PARAMS:
        if digits <= row[0]:
            return row[1], row[2]
    row = QS_PARAMS[len(QS_PARAMS) - 1]
    return row[1], row[2]

def build_factor_base(n, primes):
    #  keep only the primes for which n is a quadratic residue - no other prime
    #  can ever divide a value of the polynomial. 2 is kept but handled apart
    fb_primes = []
    fb_roots = []
    fb_logs = []
    for p in primes:
        if p == 2:
            fb_primes.append(2)
            fb_roots.append(1)
            fb_logs.append(1)
        elif legendre(n, p) == 1:
            fb_primes.append(p)
            fb_roots.append(sqrt_mod(n, p))
            fb_logs.append(int(log(p, 2) + 0.5))
    return fb_primes, fb_roots, fb_logs

def next_prime_q(n, start):
    #  the next prime q > start with q % 4 == 3 and n a square mod q. The first
    #  condition gives sqrt_mod a closed form, the second is what makes a root
    #  exist at all
    q = start
    if q % 2 == 0:
        q = q + 1
    while 1:
        q = q + 2
        if q % 4 != 3:
            continue
        if not is_small_prime(q):
            continue
        if legendre(n, q) == 1:
            return q

def polynomial(n, q):
    #  v(x) = a*x*x + 2*b*x + c with a = q*q and b*b - a*c = n, so that
    #  a * v(x) == (a*x + b)**2 - n. Sieving v rather than (a*x+b)**2 - n is the
    #  whole point: v is smaller by a factor of a
    a = q * q
    r = pow(n % q, (q + 1) // 4, q)
    #  Hensel: lift the root r mod q to a root b mod q*q
    t = ((n - r * r) // q) % q
    t = t * pow(2 * r, -1, q) % q
    b = (r + q * t) % a
    c = (b * b - n) // a
    return a, b, c

def sieve_polynomial(a, b, c, m, fb_primes, fb_roots, fb_logs, thresh):
    #  add each prime's log at the positions it divides, then return the x whose
    #  total reached thresh. Those are the candidates worth trial dividing
    size = 2 * m + 1
    sieve = bytearray(size)
    count = len(fb_primes)
    j = 0
    while j < count:
        p = fb_primes[j]
        lg = fb_logs[j]
        if p == 2:
            #  n is odd and a is odd, so v(x) == x + c (mod 2): one parity class
            i = (c + m) % 2
            while i < size:
                sieve[i] = sieve[i] + lg
                i = i + 2
            j = j + 1
            continue
        if a % p == 0:
            j = j + 1
            continue
        #  a*v(x) == (a*x + b)**2 - n, so v(x) == 0 (mod p) exactly when
        #  a*x + b == +-root (mod p)
        #  b is reduced first: it runs to about q*q, and multiplying that by ainv
        #  once per prime per polynomial is the difference between small-int and
        #  big-int arithmetic in the hottest loop in the file
        ainv = pow(a % p, -1, p)
        r = fb_roots[j]
        bp = b % p
        x1 = ainv * (r - bp) % p
        x2 = ainv * (-r - bp) % p
        i = (x1 + m) % p
        while i < size:
            sieve[i] = sieve[i] + lg
            i = i + p
        if x2 != x1:
            i = (x2 + m) % p
            while i < size:
                sieve[i] = sieve[i] + lg
                i = i + p
        j = j + 1
    #  scanning 2*m+1 bytes in a Python loop would cost more than the sieve fill
    #  itself, so map "at or above thresh" to 1 and let bytes.find do the walking
    table = bytearray(256)
    i = thresh
    while i < 256:
        table[i] = 1
        i = i + 1
    marks = sieve.translate(bytes(table))
    hits = []
    pos = marks.find(1)
    while pos >= 0:
        hits.append(pos - m)
        pos = marks.find(1, pos + 1)
    return hits

def factor_smooth(v, fb_primes):
    #  divide v by the factor base; returns the exponents and whatever is left,
    #  which is 1 exactly when v was smooth
    exps = []
    if v < 0:
        v = -v
    count = len(fb_primes)
    j = 0
    while j < count and v > 1:
        p = fb_primes[j]
        if v % p == 0:
            e = 0
            while v % p == 0:
                v = v // p
                e = e + 1
            exps.append((j, e))
        j = j + 1
    return exps, v

def relation_mask(exps, sign):
    #  only the parity of each exponent matters to the matrix. Bit 0 is the sign,
    #  so a dependency is forced to have a positive product too
    mask = sign
    for pair in exps:
        if pair[1] % 2:
            mask = mask ^ (1 << (pair[0] + 1))
    return mask

def merge_exps(e1, e2):
    #  add two exponent lists, both in increasing index order
    total = {}
    for pair in e1:
        total[pair[0]] = pair[1]
    for pair in e2:
        if pair[0] in total:
            total[pair[0]] = total[pair[0]] + pair[1]
        else:
            total[pair[0]] = pair[1]
    out = []
    for j in sorted(total):
        out.append((j, total[j]))
    return out

def try_dependency(n, mask, rel_u, rel_q, rel_exp, fb_primes):
    #  the relations named by mask multiply to x*x == (q_product * y)**2 (mod n).
    #  Exponents are accumulated and halved rather than forming the product of
    #  the values and taking isqrt, which would be enormous
    x = 1
    y = 1
    total = {}
    i = 0
    while mask:
        if mask & 1:
            x = x * rel_u[i] % n
            y = y * rel_q[i] % n
            for pair in rel_exp[i]:
                if pair[0] in total:
                    total[pair[0]] = total[pair[0]] + pair[1]
                else:
                    total[pair[0]] = pair[1]
        mask = mask >> 1
        i = i + 1
    for j in total:
        y = y * pow(fb_primes[j], total[j] // 2, n) % n
    d = gcd((x - y) % n, n)
    if d > 1 and d < n:
        return d
    d = gcd((x + y) % n, n)
    if d > 1 and d < n:
        return d
    return 0

def qs_factor(n):
    #  Returns a non-trivial divisor of composite n, or 0 if the sieve gives up.
    #  The divisor is not necessarily prime - the same contract find_factor has,
    #  so collect_factors keeps recursing on both halves
    #
    #  Seconds on balanced semiprimes, this machine, medians for the small sizes
    #  and single runs for the large ones:
    #    digits    20     23     27     32     43     48     52     56
    #    seconds  0.03   0.04   0.21   0.53   10.7   32.3     94    292
    #  60 digits was still running after 21 minutes of CPU and was abandoned, so
    #  treat the upper 50s as the practical ceiling for pure Python
    #
    #  The large prime variation is worth 1.16x at 43 digits and 1.63x at 48 on
    #  identical inputs, but only at LP_TOLERANCE_BITS = 8. Lowering the sieve
    #  threshold by the large prime's full width instead catches more partials
    #  and is 2.2x SLOWER at 43 digits - the extra candidates cost more trial
    #  division than the extra relations return. Re-tune by measurement, and on
    #  the same inputs, before changing either constant
    if n < 2:
        return 0
    if n % 2 == 0:
        return 2
    b = perfect_power(n)
    if b:
        return b
    if n.bit_length() < 60:
        #  the parameters below are meaningless this small; rho owns this range
        return 0
    bound, m = qs_parameters(n)
    primes = sieve_primes(bound)
    for p in primes:
        if n % p == 0:
            return p
    fb_primes, fb_roots, fb_logs = build_factor_base(n, primes)
    count = len(fb_primes)
    #  |v(x)| stays near m * sqrt(n/2) across the interval, which is what the
    #  choice of a below arranges. Anything much smaller than that has to be
    #  carrying extra small factors, which is what smooth means
    target = int(log(m, 2) + 0.5 * log(n, 2) - 0.5)
    if target > 250:
        #  the accumulators are single bytes, so a number this size would
        #  overflow them. Nothing this large is going to be sieved in Python
        #  anyway - decline and let rho have it
        return 0
    #  the threshold also has to allow for the large prime: a partial's sieve
    #  total is short by exactly that prime's worth of logs, since it is not in
    #  the factor base and so was never added anywhere
    lp_bound = bound * LARGE_PRIME_MULT
    thresh = target - SIEVE_ERROR_BITS - LP_TOLERANCE_BITS
    if thresh < 1:
        thresh = 1
    needed = count + RELATION_MARGIN
    rel_u = []
    rel_q = []
    rel_exp = []
    rel_mask = []
    #  a == sqrt(2n)/m keeps |v| smallest over x in [-m, m], and a == q*q
    q = isqrt(isqrt(2 * n) // m)
    if q < 3:
        q = 3
    partials = {}
    polys = 0
    while len(rel_mask) < needed and polys < MAX_POLYNOMIALS:
        q = next_prime_q(n, q)
        polys = polys + 1
        a, b, c = polynomial(n, q)
        hits = sieve_polynomial(a, b, c, m, fb_primes, fb_roots, fb_logs, thresh)
        for x in hits:
            v = (a * x + 2 * b) * x + c
            if v == 0:
                continue
            exps, residue = factor_smooth(v, fb_primes)
            sign = 0
            if v < 0:
                sign = 1
            u = (a * x + b) % n
            if residue == 1:
                rel_u.append(u)
                rel_q.append(q)
                rel_exp.append(exps)
                rel_mask.append(relation_mask(exps, sign))
            elif residue < lp_bound:
                #  v is smooth times one prime too big for the base. Useless
                #  alone, but paired with another value carrying the SAME prime
                #  the two multiply to smooth * residue**2 - the square drops
                #  out of the matrix and only the residue is carried into y.
                #  The first sighting is kept and paired with every later one
                if residue in partials:
                    u2, exps2, q2, sign2 = partials[residue]
                    merged = merge_exps(exps, exps2)
                    rel_u.append(u * u2 % n)
                    rel_q.append(q * q2 % n * residue % n)
                    rel_exp.append(merged)
                    rel_mask.append(relation_mask(merged, sign ^ sign2))
                else:
                    partials[residue] = (u, exps, q, sign)
    if len(rel_mask) < count + 1:
        return 0
    #  Gaussian elimination over GF(2). Each row is one int of exponent parities
    #  and one int recording which relations were combined into it
    rows = []
    i = 0
    while i < len(rel_mask):
        rows.append([rel_mask[i], 1 << i])
        i = i + 1
    total_rows = len(rows)
    used = bytearray(total_rows)
    col = 0
    ncols = count + 1
    while col < ncols:
        bit = 1 << col
        pivot = -1
        i = 0
        while i < total_rows:
            if used[i] == 0 and (rows[i][0] & bit):
                pivot = i
                break
            i = i + 1
        if pivot >= 0:
            used[pivot] = 1
            prow = rows[pivot]
            i = 0
            while i < total_rows:
                if i != pivot and (rows[i][0] & bit):
                    rows[i][0] = rows[i][0] ^ prow[0]
                    rows[i][1] = rows[i][1] ^ prow[1]
                i = i + 1
        col = col + 1
    #  a row that reduced to zero is a dependency: its relations multiply to a
    #  perfect square
    i = 0
    while i < total_rows:
        if rows[i][0] == 0:
            d = try_dependency(n, rows[i][1], rel_u, rel_q, rel_exp, fb_primes)
            if d:
                return d
        i = i + 1
    return 0
