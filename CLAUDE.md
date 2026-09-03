# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Two-file Python console program wrapped in a Visual Studio Python project
(`PrimeFactors.pyproj` / `PrimeFactors.sln`). No dependencies beyond the stdlib, no test
suite, no package manifest.

- `PrimeFactors.py` — expression parser, REPL, and the factorization pipeline.
- `quadratic_sieve.py` — the multiple-polynomial quadratic sieve.

The dependency runs one way only: `PrimeFactors.py` imports `quadratic_sieve`, never the
reverse. That is why the sieve carries its own small Miller-Rabin (`is_small_prime`) instead
of reusing `is_prime` — importing back would be a cycle. Keep it that way.

Both files are listed as `<Compile Include=...>` in the `.pyproj`. A new file has to be added
there or Visual Studio will not show it, even though Python imports it fine at runtime.

## Commands

```powershell
python PrimeFactors.py     # run the interactive REPL loop
```

The `.pyproj` has `CoreCompile` commented out, so **Build in Visual Studio is a no-op** —
F5 (StartupFile = `PrimeFactors.py`) is the only thing the project file does.

There is no test runner configured. Exercise functions ad hoc:

```powershell
python -c "import PrimeFactors as p; print(p.prime_factors(360))"
python -c "import quadratic_sieve as q; print(q.qs_factor(2**128+1))"
```

## Architecture

Two layers: the factorizer (`PrimeFactors.py` plus `quadratic_sieve.py`) and the expression
evaluator (`PrimeFactors.py` only).

**Factorization** — a three-stage pipeline, all integer arithmetic:

1. `prime_factors(n)` trial-divides by primes below 1000. Cheap, and it clears the common
   cases before the expensive machinery starts.
2. `collect_factors(n, factors)` recurses on what's left: if `is_prime` says prime, keep it;
   otherwise `find_factor` splits it and both halves recurse.
3. `prime_factors` sorts the flat prime list and collapses runs into `(p)` / `(p**e)` groups.

The output contract is unchanged from the original trial-division version: a *string* of
parenthesized groups in ascending order, e.g. `360` → `(2**3)(3**2)(5)`, exponent `1` printed
without `**`, and `""` for `n < 2`.

`find_factor` is a **dispatcher over two algorithms that fail on opposite inputs**. It returns
*some* non-trivial divisor, not necessarily a prime one — that's why stage 2 recurses on both
halves instead of treating the result as final. Both `rho_factor` and `qs_factor` honour that
same contract.

- `rho_factor(n, budget)` is **Pollard's rho with Brent's cycle detection**. Its cost is about
  √p where p is the *smallest* factor, so it is superb whenever `n` has any modest factor and
  hopeless on a balanced semiprime. `budget` caps inner iterations; `0` means uncapped, which
  is the original behaviour. A capped call that runs out returns `0`.
- `qs_factor(n)` in `quadratic_sieve.py` is **MPQS**. Its cost depends only on the *size* of
  `n`, not on its factors — steady exactly where rho is hopeless, and pointless where rho is
  quick. It returns `0` rather than raising when it declines or gives up.

Above `QS_BITS`, `find_factor` runs a **budgeted rho pass first**, then the sieve, then
uncapped rho as a last resort. Keep that order. The pre-pass is not a formality: on a 49-digit
number whose smallest factor is 1000003, rho returns in **0.002s** where the sieve grinds for
**77s** — a 47000× difference. And the final uncapped call means the worst case here is the
old worst case, never worse.

Measured medians on balanced semiprimes (the case that motivated the sieve):

| digits | 20 | 23 | 27 | 32 | 43 | 48 | 52 | 56 |
|---|---|---|---|---|---|---|---|---|
| rho | 0.08s | 0.16s | 5.3s | 43.2s | hours+ | — | — | — |
| sieve | 0.03s | 0.04s | 0.21s | 0.53s | 10.7s | 32.3s | 94s | 292s |

60 digits was still running after 21 minutes of CPU and was abandoned, so treat the upper 50s
as the practical ceiling in pure Python.

The sieve is ahead everywhere it is defined, and the gap widens fast. `2^128+1` went from
**60-577s to 2.9s**; note the budgeted rho pre-pass does *not* find its 17-digit factor
(it returns `0` after 0.32s) — the sieve does the work there. Below 60 bits (~18 digits) the
sieve declines outright and rho handles those in well under a millisecond, so nothing is lost
by the floor.

Rho is randomized, and the runtime spread on one input is **large** — measured 0.6s to 5.3s
over 30 runs of `2^101-1`, and 60s to 577s on `2^128+1` before the sieve existed. Results stay
exact; only the wait varies. Benchmark over many runs and quote a median, never a single
timing, and don't treat a slow call as a hang. The sieve is far more predictable, which is a
second reason to prefer it where it applies.

### The sieve, in one paragraph

Collect many `u` for which `u² - n` is smooth over a fixed factor base, then combine them so
the product of the smooth parts is a perfect square; that yields `x² ≡ y² (mod n)` and
`gcd(x - y, n)` splits `n`. Sieving is what makes smooth values cheap to find: rather than
test each candidate, add `log p` at every position `p` divides and look only where the running
total approaches the expected size. The "multiple polynomial" part keeps the sieve interval
small by switching to a fresh `v(x) = ax² + 2bx + c` (with `a = q²`, `b² - ac = n`) whenever
one is exhausted, instead of widening a single interval until the values get too big to be
smooth.

Three things that are easy to break in a refactor:

- **`a = q²` is a perfect square on purpose.** Its exponent vector is all-even, so it drops
  out of the GF(2) matrix entirely and only `q` itself has to be carried into the final square
  root. Choosing `a` as a product of factor base primes (that would be SIQS) breaks this.
- **The large-prime tolerance is a trap.** `LP_TOLERANCE_BITS` lowers the sieve threshold so
  values that are smooth apart from one oversized prime still get looked at. Lowering it by
  that prime's *full* width is the obvious choice and is **2.2× slower** at 43 digits — the
  extra candidates cost more trial division than the extra relations return. At the tuned
  value of 8 the variation is worth 1.16× at 43 digits and 1.63× at 48. Re-tune by
  measurement on identical inputs before touching it.
- **`perfect_power` must run before the sieve.** This is correctness, not optimization: on
  `n = pᵏ` the sieve produces only trivial gcds and never terminates. `qs_factor` checks first.

The sieve's byte accumulators cap the size it can handle; past roughly 250 bits of expected
value it declines with `0` rather than overflowing, and rho takes over.

`is_prime` is Miller-Rabin, exact below `MR_PROVEN_BELOW` (3.3e24) with the fixed
`SMALL_PRIMES` witnesses, plus 20 random witnesses above it. Correctness of the whole
factorization rests on it — if it wrongly called a composite prime, that composite would be
emitted as a factor.

**Do not reintroduce `math.sqrt`.** An earlier version bounded the trial-division loop with
`while j <= math.sqrt(n)`, which broke two ways: it raised `OverflowError: int too large to
convert to float` above ~1.8e308 (`3**700` crashed outright), and past 2^53 the rounding could
fall *below* the true integer square root, silently skipping a divisor and reporting a
composite as prime. Integer forms like `d * d <= n` are exact at any size and faster besides.

The ban is on the **float** `sqrt` specifically, and does not extend to everything in `math`.
Verified in this environment on `10**400`:

| call | behaviour | verdict |
|---|---|---|
| `math.sqrt(n)` | `OverflowError: int too large to convert to float` | banned |
| `math.isqrt(n)` | exact, returns `10**200` | fine, and required by the sieve |
| `math.log(n)` | `921.034...`, no overflow — CPython reads the exponent rather than converting to a double | fine, and only ever feeds a heuristic |

So `PrimeFactors.py` imports `gcd` only, and `quadratic_sieve.py` imports `gcd`, `isqrt` and
`log`. `isqrt` is integer-exact at any size and is not a way around the rule; treating it as
one and "fixing" it back to `d * d` comparisons would just make the sieve slower.

**Expression evaluator** — `__main__` reads a line, runs it through `tokenize`, and reduces
the token list in place:

1. resolve all `^` (repeated multiplication, left to right),
2. then resolve `+` / `-`, leftmost operator first,
3. when one token remains, `int()` it, run the size edits, and factor it.

`tokenize(s)` is a character scan that emits `^`, `+`, `-` as their own tokens and treats
whitespace (plus `\r` and a BOM) only as a separator, so **spaces around operators are
optional**: `2^10-24`, `2 ^ 10 - 24` and `2^ 10 -24` all yield the same five tokens.
`is_valid` then rejects any token that is not an operator or a digit run, so `2*3` reports
"Cannot read" instead of dying in `int()`.

`exp_loc` / `sum_loc` / `dif_loc` are the same function over `^`, `+`, `-`: they scan a
*token list* (despite the `s`/`len_str` naming suggesting a string) and return the index of
the first matching operator, or `0` for "not found". **Index 0 doubles as the sentinel**, so
an operator in position 0 is invisible to the reducer — a well-formed infix expression never
puts one there, but this is why the `len_str >= 3` guard exists and why refactors that change
the return convention break the reduce loops. Each reduction rewrites `s1[loc-1]` with the
result and deletes the operator and right operand.

The `+`/`-` reducer applies **one** operator per pass and re-scans afterward, picking the
smaller of `sum_loc`/`dif_loc` so evaluation runs left to right. Keep that shape — an earlier
version applied both branches in a single pass and re-scanned only inside the `-` branch,
which made any input containing `+` read past the end of the shortened list.

Negative *intermediate* values are fine (`2-10+20` → `-8` → `12`): a reduction writes the
single token `'-8'`, the list is never re-tokenized, and the `*_loc` functions compare whole
tokens, so `'-8'` is not mistaken for the operator `'-'`.

### The size edits

Two separate limits, for two different failure modes. Neither is cosmetic.

`MAX_DIGITS` (50) is about *patience*. Above it the REPL prints the digit count and asks
before starting; anything else re-prompts. Nothing is refused — answer `y` and it runs. 50 is
where the measured curve turns: 48 digits is ~32s, 52 is ~94s, and 60 was still going after 21
CPU-minutes.

`MAX_EXP_DIGITS` (4000) is about *not crashing*. Python 3.11+ refuses `int` → `str` past 4300
digits, and the `^` reducer calls `str(c)` on every result, so before this edit **`2^20000`
died with a `ValueError` traceback** and took the program down with it. Three places enforce
it, and all three are needed:

- `too_long(s1)` rejects a typed-in digit run over the limit — this is what catches a 5000-digit
  literal, including one inside an expression like `999…9+1`, where the `+` reducer would
  otherwise `str()` the result.
- the `^` reducer breaks out of its multiply loop the moment `c.bit_length()` passes
  `MAX_EXP_BITS`. Checking *during* the loop rather than after is the point: `2^100000000`
  now reports in 0.2s instead of trying to build a 30-million-digit number first.
- `digit_count(n)` uses `bit_length()` and only falls back to `len(str(n))` below the limit,
  where that call is both safe and exact. Never reach for `len(str(n))` on an unbounded int
  in this codebase.

Together these keep every value the program prints or factors under ~4001 digits, comfortably
inside Python's ceiling.

Remaining limits, all deliberate: a **leading** `-` is unsupported (`-5+8` reports "Cannot
read") because index `0` is the not-found sentinel; there are no parentheses, `*`, or `/`;
`^` binds left to right, so `2^3^2` is 64, not 512; a result below 2 is reported as having no
prime factorization rather than passed to `prime_factors`, which would hang on 0.

## Conventions

The code is deliberately plain, imperative, and stdlib-only — explicit `while` loops, index
counters rather than `range`/`enumerate`, `0`/`1` in place of `False`/`True`, and string
concatenation for output. Match that style rather than introducing comprehensions, regex, or
helper modules.

## Verifying a change to the factorizer

There is no test suite. Check any change to `prime_factors` and friends against an independent
trial-division reference — exhaustively over `2..50000`, then on random values up to 1e11.
For large inputs where a reference is too slow, verify two invariants instead: the emitted
factors multiply back to `n`, and each one satisfies `is_prime`.

For a change to the sieve, that end-to-end check is necessary but far too coarse to localize a
bug — a wrong `sqrt_mod` just means no relations are ever found, which looks identical to
"slow". Test the pieces separately, each against a property rather than a golden value:

- `sqrt_mod(a, p)` → `r * r % p == a % p`. **Cover `p % 8 == 1` specifically**: everything else
  takes the `p % 4 == 3` closed form and leaves the whole Tonelli-Shanks loop unexecuted.
- `int_root(n, k)` → `r**k <= n < (r+1)**k`.
- `perfect_power(n)` → returns *some* base, not the smallest: `perfect_power(3**100)` is
  `3**50`, because `k = 2` is tried first. That is correct and useful, since any non-trivial
  divisor is enough for `collect_factors`. Assert the contract, not a particular base.
- `polynomial(n, q)` → `a == q*q`, `b*b - a*c == n`, and the identity the whole method rests
  on, `a * v(x) == (a*x + b)**2 - n` at several `x`.
- `qs_factor(n)` on balanced semiprimes built from known primes → `1 < d < n` and `n % d == 0`.

Then time it. Rho's 6-9× spread on one input means a single timing proves nothing; take a
median over several runs *and* several distinct semiprimes of the same size. Measure rho and
the sieve in the same conditions — a benchmark run while something else is loading the machine
will silently favour whichever one you ran second.
