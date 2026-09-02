# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Single-file Python console program (`PrimeFactors.py`) wrapped in a Visual Studio Python
project (`PrimeFactors.pyproj` / `PrimeFactors.sln`). No dependencies beyond the stdlib,
no test suite, no package manifest.

## Commands

```powershell
python PrimeFactors.py     # run the interactive REPL loop
```

The `.pyproj` has `CoreCompile` commented out, so **Build in Visual Studio is a no-op** —
F5 (StartupFile = `PrimeFactors.py`) is the only thing the project file does.

There is no test runner configured. Exercise functions ad hoc:

```powershell
python -c "import PrimeFactors as p; print(p.prime_factors(360))"
```

## Architecture

Two layers, both in `PrimeFactors.py`:

**Factorization** — a three-stage pipeline, all integer arithmetic:

1. `prime_factors(n)` trial-divides by primes below 1000. Cheap, and it clears the common
   cases before the expensive machinery starts.
2. `collect_factors(n, factors)` recurses on what's left: if `is_prime` says prime, keep it;
   otherwise `find_factor` splits it and both halves recurse.
3. `prime_factors` sorts the flat prime list and collapses runs into `(p)` / `(p**e)` groups.

The output contract is unchanged from the original trial-division version: a *string* of
parenthesized groups in ascending order, e.g. `360` → `(2**3)(3**2)(5)`, exponent `1` printed
without `**`, and `""` for `n < 2`.

`find_factor` is **Pollard's rho with Brent's cycle detection** — expected O(n^¼) rather than
O(√n). It returns *some* non-trivial divisor, not necessarily a prime one; that's why stage 2
recurses on both halves instead of treating the result as final. The batched-gcd inner loop is
the fast path; the `g == n` backtrack walks one step at a time to recover the factor the batch
overshot.

It is randomized, and the runtime spread on one input is **large** — measured 0.6s to 5.3s
over 30 runs of `2^101-1`, and 60s to 577s on `2^128+1`. Results stay exact; only the wait
varies. Benchmark it over many runs and quote a median, never a single timing, and don't treat
a slow call as a hang. The measured figures are recorded in the comment on `find_factor`.

`is_prime` is Miller-Rabin, exact below `MR_PROVEN_BELOW` (3.3e24) with the fixed
`SMALL_PRIMES` witnesses, plus 20 random witnesses above it. Correctness of the whole
factorization rests on it — if it wrongly called a composite prime, that composite would be
emitted as a factor.

**Do not reintroduce `math.sqrt`.** Only `gcd` is imported from `math`. An earlier version
bounded the trial-division loop with `while j <= math.sqrt(n)`, which broke two ways: it
raised `OverflowError: int too large to convert to float` above ~1.8e308 (`3**700` crashed
outright), and past 2^53 the rounding could fall *below* the true integer square root,
silently skipping a divisor and reporting a composite as prime. Integer forms like `d * d <= n`
are exact at any size and faster besides.

**Expression evaluator** — `__main__` reads a line, runs it through `tokenize`, and reduces
the token list in place:

1. resolve all `^` (repeated multiplication, left to right),
2. then resolve `+` / `-`, leftmost operator first,
3. when one token remains, `int()` it and factor it.

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
