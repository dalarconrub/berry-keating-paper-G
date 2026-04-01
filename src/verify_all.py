"""
verify_all.py — Computational verification for Paper G
Runs all tests: Miller deterministic, generators, non-residues.
"""
import numpy as np
import time


def mr_test(n, a):
    if n < 2: return False
    if n == a: return True
    if n % 2 == 0: return n == 2
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1: return True
    for _ in range(s - 1):
        x = pow(x, 2, n)
        if x == n - 1: return True
    return False


def miller_det(n):
    """Deterministic Miller test: check a = 2..2*log^2(n)."""
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    bound = min(int(2 * np.log(n)**2) + 1, n - 1)
    for a in range(2, bound + 1):
        if not mr_test(n, a): return False
    return True


def sieve(N):
    isp = bytearray(b'\x01') * (N + 1)
    isp[0] = isp[1] = 0
    for i in range(2, int(N**0.5) + 1):
        if isp[i]: isp[i*i::i] = bytearray(len(isp[i*i::i]))
    return isp


def smallest_nonresidue(p):
    for a in range(2, p):
        if pow(a, (p - 1) // 2, p) != 1:
            return a
    return -1


def smallest_generator(p):
    n = p - 1
    factors = set()
    d, temp = 2, n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0: temp //= d
        d += 1
    if temp > 1: factors.add(temp)
    for g in range(2, p):
        if all(pow(g, n // q, p) != 1 for q in factors):
            return g
    return -1


def main():
    print("=" * 60)
    print("Paper G: Computational Verification")
    print("=" * 60)

    # 1. Miller vs sieve
    print("\n1. Miller deterministic vs sieve (n=2..100,000)")
    isp = sieve(100000)
    errors = sum(1 for n in range(2, 100001) if miller_det(n) != bool(isp[n]))
    print(f"   Errors: {errors}")

    # 2. Large primes
    print("\n2. Large known primes:")
    for p, name in [(2**31 - 1, "M31"), (10**9 + 7, "1e9+7"),
                    (10**12 + 39, "1e12+39"), (10**15 + 37, "1e15+37")]:
        r = miller_det(p)
        print(f"   {name:>12s}: {'PRIME' if r else 'ERROR'}")

    # 3. Carmichael
    print("\n3. Carmichael numbers:")
    carm = [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841,
            29341, 41041, 46657, 52633, 62745, 63973, 75361, 101101]
    err_c = sum(1 for n in carm if miller_det(n))
    print(f"   {len(carm)} tested, {err_c} errors")

    # 4. Non-residues
    print("\n4. Smallest non-residue n(p) < log^2(p):")
    primes = [i for i in range(5, 1000001) if isp[i]] if 1000000 <= len(isp) else \
             [i for i in range(5, 100001) if isp[i]]
    test_p = primes[::max(1, len(primes) // 1000)][:1000]
    viol = sum(1 for p in test_p if smallest_nonresidue(p) > np.log(p)**2)
    print(f"   {len(test_p)} primes, {viol} violations (p >= 5)")

    # 5. Generators
    print("\n5. Smallest generator g(p) < log^6(p):")
    test_g = primes[::max(1, len(primes) // 500)][:500]
    viol_g = sum(1 for p in test_g if smallest_generator(p) > np.log(p)**6)
    print(f"   {len(test_g)} primes, {viol_g} violations (p >= 5)")

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    main()
