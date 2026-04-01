"""
generate_figures.py — All figures for Paper G
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.size'] = 11
rcParams['figure.dpi'] = 150


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


def sieve(N):
    isp = bytearray(b'\x01') * (N + 1)
    isp[0] = isp[1] = 0
    for i in range(2, int(N**0.5) + 1):
        if isp[i]: isp[i*i::i] = bytearray(len(isp[i*i::i]))
    return isp


# ============================================================
# Fig 1: Chain RH -> GRH -> Computational consequences
# ============================================================
def fig1_chain():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis('off')

    # Boxes
    boxes = [
        (0.05, 0.7, 'Paper C+D\nRiemann Hypothesis\nfor ζ(s)', '#4CAF50'),
        (0.35, 0.7, 'Paper F\nGeneralised RH\nfor all L(s,χ)', '#2196F3'),
        (0.7, 0.9, 'Miller det.\nO(log⁴n)', '#FF9800'),
        (0.7, 0.65, 'Generators\nO(log⁶q)', '#FF9800'),
        (0.7, 0.4, 'Non-residues\nO(log²p)', '#FF9800'),
        (0.7, 0.15, 'Goldbach\nE(N)=O(N^½⁺ᵋ)', '#FF9800'),
    ]

    for x, y, text, color in boxes:
        bbox = dict(boxstyle='round,pad=0.4', facecolor=color, alpha=0.3, edgecolor=color)
        ax.text(x, y, text, fontsize=11, ha='left', va='center', bbox=bbox)

    # Arrows
    ax.annotate('', xy=(0.34, 0.7), xytext=(0.22, 0.7),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    for y in [0.9, 0.65, 0.4, 0.15]:
        ax.annotate('', xy=(0.69, y), xytext=(0.58, 0.7),
                    arrowprops=dict(arrowstyle='->', lw=1.5, color='gray'))

    ax.text(0.5, 0.02, 'RSA/DH/ECC security: UNAFFECTED', fontsize=12,
            ha='center', color='red', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#FFEBEE', edgecolor='red'))

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title('From RH to computational consequences', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('figures/fig1_chain.pdf', bbox_inches='tight')
    plt.savefig('figures/fig1_chain.png', bbox_inches='tight')
    print('Saved fig1_chain')


# ============================================================
# Fig 2: Miller witnesses vs n (log scale)
# ============================================================
def fig2_witnesses():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: number of witnesses
    exponents = np.arange(3, 310, 1)
    n_vals = 10.0**exponents
    witnesses = 2 * (exponents * np.log(10))**2

    axes[0].semilogy(exponents, witnesses, 'b-', lw=2, label='2 log²n (Miller bound)')
    axes[0].axvline(617, color='red', ls='--', alpha=0.7, label='RSA-2048')
    axes[0].axvline(100, color='orange', ls='--', alpha=0.7, label='100-digit prime')
    axes[0].set_xlabel('digits of n', fontsize=12)
    axes[0].set_ylabel('witnesses needed', fontsize=12)
    axes[0].set_title('(a) Miller deterministic: witnesses vs size', fontsize=12)
    axes[0].legend(fontsize=9)
    axes[0].set_xlim(0, 310)

    # Panel B: comparison AKS vs Miller
    digits = np.array([10, 50, 100, 200, 500, 1000])
    bits = digits * np.log2(10)
    aks_ops = bits**7.5
    miller_ops = bits**4

    axes[1].semilogy(digits, aks_ops, 'r-o', lw=2, label='AKS: O(log^{7.5} n)')
    axes[1].semilogy(digits, miller_ops, 'b-s', lw=2, label='Miller det.: O(log⁴ n)')
    axes[1].set_xlabel('digits of n', fontsize=12)
    axes[1].set_ylabel('bit operations (relative)', fontsize=12)
    axes[1].set_title('(b) AKS vs Miller deterministic', fontsize=12)
    axes[1].legend(fontsize=10)
    axes[1].fill_between(digits, miller_ops, aks_ops, alpha=0.1, color='green')
    axes[1].text(300, 1e15, 'speedup\nlog^{3.5}n', fontsize=10, ha='center', color='green')

    plt.tight_layout()
    plt.savefig('figures/fig2_witnesses.pdf', bbox_inches='tight')
    plt.savefig('figures/fig2_witnesses.png', bbox_inches='tight')
    print('Saved fig2_witnesses')


# ============================================================
# Fig 3: Non-residues and generators vs bounds
# ============================================================
def fig3_bounds():
    isp = sieve(200000)
    primes = [p for p in range(5, 200001) if isp[p]]

    # Non-residues
    nr_list = []
    for p in primes:
        for a in range(2, p):
            if pow(a, (p-1)//2, p) != 1:
                nr_list.append((p, a))
                break

    # Generators (subset for speed)
    gen_list = []
    gen_primes = primes[::20][:500]
    for p in gen_primes:
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
                gen_list.append((p, g))
                break

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: non-residues
    ps_nr = np.array([x[0] for x in nr_list])
    nrs = np.array([x[1] for x in nr_list])
    bound_nr = np.log(ps_nr)**2

    axes[0].scatter(ps_nr[::5], nrs[::5], s=1, alpha=0.3, color='blue', label='n(p) measured')
    axes[0].plot(sorted(ps_nr), np.log(np.array(sorted(ps_nr)))**2, 'r-', lw=2, label='log²p (GRH bound)')
    axes[0].set_xlabel('prime p', fontsize=12)
    axes[0].set_ylabel('smallest non-residue n(p)', fontsize=12)
    axes[0].set_title('(a) Quadratic non-residues', fontsize=12)
    axes[0].legend(fontsize=10)
    axes[0].set_ylim(0, 50)

    # Panel B: generators
    ps_gen = np.array([x[0] for x in gen_list])
    gens = np.array([x[1] for x in gen_list])

    axes[1].scatter(ps_gen, gens, s=8, alpha=0.5, color='blue', label='g(p) measured')
    xs = np.linspace(5, ps_gen.max(), 200)
    axes[1].plot(xs, np.log(xs)**6 / 1000, 'r-', lw=2, label='log⁶p / 1000')
    axes[1].set_xlabel('prime p', fontsize=12)
    axes[1].set_ylabel('smallest generator g(p)', fontsize=12)
    axes[1].set_title('(b) Primitive roots', fontsize=12)
    axes[1].legend(fontsize=10)
    axes[1].set_ylim(0, 60)

    plt.tight_layout()
    plt.savefig('figures/fig3_bounds.pdf', bbox_inches='tight')
    plt.savefig('figures/fig3_bounds.png', bbox_inches='tight')
    print('Saved fig3_bounds')


# ============================================================
# Fig 4: Before/After table as visual
# ============================================================
def fig4_before_after():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis('off')

    headers = ['Problem', 'Before Paper F\n(conditional)', 'After Paper F\n(unconditional)']
    rows = [
        ['Primality test\n(deterministic)', 'AKS: O(log^{7.5} n)\nMiller: O(log⁴n) IF GRH',
         'Miller: O(log⁴n)\nPROVEN'],
        ['Find generator\nof (Z/qZ)*', 'g ≤ q^{1/4+ε}\n(Burgess, weak)',
         'g ≤ log⁶q\n(Shoup, PROVEN)'],
        ['Find non-residue\nmod p', 'n ≤ p^{0.15}\n(Burgess, weak)',
         'n ≤ log²p\n(Bach, PROVEN)'],
        ['RSA security', 'SECURE', 'SECURE\n(unchanged)'],
        ['DH/ECC security', 'SECURE', 'SECURE\n(unchanged)'],
    ]

    colors_before = ['#FFF3E0', '#FFF3E0', '#FFF3E0', '#E8F5E9', '#E8F5E9']
    colors_after = ['#E3F2FD', '#E3F2FD', '#E3F2FD', '#E8F5E9', '#E8F5E9']

    table = ax.table(cellText=rows, colLabels=headers, loc='center',
                     cellLoc='center', colColours=['#ECEFF1']*3)
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2.2)

    for i, row in enumerate(rows):
        table[i+1, 1].set_facecolor(colors_before[i])
        table[i+1, 2].set_facecolor(colors_after[i])

    ax.set_title('What changes and what does not', fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig('figures/fig4_before_after.pdf', bbox_inches='tight')
    plt.savefig('figures/fig4_before_after.png', bbox_inches='tight')
    print('Saved fig4_before_after')


if __name__ == '__main__':
    fig1_chain()
    fig2_witnesses()
    fig3_bounds()
    fig4_before_after()
    print('\nAll figures generated.')
