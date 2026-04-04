#!/usr/bin/env python3
"""
SCALE SEPARATION vs ERROR RATE SWEEP
======================================

At each of the 14 error rates, measure independently:
  - B31 baseline (static, no rotation)
  - T75 Scale 2 only (dynamic detection, lattice correction)
  - Scale 1 analytical (per-unit torus, S_torus = exp(3))
  - Combined: S_torus × S_L2 and S_torus × S_L1+L2

The question: does the multiplicative two-scale structure hold
uniformly across error rates, or does it break down at high/low ε?

Also: where does the rotation gap (T75 det - B31 det) peak?
Does the gap widen or narrow with error rate?

Usage: python scale_separation_sweep.py
Output: scale_separation_sweep_output.txt
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys
import os

# ============================================================================
# CONSTANTS
# ============================================================================

GATES = ['R', 'T', 'P', 'F', 'S']
NUM_GATES = 5
RANDOM_SEED = 42
TAU_STATIC = 1
TAU_DYNAMIC = 5

B31_ACTIVE_GATES = {0, 1, 4}  # R, T, S

EPSILON_RAW_VALUES = [
    3e-1, 2e-1, 1e-1,
    5e-2, 3e-2, 2e-2, 1e-2,
    5e-3, 3e-3, 2e-3, 1e-3,
    5e-4, 2e-4, 1e-4,
]

CELL_SIZES = [7, 19, 37]
MC_TRIALS_BASE = 50_000
MC_ASSIGNMENTS = 20

S_TORUS_ANALYTICAL = np.exp(3)  # ~20.09

# ============================================================================
# EISENSTEIN CELL (inline)
# ============================================================================

class EisensteinCell:
    UNIT_VECTORS = [(1,0),(-1,0),(0,1),(0,-1),(-1,-1),(1,1)]
    def __init__(self, radius):
        self.radius = radius
        self.r_sq = radius * radius
        self.nodes = []
        for a in range(-radius-1, radius+2):
            for b in range(-radius-1, radius+2):
                if a*a - a*b + b*b <= self.r_sq:
                    self.nodes.append((a,b))
        self.num_nodes = len(self.nodes)
        self.node_index = {n: i for i, n in enumerate(self.nodes)}
        self.edges = []
        self.neighbours = defaultdict(list)
        node_set = set(self.nodes)
        for i, (a1,b1) in enumerate(self.nodes):
            for da,db in self.UNIT_VECTORS:
                nb = (a1+da, b1+db)
                if nb in node_set:
                    j = self.node_index[nb]
                    if j > i: self.edges.append((i,j))
                    self.neighbours[i].append(j)
        self.sublattice = [(a+b)%3 for (a,b) in self.nodes]
        self.chirality = [0 if s==0 else (+1 if s==1 else -1) for s in self.sublattice]
        self.coordination = [len(self.neighbours[i]) for i in range(self.num_nodes)]

# ============================================================================
# DETECTION FUNCTIONS
# ============================================================================

def absent_gate_static(base, t):
    return base % NUM_GATES

def absent_gate_dynamic(base, chirality, t):
    return (base + chirality * t) % NUM_GATES

def detect_and_correct_b31(cell, assignment, error_node, error_gate):
    """B31: static detection, tau=1, chirality=0."""
    if error_gate not in B31_ACTIVE_GATES:
        return False, False
    node_votes = Counter()
    node_edges = defaultdict(set)
    gate_votes = Counter()
    detected = False
    for nbr in cell.neighbours[error_node]:
        ai = assignment[error_node] % NUM_GATES
        an = assignment[nbr] % NUM_GATES
        if ai == an: continue
        if an == error_gate:
            detected = True
            edge = (min(error_node, nbr), max(error_node, nbr))
            node_votes[error_node] += 1
            node_votes[nbr] += 1
            node_edges[error_node].add(edge)
            node_edges[nbr].add(edge)
            gate_votes[error_gate] += 1
    if not detected: return False, False
    pg = gate_votes.most_common(1)[0][0]
    cands = list(node_votes.keys())
    if len(cands) == 1: pn = cands[0]
    else:
        ranked = sorted(cands, key=lambda n: (-len(node_edges[n]),-node_votes[n],cell.coordination[n],n))
        pn = ranked[0]
    if pn == error_node and pg == error_gate:
        for nbr in cell.neighbours[pn]:
            an = assignment[nbr] % NUM_GATES
            if an != pg: return True, True
        return True, False
    return True, False

def detect_and_correct_t75(cell, assignment, error_node, error_gate, tau=TAU_DYNAMIC):
    """T75: dynamic detection, tau=5, real chirality."""
    node_votes = Counter()
    node_edges = defaultdict(set)
    gate_votes = Counter()
    detected = False
    for t in range(tau):
        for nbr in cell.neighbours[error_node]:
            ai = absent_gate_dynamic(assignment[error_node], cell.chirality[error_node], t)
            an = absent_gate_dynamic(assignment[nbr], cell.chirality[nbr], t)
            if ai == an: continue
            if an == error_gate:
                detected = True
                edge = (min(error_node, nbr), max(error_node, nbr))
                node_votes[error_node] += 1
                node_votes[nbr] += 1
                node_edges[error_node].add(edge)
                node_edges[nbr].add(edge)
                gate_votes[error_gate] += 1
    if not detected: return False, False
    pg = gate_votes.most_common(1)[0][0]
    cands = list(node_votes.keys())
    if len(cands) == 1: pn = cands[0]
    else:
        ranked = sorted(cands, key=lambda n: (-len(node_edges[n]),-node_votes[n],cell.coordination[n],n))
        pn = ranked[0]
    if pn == error_node and pg == error_gate:
        for t in range(tau):
            for nbr in cell.neighbours[pn]:
                an = absent_gate_dynamic(assignment[nbr], cell.chirality[nbr], t)
                if an != pg: return True, True
        return True, False
    return True, False

# ============================================================================
# GREEDY ASSIGNMENT FINDER
# ============================================================================

def find_valid_assignments(cell, rng, count, max_attempts=100_000):
    valid = []
    attempts = 0
    while len(valid) < count and attempts < max_attempts:
        attempts += 1
        n = cell.num_nodes
        assignment = [-1] * n
        order = rng.permutation(n)
        ok = True
        for idx in order:
            forbidden = set()
            for nbr in cell.neighbours[idx]:
                if assignment[nbr] >= 0:
                    forbidden.add(assignment[nbr] % NUM_GATES)
            available = [g for g in range(NUM_GATES) if g not in forbidden]
            if not available: ok = False; break
            assignment[idx] = int(rng.choice(available))
        if not ok: continue
        # verify
        good = True
        for (i,j) in cell.edges:
            if assignment[i] % NUM_GATES == assignment[j] % NUM_GATES:
                good = False; break
        if good:
            valid.append(tuple(assignment))
    return valid

# ============================================================================
# MC ENGINE
# ============================================================================

def mc_sweep_point(cell, assignments, eps_raw, seed):
    """Run both B31 and T75 at a single (cell, eps_raw) point."""
    rng = np.random.default_rng(seed)
    b31_gates = list(B31_ACTIVE_GATES)

    if eps_raw >= 0.01: num_trials = MC_TRIALS_BASE
    elif eps_raw >= 0.001: num_trials = MC_TRIALS_BASE * 2
    else: num_trials = MC_TRIALS_BASE * 4

    trials_per = max(1, num_trials // len(assignments))

    # B31 counters
    b31_total = 0; b31_inj = 0; b31_det = 0; b31_corr = 0; b31_unc = 0
    # T75 counters
    t75_total = 0; t75_inj = 0; t75_det = 0; t75_corr = 0; t75_unc = 0

    for assignment in assignments:
        for trial in range(trials_per):
            for node in range(cell.num_nodes):
                # ── B31 error model ──
                b31_total += 1
                if rng.random() < eps_raw:
                    absent = assignment[node] % NUM_GATES
                    possible = [g for g in b31_gates if g != absent]
                    if not possible: continue
                    g_err = int(rng.choice(possible))
                    b31_inj += 1
                    det, corr = detect_and_correct_b31(cell, assignment, node, g_err)
                    if det: b31_det += 1
                    if corr: b31_corr += 1
                    else: b31_unc += 1

                # ── T75 error model (separate random draw, same rate) ──
                t75_total += 1
                if rng.random() < eps_raw:
                    possible = [g for g in range(NUM_GATES) if g != assignment[node]]
                    g_err = int(rng.choice(possible))
                    t75_inj += 1
                    det, corr = detect_and_correct_t75(cell, assignment, node, g_err)
                    if det: t75_det += 1
                    if corr: t75_corr += 1
                    else: t75_unc += 1

    b31_det_rate = b31_det / b31_inj if b31_inj > 0 else 0
    b31_corr_rate = b31_corr / b31_inj if b31_inj > 0 else 0
    b31_logical = b31_unc / b31_total if b31_total > 0 else eps_raw
    b31_supp = eps_raw / b31_logical if b31_logical > 0 else float('inf')

    t75_det_rate = t75_det / t75_inj if t75_inj > 0 else 0
    t75_corr_rate = t75_corr / t75_inj if t75_inj > 0 else 0
    t75_logical = t75_unc / t75_total if t75_total > 0 else eps_raw
    t75_supp = eps_raw / t75_logical if t75_logical > 0 else float('inf')

    return {
        'eps_raw': eps_raw,
        'b31_det': b31_det_rate, 'b31_corr': b31_corr_rate,
        'b31_logical': b31_logical, 'b31_supp': b31_supp, 'b31_inj': b31_inj,
        't75_det': t75_det_rate, 't75_corr': t75_corr_rate,
        't75_logical': t75_logical, 't75_supp': t75_supp, 't75_inj': t75_inj,
        'rotation_gap_det': t75_det_rate - b31_det_rate,
        'rotation_gap_corr': t75_corr_rate - b31_corr_rate,
        'gap_supp': t75_supp / b31_supp if b31_supp > 0 else float('inf'),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():
    t_start = time.time()
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'scale_separation_sweep_output.txt')

    class Tee:
        def __init__(self, fp):
            self.terminal = sys.stdout
            self.file = open(fp, 'w', encoding='utf-8')
        def write(self, msg):
            try: self.terminal.write(msg)
            except UnicodeEncodeError:
                self.terminal.write(msg.encode('ascii','replace').decode('ascii'))
            self.file.write(msg)
        def flush(self):
            self.terminal.flush(); self.file.flush()

    tee = Tee(output_path)
    sys.stdout = tee

    print("=" * 90)
    print("  SCALE SEPARATION vs ERROR RATE SWEEP")
    print("  Does the two-scale multiplicative structure hold across all error rates?")
    print("=" * 90)
    print()

    # Build cells
    cells = {}
    cell_assignments = {}
    for radius, n in [(1,7),(2,19),(3,37)]:
        cell = EisensteinCell(radius)
        assert cell.num_nodes == n
        cells[n] = cell
        rng = np.random.default_rng(RANDOM_SEED + n)
        assignments = find_valid_assignments(cell, rng, MC_ASSIGNMENTS)
        cell_assignments[n] = assignments
        print(f"  {n}-node: {len(assignments)} assignments found")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 1: FULL SWEEP — B31 vs T75 at all 14 error rates
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 1: B31 vs T75 SWEEP (all 14 error rates)")
    print("=" * 90)
    print()

    all_results = {}

    for n_nodes in CELL_SIZES:
        cell = cells[n_nodes]
        assignments = cell_assignments[n_nodes]

        print(f"  {n_nodes}-node cell:")
        print(f"  {'eps_raw':>10} | {'B31_det':>8} {'B31_corr':>9} {'B31_supp':>9} | "
              f"{'T75_det':>8} {'T75_corr':>9} {'T75_supp':>9} | "
              f"{'det_gap':>8} {'corr_gap':>9} {'supp_ratio':>10}")
        print("  " + "-" * 100)

        for eps_raw in EPSILON_RAW_VALUES:
            seed = RANDOM_SEED + hash((n_nodes, eps_raw)) % 100000

            r = mc_sweep_point(cell, assignments, eps_raw, seed)
            all_results[(n_nodes, eps_raw)] = r

            print(f"  {eps_raw:>10.0e} | "
                  f"{r['b31_det']*100:>6.1f}% {r['b31_corr']*100:>7.1f}% {r['b31_supp']:>7.1f}x | "
                  f"{r['t75_det']*100:>6.1f}% {r['t75_corr']*100:>7.1f}% {r['t75_supp']:>7.1f}x | "
                  f"{r['rotation_gap_det']*100:>6.1f}pp {r['rotation_gap_corr']*100:>7.1f}pp "
                  f"{r['gap_supp']:>8.1f}x")

        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 2: SCALE COMPOSITION TABLE
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 2: SCALE COMPOSITION (S_torus x S_L2 x S_L1)")
    print(f"  S_torus = {S_TORUS_ANALYTICAL:.2f}x (analytical)")
    print("=" * 90)
    print()

    for n_nodes in CELL_SIZES:
        print(f"  {n_nodes}-node cell:")
        print(f"  {'eps_raw':>10} | {'S_L2(MC)':>10} | {'S_L1+L2':>10} | "
              f"{'S_torus*L2':>12} | {'S_torus*L1L2':>14} | "
              f"{'St=10*L1L2':>12}")
        print("  " + "-" * 80)

        for eps_raw in EPSILON_RAW_VALUES:
            r = all_results.get((n_nodes, eps_raw))
            if r is None: continue

            s_L2 = r['t75_supp']
            s_L1L2_cons = s_L2 / (1 - 0.5)  # f_sym = 0.5
            s_combined_L2 = S_TORUS_ANALYTICAL * s_L2
            s_combined_L1L2 = S_TORUS_ANALYTICAL * s_L1L2_cons
            s_conservative = 10.0 * s_L1L2_cons

            print(f"  {eps_raw:>10.0e} | {s_L2:>8.1f}x | {s_L1L2_cons:>8.1f}x | "
                  f"{s_combined_L2:>10.0f}x | {s_combined_L1L2:>12.0f}x | "
                  f"{s_conservative:>10.0f}x")

        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 3: ROTATION GAP vs ERROR RATE
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 3: ROTATION GAP vs ERROR RATE")
    print("  Does the gap widen, narrow, or stay constant?")
    print("=" * 90)
    print()

    for n_nodes in CELL_SIZES:
        print(f"  {n_nodes}-node cell:")
        print(f"  {'eps_raw':>10} | {'det gap (pp)':>12} | {'corr gap (pp)':>13} | "
              f"{'supp ratio':>10} | {'B31 det':>8} | {'T75 det':>8}")
        print("  " + "-" * 75)

        gaps = []
        for eps_raw in EPSILON_RAW_VALUES:
            r = all_results.get((n_nodes, eps_raw))
            if r is None: continue
            gaps.append(r['rotation_gap_det'] * 100)
            print(f"  {eps_raw:>10.0e} | {r['rotation_gap_det']*100:>10.1f}pp | "
                  f"{r['rotation_gap_corr']*100:>11.1f}pp | "
                  f"{r['gap_supp']:>8.1f}x | "
                  f"{r['b31_det']*100:>6.1f}% | {r['t75_det']*100:>6.1f}%")

        if gaps:
            print(f"\n    Gap statistics: mean={np.mean(gaps):.1f}pp, "
                  f"std={np.std(gaps):.1f}pp, "
                  f"min={np.min(gaps):.1f}pp, max={np.max(gaps):.1f}pp")
            # Trend
            if len(gaps) >= 3:
                x = np.arange(len(gaps))
                slope = np.polyfit(x, gaps, 1)[0]
                trend = "WIDENING" if slope > 0.5 else ("NARROWING" if slope < -0.5 else "STABLE")
                print(f"    Trend (high→low eps): {trend} (slope={slope:.2f}pp/step)")
        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 4: DOES MULTIPLICATIVITY HOLD?
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 4: MULTIPLICATIVITY CHECK")
    print("  If scales are independent, S_combined / (S_torus * S_L2) ≈ 1.0")
    print("  Deviation from 1.0 indicates scale coupling or breakdown")
    print("=" * 90)
    print()

    # We can check this indirectly: does the T75 suppression ratio (T75/B31)
    # stay constant across error rates? If the rotation contribution is
    # truly independent of error rate, the ratio should be flat.

    for n_nodes in CELL_SIZES:
        print(f"  {n_nodes}-node cell:")
        print(f"  {'eps_raw':>10} | {'T75/B31':>8} | {'consistent?':>12}")
        print("  " + "-" * 40)

        ratios = []
        for eps_raw in EPSILON_RAW_VALUES:
            r = all_results.get((n_nodes, eps_raw))
            if r is None: continue
            ratio = r['gap_supp']
            ratios.append(ratio)
            # Check if ratio is within 2x of mean so far
            mean_so_far = np.mean(ratios)
            consistent = "YES" if abs(ratio - mean_so_far) / mean_so_far < 0.5 else "DEVIATION"
            print(f"  {eps_raw:>10.0e} | {ratio:>6.1f}x | {consistent:>12}")

        if ratios:
            cv = np.std(ratios) / np.mean(ratios) if np.mean(ratios) > 0 else 0
            print(f"\n    Mean ratio: {np.mean(ratios):.1f}x, "
                  f"CV: {cv:.2f} ({'STABLE' if cv < 0.3 else 'VARIABLE'})")
        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 5: SUMMARY TABLE
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 5: SUMMARY — THE SCALE SEPARATION LANDSCAPE")
    print("=" * 90)
    print()

    print(f"  {'':>10} |  {'7-node':^30} | {'19-node':^30} | {'37-node':^30}")
    print(f"  {'eps_raw':>10} | {'B31':>6} {'T75':>6} {'S1*S2':>8} {'gap':>6} | "
          f"{'B31':>6} {'T75':>6} {'S1*S2':>8} {'gap':>6} | "
          f"{'B31':>6} {'T75':>6} {'S1*S2':>8} {'gap':>6}")
    print("  " + "-" * 100)

    for eps_raw in EPSILON_RAW_VALUES:
        row = f"  {eps_raw:>10.0e} |"
        for n_nodes in CELL_SIZES:
            r = all_results.get((n_nodes, eps_raw))
            if r:
                s_comb = S_TORUS_ANALYTICAL * r['t75_supp'] / (1 - 0.5)
                row += (f" {r['b31_supp']:>5.1f} {r['t75_supp']:>5.1f} "
                        f"{s_comb:>7.0f} {r['rotation_gap_det']*100:>5.1f} |")
            else:
                row += f" {'---':>5} {'---':>5} {'---':>7} {'---':>5} |"
        print(row)

    print()
    print("  Columns: B31=B31 suppression, T75=T75 L2 suppression,")
    print("           S1*S2=S_torus(20x)*S_L1L2(f=0.5), gap=detection gap (pp)")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # FINDINGS
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  FINDINGS")
    print("=" * 90)
    print()

    # Compute overall statistics
    all_gaps = []
    all_ratios = []
    for key, r in all_results.items():
        all_gaps.append(r['rotation_gap_det'] * 100)
        all_ratios.append(r['gap_supp'])

    print(f"  1. ROTATION GAP ACROSS ALL ERROR RATES:")
    print(f"     Mean: {np.mean(all_gaps):.1f} pp")
    print(f"     Range: {np.min(all_gaps):.1f} – {np.max(all_gaps):.1f} pp")
    print(f"     Std: {np.std(all_gaps):.1f} pp")
    print()

    print(f"  2. SUPPRESSION RATIO (T75/B31) ACROSS ALL ERROR RATES:")
    print(f"     Mean: {np.mean(all_ratios):.1f}x")
    print(f"     Range: {np.min(all_ratios):.1f} – {np.max(all_ratios):.1f}x")
    print(f"     CV: {np.std(all_ratios)/np.mean(all_ratios):.2f}")
    print()

    # Two-scale composition range
    for n_nodes in CELL_SIZES:
        supps = []
        for eps_raw in EPSILON_RAW_VALUES:
            r = all_results.get((n_nodes, eps_raw))
            if r:
                supps.append(S_TORUS_ANALYTICAL * r['t75_supp'] / (1 - 0.5))
        if supps:
            print(f"  3. COMBINED SUPPRESSION (S_torus=20x, f=0.5), {n_nodes}-node:")
            print(f"     Range: {np.min(supps):.0f}x – {np.max(supps):.0f}x")
            print(f"     Mean: {np.mean(supps):.0f}x")
            print(f"     At eps=1e-3: {supps[EPSILON_RAW_VALUES.index(1e-3)]:.0f}x")
            print()

    elapsed = time.time() - t_start
    print(f"  Runtime: {elapsed:.1f}s")
    print()
    print("=" * 90)

    sys.stdout = tee.terminal
    tee.file.close()
    print(f"\n  Output saved to: {output_path}")


if __name__ == '__main__':
    main()
