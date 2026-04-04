#!/usr/bin/env python3
"""
ASYMPTOTE MAPPING: 7 → 19 → 37 → 61 → 91 nodes
==================================================

Maps the suppression asymptote across all five Eisenstein cell sizes.
Key questions:
  1. Does T75 suppression saturate or keep climbing?
  2. Does the rotation gap (T75 - B31) widen, narrow, or stay flat?
  3. Where does the asymptote sit?

Cell sizes: r=1(7), r=2(19), r=3(37), r=4(61), r=5(91)
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys
import os

GATES = ['R', 'T', 'P', 'F', 'S']
NUM_GATES = 5
RANDOM_SEED = 42
B31_ACTIVE_GATES = {0, 1, 4}  # R, T, S

# Focus on key error rates — enough to see the structure
EPSILON_VALUES = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
CELL_CONFIGS = [(1,7), (2,19), (3,37), (4,61), (5,91)]
MC_ASSIGNMENTS = 15
MC_TRIALS_BASE = 40_000


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
        # Boundary vs interior
        self.boundary_nodes = [i for i in range(self.num_nodes) if self.coordination[i] < 6]
        self.interior_nodes = [i for i in range(self.num_nodes) if self.coordination[i] == 6]
        self.boundary_frac = len(self.boundary_nodes) / self.num_nodes


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


def detect_and_correct_t75(cell, assignment, error_node, error_gate, tau=5):
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


def find_valid_assignments(cell, rng, count, max_attempts=200_000):
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
        good = True
        for (i,j) in cell.edges:
            if assignment[i] % NUM_GATES == assignment[j] % NUM_GATES:
                good = False; break
        if good:
            valid.append(tuple(assignment))
    return valid


def mc_sweep_point(cell, assignments, eps_raw, seed):
    rng = np.random.default_rng(seed)
    b31_gates = list(B31_ACTIVE_GATES)

    if eps_raw >= 0.01: num_trials = MC_TRIALS_BASE
    elif eps_raw >= 0.001: num_trials = MC_TRIALS_BASE * 2
    else: num_trials = MC_TRIALS_BASE * 3

    trials_per = max(1, num_trials // len(assignments))

    b31_total = 0; b31_inj = 0; b31_det = 0; b31_corr = 0; b31_unc = 0
    t75_total = 0; t75_inj = 0; t75_det = 0; t75_corr = 0; t75_unc = 0

    for assignment in assignments:
        for trial in range(trials_per):
            for node in range(cell.num_nodes):
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


def main():
    t_start = time.time()
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'asymptote_mapping_output.txt')

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
    print("  ASYMPTOTE MAPPING: 7 -> 19 -> 37 -> 61 -> 91 nodes")
    print("  Does suppression saturate? Where is the asymptote?")
    print("=" * 90)
    print()

    # ── Build all cells ──
    cells = {}
    cell_assignments = {}
    for radius, expected_n in CELL_CONFIGS:
        cell = EisensteinCell(radius)
        assert cell.num_nodes == expected_n, f"r={radius}: got {cell.num_nodes}, expected {expected_n}"
        cells[expected_n] = cell
        rng = np.random.default_rng(RANDOM_SEED + expected_n)
        assignments = find_valid_assignments(cell, rng, MC_ASSIGNMENTS)
        cell_assignments[expected_n] = assignments
        print(f"  {expected_n:>3}-node (r={radius}): {len(assignments):>3} assignments, "
              f"{cell.num_nodes - len(cell.boundary_nodes)} interior, "
              f"{len(cell.boundary_nodes)} boundary ({cell.boundary_frac*100:.0f}%), "
              f"mean coord={np.mean(cell.coordination):.2f}")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 1: STRUCTURAL PROPERTIES vs CELL SIZE
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 1: STRUCTURAL PROPERTIES")
    print("=" * 90)
    print()
    print(f"  {'Nodes':>5} {'r':>3} {'Edges':>6} {'Boundary%':>10} {'Interior%':>10} "
          f"{'MeanCoord':>10} {'MinCoord':>9} {'Chiral+':>8} {'Chiral0':>8} {'Chiral-':>8}")
    print("  " + "-" * 90)

    for radius, n in CELL_CONFIGS:
        cell = cells[n]
        chiral_counts = Counter(cell.chirality)
        print(f"  {n:>5} {radius:>3} {len(cell.edges):>6} "
              f"{cell.boundary_frac*100:>8.1f}% {(1-cell.boundary_frac)*100:>8.1f}% "
              f"{np.mean(cell.coordination):>10.2f} {min(cell.coordination):>9} "
              f"{chiral_counts.get(1,0):>8} {chiral_counts.get(0,0):>8} {chiral_counts.get(-1,0):>8}")

    print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 2: FULL SWEEP
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 2: B31 vs T75 ACROSS ALL CELL SIZES")
    print("=" * 90)
    print()

    all_results = {}

    for radius, n_nodes in CELL_CONFIGS:
        cell = cells[n_nodes]
        assignments = cell_assignments[n_nodes]

        print(f"  {n_nodes}-node cell (r={radius}):")
        print(f"  {'eps':>10} | {'B31_det':>7} {'B31_corr':>8} {'B31_supp':>8} | "
              f"{'T75_det':>7} {'T75_corr':>8} {'T75_supp':>8} | "
              f"{'det_gap':>7} {'supp_rat':>8}")
        print("  " + "-" * 88)

        for eps_raw in EPSILON_VALUES:
            seed = RANDOM_SEED + hash((n_nodes, eps_raw)) % 100000
            t0 = time.time()
            r = mc_sweep_point(cell, assignments, eps_raw, seed)
            dt = time.time() - t0
            all_results[(n_nodes, eps_raw)] = r

            print(f"  {eps_raw:>10.0e} | "
                  f"{r['b31_det']*100:>5.1f}% {r['b31_corr']*100:>6.1f}% {r['b31_supp']:>6.1f}x | "
                  f"{r['t75_det']*100:>5.1f}% {r['t75_corr']*100:>6.1f}% {r['t75_supp']:>6.1f}x | "
                  f"{r['rotation_gap_det']*100:>5.1f}pp {r['gap_supp']:>6.1f}x  ({dt:.1f}s)")

        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 3: ASYMPTOTE ANALYSIS — SUPPRESSION vs CELL SIZE
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 3: SUPPRESSION vs CELL SIZE (the asymptote)")
    print("=" * 90)
    print()

    node_sizes = [n for _, n in CELL_CONFIGS]

    for eps_raw in EPSILON_VALUES:
        b31_supps = []
        t75_supps = []
        gaps = []
        ratios = []

        for n in node_sizes:
            r = all_results.get((n, eps_raw))
            if r:
                b31_supps.append(r['b31_supp'])
                t75_supps.append(r['t75_supp'])
                gaps.append(r['rotation_gap_det'] * 100)
                ratios.append(r['gap_supp'])

        print(f"  eps = {eps_raw:.0e}:")
        print(f"    {'Nodes':>5} | {'B31_supp':>9} | {'T75_supp':>9} | {'T75/B31':>8} | {'det_gap':>8}")
        for i, n in enumerate(node_sizes):
            print(f"    {n:>5} | {b31_supps[i]:>7.1f}x | {t75_supps[i]:>7.1f}x | "
                  f"{ratios[i]:>6.1f}x | {gaps[i]:>6.1f}pp")

        # Fit: does suppression look like a + b/n (saturation) or a*log(n)?
        if len(t75_supps) >= 3:
            ns = np.array(node_sizes, dtype=float)
            s_b31 = np.array(b31_supps)
            s_t75 = np.array(t75_supps)

            # Linear fit in 1/n (saturation model: S = a + b/n)
            x_inv = 1.0 / ns
            if np.std(s_t75) > 0:
                p_t75 = np.polyfit(x_inv, s_t75, 1)
                t75_asymptote = p_t75[1]  # value at 1/n -> 0

                p_b31 = np.polyfit(x_inv, s_b31, 1)
                b31_asymptote = p_b31[1]

                # Log fit (S = a + b*log(n))
                log_ns = np.log(ns)
                p_t75_log = np.polyfit(log_ns, s_t75, 1)
                p_b31_log = np.polyfit(log_ns, s_b31, 1)

                # R^2 for both models
                def r_squared(y, y_pred):
                    ss_res = np.sum((y - y_pred)**2)
                    ss_tot = np.sum((y - np.mean(y))**2)
                    return 1 - ss_res / ss_tot if ss_tot > 0 else 0

                r2_sat_t75 = r_squared(s_t75, np.polyval(p_t75, x_inv))
                r2_log_t75 = r_squared(s_t75, np.polyval(p_t75_log, log_ns))
                r2_sat_b31 = r_squared(s_b31, np.polyval(p_b31, x_inv))
                r2_log_b31 = r_squared(s_b31, np.polyval(p_b31_log, log_ns))

                print(f"\n    T75 saturation model (a + b/n): asymptote = {t75_asymptote:.1f}x, R2 = {r2_sat_t75:.3f}")
                print(f"    T75 log model (a + b*ln n):     slope = {p_t75_log[0]:.2f}x/ln(n), R2 = {r2_log_t75:.3f}")
                print(f"    B31 saturation model (a + b/n): asymptote = {b31_asymptote:.1f}x, R2 = {r2_sat_b31:.3f}")
                print(f"    B31 log model (a + b*ln n):     slope = {p_b31_log[0]:.2f}x/ln(n), R2 = {r2_log_b31:.3f}")

                if r2_sat_t75 > r2_log_t75:
                    print(f"    --> SATURATION model fits T75 better (asymptote ~ {t75_asymptote:.0f}x)")
                else:
                    print(f"    --> LOG model fits T75 better (still growing)")

        # Gap trend
        if len(gaps) >= 3:
            gap_slope = np.polyfit(np.arange(len(gaps)), gaps, 1)[0]
            print(f"    Detection gap trend: {'WIDENING' if gap_slope > 0.5 else 'NARROWING' if gap_slope < -0.5 else 'FLAT'} "
                  f"(slope={gap_slope:.2f}pp/step)")
        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 4: COMBINED TWO-SCALE COMPOSITION AT ALL SIZES
    # ══════════════════════════════════════════════════════════════════════
    S_TORUS = 35.0  # From eisenstein_torus_simulation.py direct measurement

    print("=" * 90)
    print(f"  PART 4: TWO-SCALE COMPOSITION (S_torus = {S_TORUS:.0f}x from direct simulation)")
    print("=" * 90)
    print()

    print(f"  {'Nodes':>5} {'eps':>10} | {'S_L2(MC)':>9} | {'S_L1+L2':>9} | {'S_torus*L2':>11} | {'S_torus*L1L2':>13}")
    print("  " + "-" * 70)

    for eps_raw in [1e-2, 1e-3, 1e-4]:
        for n in node_sizes:
            r = all_results.get((n, eps_raw))
            if r:
                s_L2 = r['t75_supp']
                s_L1L2 = s_L2 / (1 - 0.5)  # f_sym=0.5 conservative
                s_comb_L2 = S_TORUS * s_L2
                s_comb_L1L2 = S_TORUS * s_L1L2
                print(f"  {n:>5} {eps_raw:>10.0e} | {s_L2:>7.1f}x | {s_L1L2:>7.1f}x | "
                      f"{s_comb_L2:>9.0f}x | {s_comb_L1L2:>11.0f}x")
        print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 5: ROTATION GAP STABILITY ACROSS ALL SIZES
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 5: IS THE ROTATION GAP FLAT ACROSS CELL SIZES?")
    print("=" * 90)
    print()

    print(f"  {'eps':>10} |", end="")
    for n in node_sizes:
        print(f" {n:>5}n", end="")
    print(" | mean    std     CV")
    print("  " + "-" * 70)

    for eps_raw in EPSILON_VALUES:
        print(f"  {eps_raw:>10.0e} |", end="")
        gaps = []
        for n in node_sizes:
            r = all_results.get((n, eps_raw))
            if r:
                gap = r['rotation_gap_det'] * 100
                gaps.append(gap)
                print(f" {gap:>4.1f}pp", end="")
            else:
                print(f"  {'---':>5}", end="")
        if gaps:
            mean_g = np.mean(gaps)
            std_g = np.std(gaps)
            cv = std_g / mean_g if mean_g > 0 else 0
            print(f" | {mean_g:>5.1f}  {std_g:>5.2f}  {cv:>6.3f}")
        else:
            print()

    print()

    # ══════════════════════════════════════════════════════════════════════
    # PART 6: BOUNDARY FRACTION vs SUPPRESSION
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  PART 6: BOUNDARY FRACTION vs SUPPRESSION")
    print("  Does the boundary fraction explain the suppression curve?")
    print("=" * 90)
    print()

    for eps_raw in [1e-2, 1e-3]:
        print(f"  eps = {eps_raw:.0e}:")
        print(f"  {'Nodes':>5} {'Bound%':>7} {'Inter%':>7} | {'B31_supp':>8} {'T75_supp':>8} {'T75/B31':>8}")
        print("  " + "-" * 55)
        for n in node_sizes:
            cell = cells[n]
            r = all_results.get((n, eps_raw))
            if r:
                print(f"  {n:>5} {cell.boundary_frac*100:>5.1f}% {(1-cell.boundary_frac)*100:>5.1f}% | "
                      f"{r['b31_supp']:>6.1f}x {r['t75_supp']:>6.1f}x {r['gap_supp']:>6.1f}x")

        # Correlation: boundary_frac vs suppression
        bfracs = np.array([cells[n].boundary_frac for n in node_sizes])
        t75s = np.array([all_results[(n, eps_raw)]['t75_supp'] for n in node_sizes])
        b31s = np.array([all_results[(n, eps_raw)]['b31_supp'] for n in node_sizes])
        corr_t75 = np.corrcoef(bfracs, t75s)[0,1]
        corr_b31 = np.corrcoef(bfracs, b31s)[0,1]
        print(f"\n    Correlation(boundary%, T75_supp) = {corr_t75:.3f}")
        print(f"    Correlation(boundary%, B31_supp) = {corr_b31:.3f}")
        print(f"    {'Strong negative = boundary fraction drives suppression' if corr_t75 < -0.8 else 'Weak correlation = other factors dominate'}")
        print()

    # ══════════════════════════════════════════════════════════════════════
    # SUMMARY
    # ══════════════════════════════════════════════════════════════════════
    print("=" * 90)
    print("  SUMMARY: THE ASYMPTOTE")
    print("=" * 90)
    print()

    # At eps=1e-3, what are the numbers?
    eps_ref = 1e-3
    print(f"  At eps = {eps_ref:.0e}:")
    print()

    b31_at_ref = [all_results[(n, eps_ref)]['b31_supp'] for n in node_sizes]
    t75_at_ref = [all_results[(n, eps_ref)]['t75_supp'] for n in node_sizes]
    gap_at_ref = [all_results[(n, eps_ref)]['rotation_gap_det']*100 for n in node_sizes]

    print(f"    B31 suppression:  {' -> '.join(f'{s:.1f}x' for s in b31_at_ref)}")
    print(f"    T75 suppression:  {' -> '.join(f'{s:.1f}x' for s in t75_at_ref)}")
    print(f"    Detection gap:    {' -> '.join(f'{g:.1f}pp' for g in gap_at_ref)}")
    print()

    # Asymptote estimates
    ns = np.array(node_sizes, dtype=float)
    x_inv = 1.0 / ns
    p_t75 = np.polyfit(x_inv, t75_at_ref, 1)
    p_b31 = np.polyfit(x_inv, b31_at_ref, 1)

    print(f"    T75 asymptote (n->inf): {p_t75[1]:.1f}x")
    print(f"    B31 asymptote (n->inf): {p_b31[1]:.1f}x")
    print(f"    Gap asymptote: {p_t75[1] - p_b31[1]:.1f}x (ratio)")
    print()

    # Two-scale at asymptote
    print(f"    Two-scale at asymptote (S_torus=35x, f_sym=0.5):")
    s_L2_asym = p_t75[1]
    s_L1L2_asym = s_L2_asym / (1 - 0.5)
    s_combined = S_TORUS * s_L1L2_asym
    print(f"      S_L2 = {s_L2_asym:.1f}x")
    print(f"      S_L1+L2 = {s_L1L2_asym:.1f}x")
    print(f"      S_total = S_torus * S_L1L2 = {S_TORUS:.0f} * {s_L1L2_asym:.1f} = {s_combined:.0f}x")
    print()

    elapsed = time.time() - t_start
    print(f"  Runtime: {elapsed:.1f}s ({elapsed/60:.1f}m)")
    print()
    print("=" * 90)

    sys.stdout = tee.terminal
    tee.file.close()
    print(f"\n  Output saved to: {output_path}")


if __name__ == '__main__':
    main()
