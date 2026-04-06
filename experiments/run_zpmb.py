#!/usr/bin/env python3
"""
Zero-Point Merkabit Benchmark (ZPMB) — ZP-ORF on IBM Hardware
==============================================================

The P gate breakthrough (April 6 2026):
  P^(0)(phi) = Rz(-phi) on q+ (forward spinor)
               Rz(+phi) on q- (inverse spinor)
compiles to two native IBM Rz gates. Zero entangling gates in U_0.

Runs two experiments:
  1. ZP-ORF (paired):   |00> -> U_0 -> U_0† -> measure   — real merkabit
  2. ZP-ORF (unpaired): same but Rz(-phi) on q- removed  — control

ZP-ORF = Pr(q+=0, q-=0).
Suppression ratio = paired / unpaired → Level 1 noise cancellation.

Authors: Stenberg & Hetland, April 2026
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "zpmb"

# ─── Gate angles (E₆ geometry, zero free parameters) ──────────────────────────
# Reproduced verbatim from merkabit_verification.py (Stenberg & Claude, Feb 2026)

T_CYCLE = 12
STEP_PHASE = 2 * np.pi / T_CYCLE   # π/6
GATES = ['S', 'R', 'T', 'F', 'P']


def get_gate_angles(k: int) -> tuple[float, float, float]:
    """Return (p_angle, rz_angle, rx_angle) for ouroboros step k."""
    absent     = k % 5
    gate_label = GATES[absent]
    p_angle    = STEP_PHASE
    sym_base   = STEP_PHASE / 3
    omega_k    = 2 * np.pi * k / T_CYCLE
    rx_angle   = sym_base * (1.0 + 0.5 * np.cos(omega_k))
    rz_angle   = sym_base * (1.0 + 0.5 * np.cos(omega_k + 2 * np.pi / 3))
    if gate_label == 'S':
        rz_angle *= 0.4;  rx_angle *= 1.3
    elif gate_label == 'R':
        rx_angle *= 0.4;  rz_angle *= 1.3
    elif gate_label == 'T':
        rx_angle *= 0.7;  rz_angle *= 0.7
    elif gate_label == 'P':
        p_angle  *= 0.6;  rx_angle *= 1.8;  rz_angle *= 1.5
    return p_angle, rz_angle, rx_angle


# ─── Circuit builders ──────────────────────────────────────────────────────────
#
# Qiskit Rz(λ) = diag(exp(-iλ/2), exp(+iλ/2))
# From merkabit_verification.py step_unitary():
#   Forward spinor P:  Pf = diag(exp(+ip/2), exp(-ip/2)) = Rz(-p)
#   Inverse spinor P:  Pi = diag(exp(-ip/2), exp(+ip/2)) = Rz(+p)
#   Symmetric Rz:      Rz(rz_angle) on both
# Merged per qubit (Rz gates commute):
#   q_fwd: Rz(-p_angle) then Rz(rz_angle) = Rz(rz_angle - p_angle)
#   q_inv: Rz(+p_angle) then Rz(rz_angle) = Rz(rz_angle + p_angle)

def _append_u0(qc: QuantumCircuit, q_fwd: int, q_inv: int, paired: bool) -> None:
    """Append U_0 (12 ouroboros steps) in-place. Steps k=0..11."""
    for k in range(T_CYCLE):
        p, rz, rx = get_gate_angles(k)
        qc.rz(rz - p, q_fwd)
        qc.rz(rz + p if paired else rz, q_inv)
        qc.rx(rx, q_fwd)
        qc.rx(rx, q_inv)


def _append_u0_dagger(qc: QuantumCircuit, q_fwd: int, q_inv: int, paired: bool) -> None:
    """Append U_0† (reverse steps, inverted gate order). Steps k=11..0."""
    for k in range(T_CYCLE - 1, -1, -1):
        p, rz, rx = get_gate_angles(k)
        # Step dagger: Rx†, then merged Rz†/P†
        # Rx†(theta) = Rx(-theta); Rz†(-p)⊗Rz†(rz) = Rz(p)⊗Rz(-rz)
        # Merged q_fwd: Rz(p - rz);  q_inv: Rz(-p - rz)
        qc.rx(-rx, q_fwd)
        qc.rx(-rx, q_inv)
        qc.rz(p - rz, q_fwd)
        qc.rz((-p - rz) if paired else (-rz), q_inv)


def build_zporf_circuit(paired: bool) -> QuantumCircuit:
    """
    ZP-ORF: |00> -> U_0 -> U_0† -> measure.

    Virtual qubit 0 = forward spinor (q+)
    Virtual qubit 1 = inverse spinor (q-)
    Physical layout set at transpilation time via initial_layout.
    """
    qr = QuantumRegister(2, 'q')
    cr = ClassicalRegister(2, 'c')
    qc = QuantumCircuit(qr, cr)
    _append_u0(qc, 0, 1, paired)
    _append_u0_dagger(qc, 0, 1, paired)
    qc.measure(qr[0], cr[0])
    qc.measure(qr[1], cr[1])
    return qc


# ─── Simulator baseline (sanity check) ────────────────────────────────────────

def simulate_zporf(paired: bool) -> float:
    """Ideal (noiseless) ZP-ORF via statevector simulation."""
    from qiskit_aer import AerSimulator
    qc = build_zporf_circuit(paired)
    qc_sv = qc.copy()
    # Use statevector via AerSimulator
    sim = AerSimulator(method='statevector')
    from qiskit import transpile
    t = transpile(qc_sv, sim)
    result = sim.run(t, shots=1).result()
    # For noiseless: measure manually via statevector
    from qiskit.quantum_info import Statevector
    qc_no_measure = build_zporf_circuit(paired).remove_final_measurements(inplace=False)
    sv = Statevector(qc_no_measure)
    # P(|00>) = |<00|sv>|^2  — qubit 0 in cr[0], qubit 1 in cr[1]
    probs = sv.probabilities_dict()
    return probs.get('00', 0.0)


# ─── Hardware run ──────────────────────────────────────────────────────────────

def run_zporf_hardware(backend, shots: int, q_fwd: int, q_inv: int,
                       paired: bool, label: str) -> dict:
    qc = build_zporf_circuit(paired)
    print(f"\n── {label} ──────────────────────────────────────")
    print(f"   Abstract depth : {qc.depth()}   gates : {qc.size()}")

    pm = generate_preset_pass_manager(
        optimization_level=1,
        backend=backend,
        initial_layout=[q_fwd, q_inv],
    )
    transpiled = pm.run(qc)
    print(f"   Transpiled depth: {transpiled.depth()}")

    sampler = Sampler(backend)
    job = sampler.run([transpiled], shots=shots)
    print(f"   Job ID: {job.job_id()}  — waiting …")
    result = job.result()

    counts = result[0].data.c.get_counts()
    total  = sum(counts.values())
    p00    = counts.get('00', 0) / total

    print(f"   ZP-ORF = {p00:.4f}   ({total} shots)")
    print(f"   Counts: {dict(sorted(counts.items()))}")

    return {
        "label":            label,
        "paired":           paired,
        "q_fwd":            q_fwd,
        "q_inv":            q_inv,
        "shots":            total,
        "zp_orf":           p00,
        "counts":           counts,
        "job_id":           job.job_id(),
        "abstract_depth":   qc.depth(),
        "transpiled_depth": transpiled.depth(),
        "backend":          backend.name,
        "timestamp":        datetime.now().isoformat(),
    }


# ─── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="ZPMB ZP-ORF — merkabit hardware run")
    parser.add_argument("--backend",  default="ibm_strasbourg")
    parser.add_argument("--token",    default=None)
    parser.add_argument("--shots",    type=int, default=8192)
    parser.add_argument("--q-fwd",   type=int, default=62,
                        help="Physical qubit: forward spinor q+ (default 62)")
    parser.add_argument("--q-inv",   type=int, default=81,
                        help="Physical qubit: inverse spinor q- (default 81)")
    parser.add_argument("--sim-only", action="store_true",
                        help="Run ideal simulation only (no hardware)")
    args = parser.parse_args()

    # ── Simulation sanity check ───────────────────────────────────────────────
    try:
        p00_paired   = simulate_zporf(paired=True)
        p00_unpaired = simulate_zporf(paired=False)
        print(f"\n── Ideal simulation (noiseless) ──────────────────────")
        print(f"   ZP-ORF paired   = {p00_paired:.6f}  (expect ~1.0)")
        print(f"   ZP-ORF unpaired = {p00_unpaired:.6f}  (expect ~1.0)")
    except Exception as e:
        print(f"   Simulation skipped ({e})")

    if args.sim_only:
        return

    # ── Hardware ──────────────────────────────────────────────────────────────
    token = args.token or os.environ.get("IBM_QUANTUM_TOKEN")
    if token:
        service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    else:
        service = QiskitRuntimeService(channel="ibm_quantum_platform")

    backend = service.backend(args.backend)
    print(f"\nBackend: {backend.name}  ({backend.num_qubits} qubits)")
    print(f"q+ = {args.q_fwd}  (forward spinor)")
    print(f"q- = {args.q_inv}  (inverse spinor)")
    print(f"Shots: {args.shots}")

    results = []

    r_paired = run_zporf_hardware(
        backend, args.shots, args.q_fwd, args.q_inv,
        paired=True,  label="ZP-ORF — paired (merkabit P gate active)"
    )
    results.append(r_paired)

    r_control = run_zporf_hardware(
        backend, args.shots, args.q_fwd, args.q_inv,
        paired=False, label="ZP-ORF — unpaired (control, P gate removed)"
    )
    results.append(r_control)

    # ── Summary ───────────────────────────────────────────────────────────────
    zorf_p = r_paired["zp_orf"]
    zorf_u = r_control["zp_orf"]
    ratio  = zorf_p / zorf_u if zorf_u > 0 else float("inf")

    print(f"\n══ ZPMB Summary ═══════════════════════════════════════")
    print(f"  ZP-ORF paired   : {zorf_p:.4f}   {'PASS ✓' if zorf_p > 0.90 else 'BELOW 0.90'}")
    print(f"  ZP-ORF unpaired : {zorf_u:.4f}")
    print(f"  Ratio (Level 1) : {ratio:.3f}x")
    print(f"  Paper P1 Berry phase separation: 7.13 rad over 12 steps")
    print(f"  P gate: Rz(-p) on q+, Rz(+p) on q-  →  2 native Rz gates, 0 ECR")

    # ── Save ──────────────────────────────────────────────────────────────────
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts       = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = RESULTS_DIR / f"zpmb_zporf_{backend.name}_{ts}.json"
    with open(out_path, "w") as f:
        json.dump({"results": results, "ratio": ratio, "backend": backend.name}, f, indent=2)
    print(f"\nResults → {out_path}")


if __name__ == "__main__":
    main()
