#!/usr/bin/env python3
"""
Merkabit U₀ Forward Pass + ZP-PPW on IBM Hardware
==================================================

Runs U₀ (12 ouroboros steps) without the dagger — the circuit cannot
be optimized away. Measures two things:

  1. Direct readout: measure q+ and q- after U₀
     Shows the output state distribution. Non-trivial superposition expected.

  2. ZP-PPW (π-Lock Parity Witness):
     Append ancilla parity map → CX(q+, anc), CX(q-, anc) → measure anc.
     ZP-PPW = Pr(anc=0) - Pr(anc=1)
     Tests whether forward/inverse spinors stay phase-locked after U₀.

Authors: Stenberg & Hetland, April 2026
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "zpmb"

# ─── Gate angles (E₆ geometry, zero free parameters) ──────────────────────────

T_CYCLE    = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES      = ['S', 'R', 'T', 'F', 'P']


def get_gate_angles(k: int) -> tuple[float, float, float]:
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


def _append_u0(qc: QuantumCircuit, q_fwd: int, q_inv: int) -> None:
    """Append U₀ (12 ouroboros steps) in-place on virtual qubits q_fwd, q_inv."""
    for k in range(T_CYCLE):
        p, rz, rx = get_gate_angles(k)
        # Merged P + symmetric Rz per qubit
        # q_fwd: Rz(-p) then Rz(rz) → Rz(rz - p)
        # q_inv: Rz(+p) then Rz(rz) → Rz(rz + p)
        qc.rz(rz - p, q_fwd)
        qc.rz(rz + p, q_inv)
        qc.rx(rx, q_fwd)
        qc.rx(rx, q_inv)


# ─── Ideal simulation ──────────────────────────────────────────────────────────

def simulate_u0():
    """Compute ideal output state after U₀ via statevector."""
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(2)
    _append_u0(qc, 0, 1)
    sv = Statevector(qc)
    probs = sv.probabilities_dict()
    # Also compute ZZ correlation ⟨Z₀ Z₁⟩ and ZP-PPW (parity)
    p00 = probs.get('00', 0.0)
    p01 = probs.get('01', 0.0)
    p10 = probs.get('10', 0.0)
    p11 = probs.get('11', 0.0)
    zz  = p00 - p01 - p10 + p11          # ⟨Z₀Z₁⟩
    parity_0 = p00 + p11                  # both same → ancilla stays 0
    parity_1 = p01 + p10                  # different → ancilla flips to 1
    zppw = parity_0 - parity_1
    return {
        "probs": probs,
        "zz_corr": zz,
        "zp_ppw": zppw,
        "state_vector": sv.data.tolist(),
    }


# ─── Circuit builders ──────────────────────────────────────────────────────────

def build_direct_circuit() -> QuantumCircuit:
    """U₀ on q[0], q[1] then measure both."""
    qr = QuantumRegister(2, 'q')
    cr = ClassicalRegister(2, 'c')
    qc = QuantumCircuit(qr, cr)
    _append_u0(qc, 0, 1)
    qc.measure(qr[0], cr[0])
    qc.measure(qr[1], cr[1])
    return qc


def build_zppw_circuit() -> QuantumCircuit:
    """
    U₀ on q[0], q[1], then ZZ parity witness via Hadamard test on ancilla q[2].

    Uses NATIVE gate directions (anc=72 controls q+=62, q-=81 on ibm_strasbourg):

      q[0] (q+):  U₀ ────────── CX ──────── CX ──────── measure
      q[1] (q-):  U₀ ──────────────── CX ── CX ──────── measure
      q[2] (anc): |0⟩ ── H ─── CX ── CX ── H ─────────── measure

    Measures ⟨Z_q+ Z_q-⟩ = Pr(same parity) - Pr(different parity) = ZP-PPW.
    Identical to CX(q+→anc); CX(q-→anc) but uses native 72→62, 72→81 directions.
    """
    qr = QuantumRegister(3, 'q')
    cr = ClassicalRegister(3, 'c')
    qc = QuantumCircuit(qr, cr)
    _append_u0(qc, 0, 1)
    qc.h(qr[2])            # anc → |+⟩
    qc.cx(qr[2], qr[0])   # native: anc controls q+  (72→62)
    qc.cx(qr[2], qr[1])   # native: anc controls q-  (72→81)
    qc.h(qr[2])            # close interference
    qc.measure(qr[0], cr[0])
    qc.measure(qr[1], cr[1])
    qc.measure(qr[2], cr[2])
    return qc


# ─── Hardware run ──────────────────────────────────────────────────────────────

def run_circuit(backend, qc: QuantumCircuit, layout: list,
                shots: int, label: str) -> dict:
    print(f"\n── {label} ─────────────────────────────────")
    print(f"   Abstract depth : {qc.depth()}   gates : {qc.size()}")

    pm = generate_preset_pass_manager(
        optimization_level=1,
        backend=backend,
        initial_layout=layout,
    )
    transpiled = pm.run(qc)
    print(f"   Transpiled depth: {transpiled.depth()}")

    sampler = Sampler(backend)
    job     = sampler.run([transpiled], shots=shots)
    print(f"   Job ID: {job.job_id()}  — waiting …")
    result  = job.result()

    counts = result[0].data.c.get_counts()
    total  = sum(counts.values())
    print(f"   Shots: {total}")
    print(f"   Counts: {dict(sorted(counts.items()))}")

    return {
        "label":            label,
        "shots":            total,
        "counts":           counts,
        "job_id":           job.job_id(),
        "abstract_depth":   qc.depth(),
        "transpiled_depth": transpiled.depth(),
        "backend":          backend.name,
        "layout":           layout,
        "timestamp":        datetime.now().isoformat(),
    }


# ─── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Merkabit U₀ forward pass + ZP-PPW")
    parser.add_argument("--backend",  default="ibm_strasbourg")
    parser.add_argument("--token",    default=None)
    parser.add_argument("--shots",    type=int, default=8192)
    parser.add_argument("--q-fwd",   type=int, default=62,
                        help="Physical qubit: forward spinor q+ (default 62)")
    parser.add_argument("--q-inv",   type=int, default=81,
                        help="Physical qubit: inverse spinor q- (default 81)")
    parser.add_argument("--q-anc",   type=int, default=72,
                        help="Physical qubit: ancilla for ZP-PPW (default 72)")
    parser.add_argument("--sim-only", action="store_true")
    args = parser.parse_args()

    # ── Ideal simulation ──────────────────────────────────────────────────────
    sim = simulate_u0()
    print(f"\n── Ideal U₀ output (noiseless statevector) ────────────────")
    for outcome, p in sorted(sim["probs"].items(), key=lambda x: -x[1]):
        if p > 1e-6:
            print(f"   |{outcome}⟩  p = {p:.6f}")
    print(f"   ⟨Z₀Z₁⟩    = {sim['zz_corr']:+.6f}")
    print(f"   ZP-PPW     = {sim['zp_ppw']:+.6f}  (ideal; +1 = fully locked, -1 = anti-locked)")

    if args.sim_only:
        return

    # ── Hardware ──────────────────────────────────────────────────────────────
    token = args.token or os.environ.get("IBM_QUANTUM_TOKEN")
    service = QiskitRuntimeService(
        channel="ibm_quantum_platform", token=token
    ) if token else QiskitRuntimeService(channel="ibm_quantum_platform")

    backend = service.backend(args.backend)
    print(f"\nBackend : {backend.name}  ({backend.num_qubits} qubits)")
    print(f"q+ = {args.q_fwd}   q- = {args.q_inv}   anc = {args.q_anc}")
    print(f"Shots   : {args.shots}")

    results = []

    # Run 1: direct U₀ readout
    qc_direct = build_direct_circuit()
    r_direct  = run_circuit(
        backend, qc_direct, [args.q_fwd, args.q_inv],
        args.shots, "U₀ direct (measure q+, q-)"
    )
    results.append(r_direct)

    # Run 2: ZP-PPW (parity witness)
    qc_zppw = build_zppw_circuit()
    r_zppw  = run_circuit(
        backend, qc_zppw, [args.q_fwd, args.q_inv, args.q_anc],
        args.shots, "ZP-PPW (parity witness, ancilla=72)"
    )
    results.append(r_zppw)

    # ── Analysis ──────────────────────────────────────────────────────────────
    # Direct: state distribution
    direct_counts = r_direct["counts"]
    total         = r_direct["shots"]
    probs_hw      = {k: v / total for k, v in direct_counts.items()}

    # ZP-PPW from hardware (ancilla = bit 0 in 3-bit string)
    zppw_counts = r_zppw["counts"]
    total_z     = r_zppw["shots"]
    p_anc0 = sum(v for k, v in zppw_counts.items() if k[0] == '0') / total_z
    p_anc1 = sum(v for k, v in zppw_counts.items() if k[0] == '1') / total_z
    zppw_hw = p_anc0 - p_anc1

    print(f"\n══ Results ═════════════════════════════════════════════")
    print(f"  Ideal ZP-PPW  : {sim['zp_ppw']:+.4f}")
    print(f"  Hardware ZP-PPW: {zppw_hw:+.4f}")
    print(f"  Ideal ⟨Z₀Z₁⟩  : {sim['zz_corr']:+.4f}")
    print()
    print("  Direct U₀ output state (hardware):")
    for k, p in sorted(probs_hw.items(), key=lambda x: -x[1]):
        ideal_p = sim["probs"].get(k, 0.0)
        print(f"    |{k}⟩  hw={p:.4f}  ideal={ideal_p:.4f}")

    # ── Save ──────────────────────────────────────────────────────────────────
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts       = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = RESULTS_DIR / f"u0_zppw_{backend.name}_{ts}.json"
    output   = {
        "results":  results,
        "analysis": {
            "ideal_zppw":   sim["zp_ppw"],
            "hardware_zppw": zppw_hw,
            "ideal_zz_corr": sim["zz_corr"],
            "ideal_probs":  sim["probs"],
            "hardware_probs": probs_hw,
        },
        "backend": backend.name,
    }
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults → {out_path}")


if __name__ == "__main__":
    main()
