#!/usr/bin/env python3
"""
ZP-GPW: Zero-Point Geometric Phase Witness on IBM Hardware
===========================================================

Measures the geometric phase δ accumulated by the merkabit ouroboros
cycle via a Hadamard test on the ancilla.

Circuit (X-basis run):
  anc (q=72): |0⟩ ── H ── [ctrl-U₀_n] ── H ── measure
  q+  (q=62): |0⟩ ── [U₀_n when anc=1] ─────── (no measure needed)
  q-  (q=81): |0⟩ ── [U₀_n when anc=1] ─────── (no measure needed)

  ⟨X_anc⟩ = Re(⟨00|U₀_n|00⟩) = |M₀₀| cos(δ_n)

Circuit (Y-basis run):
  Replace final H(anc) with Sdg(anc) → H(anc).
  ⟨Y_anc⟩ = Im(⟨00|U₀_n|00⟩) = |M₀₀| sin(δ_n)

  δ_n = atan2(⟨Y⟩, ⟨X⟩)

All CX gates use native ibm_strasbourg directions (anc=72 → q+=62, q-=81).

Controlled-Rz(θ) decomposition (native CX(anc,target)):
  Rz(θ/2, target); CX(anc, target); Rz(-θ/2, target); CX(anc, target)

Controlled-Rx(θ) decomposition (H-conjugated):
  H(target); Rz(θ/2, target); CX(anc, target);
  Rz(-θ/2, target); CX(anc, target); H(target)

Expected values (noiseless simulation):
  n=4: ⟨X⟩=+0.856  ⟨Y⟩=−0.310  δ=−19.9°  (32 CX)
  n=6: ⟨X⟩=+0.719  ⟨Y⟩=−0.497  δ=−34.6°  (48 CX)  ← default
  n=8: ⟨X⟩=+0.446  ⟨Y⟩=−0.762  δ=−59.7°  (64 CX)

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

T_CYCLE    = 12
STEP_PHASE = 2 * np.pi / T_CYCLE
GATES      = ['S', 'R', 'T', 'F', 'P']


def get_gate_angles(k: int) -> tuple[float, float, float]:
    absent     = k % 5
    gate_label = GATES[absent]
    p          = STEP_PHASE
    sym        = STEP_PHASE / 3
    w          = 2 * np.pi * k / T_CYCLE
    rx         = sym * (1.0 + 0.5 * np.cos(w))
    rz         = sym * (1.0 + 0.5 * np.cos(w + 2 * np.pi / 3))
    if gate_label == 'S': rz *= 0.4;  rx *= 1.3
    elif gate_label == 'R': rx *= 0.4; rz *= 1.3
    elif gate_label == 'T': rx *= 0.7; rz *= 0.7
    elif gate_label == 'P': p  *= 0.6; rx *= 1.8; rz *= 1.5
    return p, rz, rx


# ─── Ideal signal computation ──────────────────────────────────────────────────

def expected_signal(n_steps: int) -> dict:
    """Compute ⟨00|U₀_n|00⟩ from statevector simulation."""
    state0 = np.array([1, 0, 0, 0], dtype=complex)
    U = np.eye(4, dtype=complex)
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k)
        Rz_f = np.diag([np.exp(-1j*(rz-p)/2), np.exp(1j*(rz-p)/2)])
        Rz_i = np.diag([np.exp(-1j*(rz+p)/2), np.exp(1j*(rz+p)/2)])
        c = np.cos(rx/2); s = -1j*np.sin(rx/2)
        Rx = np.array([[c, s], [s, c]])
        U = np.kron(Rx @ Rz_f, Rx @ Rz_i) @ U
    m00     = state0.conj() @ U @ state0
    delta   = float(np.angle(m00))
    return {
        "n_steps":   n_steps,
        "m00_re":    float(m00.real),
        "m00_im":    float(m00.imag),
        "magnitude": float(abs(m00)),
        "delta_rad": delta,
        "delta_deg": float(np.degrees(delta)),
    }


# ─── Controlled-U₀ circuit builder ────────────────────────────────────────────

def _ctrl_rz(qc: QuantumCircuit, theta: float, anc: int, target: int) -> None:
    """Controlled-Rz(theta): applies Rz(theta) on target iff anc=|1⟩.
    Uses native CX(anc→target) direction.
    """
    qc.rz(theta / 2, target)
    qc.cx(anc, target)
    qc.rz(-theta / 2, target)
    qc.cx(anc, target)


def _ctrl_rx(qc: QuantumCircuit, theta: float, anc: int, target: int) -> None:
    """Controlled-Rx(theta): applies Rx(theta) on target iff anc=|1⟩.
    Uses H-conjugation to convert CX(anc→target) to controlled-Rx.
    """
    qc.h(target)
    qc.rz(theta / 2, target)
    qc.cx(anc, target)
    qc.rz(-theta / 2, target)
    qc.cx(anc, target)
    qc.h(target)


def _append_ctrl_u0(qc: QuantumCircuit, anc: int, q_fwd: int, q_inv: int,
                    n_steps: int) -> None:
    """Append n_steps of controlled-U₀ using native CX(anc→q_fwd/q_inv)."""
    for k in range(n_steps):
        p, rz, rx = get_gate_angles(k)
        # q_fwd: controlled-Rz(rz-p) then controlled-Rx(rx)
        _ctrl_rz(qc, rz - p, anc, q_fwd)
        _ctrl_rx(qc, rx,     anc, q_fwd)
        # q_inv: controlled-Rz(rz+p) then controlled-Rx(rx)
        _ctrl_rz(qc, rz + p, anc, q_inv)
        _ctrl_rx(qc, rx,     anc, q_inv)


def build_zpgpw_circuit(n_steps: int, basis: str) -> QuantumCircuit:
    """
    ZP-GPW Hadamard test circuit.

    basis='X': H(anc) → ctrl-U₀_n → H(anc) → measure
               ⟨Z_anc⟩ = ⟨X_anc_before⟩ = Re(⟨00|U₀_n|00⟩)

    basis='Y': H(anc) → ctrl-U₀_n → Sdg(anc) → H(anc) → measure
               ⟨Z_anc⟩ = ⟨Y_anc_before⟩ = Im(⟨00|U₀_n|00⟩)

    Virtual qubits: 0=anc, 1=q_fwd, 2=q_inv
    Physical layout set via initial_layout at transpile time.
    """
    assert basis in ('X', 'Y')
    qr = QuantumRegister(3, 'q')
    cr = ClassicalRegister(1, 'c')   # measure ancilla only
    qc = QuantumCircuit(qr, cr)

    qc.h(qr[0])                             # anc → |+⟩
    _append_ctrl_u0(qc, 0, 1, 2, n_steps)   # controlled-U₀_n
    if basis == 'Y':
        qc.sdg(qr[0])                       # rotate to Y basis
    qc.h(qr[0])                             # close Hadamard test
    qc.measure(qr[0], cr[0])               # measure ancilla
    return qc


# ─── Hardware run ──────────────────────────────────────────────────────────────

def run_zpgpw(backend, n_steps: int, shots: int,
              q_anc: int, q_fwd: int, q_inv: int) -> dict:
    layout = [q_anc, q_fwd, q_inv]

    results = {}
    for basis in ('X', 'Y'):
        qc = build_zpgpw_circuit(n_steps, basis)
        print(f"\n── ZP-GPW {basis}-basis  (n={n_steps} steps) ─────────────────────")
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

        counts  = result[0].data.c.get_counts()
        total   = sum(counts.values())
        p0      = counts.get('0', 0) / total
        p1      = counts.get('1', 0) / total
        exp_val = p0 - p1   # ⟨Z⟩ after basis rotation = ⟨X or Y⟩ before

        print(f"   Counts: {counts}   ⟨{basis}⟩ = {exp_val:+.4f}")
        results[basis] = {
            "counts":           counts,
            "shots":            total,
            "exp_val":          exp_val,
            "transpiled_depth": transpiled.depth(),
            "abstract_depth":   qc.depth(),
            "job_id":           job.job_id(),
        }

    x_hw = results['X']['exp_val']
    y_hw = results['Y']['exp_val']
    delta_hw = float(np.arctan2(y_hw, x_hw))

    return {
        "n_steps":   n_steps,
        "x_hw":      x_hw,
        "y_hw":      y_hw,
        "delta_hw":  delta_hw,
        "delta_deg": float(np.degrees(delta_hw)),
        "runs":      results,
        "backend":   backend.name,
        "timestamp": datetime.now().isoformat(),
    }


# ─── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="ZP-GPW: geometric phase on IBM hardware")
    parser.add_argument("--backend",  default="ibm_strasbourg")
    parser.add_argument("--token",    default=None)
    parser.add_argument("--shots",    type=int, default=8192)
    parser.add_argument("--steps",    type=int, default=6,
                        help="Ouroboros steps (default 6; try 4, 8, 12)")
    parser.add_argument("--q-anc",   type=int, default=72)
    parser.add_argument("--q-fwd",   type=int, default=62)
    parser.add_argument("--q-inv",   type=int, default=81)
    parser.add_argument("--sim-only", action="store_true")
    args = parser.parse_args()

    # ── Simulation ────────────────────────────────────────────────────────────
    print("\n── Ideal signal (noiseless) ───────────────────────────")
    for n in [4, 6, 8, 12]:
        s = expected_signal(n)
        print(f"   n={n:2d}:  ⟨X⟩={s['m00_re']:+.4f}  ⟨Y⟩={s['m00_im']:+.4f}"
              f"  δ={s['delta_deg']:+.2f}°  |M₀₀|={s['magnitude']:.4f}")

    target = expected_signal(args.steps)
    print(f"\n   Target (n={args.steps}):  ⟨X⟩={target['m00_re']:+.4f}  "
          f"⟨Y⟩={target['m00_im']:+.4f}  δ={target['delta_deg']:+.2f}°")

    if args.sim_only:
        return

    # ── Hardware ──────────────────────────────────────────────────────────────
    token = args.token or os.environ.get("IBM_QUANTUM_TOKEN")
    service = QiskitRuntimeService(
        channel="ibm_quantum_platform", token=token
    ) if token else QiskitRuntimeService(channel="ibm_quantum_platform")

    backend = service.backend(args.backend)
    print(f"\nBackend : {backend.name}  ({backend.num_qubits} qubits)")
    print(f"anc={args.q_anc}  q+={args.q_fwd}  q-={args.q_inv}")
    print(f"Steps: {args.steps}   Shots: {args.shots}")

    result = run_zpgpw(backend, args.steps, args.shots,
                       args.q_anc, args.q_fwd, args.q_inv)

    print(f"\n══ ZP-GPW Summary (n={args.steps}) ═══════════════════════════")
    print(f"  Ideal  :  ⟨X⟩={target['m00_re']:+.4f}  ⟨Y⟩={target['m00_im']:+.4f}  "
          f"δ={target['delta_deg']:+.2f}°")
    print(f"  Hardware: ⟨X⟩={result['x_hw']:+.4f}  ⟨Y⟩={result['y_hw']:+.4f}  "
          f"δ={result['delta_deg']:+.2f}°")
    print(f"  Δδ = {result['delta_deg'] - target['delta_deg']:+.2f}°  "
          f"({abs(result['delta_deg'] - target['delta_deg']):.1f}° error)")

    out = {"ideal": target, "hardware": result}
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts  = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = RESULTS_DIR / f"zpgpw_n{args.steps}_{backend.name}_{ts}.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults → {path}")


if __name__ == "__main__":
    main()
