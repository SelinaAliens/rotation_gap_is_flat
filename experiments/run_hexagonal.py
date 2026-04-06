"""
Run the 2D hexagonal plaquette syndrome experiment on ibm_fez.

Tests whether the sub-Poissonian Fano factor (F < 1) reported in
Paper 3 appears when using the correct 2D geometry — a native
12-qubit hexagonal plaquette on ibm_fez's heavy-hex coupling map.

The 1D repetition code experiments showed F >> 1 (d=3: 9.6, d=5: 19.0,
d=7: 17.5). This experiment tests whether the 2D hexagonal geometry
changes the result, as the paper's claims are rooted in 2D Eisenstein
lattice physics.

Circuit is native to ibm_fez — every CNOT uses a physical coupling edge.
We constrain the transpiler to the specific physical qubits via
initial_layout, with optimization_level=0 to prevent remapping.
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# Add experiments dir to path
sys.path.insert(0, str(Path(__file__).parent))
from hexagonal_cell import (
    build_hexagonal_syndrome_circuit,
    parse_hexagonal_counts,
    DATA_QUBITS,
    ANC_QUBITS,
)

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "hardware"


def get_service(token: str | None = None) -> QiskitRuntimeService:
    token = token or os.environ.get("IBM_QUANTUM_TOKEN")
    if token:
        return QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    try:
        return QiskitRuntimeService(channel="ibm_quantum_platform")
    except Exception:
        print("ERROR: No IBM Quantum token found.")
        print("  Set IBM_QUANTUM_TOKEN env var, or save credentials with:")
        print("  QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token='<TOKEN>')")
        sys.exit(1)


def run_hexagonal(backend, T: int, shots: int,
                  data_qubits=None, anc_qubits=None) -> tuple:
    if data_qubits is None:
        data_qubits = DATA_QUBITS
    if anc_qubits is None:
        anc_qubits = ANC_QUBITS

    all_q = sorted(set(data_qubits + anc_qubits))
    qc, q_idx = build_hexagonal_syndrome_circuit(T, data_qubits, anc_qubits)

    print(f"  Circuit: T={T}, depth={qc.depth()}, qubits={qc.num_qubits}")
    print(f"  Physical qubits: {all_q}")

    # initial_layout as list: virtual qubit i → physical qubit all_q[i]
    # optimization_level=0: assign layout only, no remapping or optimisation
    pm = generate_preset_pass_manager(
        optimization_level=0,
        backend=backend,
        initial_layout=all_q,
    )
    transpiled = pm.run(qc)
    print(f"  Transpiled depth: {transpiled.depth()}")

    sampler = Sampler(backend)
    job = sampler.run([transpiled], shots=shots)
    print(f"  Job ID: {job.job_id()} — waiting for results …")
    result = job.result()

    pub_result = result[0]
    data = pub_result.data
    syn_strings = data.s.get_bitstrings()   # list[str], each len = T*6
    fin_strings = data.f.get_bitstrings()   # list[str], each len = 6

    shots_data = parse_hexagonal_counts(syn_strings, fin_strings, T)
    return shots_data, job.job_id()


def hexagonal_fano(shots_data: list, T: int,
                   all_q=None, anc_qubits=None) -> dict:
    """
    Compute Fano factor on total syndrome weight per shot.

    shots_data: list[shots] of list[T rounds] of list[6 ancilla bits]

    Fano factor F = Var(weight) / Mean(weight).
    F < 1 → sub-Poissonian (paper's claim).
    F > 1 → super-Poissonian (what 1D repetition code showed).
    """
    # Total syndrome weight per shot (summed across all rounds)
    total_weights = np.array(
        [sum(sum(r) for r in shot) for shot in shots_data], dtype=float
    )

    mean_w = total_weights.mean()
    var_w  = total_weights.var(ddof=1)
    F      = var_w / mean_w if mean_w > 0 else float("nan")

    # Per-round syndrome weights: shape (shots, T)
    per_round = np.array(
        [[sum(r) for r in shot] for shot in shots_data], dtype=float
    )
    round_means = per_round.mean(axis=0)   # shape (T,)
    round_vars  = per_round.var(axis=0, ddof=1)
    per_round_F = np.where(round_means > 0, round_vars / round_means, np.nan)

    return {
        "geometry":           "hexagonal_2d",
        "physical_qubits":    all_q or sorted(set(DATA_QUBITS + ANC_QUBITS)),
        "n_shots":            len(shots_data),
        "T_rounds":           T,
        "mean_weight_total":  float(mean_w),
        "var_weight_total":   float(var_w),
        "fano_factor":        float(F),
        "sub_poissonian":     bool(float(F) < 1.0),
        "per_round_F_mean":   float(np.nanmean(per_round_F)),
        "per_round_F_std":    float(np.nanstd(per_round_F)),
    }


def per_ancilla_fano(shots_data: list, T: int, anc_qubits=None) -> list[dict]:
    """
    Per-ancilla Fano factor: for each of the 6 ancillas,
    count how many times it fires (=1) across T rounds per shot.
    """
    if anc_qubits is None:
        anc_qubits = ANC_QUBITS
    results = []
    for anc_i in range(6):
        counts = np.array(
            [sum(shot[t][anc_i] for t in range(T)) for shot in shots_data],
            dtype=float,
        )
        mean_c = counts.mean()
        var_c  = counts.var(ddof=1)
        F      = var_c / mean_c if mean_c > 0 else float("nan")
        results.append({
            "ancilla_idx": anc_i,
            "phys_qubit":  anc_qubits[anc_i],
            "fire_rate":   float(mean_c / T),
            "fano_factor": float(F),
        })
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Hexagonal plaquette syndrome experiment on ibm_fez"
    )
    parser.add_argument("--backend", default="ibm_fez", help="IBM backend name")
    parser.add_argument("--token",   default=None,      help="IBM Quantum API token")
    parser.add_argument("--shots",   type=int, default=4000)
    parser.add_argument("--T",       type=int, default=20, help="Syndrome rounds")
    parser.add_argument("--cell",    default=None,
                        help="Path to cells JSON (e.g. outputs/ibm_strasbourg_cells.json)")
    args = parser.parse_args()

    # Load qubit assignments from cell JSON if provided
    data_qubits, anc_qubits = None, None
    if args.cell:
        with open(args.cell) as f:
            cell = json.load(f)
        data_qubits = cell["hexagon_vertices"]
        anc_qubits  = cell["edge_qubits"]
        print(f"Cell loaded from {args.cell}")
        print(f"  Data qubits: {data_qubits}")
        print(f"  Anc  qubits: {anc_qubits}")

    print("\n" + "=" * 55)
    print("  Merkabit — Hexagonal Plaquette Experiment")
    print("  Testing 2D geometry on ibm_fez heavy-hex")
    print("  Paper 3: The Rotation Gap Is Not An Error")
    print("=" * 55)

    service = get_service(args.token)
    backend = service.backend(args.backend)
    print(f"Backend: {backend.name}  (qubits: {backend.num_qubits})")

    shots_data, job_id = run_hexagonal(backend, T=args.T, shots=args.shots,
                                       data_qubits=data_qubits, anc_qubits=anc_qubits)

    all_q = sorted(set((data_qubits or DATA_QUBITS) + (anc_qubits or ANC_QUBITS)))
    fano = hexagonal_fano(shots_data, args.T, all_q=all_q, anc_qubits=anc_qubits)
    fano["job_id"]  = job_id
    fano["backend"] = backend.name

    anc_fano = per_ancilla_fano(shots_data, args.T, anc_qubits=anc_qubits)

    print(f"\n── Hexagonal Plaquette Results ──────────────────────")
    print(f"  Job ID:            {job_id}")
    print(f"  Shots:             {fano['n_shots']}")
    print(f"  T rounds:          {fano['T_rounds']}")
    print(f"  Mean weight/shot:  {fano['mean_weight_total']:.3f}")
    print(f"  Fano factor F:     {fano['fano_factor']:.4f}")
    print(f"  Sub-Poissonian:    {fano['sub_poissonian']}")
    print(f"  Per-round F (avg): {fano['per_round_F_mean']:.4f} ± {fano['per_round_F_std']:.4f}")

    print(f"\n── Per-Ancilla Breakdown ────────────────────────────")
    for a in anc_fano:
        print(
            f"  Ancilla {a['ancilla_idx']} (q{a['phys_qubit']:3d}): "
            f"fire_rate={a['fire_rate']:.4f}  F={a['fano_factor']:.3f}"
        )

    print(f"\n── Comparison with Paper 3 ──────────────────────────")
    print(f"  Paper reports F = 0.856 ± 0.03  (sub-Poissonian, t = −131)")
    print(f"  1D rep code:   F ≈ 9.6–19.0     (super-Poissonian)")
    print(f"  2D hex cell:   F = {fano['fano_factor']:.4f}              ", end="")
    if fano["sub_poissonian"]:
        print("← ✓ CONSISTENT with paper")
    else:
        print("← ✗ INCONSISTENT (super-Poissonian)")

    # Save full results
    output = {
        "summary":    fano,
        "per_ancilla": anc_fano,
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = RESULTS_DIR / f"hexagonal_{backend.name}_{ts}.json"
    with open(path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved → {path}")


if __name__ == "__main__":
    main()
