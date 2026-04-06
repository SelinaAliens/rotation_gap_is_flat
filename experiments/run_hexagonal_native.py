"""
Run the native-direction hexagonal plaquette syndrome experiment on ibm_strasbourg.

Uses build_native_syndrome_circuit() which selects Pattern A or B for each
check based on which CX directions are physically available, avoiding all
SWAP routing.

Expected improvement over run_hexagonal.py:
  - Transpiled depth: ~42 per round (4-6 gates per check × 6 checks + overhead)
    vs 2045 total (102/round) in the routed version
  - Fire rates: should drop from 50-65% to ~5-15%
  - Fano factor: should approach paper's F = 0.856 ± 0.03
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

sys.path.insert(0, str(Path(__file__).parent))
from hexagonal_cell_native import (
    build_native_syndrome_circuit,
    parse_hexagonal_counts,
)
from run_hexagonal import hexagonal_fano, per_ancilla_fano

RESULTS_DIR = Path(__file__).parent.parent / "outputs" / "hardware"


def get_service(token: str | None = None) -> QiskitRuntimeService:
    token = token or os.environ.get("IBM_QUANTUM_TOKEN")
    if token:
        return QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    try:
        return QiskitRuntimeService(channel="ibm_quantum_platform")
    except Exception:
        print("ERROR: No IBM Quantum token found.")
        print("  Set IBM_QUANTUM_TOKEN env var or pass --token")
        sys.exit(1)


def run_native(backend, T: int, shots: int,
               data_qubits: list, anc_qubits: list) -> tuple:
    all_q = sorted(set(data_qubits + anc_qubits))
    native_edges = set(backend.coupling_map.get_edges())

    qc, q_idx = build_native_syndrome_circuit(
        T, data_qubits, anc_qubits, native_edges
    )
    print(f"  Circuit: T={T}, depth={qc.depth()}, qubits={qc.num_qubits}")
    print(f"  Physical qubits: {all_q}")

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
    syn_strings = data.s.get_bitstrings()
    fin_strings = data.f.get_bitstrings()

    shots_data = parse_hexagonal_counts(syn_strings, fin_strings, T)
    return shots_data, job.job_id()


def main():
    parser = argparse.ArgumentParser(
        description="Native hexagonal plaquette syndrome experiment"
    )
    parser.add_argument("--backend", default="ibm_strasbourg")
    parser.add_argument("--token",   default=None)
    parser.add_argument("--shots",   type=int, default=4000)
    parser.add_argument("--T",       type=int, default=20)
    parser.add_argument("--cell",    default="outputs/ibm_strasbourg_cells.json",
                        help="Path to cell JSON (default: outputs/ibm_strasbourg_cells.json)")
    args = parser.parse_args()

    # Load qubit assignments
    cell_path = Path(args.cell)
    if not cell_path.is_absolute():
        cell_path = Path(__file__).parent.parent / cell_path
    with open(cell_path) as f:
        cell = json.load(f)
    data_qubits = cell["hexagon_vertices"]
    anc_qubits  = cell["edge_qubits"]
    print(f"Cell loaded from {cell_path}")
    print(f"  Data qubits: {data_qubits}")
    print(f"  Anc  qubits: {anc_qubits}")

    print("\n" + "=" * 58)
    print("  Merkabit — Native Hexagonal Plaquette Experiment")
    print("  Native-direction CX, zero SWAP routing")
    print("  Paper 3: The Rotation Gap Is Not An Error")
    print("=" * 58)

    service = get_service(args.token)
    backend = service.backend(args.backend)
    print(f"Backend: {backend.name}  (qubits: {backend.num_qubits})")

    shots_data, job_id = run_native(
        backend, T=args.T, shots=args.shots,
        data_qubits=data_qubits, anc_qubits=anc_qubits
    )

    all_q = sorted(set(data_qubits + anc_qubits))
    fano = hexagonal_fano(shots_data, args.T, all_q=all_q, anc_qubits=anc_qubits)
    fano["job_id"]    = job_id
    fano["backend"]   = backend.name
    fano["circuit"]   = "native"

    anc_fano = per_ancilla_fano(shots_data, args.T, anc_qubits=anc_qubits)

    print(f"\n── Native Hexagonal Plaquette Results ──────────────────")
    print(f"  Job ID:            {job_id}")
    print(f"  Shots:             {fano['n_shots']}")
    print(f"  T rounds:          {fano['T_rounds']}")
    print(f"  Mean weight/shot:  {fano['mean_weight_total']:.3f}")
    print(f"  Fano factor F:     {fano['fano_factor']:.4f}")
    print(f"  Sub-Poissonian:    {fano['sub_poissonian']}")
    print(f"  Per-round F (avg): {fano['per_round_F_mean']:.4f} ± {fano['per_round_F_std']:.4f}")

    print(f"\n── Per-Ancilla Breakdown ────────────────────────────────")
    for a in anc_fano:
        print(
            f"  Ancilla {a['ancilla_idx']} (q{a['phys_qubit']:3d}): "
            f"fire_rate={a['fire_rate']:.4f}  F={a['fano_factor']:.3f}"
        )

    print(f"\n── Comparison with Paper 3 ──────────────────────────────")
    print(f"  Paper reports F = 0.856 ± 0.03  (sub-Poissonian, t = −131)")
    print(f"  Routed circuit:    F = 1.207    (ibm_strasbourg, per-round 0.469)")
    print(f"  Native circuit:    F = {fano['fano_factor']:.4f}              ", end="")
    if fano["sub_poissonian"]:
        print("← ✓ CONSISTENT with paper")
    else:
        print("← ✗ INCONSISTENT (super-Poissonian)")

    output = {"summary": fano, "per_ancilla": anc_fano}
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = RESULTS_DIR / f"hexagonal_native_{backend.name}_{ts}.json"
    with open(path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved → {path}")


if __name__ == "__main__":
    main()
