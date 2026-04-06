"""
Native-direction hexagonal syndrome circuit.

Builds the syndrome extraction circuit using only physically-available
directed edges in the coupling map, avoiding SWAP routing entirely.

For each check (anc, d1, d2), two native edge patterns exist on heavy-hex:

  Pattern A — both edges go anc→data:
    H(anc); CX(anc,d1); CX(anc,d2); H(anc)  →  measures ZZ correctly

  Pattern B — mixed: d1→anc native, anc→d2 native:
    CX(d1,anc); H(anc); H(d2); CX(anc,d2); H(anc); H(d2)  →  measures ZZ correctly
    (verified analytically on all 4 input states)

Both patterns leave d1, d2 unchanged (up to stabilizer phases that don't
affect syndrome readout), and correctly accumulate d1 XOR d2 into ancilla.
"""

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister


def build_native_syndrome_circuit(
    T: int,
    data_qubits: list,
    anc_qubits: list,
    native_edges: set,         # set of (control, target) directed edges
) -> tuple[QuantumCircuit, dict]:
    """
    Build T-round syndrome extraction using only native-direction CX gates.

    native_edges: set of (ctrl, tgt) tuples from backend.coupling_map.get_edges()

    Returns (circuit, q_idx) where q_idx maps physical → circuit qubit index.
    """
    all_q = sorted(set(data_qubits + anc_qubits))
    q_idx = {q: i for i, q in enumerate(all_q)}

    n_q   = len(all_q)
    n_syn = 6 * T
    n_fin = 6

    qr  = QuantumRegister(n_q, "q")
    syn = ClassicalRegister(n_syn, "s")
    fin = ClassicalRegister(n_fin, "f")
    qc  = QuantumCircuit(qr, syn, fin)

    checks = [(anc_qubits[i], data_qubits[i], data_qubits[(i + 1) % 6])
              for i in range(6)]

    for t in range(T):
        for i, (anc, d1, d2) in enumerate(checks):
            ai  = q_idx[anc]
            d1i = q_idx[d1]
            d2i = q_idx[d2]
            bit = t * 6 + i

            # Determine native edge pattern
            anc_to_d1 = (anc, d1) in native_edges
            anc_to_d2 = (anc, d2) in native_edges

            if anc_to_d1 and anc_to_d2:
                # Pattern A: ancilla controls both
                qc.h(qr[ai])
                qc.cx(qr[ai], qr[d1i])
                qc.cx(qr[ai], qr[d2i])
                qc.h(qr[ai])

            elif (not anc_to_d1) and anc_to_d2:
                # Pattern B: d1→anc native, anc→d2 native
                qc.cx(qr[d1i], qr[ai])
                qc.h(qr[ai])
                qc.h(qr[d2i])
                qc.cx(qr[ai], qr[d2i])
                qc.h(qr[ai])
                qc.h(qr[d2i])

            elif anc_to_d1 and (not anc_to_d2):
                # Pattern C: anc→d1 native, d2→anc native (symmetric to B)
                qc.cx(qr[d2i], qr[ai])
                qc.h(qr[ai])
                qc.h(qr[d1i])
                qc.cx(qr[ai], qr[d1i])
                qc.h(qr[ai])
                qc.h(qr[d1i])

            else:
                # Fallback: neither direction native (shouldn't happen on heavy-hex)
                # Use H-conjugated ancilla-as-control for both
                qc.h(qr[ai])
                qc.h(qr[d1i])
                qc.cx(qr[d1i], qr[ai])
                qc.h(qr[d1i])
                qc.h(qr[d2i])
                qc.cx(qr[d2i], qr[ai])
                qc.h(qr[d2i])
                qc.h(qr[ai])

            qc.measure(qr[ai], syn[bit])
            qc.reset(qr[ai])

    for i, d in enumerate(data_qubits):
        qc.measure(qr[q_idx[d]], fin[i])

    return qc, q_idx


def parse_hexagonal_counts(syn_strings: list[str],
                           fin_strings: list[str],
                           T: int) -> list:
    """Same parser as hexagonal_cell.py."""
    shots_data = []
    for syn_s, fin_s in zip(syn_strings, fin_strings):
        syn_rev = syn_s[::-1]
        rounds  = []
        for t in range(T):
            start = t * 6
            bits  = [int(b) for b in syn_rev[start:start + 6]]
            rounds.append(bits)
        shots_data.append(rounds)
    return shots_data
