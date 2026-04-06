"""
Hexagonal plaquette syndrome extraction on ibm_fez.

Uses the first 12-qubit hexagonal plaquette found in the heavy-hex
coupling map:
  - 6 data qubits  (degree-3 / lattice vertices): [21, 23, 25, 41, 43, 45]
  - 6 ancilla qubits (degree-2 / edge qubits):    [22, 24, 36, 37, 42, 44]

Each ancilla sits between two data qubits on the hexagon boundary and
measures their Z-parity. This is a native circuit — every CNOT uses
a physical coupling edge; no transpilation needed for the syndrome
extraction within the plaquette.

Topology (degree-3 = data D, degree-2 = ancilla A):
         D21 - A22 - D23
        /                \
      A36                A24
        \                /
         D41 - A42 - D43
        (mirrored on far side)
         D41 - A42 - D43
        \                /
      ...      ...      ...
  Full ring: D21-A22-D23-A24-D25-A37-D45-A44-D43-A42-D41-A36-D21

The sublattice (Z3 chirality) for the 6 data qubits is assigned
by position in the hexagonal ring: {0,2,4} → sublattice A,
{1,3,5} → alternating B/C. In the Eisenstein Z[ω] embedding:
  (a+b) % 3 gives chirality ∈ {-1, 0, +1}.
"""

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

# ── Physical qubit assignments ────────────────────────────────────────────────
DATA_QUBITS = [21, 23, 25, 45, 43, 41]   # hexagon vertices, clockwise
ANC_QUBITS  = [22, 24, 37, 44, 42, 36]   # edge qubits between adjacent data

# Parity checks: anc[i] measures data[i] XOR data[(i+1)%6]
CHECKS = [(ANC_QUBITS[i], DATA_QUBITS[i], DATA_QUBITS[(i + 1) % 6])
          for i in range(6)]

# Z3 sublattice assignment (position in hexagon ring mod 3)
SUBLATTICE = [(i % 3) for i in range(6)]        # 0,1,2,0,1,2
CHIRALITY  = [0 if s == 0 else (1 if s == 1 else -1) for s in SUBLATTICE]


def build_hexagonal_syndrome_circuit(
    T: int,
    data_qubits: list | None = None,
    anc_qubits: list | None = None,
) -> QuantumCircuit:
    """
    Build T-round syndrome extraction for the hexagonal plaquette.

    Classical layout:
      s[t*6 + i] = ancilla i, round t   (6 ancillas × T rounds)
      f[i]       = final data qubit i   (6 data qubits)

    data_qubits / anc_qubits: override module-level defaults for
    use on backends other than ibm_fez (e.g. ibm_strasbourg).
    """
    if data_qubits is None:
        data_qubits = DATA_QUBITS
    if anc_qubits is None:
        anc_qubits = ANC_QUBITS

    checks = [(anc_qubits[i], data_qubits[i], data_qubits[(i + 1) % 6])
              for i in range(6)]

    all_q = sorted(set(data_qubits + anc_qubits))
    q_idx = {q: i for i, q in enumerate(all_q)}   # physical → circuit index

    n_q   = len(all_q)     # 12 qubits
    n_syn = 6 * T
    n_fin = 6

    qr  = QuantumRegister(n_q, "q")
    syn = ClassicalRegister(n_syn, "s")
    fin = ClassicalRegister(n_fin, "f")
    qc  = QuantumCircuit(qr, syn, fin)

    for t in range(T):
        for i, (anc, d1, d2) in enumerate(checks):
            ai  = q_idx[anc]
            d1i = q_idx[d1]
            d2i = q_idx[d2]
            bit = t * 6 + i

            qc.cx(qr[d1i], qr[ai])
            qc.cx(qr[d2i], qr[ai])
            qc.measure(qr[ai], syn[bit])
            qc.reset(qr[ai])

    for i, d in enumerate(data_qubits):
        qc.measure(qr[q_idx[d]], fin[i])

    return qc, q_idx


def parse_hexagonal_counts(syn_strings: list[str],
                            fin_strings: list[str],
                            T: int) -> list:
    """
    Parse BitArray strings into structured syndrome data.

    Returns list[shots] of list[T rounds] of list[6 ancilla bits].
    Qiskit get_bitstrings(): bit 0 of register is rightmost char.
    """
    shots_data = []
    for syn_s, fin_s in zip(syn_strings, fin_strings):
        # Reverse so index 0 = bit 0 = ancilla 0, round 0
        syn_rev = syn_s[::-1]
        rounds  = []
        for t in range(T):
            start = t * 6
            bits  = [int(b) for b in syn_rev[start:start + 6]]
            rounds.append(bits)
        shots_data.append(rounds)
    return shots_data
