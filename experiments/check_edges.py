"""Check directed coupling edges for the ibm_strasbourg hexagonal plaquette."""
import json
from qiskit_ibm_runtime import QiskitRuntimeService

TOKEN = "nyqb_-xSU6pORdVplRs3GDeyeJQrFENyR8RDUY7W-nR6"

svc = QiskitRuntimeService(channel="ibm_quantum_platform", token=TOKEN)
b = svc.backend("ibm_strasbourg")
edges = set(b.coupling_map.get_edges())

with open("outputs/ibm_strasbourg_cells.json") as f:
    cell = json.load(f)

data = cell["hexagon_vertices"]   # [62, 81, 79, 77, 58, 60]
anc  = cell["edge_qubits"]        # [72, 80, 78, 71, 59, 61]

print(f"{'Check':6} {'anc':4} {'d1':4} {'d2':4}  "
      f"{'anc→d1':8} {'d1→anc':8} {'anc→d2':8} {'d2→anc':8}")
print("-" * 68)
for i in range(6):
    a  = anc[i]
    d1 = data[i]
    d2 = data[(i + 1) % 6]
    print(f"{i:6} {a:4} {d1:4} {d2:4}  "
          f"{str((a,d1) in edges):8} {str((d1,a) in edges):8} "
          f"{str((a,d2) in edges):8} {str((d2,a) in edges):8}")
