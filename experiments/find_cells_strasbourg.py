"""Find hexagonal cells in ibm_strasbourg (Eagle r3) coupling map."""
import json
import networkx as nx
from qiskit_ibm_runtime import QiskitRuntimeService

TOKEN = "nyqb_-xSU6pORdVplRs3GDeyeJQrFENyR8RDUY7W-nR6"

svc = QiskitRuntimeService(channel="ibm_quantum_platform", token=TOKEN)
b = svc.backend("ibm_strasbourg")
edges = list(set(tuple(sorted(e)) for e in b.coupling_map.get_edges()))

G = nx.Graph()
G.add_edges_from(edges)
deg = dict(G.degree())
d3 = [n for n, d in deg.items() if d == 3]
d2 = [n for n, d in deg.items() if d == 2]
print(f"Degree-3: {len(d3)}, Degree-2: {len(d2)}")

# Build vertex graph: degree-3 nodes adjacent if linked via a degree-2 node
Gv = nx.Graph()
Gv.add_nodes_from(d3)
for v in d3:
    for nb in G.neighbors(v):
        if deg[nb] == 2:
            for nb2 in G.neighbors(nb):
                if nb2 != v and deg[nb2] == 3:
                    Gv.add_edge(v, nb2)

hex6 = [c for c in nx.simple_cycles(Gv) if len(c) == 6]
print(f"Hexagonal 6-cycles: {len(hex6)}")

if hex6:
    h = hex6[0]
    print(f"First hexagon vertices: {h}")
    edge_q = []
    for i in range(6):
        a, bb = h[i], h[(i + 1) % 6]
        for nb in G.neighbors(a):
            if deg[nb] == 2 and bb in G.neighbors(nb):
                edge_q.append(nb)
                break
    print(f"Edge qubits:           {edge_q}")
    plaquette = sorted(set(h) | set(edge_q))
    print(f"Full 12-qubit plaquette: {plaquette}")

    # Save for circuit construction
    cell_info = {
        "backend": "ibm_strasbourg",
        "hexagon_vertices": h,
        "edge_qubits": edge_q,
        "plaquette_12": plaquette,
        "all_hexagons": hex6[:5],
    }
    with open("outputs/ibm_strasbourg_cells.json", "w") as f:
        json.dump(cell_info, f, indent=2)
    print("Saved → outputs/ibm_strasbourg_cells.json")
