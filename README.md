# Paper 15: The Binary-Ternary Hybrid Architecture

**Two-Scale Error Correction, the Rotation Gap, and a Migration Path to Fault-Tolerant Quantum Computation**

Selina Stenberg with Claude Anthropic | March 2026

Part of the [Merkabit Research Series](https://github.com/SelinaAliens/The_Merkabit)

## What This Paper Shows

Every existing qubit is a merkabit restricted to the binary sector B31 of PSL(2,7). This paper measures exactly what that restriction costs.

### Three Results

1. **The Rotation Gap (12.6 pp)**: The detection rate difference between B31 nodes (static, chirality=0, 86.1%) and T75 nodes (dynamic, chirality +/-1, 98.7%) is exactly 12.6 percentage points at physical error rates. This gap is **flat across four orders of magnitude of error rate** (CV < 3%), confirming it is an architectural property, not a noise-dependent one.

2. **Two-Scale Composition (931-1862x suppression)**: The merkabit corrects errors at two independent scales:
   - Scale 1: Intra-unit toroidal correction (S_torus = 35x, directly measured)
   - Scale 2: Inter-unit torsion channel correction (S_L2 = 26.6x at 37 nodes)
   - Combined: 35 x 26.6 = 931x (without pi-lock) to 35 x 53.2 = 1862x (with pi-lock)
   - Zero qubit overhead

3. **Asymmetric Optimal Coupling**: The hybrid architecture's optimal coupling is kappa*_intra = 0, kappa*_inter = 1.0, confirming the two scales are physically distinct with different cost-benefit profiles.

### Central Comparison

| Architecture | Threshold | Suppression @ 1e-3 | Overhead |
|---|---|---|---|
| Surface code | ~1% | 10^6-10^10x | 1000 qubits/logical |
| Pure binary (kappa=0) | ~1% | 2.9x | 0 |
| Hybrid (kappa*) | ~22% | 11x | 0 |
| Pure ternary (both scales) | ~30%+ | 931-1862x | 0 |

## Repository Structure

```
simulations/          # All simulation code (NumPy only, seed=42)
  binary_ternary_hybrid_simulation.py   # Main hybrid architecture simulation
  scale_separation_sweep.py             # B31 vs T75 across 14 error rates
  asymptote_mapping.py                  # Extension to 61/91 nodes, saturation
  eisenstein_torus_simulation.py        # Direct S_torus measurement on torus
  temporal_spatial_bridge.py            # Ouroboros cycle temporal-spatial equivalence

outputs/              # Deterministic simulation results
  binary_ternary_hybrid_output.txt
  scale_separation_sweep_output.txt
  asymptote_mapping_output.txt
  eisenstein_torus_output.txt
  temporal_spatial_bridge_output.txt

paper/                # Paper 15 manuscript
  Paper_15_Binary_Ternary_Hybrid_Architecture_v4.docx
```

## Running the Simulations

All simulations require only NumPy. No external dependencies.

```bash
cd simulations
python binary_ternary_hybrid_simulation.py      # ~45s
python scale_separation_sweep.py                # ~160s
python asymptote_mapping.py                     # ~390s
python eisenstein_torus_simulation.py            # ~75s
python temporal_spatial_bridge.py                # ~8s
```

All results are deterministic (seed=42) and reproducible.

## Key Insight

The R gate sits inert in every qubit on every quantum computer running right now. It is the only gate in B31 = {R, S, T} that was ever a boundary crossing in the full ouroboros cycle. When the Z62 coupling activates, R becomes a live channel: errors flow through R into the ternary substrate, corrections flow back. The interface cost is small (5% of raw error rate) because you are not building a bridge. You are opening a door that was always in the wall.

## Related Work

- **Base paper**: [The Merkabit](https://doi.org/10.5281/zenodo.18925475) (Zenodo, v4, March 2026)
- **Paper 13**: PSL(2,7) three-stratum decomposition (10.5281/zenodo.19093820)
- **Paper 14**: Stable matter is self-correcting geometry (10.5281/zenodo.19153215)

## License

This work is provided for scientific review and verification. All simulation code may be freely run and inspected. Please cite the paper if you build on these results.
