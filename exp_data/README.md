# Single-Cell Tracking Dataset of HeLa Cells

Single-cell tracking data of **HeLa cervical cancer cells** acquired by long-term time-lapse microscopy using a custom shear-free, diffusive microfluidic platform.

The dataset was used for parameter estimation and validation of a cellular automata (CA) model describing cell proliferation, cell-cycle variability, and inheritance across successive cell generations.

---

## Dataset

The repository contains **15 CSV files**, each corresponding to one independent microfluidic chamber.

**Files**

- 1.csv
- 2.csv
- ...
- 15.csv

Each file contains the complete tracking information for all cells observed within a single chamber.

---

## Variables

| Variable | Description |
|----------|-------------|
| **Position X [μm]** | Horizontal coordinate of the cell centroid. |
| **Position Y [μm]** | Vertical coordinate of the cell centroid. |
| **Time [h]** | Time elapsed from the beginning of the experiment. |
| **Lineage ID** | Hierarchical identifier describing the cell pedigree. |

---

## Cell Lineage

Founder cells are assigned single-digit identifiers.

Example:

**Generation 1**

- 1

↓

**Generation 2**

- 11
- 12

↓

**Generation 3**

- 111
- 112
- 121
- 122

The number of digits in the **Lineage ID** corresponds directly to the generation number.

---

## Training and Test Sets

| Dataset | Chambers |
|---------|----------|
| **Training** | 9, 10, 11 |
| **Test** | 1–8, 12–15 |

---

## Repository Structure

- README.md
- 1.csv
- 2.csv
- ...
- 15.csv

---

## Applications

This dataset can be used for:

- Cellular automata modeling
- Cell lineage reconstruction
- Cell-cycle analysis
- Time-lapse microscopy
- Agent-based modeling
- Parameter estimation
- Model validation

---

## Citation

If you use this dataset, please cite the associated publication.

*The citation will be added after publication.*

---

## Contact

For questions regarding the dataset, please contact the corresponding author listed in the associated publication.
