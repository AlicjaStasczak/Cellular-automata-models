```markdown
# Single-Cell Tracking Dataset of HeLa Cells

## Overview

This repository contains single-cell tracking data of **HeLa cervical cancer cells** acquired by long-term time-lapse microscopy using a custom **shear-free, diffusive microfluidic platform**.

The dataset was generated to estimate model parameters and validate a **cellular automata (CA) model** describing cell proliferation, cell-cycle variability, and the inheritance of proliferative behavior across successive cell generations.

---

## Dataset Structure

The repository contains **15 CSV files**, each representing one independent microfluidic chamber analyzed in the study.

```

1.csv
2.csv
...
15.csv

```

Each file contains complete tracking information for all cells observed within a single chamber.

---

## Data Format

Each CSV file contains the following variables:

| Column | Description |
|---------|-------------|
| **Position X [μm]** | Horizontal (X) coordinate of the cell centroid, expressed in micrometers. |
| **Position Y [μm]** | Vertical (Y) coordinate of the cell centroid, expressed in micrometers. |
| **Time [h]** | Time elapsed since the beginning of the experiment (hours). |
| **Lineage ID** | Hierarchical identifier encoding the complete cell pedigree. |

---

## Cell Lineage Encoding

Each cell is assigned a unique hierarchical identifier representing its complete division history.

### Founder cells

Cells present at the beginning of the experiment receive single-digit identifiers:

```

1
2
3
...

```

### Cell division

After mitosis, daughter cells inherit the mother's identifier with an appended **1** or **2**:

```

1
├── 11
└── 12

```

Subsequent divisions continue recursively:

```

1
├── 11
│   ├── 111
│   └── 112
└── 12
├── 121
└── 122

```

The number of digits in the **Lineage ID** directly corresponds to the **generation number** of the cell.

---

## Dataset Partitioning

For parameter optimization and independent model evaluation, the dataset was divided into separate training and testing subsets.

### Training set

- Chamber 9
- Chamber 10
- Chamber 11

### Test set

- Chambers 1–8
- Chambers 12–15

---

## Repository Contents

```

.
├── 1.csv
├── 2.csv
├── ...
├── 15.csv
└── README.md

```

---

## Applications

This dataset is suitable for research involving:

- Cellular automata models
- Agent-based modeling
- Cell-cycle analysis
- Single-cell lineage reconstruction
- Time-lapse microscopy
- Cell proliferation studies
- Parameter estimation
- Mathematical modeling of tumor growth

---

## Citation

If you use this dataset in your research, please cite the associated publication.

**Citation will be added upon publication.**

---

## Contact

For questions regarding the dataset or the accompanying cellular automata model, please contact the corresponding author listed in the associated publication.
```
