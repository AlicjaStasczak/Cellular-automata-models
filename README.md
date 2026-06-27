# Cellular Automata Model for Cancer Cell Growth

This repository contains an R implementation of a stochastic cellular automaton (CA) model developed to simulate cancer cell proliferation in a microfluidic environment.

The model reproduces experimentally observed cell population dynamics by incorporating:

- stochastic cell movement,
- cell division and death,
- experimentally estimated cell cycle distributions,
- neighborhood-dependent regulation of the cell cycle,
- generation-dependent regulation of the cell cycle,
- an iterative parameter fitting procedure.

The implementation was developed using experimental data obtained from long-term live-cell imaging of HeLa cells cultured in a microfluidic platform.

## Model Variants

Three biological hypotheses are implemented:

1. **Neighbourhood model** – cell cycle length depends on local cell density.
2. **Generation model** – cell cycle length depends on the number of previous cell divisions.
3. **Combined model** – integrates both neighbourhood- and generation-dependent effects.

## Input

The model requires:

- experimental Excel files (`.xlsx`) containing cell growth data,
- a reference `.rda` file with fitted probability distributions (cell cycle length, migration distance, and migration angle).

## Output

For each simulation the model records:

- number of cells over time,
- simulation iteration number,
- fitted parameter value,
- final simulation statistics.

## Requirements

- R (≥4.1)
- Required R packages are installed automatically when missing.

## Citation

If you use this code in your research, please cite:

**Student S., Staśczak A.**
*Modeling of living cancer cells in a microfluidic in vitro platform using cellular automata models.*
