# Parameter Sensitivity Outputs

This directory separates two tasks.

- `schakel2010_strict_reproduction/`: standalone Schakel & Smeulders (2010) reproduction using Table I parameters, Fig. 2 potential coefficients, and Table II/III energy-flux validation from the paper's definitions.
- `research_data/`: Sensitivity analysis for `global_evolution.xlsx`. One-at-a-time curves isolate model-input groups along the dissolution path. The waveform contribution bar uses signed changes in `log10(Amax)` between the first valid and last valid poroelastic snapshots; the residual is the nonlinear interaction not assigned to one parameter.
