# Low-Frequency Limit Diagnostic

This directory evaluates the model behavior at the zero-frequency / Darcy-flow
limit and compares it with near-zero positive-frequency Schakel interface
solutions.

## Main interpretation

- The exact `frequency_Hz = 0` rows are analytic Darcy-limit diagnostics. They
  report `k(0) = k0`, `L(0)`, `sigma(0)`, and the steady electrokinetic response
  index `L(0)/sigma(0)`.
- The full Schakel/Pride wave boundary-value problem is not evaluated at
  `omega = 0`, because the equations contain conductive terms proportional to
  `1/omega` and poro-electromagnetic terms proportional to `1/omega^2`.
- Near-zero positive frequencies are retained as convergence diagnostics for
  the interface solver at theta = 45 degrees.
- Across the dissolution sequence, the final normalized Darcy response index is
  0.000837126 relative to the initial value.
- The near-zero interface matrix condition number ranges from 4.226e+11 to
  3.877e+31, so the interface coefficients near zero frequency should be
  interpreted as numerical-limit diagnostics rather than a time-domain waveform.

## Files

- `low_frequency_dynamic_and_interface.csv`: all frequencies and dissolution
  times, including the `frequency_Hz = 0` Darcy-limit rows.
- `darcy_limit_vs_dissolution_time.csv`: only the exact zero-frequency Darcy
  rows.
- `low_frequency_run_summary.csv`: run settings and output overview.
- `darcy_response_index_vs_dissolution_time.png`: evolution of the steady
  Darcy electrokinetic response index.
- `low_frequency_dynamic_convergence.png`: near-zero convergence of dynamic
  properties to their Darcy limits.
- `near_zero_interface_coefficients.png`: interface coefficient and matrix
  conditioning diagnostics for positive low frequencies.

The sampled frequencies are: 0, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000 Hz.
