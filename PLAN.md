# Plan

## Cut the OMK corona at an outer radius xi_cut (tau_r_corona, tau_perp_corona)

Status: planned (2026-11-11), not started.

### Problem
`tau_r_corona` (cpp/src/freddi_state.cpp) integrates the static-corona base density
n_e,b ~ xi'^-2 along the ray from the source to the disc surface at xi = R/R_IC, from
xi_min = (z_disk/R)^2 up to xi, with no outer limit (OMK91 Eqs. 39-40, exact for h_x = 0):

    tau = (L/Ledd)/Xi_tilde * max(0, 1/xi_min - 1/xi)

So for large xi it counts "static corona" where there is none. OMK91 Sect. 2 (after
Begelman, McKee & Shields 1983, Paper I): static corona for xi <~ 0.1; quasi-static with
substantial mass loss for 0.1 <~ xi <~ 1; wind for xi >~ 1 ("it becomes the wind"). The wind
has its own optical depths (tau_r_wind, tau_perp_wind) with different formulas.

### Change
Integrate the corona only up to xi_cut:

    tau = (L/Ledd)/Xi_tilde * max(0, 1/xi_min - 1/min(xi, xi_cut))

Two options to try (and compare):
- (b) xi_cut = 0.1 -- strictly the static corona;
- (c) xi_cut = 1   -- OMK91 "corona exists for R < R_IC".

Implementation notes:
- Make xi_cut a command-line option (e.g. `--corona_xi_cut`, default to be decided) so both
  options can be run without recompiling; document it in options.cpp help.
- `tau_r_corona()`: apply min(xi, xi_cut); zero if xi_min >= xi_cut.
- `tau_perp_corona()`: zero for xi > xi_cut (no static corona above the disc there: that
  column is wind, counted by tau_perp_wind).
- `tau_r_corona_at_R_IC()` (xi = 1): unchanged for (c); for (b) it becomes
  (L/Ledd)/Xi_tilde * max(0, 1/xi_min - 1/xi_cut) (corona column to R_IC).
- Scalar columns tau_r_corona_Rhot / _Fmax follow from the grid values automatically.
- Update comments and output.cpp help texts (keep the "G.2" warning).

### Open issue: overlap with the wind
f_dir = exp(-(tau_r_corona + tau_r_wind)) in tau_perp_corona counts the region
windR_launch_factor * R_IC < R < xi_cut * R_IC twice (corona and wind). In
test_Cirr_scatter_Tic_var windR_launch_factor = 1e-3, so the overlap is large for both (b)
and (c). Decide later whether to start the wind integral at xi_cut, or accept the overlap.

### Verification
- test_Cirr_scatter_Tic_var at t = 54.25 d is NOT a test: there the corona already ends at
  xi ~ 0.07 (Shadow() at 1.14e11 cm, R_IC = 1.7e12 cm). Need a model where the unshadowed
  corona reaches xi > 0.1 (e.g. higher T_ic -> smaller R_IC, or a flaring outer disc).
- Check that tau_r_corona is flat (saturates) beyond xi_cut and tau_perp_corona is zero there;
  compare (b) vs (c) with tau_r_evol / tau_perp_evol plots.
