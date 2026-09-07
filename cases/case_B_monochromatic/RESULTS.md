# Case B: quasi-monochromatic source, moving lens

Produced by `run.py` (`python3 cases/case_B_monochromatic/run.py`, ~7 s);
every number below is in `provenance/numbers.json`. Same triple as
[Case A](../case_A_chirp/RESULTS.md), same inner binary, earlier in its
inspiral: f=0.05 Hz instead of the chirp's 10-125.6 Hz.

## Regime

Over the **6-outer-period (24 d)** observation used here, the GW frequency
drifts by only **4.0%** (0.0500 -> 0.0521 Hz) while there are **60.2** outer
periods left before merger — quasi-monochromatic, confirmed in
`tests/test_system.py::check_quasi_monochromatic_regime`. Unlike Case A, the
lens genuinely **moves**: y(t) sweeps through the diffraction pattern once
per 4-day outer period.

w_B = **0.309** throughout (fixed — only y(t) varies) — solidly in the
diffraction-dominated regime (`F_hybrid` never switches away from the exact
closed form here; `w_geo_threshold=30`).

## Repeated lensing

Exactly **half** the orbit is lensed at all (`fraction_of_time_lensed
=0.4997`): the source spends the other half in front of the lens
(`D_LS<0`, `geometry.impact_parameter_of_time`) where there is no lensing
geometry and F=1 by construction, not approximation. `caseB_y_of_t.png`
shows this directly; the impact parameter formally diverges (`y->infinity`)
exactly at each front/back crossing, where `D_LS->0`.

The minimum impact parameter reached, **y=1.589**, coincides with Case A's
y_A to 4 significant figures — not a coincidence: both are the same
geometric minimum of the same (circular, e_out=0) outer orbit, found by two
independently-written search routines (`run.py::main` in each case
directory). Treated here as a cross-check between the two cases, not
asserted in advance.

`caseB_repeated_pulses.png` shows the resulting **periodic amplification
pulses**, one per 4-day outer period, `|F|²` peaking at **1.385**
(`max_abs_F²`, from `max_abs_F=1.177`) with visible diffraction ringing on
each side of the main peak — the same ring structure as Case A's Figure 3,
just swept through in time instead of shown as a static map. Averaged over
the whole 24-day observation, `<|F|²>=1.006` and the RMS amplitude
amplification is **1.003**: close to unity, because the lensed half and the
completely-unlensed half (F=1) roughly average out — most of the
observation is spent either unlensed or only mildly amplified, with sharp
but brief peaks.

## The extended pattern, and how a detector would see it

`caseB_pattern_with_orbit.png` puts Case A's ring pattern and Case B's orbit
together explicitly: the fixed |F(w_B,y)|² diffraction pattern (concentric
rings, central Einstein-ring spike) with the source's actual sky-projected
track overlaid, split into its lensed (cyan) and unlensed (dashed green)
halves — this is the "extended, in-space" pattern and the trajectory a real
triple would trace through it, in one figure.

`caseB_detector_view.png` is the idealized-detector view asked for: top
panel, the strain envelope over the full 24-day observation, showing six
clean amplification pulses riding on the slowly-drifting quasi-monochromatic
carrier; bottom panel, a dedicated fine time grid (dt=0.1 s, since the
coarse 3600-point grid used for the envelope — 576 s spacing — wildly
undersamples the 20 s carrier period; using it directly for a "raw
waveform" panel was an aliasing bug caught by eye and fixed, see
wiki/log.md) zoomed on one pulse peak, where the lensed and unlensed
waveforms visibly dephase over just a few carrier cycles.
