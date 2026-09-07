"""Checks that the chosen system parameters (src/gwlens/system.py) actually
satisfy the regime assumptions stated in wiki/conventions.md. If these fail,
the whole "static lens" / "quasi-monochromatic" story is not justified for
the numbers actually used in cases/, so this gates them.
Run with `python3 tests/test_system.py`.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gwlens import system, chirp, units


def check_hierarchy():
    """a_in << a_out at both epochs used (Case A start and Case B) -- the
    condition for treating the binary as a single point source (no tidal /
    differential-lensing effects across it)."""
    a_in_A = system.inner_separation_m(system.F_A_START_HZ, system.M1_MSUN, system.M2_MSUN)
    a_in_B = system.inner_separation_m(system.F_B_HZ, system.M1_MSUN, system.M2_MSUN)
    ratio_A = a_in_A / system.A_OUT_M
    ratio_B = a_in_B / system.A_OUT_M
    ok = ratio_A < 1e-3 and ratio_B < 1e-3
    return ok, f"a_in/a_out at Case A start = {ratio_A:.2e}, at Case B = {ratio_B:.2e} (require < 1e-3)"


def check_static_lens_regime():
    """Case A: total merger duration (from f_A_start to formal coalescence)
    must be a tiny fraction of the outer orbital period -- the lens must not
    move at all during the observed chirp."""
    tau_A = chirp.time_to_merger(system.F_A_START_HZ, system.MCHIRP_MSUN)
    frac = tau_A / system.P_OUT_S
    ok = frac < 1e-3
    return ok, f"merger duration / outer period = {frac:.2e} (require < 1e-3)"


def check_quasi_monochromatic_regime():
    """Case B: time-to-merger from f_B must span several outer periods (so
    the orbital modulation is actually observable) while f itself barely
    changes over that span (checked separately, via the frequency drift)."""
    tau_B = chirp.time_to_merger(system.F_B_HZ, system.MCHIRP_MSUN)
    n_periods_to_merger = tau_B / system.P_OUT_S
    # fractional frequency drift over the ACTUAL Case B observing baseline
    # (system.T_OBS_B_S = 8 outer periods), not over the full remaining
    # time to merger -- an earlier version of this check used tau_B/2 and
    # failed because f formally diverges as t->t_c regardless of how the
    # observing window is chosen (see wiki/log.md).
    import numpy as np
    f_start = chirp.freq_of_time(np.array([0.0]), tau_B, system.MCHIRP_MSUN)[0]
    f_end = chirp.freq_of_time(np.array([system.T_OBS_B_S]), tau_B, system.MCHIRP_MSUN)[0]
    drift = abs(f_end - f_start) / f_start
    n_obs_periods = system.T_OBS_B_S / system.P_OUT_S
    ok = n_periods_to_merger > 3 * n_obs_periods and drift < 0.05
    return ok, (f"{n_periods_to_merger:.1f} outer periods before merger "
                f"(observing {n_obs_periods:.0f} of them); fractional freq "
                f"drift over the observing baseline = {drift:.2e} "
                f"(require merger>>observing window, drift<5%)")


def check_wave_optics_regime_nontrivial():
    """Neither epoch should sit in the trivial w->0 (F~1 everywhere, no
    lensing at all worth showing) or so deep in w>>1 that F is numerically
    indistinguishable from the pure geometric-optics envelope -- i.e. w is
    not absurd at either epoch (a loose sanity bound, not a design target)."""
    from gwlens import waveoptics as wo
    w_B = wo.w_of_frequency(system.F_B_HZ, system.M_LENS_MSUN)
    w_A_start = wo.w_of_frequency(system.F_A_START_HZ, system.M_LENS_MSUN)
    ok = 1e-2 < w_B < 1e2 and 1e-2 < w_A_start < 1e4
    return ok, f"w_B={w_B:.3f}, w_A_start={w_A_start:.3f}"


CHECKS = [
    ("hierarchy", check_hierarchy),
    ("static_lens_regime", check_static_lens_regime),
    ("quasi_monochromatic_regime", check_quasi_monochromatic_regime),
    ("wave_optics_regime_nontrivial", check_wave_optics_regime_nontrivial),
]


def main():
    print(system.summary())
    print()
    results = {}
    all_ok = True
    for name, fn in CHECKS:
        try:
            ok, msg = fn()
        except Exception as exc:  # noqa: BLE001
            ok, msg = False, f"EXCEPTION: {exc!r}"
        all_ok &= ok
        results[name] = {"passed": bool(ok), "message": msg}
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    out = Path(__file__).parent / "CHECKS_system.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
