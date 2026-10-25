# Winding-factor investigation summary

Consolidated from the 15 investigation reports dated 22–23 October 2026.
The numbers are results from that earlier investigation; nothing was recalculated
here. Only this summary is included in Git; other local investigation artifacts
are excluded.

## Identified cause

The historical 300/301 evaluations labeled "5%" and "10%" contain identical
measurement, grid, and fit data for all 14 excitation-current slices. Their order
differs: the 0 A slice appears first in one and last in the other.

`Torque_LUT_Calc.m` estimates the winding factor from adjacent slices using
`abs(ΔId) / ΔIexc`. At the transition from 40 to 0 A, the denominator becomes
negative. This changes the factor `k`, the projection `Id* = Id + k·Iexc`,
and subsequently the entire merged flux map.

The earlier cross-check sorted complete slices numerically by `Ref_exc_current`
before estimating the factor. This removed the difference. Replacing only the
final factor or individual breakpoints is insufficient. Voltage reserve alone
does not explain the observed difference; a fresh 300 → 301 import placed 0 A
last for both reserves.

## Values at the diagnostic point

Physical inputs: `Id = −405 A`, `Iq = 293 A`, `Iexc = 33.2 A`.
Stationary voltage calculation at 12,000 rpm, three pole pairs, and `Rs = 0.013 Ω`.
Flux measurements were recorded at 3,000 rpm. This diagnostic point is distinct
from the reported oscillation case at 50 → 70 Nm.

| Calculation | "5%" variant | "10%" variant |
| --- | ---: | ---: |
| Old export: winding factor k | 13.35836 | 11.28958 |
| Old export: Uq | 126.57 V | 11.87 V |
| Old export: voltage magnitude | 394.45 V | 386.96 V |
| Recalculation with main at the time: Uq | 136.33 V | 15.64 V |
| main with sorting: Uq | 136.33 V | 136.33 V |

The small q-component does not imply a correspondingly small voltage magnitude.
Differences between the old export and main at the time arise from different
calculation revisions.

The later raw-data cross-check interpolates over **measured** excitation
current: the nominal 32/36 A slices are approximately 32.4025/36.4030 A at the
diagnostic point. At an actual 33.2 A, this gives approximately **159.55 V Uq**.
The older cross-check using **nominal** excitation currents gave **169.40 V**.
Both are stationary projections from 3,000 rpm data. The expected 200 V is not
supported by a matching independent measurement at 12,000 rpm.

## Axes and handover

- `Inv_Param.CtrlMaps` contains 2D tables over `Id*` and `Iq`.
- `Inv_Param.Maps4D.CtrlMaps` contains physical Id/Iq tables per slice.
  Different Id breakpoints are therefore not inherently an error.
- Axes and tables must match. Legacy and Maps4D breakpoints must not be interchanged.
- The original Maps4D slice order is in `LUT.CtrlMaps.complete.I_exc_vec`;
  this vector is missing as a separate excitation-current axis from the inspected
  `Inv_Param.Maps4D.CtrlMaps` export.
- ParamGUI copies the legacy tables; the VICE Python script transposes them for
  the C access path. Details: [field mapping](../Mapping/README.md).
- Invalid values are sometimes stored as zeros with `*_ZERO` masks.
  These zeros are not valid physical measurements.

## Open questions

The reported **1.2 Hz oscillation remains unexplained**. Synchronous time series
and active controller parameters were unavailable in the material at the time.
The statement that older inductances eliminated the oscillation came from team
communication at the time. The synthetic testbench did not reproduce 1.2 Hz.

The 3D→2D projection remains a modeling assumption even after correcting slice
order. Flux and differential inductance must be considered together;
`L_exc = k²·L_dd` is not a general EESM relationship. The inspected tables showed
cross-derivative discrepancies and, in some cases, negative symmetric eigenvalues.
These findings alone do not establish the cause of the oscillation.

LOWESS, polynomials, 3D RBF, and potential fits were compared offline. RBF performed
slightly better on average but had larger maximum errors than LOWESS.
A general change of fitting method does not resolve the proven ordering error.

Other limitations at the time: no complete ParamGUI/COM/hardware round trip;
an available SLX using linear extrapolation versus C code using boundary clamping;
contradictory historical GUI version information. Runs used R2025b;
R2018a could not be executed at the time.

## Further investigation

1. Compare winding factors and maps for different import orders; sort complete
   excitation-current slices before merging.
2. Use physical 3D slices as a reference for the reduced 2D LUT, keeping measured
   and requested excitation currents separate.
3. For the oscillation, record requested/actual currents, active L/gains, voltage
   before/after limiting, DC voltage, controller modes, and delays together.
4. Check the VICE build actually used and the map values after parameter transfer.

According to the user's statement at the time, the 200-series and 300-series
measurements additionally differ by a −24 µs current-sensor delay; no causal
relationship with the delay was established. Files 201/202 complete the 40 A
measurement points of the 200-series; 300 contains 2–40 A and 301 adds 0 A.

The earlier report records a local sorting fix on `codex/eesm-excitation-order`,
commit `39fed85e2c4bfd462f04474100f25f075aa2a868`, based on `3189ae83`.
The current integration status of that fix was not investigated here.
