# Winding Factor Analysis

Windows 64-bit application using the original MATLAB contour-factor analysis.

1. Copy `installer/WindingFactorAnalysis_Setup.exe` to the other PC and run it.
   MATLAB Runtime R2025b is included for offline installation; MATLAB and a
   MATLAB license are not required. If the installer is delivered as a ZIP,
   extract the complete ZIP first and run its setup EXE.
2. Start Winding Factor Analysis and click **Select Flux Eval**.
3. Select a FluxEval `.mat` containing `Total.<speed>_rpm.<excitation>_A.Meas`.
   The original interactive plot tabs open automatically, separately per speed.
4. **Save results...** optionally saves the imported data and contour results.
   Pointwise/global alignment remains available in the plot tabs, calculated
   when selected. Without a saved method-F reference, it uses the original
   conventional-seed fallback.

`application/WindingFactorAnalysis.exe` can be copied on its own to a PC that
already has the matching MATLAB Runtime. Closing the launcher closes its plots.

Rebuild with MATLAB R2025b + MATLAB Compiler:
`addpath('WindingFactorAnalysis'); deploymentScript`

Branding images are reused from GUI_VICE-MeasurementEvalKit.
