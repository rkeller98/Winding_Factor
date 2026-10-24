# Winding Factor – Messdatenaustausch

[analyze_winding_factors.m](analyze_winding_factors.m) untersucht die EESM→PSM-Reduktion:
jeweils `kd` und `kq` aus konstantem d-/q-Fluss, Mesh-Kreuzvergleiche,
alle Faktoren über dem Fluss sowie Originalkennfelder neben der bisherigen
Transformation. Start mit Set 2; MAT und Plotbilder landen lokal in `results/`.

[validate_winding_factor.m](validate_winding_factor.m) ergänzt lokale
Flussableitungen, gemeinsame Faktoren, Erregerstrom-Holdout und Rasterfehler für
Set 1 und Set 2. Die englische mathematische Ausarbeitung liegt in
[winding_factor_study](winding_factor_study/main.tex).

Die Datensätze unter [Resources/DataSets](Resources/DataSets/README.md) enthalten
unveränderte Rohdateien und gruppierte Messdaten als MATLAB-MAT-Dateien.
Die großen vollständigen FluxEval-Auswertungen bleiben lokal und sind ignoriert.

[analyze_flux_map_reduction.m](analyze_flux_map_reduction.m) vergleicht die
gespeicherten Einzelkennfelder und das gemeinsame 2D-Kennfeld mit ihren
Messwerten. Es benötigt die lokalen vollständigen FluxEval-Dateien von Set 1
und Set 2, untersucht zusätzlich Rastergröße und wenige konstante Faktoren
und speichert die Ergebnisse samt Punktmasken in `results/flux_map_reduction.mat`.
Die drei zugehörigen Abbildungen und die englische Herleitung sind in der
bestehenden Ausarbeitung enthalten.

Eine einzelne Auswertung extrahieren:

```matlab
addpath('C:\Git\Winding_Factor\Resources');
measFile = extract_meas_from_fluxeval('C:\Pfad\zur\FluxEval.mat');
data = load(measFile);
% Cell-Zeilen: Erregerstrom, Cell-Spalten: Drehzahl; numerisch sortiert.
id = data.Meas.PE_Id_meas_avg{1, 1};
% Originalhierarchie bleibt ebenfalls erhalten:
% data.Total.pos_3000_rpm.pos_2_A.Meas
```

[calculate_winding_factor.m](Resources/calculate_winding_factor.m) lädt den
kompakten Messdatenexport von Set 2. Die Auswertung für Set 3 steht noch aus;
dessen Rohdaten sind bereits vollständig gruppiert verfügbar.

Die Beziehung zwischen MeasEval, ParamGUI und VICE steht als kompakte
[Mapping-Referenz](Resources/Mapping/README.md) zur Verfügung. Die
vollständige ParamGUI bleibt lokal und ist vom Git ausgeschlossen.

Die frühere Untersuchung ist in einer [Zusammenfassung](Resources/investigation/README.md)
konsolidiert: Reihenfolgefehler, Vergleichswerte und offene Fragen.
