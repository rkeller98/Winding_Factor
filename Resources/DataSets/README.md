# Messdaten

| Set | Drehzahl | Erregerstrom | Auswertung |
| --- | --- | --- | --- |
| Set1 | 3000 rpm | 2–40 A | 2 Varianten vorhanden |
| Set2 | 3000 rpm | 0–40 A | 3 Varianten vorhanden |
| Set3 | 2250 rpm | 3–30 A | neue Forschungsauswertung vorhanden |

- `raw/`: unveränderte Rohdateien.
- `FluxEval/`: vollständige Auswertungen, bleiben lokal und sind im Gitignore.
- `extracted/raw/`: alle Rohkanäle nach Drehzahl und Erregerstrom gruppiert.
- `extracted/evaluation/`: alle gespeicherten `Meas`-Felder aus den Auswertungen.

`research_flux_measurements.mat` wird für alle drei Sets mit
`prepare_reduction_data.m` neu aus den Rohkanälen erzeugt. Es enthält
arithmetisch gemittelte Betriebspunkte, gemessene Ströme und stationär aus den
gemessenen dq-Spannungen berechnete Flussverkettungen. Es ist ein kompakter
Forschungsexport ohne produktive GUI- oder LUT-Zustände. Set3 ist damit in der
aktuellen Untersuchung enthalten. Die übrigen historischen Exporte bleiben erhalten.

Einzelne FluxEval-Datei als MAT extrahieren:

```matlab
addpath('C:\Git\Winding_Factor\Resources');
measFile = extract_meas_from_fluxeval('C:\Pfad\FluxEval.mat');
data = load(measFile);
id = data.Meas.PE_Id_meas_avg{1,1};
```

Optionaler zweiter Parameter: Pfad der gewünschten Ausgabe-MAT.
Bei Quellen in `FluxEval/` landet die Ausgabe unter `extracted/evaluation/`.
Sonst entsteht `<Originalname>_measurements.mat` neben der Quelle.

`Meas.<Kanal>{Erregerstromzeile,Drehzahlspalte}` enthält alle Felder im
Cell-Layout aus `calculate_winding_factor.m`. Drehzahl und Erregerstrom sind
numerisch sortiert; `rpm_ref` und `exc_ref` geben die Referenzwerte an.
Zusätzlich bleibt `Total.<rpm>.<Iexc>.Meas` erhalten, samt Matrizen und Outliers.

Alle lokalen Dateien erneut extrahieren:

```matlab
addpath('C:\Git\Winding_Factor\Resources\DataSets');
extract_measurements;
```

Rohdaten behalten alle Wiederholungen und Kanäle: 268 in Set1/2, 195 in Set3.
`SampleIndex` und `Time` verweisen auf die ursprünglichen Samples.
Es wird nichts gemittelt oder gefiltert. Auch die 284 Index-0-Samples in Set3
bleiben erhalten. Die vorhandenen Auswertungen haben je 303 `Meas`-Felder.
