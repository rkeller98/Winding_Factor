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

Alle Rohdateien und extrahierten Messdaten liegen zusätzlich als gzip-komprimiertes
UTF-8-JSON (`.json.gz`) mit gleichem Dateinamen vor. Die MAT-Dateien bleiben unverändert erhalten.
Die JSON-Objekte entsprechen den Variablen und Feldern der MAT-Dateien;
es wird nichts gefiltert oder gemittelt. Nicht endliche Werte (`NaN`, `Inf`)
werden als `null` gespeichert. Vektoren werden zu Listen, Matrizen zu
Listen von Zeilen; die ursprüngliche Zeilen-/Spaltenorientierung eines
Vektors ist in JSON nicht enthalten. `SampleIndex` bleibt MATLAB-basiert
(erstes Sample = 1).

In Python ohne Zusatzpakete lesen:

```python
import gzip
import json
from pathlib import Path

path = Path("Resources/DataSets/Set1/extracted/evaluation/research_flux_measurements.json.gz")
with gzip.open(path, "rt", encoding="utf-8") as file:
    data = json.load(file)

group = next(iter(next(iter(data["Total"].values())).values()))
measurements = group["Meas"]
```

Bei Rohdateien steht die ursprüngliche Messstruktur unter dem ursprünglichen
Variablennamen (`next(iter(data.values()))`); `Y` enthält die Kanäle mit
`Name`, `Unit` und `Data`, `X` die Zeitachse. Für gruppierte Messdaten ist
`Total` der direkte Zugriff ohne die MATLAB-Cell-Indexierung von `Meas`.

JSON-Kopien in MATLAB neu erzeugen:

```matlab
addpath('Resources/DataSets');
export_measurements_json;
```

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
