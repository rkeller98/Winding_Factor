# Ergebnis und Ausführung

Die Untersuchung bleibt auf Englisch. MATLAB-Code, Rohdatenaufbereitung,
Rückrekonstruktion und Vorhersageprüfung wurden für Set1, Set2 und Set3 ausgeführt.
Die produktive Auswertungssoftware und sämtliche Rohdateien wurden unverändert gelassen.

Set1 enthält 7.476, Set2 8.092 und die neu erzeugte Set3-Auswertung 9.022 Betriebspunkte.
Set3 ist eine andere Maschinenkonfiguration (vier statt drei Polpaare). Die nominalen
Erregerströme gruppieren die Daten; gerechnet wird mit den gemessenen Strömen.

Gemeinsamer RMSE beider Flusskomponenten bei vollständig ausgelassenen inneren
Erregerstromgruppen, in mVs:

| Set | A: bisherig | C: konstanter gemeinsamer Fit | E: q-/erregerabhängig | F: zwei konstante Achsenfaktoren | G: drei Ströme | H: Dreistromreferenz |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 5,998 | 5,518 | 5,453 | 5,490 | 5,513 | 1,787 |
| 2 | 6,309 | 5,736 | 5,470 | 5,653 | 5,443 | 1,780 |
| 3 | 6,970 | 6,996 | 6,839 | 6,529 | 6,425 | 1,480 |

Jede Zeile verwendet für alle acht Methoden dieselben unterstützten physikalischen
Prüfpunkte. Die Abdeckung beträgt 89,15 %, 89,03 % beziehungsweise 90,75 % der
inneren Gruppen. Randgruppen stehen separat in `results/boundary_groups.csv`.
Das gemeinsame Maß ist die Wurzel aus dem Mittel der quadrierten Fehler beider
Komponenten; es ist kein arithmetischer Mittelwert einzelner Gruppen-RMSEs.

Empfehlung: Wenn eine Reduktion auf zwei Kennfeldachsen erforderlich ist und beide
Achsen verschoben werden dürfen, ist F der überzeugendste einfache Kompromiss.
Zwei konstante Parameter verbessern den gemeinsamen Fehler gegenüber A in allen
drei Datensätzen um etwa 8,5 %, 10,4 % und 6,3 %. E ist für Set1/2 eine sinnvolle
Alternative bei unveränderter physikalischer q-Achse. G gewinnt nicht durchgängig;
sein geringer Zusatzgewinn gegenüber einfacheren Varianten rechtfertigt die
zusätzlichen Parameter und Invertierbarkeitsbedingungen nicht allgemein.
Bei vorrangigem Genauigkeitsziel ist H vorzuziehen: etwa 21-mal so viele
Kennfeldknoten, dafür deutlich kleinere Vorhersagefehler.

Die exakte konstante d-Reduktion setzt gleichzeitig `Lde = k Ldd` und
`Lqe = k Lqd` voraus. Die Ableitungsdiagnostik lässt Residuen zurück, und die
verschobenen Originalflächen verlangen teilweise unterschiedliche Flüsse an
derselben reduzierten Position. Daraus folgt keine isolierte physikalische
Fehleruntergrenze: Messung, Temperatur, Vorverarbeitung und Interpolation bleiben
mit enthalten. Ein wechselndes iq erklärt die Ergebnisse nicht allein.

Grenzen: Es gibt keinen anwendungsseitig begründeten Akzeptanzgrenzwert und keine
unabhängige Absolutkalibrierung der Flüsse. Set3 fehlen unter anderem
Phasen-/RMS-Metadaten; sein Rotor-Temperaturkanal ist nicht als kalibrierter Verlauf
verwendbar. Vier G-Anpassungen erreichten die endliche Suchgrenze; ein globales
Optimum wird nicht behauptet. Es wurde keine Extrapolation der Interpolanten
zugelassen. Unterstützte reduzierte Abfragen einer Randgruppe sind trotzdem von
einer Interpolation innerhalb des gemessenen Dreistrombereichs zu unterscheiden.

Die schlanken analytischen Kontrollen liefern die erwarteten Ergebnisse:
konstante und integrierte variable Reduktion stimmen bis auf Rundung überein;
das falsche Produkt aus Erregerstrom und lokalem Faktor weicht ab; der bewusst
inkompatible Fall hat unterschiedliche q-Flüsse bei gleichen reduzierten Koordinaten.
Ein Testgerüst oder Assertions wurden nicht angelegt.

Exakte serielle Wiederholung in PowerShell aus dem Projektordner:

```powershell
Set-Location C:/Git/Winding_Factor
matlab -batch "restoredefaultpath; addpath(pwd); prepare_reduction_data; analyze_winding_factors; validate_winding_factor; audit_reduction_sources; analyze_flux_map_reduction; plot_reduction_study"
Set-Location C:/Git/Winding_Factor/winding_factor_study
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
pdflatex -interaction=nonstopmode -halt-on-error atlas.tex
```

Im durchgeführten Lauf wurden die drei unabhängigen Datensätze parallel in eigenen
MATLAB-Prozessen mit `maxNumCompThreads(2)` und `analyze_flux_map_reduction(1)`,
`(2)` beziehungsweise `(3)` berechnet. Die korrigierten Läufe benötigten ungefähr
1.000, 1.226 und 943 Sekunden. Die Logs liegen in `results/reduction_set*.log`.
MATLAB-Version: R2025b. Die Literaturangaben wurden anhand der Originalquellen
beziehungsweise DOI-Metadaten geprüft.

Die Datenherkunft, Detailausgaben und Zuordnung der Review-Anforderungen stehen
in [RUN.md](RUN.md). Die vollständigen Einzelgruppen-Vergleiche stehen im Atlas;
die numerischen MAT-/CSV-Dateien werden zusätzlich als Ergebnisarchiv übergeben.
Die Achsenbeschriftungen verwenden „in A“ und „in mVs“. Die 3D-Ansichten von ψd
haben vertauschte d-/q-Achsen; die Heatmaps behalten ihre bisherige Anordnung.
Figure 1 verwendet offene Meshes, damit die einzelnen Erregerstromstufen sichtbar bleiben.
Eine gemeinsame LaTeX-Vorlage erzeugt die Atlas-Seiten per Schleife aus drei
Datensatzdefinitionen statt aus wiederholten Seitenblöcken.
Die Arbeit bleibt auf `codex/flux-reduction-study`; keine Zusammenführung und kein Deployment.
