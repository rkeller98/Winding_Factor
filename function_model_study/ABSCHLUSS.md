# Ergebnis der Funktionsmodell-Studie

Die separaten MATLAB-Auswertungen für Set 1–3 sind ausgeführt. Gauß-RBFs und
kubische Tensorprodukt-B-Splines verwenden gemeinsame Basiswerte für beide Flüsse.
Verglichen werden eine neu aus denselben Trainingsdaten aufgebaute vollständige
LUT sowie die bestehenden k-Varianten A–G. Das bisherige k-Faktor-Paper und die
produktive Regelsoftware wurden nicht verändert.

Gemeinsamer Fluss-RMSE auf vollständig zurückgehaltenen inneren Erregerstromgruppen,
in mVs:

| Set | RBF | Spline | vollständige LUT | kA | kF | gemeinsame Punkte |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3,074 | 1,450 | 1,744 | 4,521 | 4,570 | 1.034 / 1.156 |
| 2 | 3,200 | 1,552 | 1,823 | 5,726 | 5,070 | 1.551 / 1.734 |
| 3 | 1,853 | 1,210 | 1,421 | 6,115 | 5,706 | 1.654 / 1.819 |

Alle zehn Methoden werden je Set auf denselben unterstützten physischen Punkten
bewertet. Die gemeinsamen Abdeckungen liegen bei etwa 89,45 %, 89,45 % und 90,93 %.
Gesamtfehler werden aus quadrierten Punktfehlern und Punktzahlen gepoolt. Die CSVs
enthalten zusätzlich d-/q-RMSE, getrennte 95-%-Quantile und Maxima, alle Gruppen,
Rekonstruktionsfehler, eigene Abdeckung sowie Unterschiede zur frischen und zur
bereits vorhandenen LUT. Die vorhandene LUT hatte sämtliche Gruppen gesehen und
ist deshalb ausdrücklich keine unabhängige Prüfquelle.

Die Splines benötigen für beide Flüsse zusammen 972 / 972 / 490 Koeffizienten und
etwa 7,98 / 7,98 / 4,17 KiB numerischen Online-Speicher. Die vollständige LUT liegt
bei etwa 2.154 / 2.154 / 1.539 KiB. Die gewählten RBFs haben jeweils 128 Zentren,
264 Flusskoeffizienten einschließlich affinem Anteil und etwa 5,16 KiB Online-Daten.
Dokumentation, MATLAB-Strukturaufwand, Arbeitsbereich und ein zusätzlicher
geometrischer Gültigkeitswächter sind in diesen Online-Zahlen nicht enthalten.

Empfehlung: Als nächsten Genauigkeits-/Speicherkandidaten einen lokalen kubischen
Spline-Auswerter für die Zielhardware prototypisieren; die RBF als Laufzeitvergleich
mitführen. Der generische MATLAB-RBF-Auswerter ist hier schneller: ungefähr
14–27 µs pro vollständiger Abfrage gegenüber etwa 69–137 µs für die Splines.
Das Ausgabepaket besteht bei jeder Methode aus zwei Flüssen, sechs partiellen
Ableitungen und dem Gültigkeitsstatus. Das sind PC-Batchmittelwerte, keine
Einzellaufzeitgrenzen und kein 20-kHz-Nachweis.

Der untersuchte vorhandene Regelblock verarbeitet Sollströme und einen zusätzlichen
d-Stromanteil. Er benötigt zwei Flussauswertungen an unterschiedlichen Strompunkten
und vier d/q-Induktivitäten. Sein vollständiger Ersatzaufwand ist daher nicht mit
einer einzelnen Forschungsabfrage gleichzusetzen. Eine exakte Wiederverwendung
setzt unveränderte vollständige Eingänge, Parameter und Modi voraus; auch LUTs
könnten sie nutzen. Eine Wiederverwendungsquote wurde nicht gemessen.

Es fehlt ein anwendungsseitiger Grenzwert für Fluss-/Ableitungsfehler sowie ein
Zielhardware-Nachweis einschließlich freiem Taskbudget, tatsächlich benötigten
Ausgaben, Zahlenformat und Gültigkeitsprüfung. Unabhängig gemessene Ableitungen
liegen nicht vor. Die lokalen Ableitungsreferenzen stammen aus Trainingsflussdaten;
ihre Abweichungen sind Diagnostik. Temperaturen waren in Set 1–3 nicht konstant.
Die separat ausgewertete Temperatur-Testdatei mit 680 Punkten hat keine verfügbare
Rohdatei und nur zwei Temperaturbänder. Ihre 44 gemeinsamen Prüfpunkte betreffen
räumlich ausgelassene Strombereiche; daraus folgt kein physischer Temperatur- oder
4D-Nachweis. Dort liegt die RBF vorn, was ebenfalls gegen eine pauschale Rangfolge spricht.

Die obere Randgruppe wird getrennt geführt und außerhalb des Trainingsbereichs
zurückgewiesen. Auch innerhalb der Strombox kann die konservative gemeinsame
LUT-Maske Randflächen ausblenden. Leere Fehlerkarten sind fehlende gemeinsame
Unterstützung, kein Nullfehler. Der Atlas umfasst alle 376 Gruppen-/Methodenseiten
und ein Deckblatt. Das englische Paper hat 11 Seiten; beide PDFs wurden kompiliert,
gerendert und visuell geprüft.

Direkte numerische Checks: größte Abweichung der analytischen Ableitungen von
zentrierten Differenzen etwa 2,22e-11 in den jeweiligen physikalischen Einheiten;
LUT-Werte stimmen mit dem bestehenden Helfer bis etwa 5,55e-17 Vs überein.
Keine gespeicherten Trainingsindizes überlappen Prüf- oder Randpunkte. Alle
k-Optimierer meldeten in diesem Lauf einen positiven Exit-Flag. Ein globales
Optimum oder eine unabhängige physikalische Validierung folgt daraus nicht.

Start aus `C:/Git/Winding_Factor`:
```powershell
matlab -batch "restoredefaultpath; maxNumCompThreads(2); addpath(pwd); run_flux_function_study"
```
Modelle, Eingangspunkte und Tabellen sind in `model_results.zip` gesichert.
Zum Auswerten ohne Neufit nach `results/flux_fits/` entpacken und
`compare_flux_fits; compare_fits_with_k_factor; plot_flux_function_study` ausführen.
Die Dateireferenzen und vollständigen LaTeX-Befehle stehen in [RUN.md](RUN.md).
