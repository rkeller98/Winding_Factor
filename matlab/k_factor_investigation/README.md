# k-factor investigation — Arbeitsnotiz

Frage: Welche Variation entsteht durch Konturgeometrie/Interpolation, welche
durch die Flussberechnung, und was bleibt als mögliche physikalische Grenze
einer konstanten Stromverschiebung? Der konventionelle Estimator und seine
gespeicherten Ergebnisse bleiben unverändert. Keine neuen globalen Fits.

## Ausführen

In MATLAB R2025b vom Repository-Hauptordner:

```matlab
addpath(fullfile(pwd,'matlab','k_factor_investigation'));
study = run_k_factor_investigation; % alle Sets, insgesamt fünf Figures
% Alternativ gezielt: run_k_factor_investigation(1)
```

Eine Figure je Set bündelt Konturverfahren, Flusssymmetrie, Statorwiderstand,
lokale Faktoren und den Methodenvergleich in Tabs mit Diagnose-Untertabs.
Die ideale Maschine und der Vergleich über die ausgewählten Sets haben jeweils
eine eigene Figure. Wiederholte Aufrufe ersetzen nur den jeweiligen Plot-Tab;
sie öffnen keine weiteren Figures für dasselbe Set. Die Auswahl der
Erregerstromgruppen gilt gemeinsam für alle räumlichen Faktor-Tabs.
Konturverfahren und Differentialfaktoren werden nach getrennten Kriterien
verglichen; Quantilstriche beschreiben Kennfeldstreuung, keine Unsicherheit.

Die bisherigen Einzelaufrufe ergänzen dieselben Set-Figures:

```matlab
ideal = investigate_ideal_machine;
averaging = investigate_contour_averaging(1,[6 7],50);
symmetry = investigate_flux_symmetry(1,[6 7]);
resistance = investigate_resistance_sensitivity(1);
local = investigate_local_factors;  % alle drei Sets, gespeicherte Ableitungen
plot_local_factor_comparison(local{1});
plot_local_factor(local{1},'psi_d'); % auch 'psi_q', 'ls', 'kd', 'kq'
```

Die Funktionen sind unabhängig aufrufbar. Sie lesen die vorhandenen
`results/reduction_data.mat`, `winding_validation.mat` und
`flux_map_reduction.mat` im Repository-Hauptordner, ohne sie zu überschreiben.
Kein automatischer Export. Explizite Beispiele:

```matlab
investigate_ideal_machine('png');
investigate_contour_averaging(1,[6 7],50,'pdf');
investigate_flux_symmetry(1,[6 7],'png');
investigate_resistance_sensitivity(1,[-.1 0 .1],[6 7],'png');
plot_local_factor(local{1},'psi_q',[],[],'png');
```

Neue Exporte landen im lokalen, Git-ignorierten `results/` dieses Ordners.
`figures/` und `manuscript/` enthalten unveränderte historische Artefakte des
vorherigen Auftrags; hier wird kein Paper oder Atlas fortgeschrieben.

## Definitionen und Konturvergleich

Die unveränderte Referenz verwendet 100 Niveaus im mittleren 25–75-%-Bereich
der Flusswertespanne, Delaunay-Kanten in gemessenen `(id,iq)` und die bestehende
strenge Langkantenregel. Behaltene Schnittpunkte werden in `(id,iq,ie)` gemittelt:
`k=abs(Delta mean(id))/Delta mean(ie)`. Der Nenner behält sein Vorzeichen.
Der globale Faktor mittelt die gültigen Quotienten benachbarter Gruppen/Niveaus.

`contour_pair_statistics` vergleicht diese Schnittpunktmittel mit dem exakten
Längenmittelpunkt, komponentenweisem Median/10-%-Trim-Mittel,
Bogenlängen-/Segmentabtastung und einer ausgedünnten Triangulierung.
Physikalische Segmentlängen in `(id,iq)` gewichten die Segmentmittelpunkte;
doppelte Segmente werden entfernt. Isolierte Referenzschnittpunkte haben keine
Länge. Robuste Statistiken sind nicht automatisch geometrische Konturzentren.
Ausdünnen verändert auch den Interpolanten, Neuabtastung nur die Gewichte.

Direkter fixed-q-Vergleich hält **Flussniveau und physikalisches iq** fest:
`k_fixed_q=-Delta id/Delta ie`. Nur genau eine Wurzel pro Gruppe ergibt einen
Faktor; Mehrdeutigkeit und fehlende Abdeckung bleiben sichtbar. Der Nenner
verwendet interpolierte **gemessene** Erregerströme; Sollwerte benennen die Gruppen.
Die vorhandene Methode B optimiert dagegen die d-Fluss-Vorhersage einer
gemeinsamen Tabelle bei festem iq und ist kein direkter Konturwurzelvergleich.

## A/B: analytischer Nachweis und numerische Reproduktion

Sieben gezielte ideale Fälle verwenden die Koenergie
`W'=Ld*s^2/2+c*s*iq^2+Lq*iq^2/2` mit

```
s = id + k0*ie,              k0 = 12
psi_d = Ld*s + c*iq^2,       Ld = 0.001 H
psi_q = Lq*iq + 2*c*s*iq,    Lq = 0.010 H.
```

Die Reduktion ist exakt, d-Fluss gerade und q-Fluss ungerade in iq.
Der symmetrische Stator-Hessian `J=[Ld,2*c*iq;2*c*iq,Lq+2*c*s]` ist auf den
gewählten Bereichen positiv definit: kleinste untere Schranken rund `4.66e-6 H^2`
für seine Determinante und `5.59e-4 H` für seine Eigenwerte.

**Analytisch** gilt für Mittelwerte exakter Konturpunkte:
`k_signed=k0+(c/Ld)*Delta mean(iq^2)/Delta mean(ie)`.
Bei interpolierten Punkten kommt mit `eta=psi_d_analytic(intersection)-lambda`
der Term `-Delta mean(eta)/(Ld*Delta mean(ie))` hinzu. Gleiche Mittelwerte von iq
reichen also nicht: Die zweiten Momente sind ebenfalls entscheidend.

**Numerisch reproduziert**, ausgewähltes Niveau 50:

| Fall | Schnittpunktmittel k | Längenmittel k | Direktes fixed-q k |
|---|---:|---:|---:|
| Gleichmäßig, symmetrisch | 12 | 12 | 12 |
| Symmetrisch, unterschiedliche Punktdichte | 5.74983 | 12.06635 | 11.90–12.45 |
| Gleiche Punkte, andere Zelldiagonalen | 13.33925 | 12 | 12 |
| Symmetrisch, zweite Gruppe schmaler | 2.32204 | 0.81194 | 12 im gemeinsamen Bereich |

Im Dichtebeispiel ist die q-Schwerpunktdrift praktisch null, jedoch
`mean(iq^2)=31079.137/18525.761 A^2`. Die Identität schließt numerisch:
`12-6.276688+0.026522=5.749834`. Der ursprüngliche globale 100-Niveau-Estimator
liefert dort 6.074818. Analytisches fixed-q ergibt immer 12; die ausgewiesene
Streuung betrifft die lineare Interpolation. Bei c=0 verschwindet der Dichtefehler;
stärkere Krümmung verstärkt ihn.

Beim Diagonalenbeispiel sind die vier d-Flusswerte jeder Zelle koplanar.
Der Interpolant bleibt gleich, nur die Kantenschnittpunktgewichte ändern sich.
Der unveränderte Estimator mit ursprünglicher Triangulierung liefert weiterhin 12.
Längengewichtung behebt dagegen keine unterschiedliche Konturabdeckung.

**Empirisch:** Set1, 17/20 A, ursprüngliche Niveaus 25/50/75:
Schnittpunktmittel `13.13910/13.00107/13.49304`, Längenmittel
`13.30149/13.37512/13.60731 A/A`. Referenzmittel stimmen exakt überein (0 A Differenz).
Auf Niveau 50 reicht fixed-q von `10.17872` bis `15.49019 A/A`. Alle abgefragten
gemeinsamen q-Slices sind eindeutig; die je Niveau unterschiedlichen Intervalle
werden ausgegeben und gespeichert. Das ist keine Vollabdeckung des Kennfelds
und begründet noch keine Rangfolge der Methoden. Der ältere Spiegelindikator
im Konturplot nutzt die konvexe Hülle; kantenbegrenzte Spiegelabdeckung folgt in C.

## C/D: Symmetrie und Widerstand — empirische Sensitivität

Spiegelabfragen vergleichen gleiches id und entgegengesetztes iq innerhalb einer
nominalen Gruppe, nur in beidseitig unterstützten Dreiecken ohne überlange Kanten.
Mit `f+/-=f(id,+/-iq)` sind gerader/ungerader Anteil `(f+ + f-)/2` und
`(f+ - f-)/2`. Vollständige Residuen sind `psi_d+ - psi_d-` und `psi_q+ + psi_q-`.
Keine Extrapolation oder erzwungene Symmetrie.

Set1, 17/20 A: d-Residuen-RMS `4.709/4.594 mVs`, q-Residuen-RMS
`3.142/3.502 mVs`. Gemessenes ie unterscheidet sich an Spiegelpunkten aber um
`0.957/0.950 A RMS`: erhebliche Einschränkung einer Interpretation bei gleichem ie.

Die Rs-Sensitivität berechnet Flüsse aus originalen Spannungen, Strömen, omega
und punktweisem Rs plus konstantem Offset neu. Fehlende Rohkanäle führen zum
Auslassen. Set1 reproduziert beim Null-Offset die gespeicherten Flüsse exakt.
**Analytisch:** `Delta psi_d=-Delta Rs*iq/omega`,
`Delta psi_q=+Delta Rs*id/omega`. Bei konstantem omega ändert sich `Ldq-Lqd`
um `-2*Delta Rs/omega`; variables omega und Ableitungsschätzung ergänzen Terme.
Numerische Steigung `-0.002121966 H/Ohm`, gegenüber `median(-2/omega)=-0.002121968`.
Dies reproduziert Sensitivität, identifiziert keinen tatsächlichen Rs-Fehler.

Bei Set1 und `Delta Rs=[-0.0013,0,+0.0013] Ohm`:

| Größe | -10 % Rs | Original | +10 % Rs |
|---|---:|---:|---:|
| Globaler konventioneller Faktor | 13.49548 | 13.49212 | 13.48910 |
| Schnittpunktmittel, festes Niveau/Paar | 13.03357 | 13.00107 | 12.95372 |
| Fixed-q, iq=449.232 A, gleiches Niveau/Paar | 13.11409 | 13.10274 | 13.09138 |
| Median lokales LS, gemeinsame gültige Punkte | 13.17548 | 13.17895 | 13.18004 |

Global wird die ursprüngliche 100-Niveau-Regel auf jede neu berechnete Karte
angewendet; der Paarvergleich hält das ursprüngliche Niveau fest. Lokale
Ableitungen werden durch dieselbe lineare 80-Nachbar-Regression der
Flusskorrekturen ergänzt. Gedruckte/geplottete lokale Mediane verwenden gemeinsame
Masken; `study.summary` enthält dagegen pro Szenario gültige lokale Statistiken.
Bei konstantem omega gilt fixed-q:
`k_R(lambda,q)=k_original(lambda+Delta Rs*q/omega,q)`; Invarianz ist nicht garantiert.
iq=0 ist ein Sonderfall: `investigate_resistance_sensitivity(1,[],[],'',0)`.

## E/F: lokale Faktoren und interaktive Darstellung

Gespeicherte 80-Nachbar-Ableitungen liefern

```
kd(psi_d) = Lde/Ldd,       kd(psi_q) = Lqe/Lqd
kLS = (Ldd*Lde+Lqd*Lqe)/(Ldd^2+Lqd^2).
```

Komponentenquotienten brauchen unabhängig je einen endlichen Nenner ungleich null
mit `abs(Nenner)>=0.1*median(abs(Nenner))`. LS behält seine ursprüngliche Maske.
Kleine Lqd schließen `501/518/398` Punkte in Set1/2/3 aus.
Mediane auf gemeinsamer Komponenten-/LS-Abdeckung:

| Set | kd(psi_d) | kd(psi_q) | kLS |
|---|---:|---:|---:|
| 1 | 13.017 | 16.135 | 13.219 |
| 2 | 13.038 | 15.996 | 13.250 |
| 3 | 13.927 | 19.399 | 14.323 |

q-Quotienten streuen auch nach dieser einfachen Maske stark; kein Beweis einer
physikalischen Faktoränderung. Der Vergleichsplot zeigt Differenzen auf identischen
Punkten, kleine Lqd, Ldd und Gruppenmediane mit globalen Kontur-/Differentialfaktoren.
Ausgegebene Residuen bewerten Ableitungen innerhalb der Daten, keine Flussvorhersage.

Getrennte Schätzungen **einer d-Verschiebung** unterscheiden sich von zweiachsigen
Faktoren aus `J*[kd;kq]=[Lde;Lqe]`. Diese haben unabhängig `rcond(J)>=1e-3`;
1e-4/1e-2 werden ebenfalls betrachtet. Alle aktuellen Punkte bestehen diese
Schwellen, minimale rcond je Set `0.170/0.169/0.210`; keine Unsicherheitsgrenze.
Globale Methode F bleibt unverändert. Variable lokale Differentialrichtungen
dürfen nicht ungeprüft als endliche Koordinaten `id+kd*ie,iq+kq*ie` eingesetzt werden.

Interaktive Plots: physikalisches id/iq, Faktor als Höhe, gemessenes ie als Farbe.
Jede Gruppe hat eigenes Mesh und Checkbox, dazu `Alle`/`Keine`. Ursprüngliche
Gruppenpunkte werden trianguliert, dann ungültige/überlange Dreiecke entfernt;
Löcher bleiben offen und Gruppen werden nie verbunden. Keine Glättung oder
Wertbegrenzung. MATLAB-Drehung und Zoom bleiben verfügbar. Zusätzliche Maske:
`plot_local_factor(local{1},'psi_d',[],local{1}.common_scalar)`.

## Offene Ursachen und nächste Experimente

| Möglicher Einfluss | Bedeutung / zusätzlich benötigte Daten |
|---|---|
| Eisenverluste | Terminaler versus magnetisierender Strom; unabhängige Verlust-/Leistungs- und Phaseninformation nötig. |
| Winkeloffset | Kann Symmetrie um angenommene Achsen verändern; Winkelkalibrierung nötig. Konsistente konstante Drehung von Strom und Fluss erhält Reziprozität. |
| Sensorverzögerungen | Inkonsistenter Zeit-/Winkelbezug; synchronisierte Zeitreihen und genaue Verzögerungen nötig. |
| Wechselrichter | Spannungsfehler verändern berechneten Fluss; unabhängige Klemmenspannung, Totzeiten und Schaltbedingungen nötig. |
| Temperatur | Rs ändert sich; kalibriertes Wicklungs-R(T) und zugeordnete Wicklungstemperatur nötig. |
| Nichtstationarität | Stationäre Rechnung lässt Ableitungsterme weg; Zeitverläufe/Einschwingprüfung nötig. |
| Messbedingungen | Gruppen-/Kampagnenunterschiede; Wiederholungen unter identischen Bedingungen nötig. |

Spannungen, Ströme, omega und Rs liegen in allen drei kanonischen Sets vor.
Rs ist je Set konstant (`0.013/0.013/0.0366185 Ohm`). Temperaturen und teilweise
konfigurierte Sensorverzögerungen liegen vor, ersetzen aber keine Kalibrierung.
Einzelne Temperaturwerte sind unplausibel; keine direkte Temperaturkorrektur daraus.
Zeitaufgelöste Stationarität, unabhängige Winkel-/Spannungskalibrierung und
Verlusttrennung sind hier nicht nachgewiesen.

**Offene Hypothesen:** Welchen Anteil erklären Messfehler, Ableitungsfehler und
physikalische Nichtreduzierbarkeit jeweils? Keine eindeutige Zuordnung aus diesen
Daten. Sinnvolle nächste Experimente: Spiegelmessungen bei gleichem tatsächlichem
ie, unabhängige Rs-/Winkel-/Spannungskalibrierung und zusätzliche direkte fixed-q-
Vergleiche auf gemeinsamer Abdeckung. Erst nach deren Bewertung Verfahren verbessern.
