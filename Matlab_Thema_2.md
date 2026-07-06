Einfluss unterschiedlicher FFT-Fensterfunktionen auf die
HRV-Analyse
June 9, 2026
Ziel des Projekts ist die Untersuchung des Einflusses verschiedener Fensterfunktionen
aufdiespektraleAnalysederHerzfrequenzvariabilit¨at(HRV)ausEKG-Daten. Dabeisollen
unterschiedliche FFT-Fenster hinsichtlich ihrer Auswirkungen auf:
• Frequenzaufl¨osung
• Spektralleckeffekte
• Amplitudenverhalten
• HRV-Kennwerte
• sowie die 3D-Wasserfalldarstellung
analysiert und bewertet werden. Die Implementierung erfolgt vollst¨andig in MatLab.
1 Hintergrund
BeiderdiskretenFourier-TransformationendlicherSignalausschnitteentstehensogenannte
Spektralleckeffekte (“Spectral Leakage”). Fensterfunktionen dienen dazu, diese Effekte zu
reduzieren. Je nach Fenster ergeben sich unterschiedliche Eigenschaften hinsichtlich:
• Hauptkeulenbreite (Breite des zentralen Peaks im Frequenzspektrum – bestimmt die
Frequenzaufl¨osung)
• Nebenkeulenunterdru¨ckung(D¨ampfungderunerwu¨nschtenAusl¨aufernebendemHaupt-
peak – bestimmt wie stark Leckeffekte unterdru¨ckt werden)
• Frequenzaufl¨osung
• Amplitudengenauigkeit
1

Diese Unterschiede beeinflussen unmittelbar die HRV-Frequenzanalyse.
2 Aufgabenstellung
| 2.1 Import              | und | Vorverarbeitung |     | von EKG-Daten  |
| ----------------------- | --- | --------------- | --- | -------------- |
| Es sind EKG-Datens¨atze |     | einzulesen      | und | aufzubereiten. |
Mindestanforderungen
• Import eines mindestens 1-stu¨ndigen EKG-Datensatzes im EDF-Format (European
| Data Format)         |                |                |     |     |
| -------------------- | -------------- | -------------- | --- | --- |
| • Entfernung         | von            | Gleichanteilen |     |     |
| • Filterung          | von St¨orungen |                |     |     |
| – Baseline-Wanderung |                |                |     |     |
| – Netzbrummen        |                |                |     |     |
| – hochfrequentes     |                | Rauschen       |     |     |
• Darstellung
| – Rohsignal          |        |     |                        |     |
| -------------------- | ------ | --- | ---------------------- | --- |
| – gefiltertes        | Signal |     |                        |     |
| 2.2 R-Peak-Erkennung |        | und | RR-Intervallberechnung |     |
Anforderungen
| • automatische | Detektion | der | R-Zacken |     |
| -------------- | --------- | --- | -------- | --- |
•
| Berechnung  | der          | RR-Intervalle |     |     |
| ----------- | ------------ | ------------- | --- | --- |
| • Erkennung | fehlerhafter | Intervalle    |     |     |
• Plausibilit¨atspru¨fung
Darstellung
| • markierte | R-Peaks | im EKG |     |     |
| ----------- | ------- | ------ | --- | --- |
• RR-Intervall-Zeitreihe
2

| 2.3 Erzeugung | eines | gleichm¨aßig | abgetasteten | HRV-Signals |
| ------------- | ----- | ------------ | ------------ | ----------- |
Da RR-Intervalle ungleichm¨aßig verteilt sind, muss eine Interpolation erfolgen.
Anforderungen
| • Interpolation | der RR-Zeitreihe |     |     |     |
| --------------- | ---------------- | --- | --- | --- |
•
| Auswahl         | geeigneter | Interpolationsverfahren |                   |     |
| --------------- | ---------- | ----------------------- | ----------------- | --- |
| • Wahl einer    | geeigneten | Abtastrate              | (z.B. 4 Hz)       |     |
| 2.4 FFT-Analyse | mit        | verschiedenen           | Fensterfunktionen |     |
Die FFT ist mehrfach mit unterschiedlichen Fensterfunktionen durchzufu¨hren.
| Zu untersuchende |     | Fenster |     |     |
| ---------------- | --- | ------- | --- | --- |
Mindestens:
•
Rechteckfenster
• Hann-Fenster
•
Hamming-Fenster
• Blackman-Fenster
Optional:
• Kaiser-Fenster
•
Flat-Top-Fenster
Hinweis: Alle genannten Fensterfunktionen sind in MATLAB als eingebaute Funktio-
nen verfu¨gbar, z.B. hann(N), hamming(N), blackman(N), kaiser(N, beta). Die Fen-
sterl¨ange N entspricht der Anzahl der Abtastwerte im jeweiligen Analysefenster.
3

| 2.5 Untersuchung |              | der      | Fenstereigenschaften |     |                 |     |
| ---------------- | ------------ | -------- | -------------------- | --- | --------------- | --- |
| Fu¨r jedes       | Fenster sind | folgende | Eigenschaften        |     | zu analysieren: |     |
Frequenzaufl¨osung
| Bewertung | der Trennsch¨arfe |     | benachbarter |     | Frequenzen. |     |
| --------- | ----------------- | --- | ------------ | --- | ----------- | --- |
Spektralleckeffekte
Bewertung unerwu¨nschter Energieverteilungen außerhalb der eigentlichen Frequenzen.
Amplitudenverhalten
| Vergleich  | der Spektralamplituden |     |      | zwischen | den Fenstern. |     |
| ---------- | ---------------------- | --- | ---- | -------- | ------------- | --- |
| Einfluss   | auf HRV-Kennwerte      |     |      |          |               |     |
| Berechnung | und Vergleich          |     | von: |          |               |     |
•
VLF-Leistung
• LF-Leistung
• HF-Leistung
• LF/HF-Verh¨altnis
2.6 HRV-Frequenzbandanalyse
| Die Spektren | sind gem¨aß | HRV-Leitlinie |                    | aufzuteilen: |             |               |
| ------------ | ----------- | ------------- | ------------------ | ------------ | ----------- | ------------- |
|              |             | Frequenzband  |                    |              | Bereich     |               |
|              |             |               | VLF                |              | 0.0033 -    | 0.04 Hz       |
|              |             |               |                    | LF           | 0.04 - 0.15 | Hz            |
|              |             |               |                    | HF           | 0.15 - 0.40 | Hz            |
|              |             | Table         | 1: Frequenzb¨ander |              | gem¨aß      | HRV-Leitlinie |
| Zu berechnen | sind:       |               |                    |              |             |               |
| • absolute   | Leistung    |               |                    |              |             |               |
| • relative   | Leistung    |               |                    |              |             |               |
• LF/HF-Verh¨altnis
4

2.7 3D-Wasserfalldarstellung
Fu¨r jede Fensterfunktion ist eine Wasserfalldarstellung der zeitlichen Frequenzentwicklung
zu erzeugen.
Hinweis: Die Wasserfalldarstellung entsteht, indem die Messung in gleichlange, sich
u¨berlappende Zeitfenster aufgeteilt wird (siehe Abschnitt 3). U¨ber jedem Fenster wird
eine eigene FFT berechnet. Die so entstehenden Spektren werden hintereinander als 3D-
Diagramm dargestellt, wobei die x-Achse die Frequenz, die y-Achse die Zeit (Fensterindex)
und die z-Achse die spektrale Leistung zeigt. Da die Fensterfunktion das Spektrum jedes
einzelnen Fensters beeinflusst, ver¨andert sie direkt das Aussehen des gesamten Wasserfall-
diagramms.
Anforderungen
| • Darstellung  | der Spektren             | u¨ber der Zeit |
| -------------- | ------------------------ | -------------- |
| • Vergleich    | der Fensterfunktionen    |                |
| • Markierung   | der HRV-Frequenzbereiche |                |
| • einheitliche | Skalierung               |                |
| 2.8 Vergleich  | und Bewertung            |                |
Die Fensterfunktionen sind hinsichtlich ihrer Eignung fu¨r die HRV-Analyse zu bewerten.
| Zu diskutierende | Kriterien |     |
| ---------------- | --------- | --- |
• Frequenzaufl¨osung
•
| Stabilit¨at      | der Spektren |     |
| ---------------- | ------------ | --- |
| • Unterdru¨ckung | von Leakage  |     |
•
Rechenaufwand
| • Eignung | fu¨r Wasserfalldarstellungen |     |
| --------- | ---------------------------- | --- |
•
| Interpretierbarkeit | der | HRV-B¨ander |
| ------------------- | --- | ----------- |
5

| 3 Technische |     | Mindestanforderungen |     |
| ------------ | --- | -------------------- | --- |
Analysefenster
| • Fensterbreite: | 5   | Minuten |     |
| ---------------- | --- | ------- | --- |
• U¨berlappung:
|     | 50  | %   |     |
| --- | --- | --- | --- |
4 MatLab-Struktur
| Das Projekt | soll modular | aufgebaut | sein. |
| ----------- | ------------ | --------- | ----- |
Mindestfunktionen
| • load | ecg data() |     |     |
| ------ | ---------- | --- | --- |
•
| preprocess | ecg()     |     |     |
| ---------- | --------- | --- | --- |
| • detect   | r peaks() |     |     |
•
| calculate     | rr intervals()    |          |     |
| ------------- | ----------------- | -------- | --- |
| • interpolate | rr                | signal() |     |
| • apply       | window function() |          |     |
•
| calculate   | fft() |         |     |
| ----------- | ----- | ------- | --- |
| • calculate | hrv   | bands() |     |
•
| plot   | spectrum()  |     |     |
| ------ | ----------- | --- | --- |
| • plot | waterfall() |     |     |
•
| compare          | window       | functions()       |     |
| ---------------- | ------------ | ----------------- | --- |
| 5 Erwartete      | Ergebnisse   |                   |     |
| Am Ende          | des Projekts | sollen vorliegen: |     |
| • vollst¨andiger | MATLAB-Code  |                   |     |
6

| • Dokumentation      | der | Implementierung       |     |     |
| -------------------- | --- | --------------------- | --- | --- |
| • Vergleichsgrafiken |     | der Fensterfunktionen |     |     |
• FFT-Spektren
• Wasserfalldiagramme
| • Tabellen  | der HRV-Kennwerte |             |     |     |
| ----------- | ----------------- | ----------- | --- | --- |
| • fachliche | Bewertung         | der Fenster |     |     |
6 Leitfragen
| Fachliche | Fragestellungen |     |     |     |
| --------- | --------------- | --- | --- | --- |
1. Wie beeinflussen unterschiedliche Fensterfunktionen die HRV-Spektralanalyse?
| 2. Welche | Fensterfunktion | liefert | den besten Kompromiss | zwischen: |
| --------- | --------------- | ------- | --------------------- | --------- |
•
Frequenzaufl¨osung
• Leakage-Unterdru¨ckung
•
|           | Stabilit¨at der | HRV-Kennwerte? |                                  |     |
| --------- | --------------- | -------------- | -------------------------------- | --- |
| 3. Welche | Auswirkungen    | zeigen sich    | in der 3D-Wasserfalldarstellung? |     |
7
