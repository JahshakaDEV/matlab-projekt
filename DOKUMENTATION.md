# Dokumentation: Einfluss von FFT-Fensterfunktionen auf die HRV-Analyse

Diese Datei dokumentiert die Implementierung in [`hrv_fft_analyse.m`](hrv_fft_analyse.m)
und beantwortet die fachlichen Leitfragen aus der Aufgabenstellung.

Der gesamte Code ist **modular** aufgebaut: Das Skript ganz oben ruft nur die
Verarbeitungsschritte der Reihe nach auf, die eigentliche Logik steckt in den
lokalen Funktionen am Dateiende. Es werden durchgehend **eingebaute MATLAB- und
Signal-Processing-Toolbox-Funktionen** verwendet, statt Algorithmen selbst
nachzubauen.

---

## 1. Benötigte Software

- MATLAB
- **Signal Processing Toolbox** (für `edfread`, `highpass`, `lowpass`,
  `iirnotch`, `findpeaks`, `hann`, `hamming`, `blackman`, `kaiser`,
  `flattopwin`, `enbw`)

Start: Datei öffnen, EDF-Datei im Parameter `EDF_FILE` eintragen, ausführen (F5).

---

## 2. Verarbeitungskette im Überblick

```
load_ecg_data      →  EDF einlesen (Rohsignal + Abtastrate)
preprocess_ecg     →  DC entfernen, Baseline-, Netz-, Rausch-Filter
detect_r_peaks     →  R-Zacken finden
calculate_rr_intervals → RR-Intervalle + Plausibilitätsprüfung
interpolate_rr_signal  → gleichmäßig abgetastetes HRV-Signal (4 Hz)
compare_window_functions → je Fenster: Segmentierung, FFT, Mittelung
   ├─ apply_window_function → Fenster auf Segment legen
   ├─ calculate_fft         → Leistungsdichtespektrum (PSD)
   └─ calculate_hrv_bands   → VLF/LF/HF-Leistung
plot_spectrum / plot_waterfall → Darstellung
```

---

## 3. Zuordnung: Aufgabenpunkt → Umsetzung im Code

| PDF-Abschnitt | Anforderung | Umsetzung (Funktion / Zeile) | Genutzte native Funktion |
|---|---|---|---|
| 2.1 | EDF ≥ 1 h importieren | `load_ecg_data` | `edfinfo`, `edfread` |
| 2.1 | Gleichanteil entfernen | `preprocess_ecg`: `ecg - mean(ecg)` | `mean` |
| 2.1 | Baseline-Wanderung filtern | `preprocess_ecg`: `highpass(...,0.5,fs)` | `highpass` |
| 2.1 | Netzbrummen (50 Hz) filtern | `preprocess_ecg`: `iirnotch` + `filtfilt` | `iirnotch`, `filtfilt` |
| 2.1 | Hochfrequentes Rauschen filtern | `preprocess_ecg`: `lowpass(...,40,fs)` | `lowpass` |
| 2.1 | Roh- + gefiltertes Signal plotten | Skriptteil 2 | `plot`, `subplot` |
| 2.2 | R-Zacken automatisch erkennen | `detect_r_peaks` | `findpeaks` |
| 2.2 | RR-Intervalle berechnen | `calculate_rr_intervals`: `diff(t_r)` | `diff` |
| 2.2 | Fehlerhafte Intervalle erkennen / Plausibilität | `calculate_rr_intervals` (feste Grenzen + lokaler Median) | `movmedian` |
| 2.2 | R-Peaks markiert + RR-Zeitreihe plotten | Skriptteil 3 + 4 | `plot` |
| 2.3 | RR-Reihe interpolieren, ~4 Hz | `interpolate_rr_signal` (Spline) | `interp1`, `detrend` |
| 2.4 | FFT mit ≥ 4 Fenstern | `apply_window_function` + `calculate_fft` | `fft`, `hann`, `hamming`, `blackman`, `kaiser`, `flattopwin` |
| 2.5 | Fenstereigenschaften (Auflösung, Leakage, Amplitude) | `window_metrics`, `plot_window_characteristics`, `plot_overlaid_spectra` | `enbw`, `fft` |
| 2.6 | HRV-Bänder VLF/LF/HF, absolut/relativ, LF/HF | `calculate_hrv_bands` | `trapz` |
| 2.7 | 3D-Wasserfall je Fenster, HRV-Bänder markiert, einheitliche Skalierung | `plot_waterfall` (gemeinsames `zmax_all`) | `waterfall` |
| 2.8 | Vergleich / Bewertung | `print_comparison_table` + dieses Dokument (Abschnitt 6) | `fprintf` |
| 3 | Segment 5 min, 50 % Überlappung | Parameter `SEG_LEN_S=300`, `OVERLAP=0.5` | — |
| 4 | Alle 11 Mindestfunktionen | vorhanden (siehe Abschnitt 4) | — |

---

## 4. Erklärung der einzelnen Funktionen

### `load_ecg_data(filename)`
Liest die EDF-Datei. Die **Abtastrate** wird aus dem Header berechnet
(`NumSamples / DataRecordDuration`), weil EDF die Samplerate nicht direkt
speichert. `edfread` liefert die Daten als Timetable; der erste Kanal wird zu
einem durchgehenden Spaltenvektor zusammengefügt.

### `preprocess_ecg(ecg, fs)`
Vier Filterstufen, alle **nullphasig** mit `filtfilt` bzw. den nullphasigen
`highpass`/`lowpass` — dadurch werden die R-Zacken zeitlich **nicht verschoben**:
1. `ecg - mean(ecg)` → Gleichanteil (DC) weg.
2. `highpass(0.5 Hz)` → entfernt langsame Baseline-Wanderung (Atmung, Bewegung).
3. `iirnotch(50 Hz)` → schmalbandiger Kerbfilter gegen Netzbrummen.
4. `lowpass(40 Hz)` → dämpft hochfrequentes Rauschen (z. B. Muskelartefakte).

### `detect_r_peaks(ecg, fs)`
`findpeaks` mit zwei Kriterien:
- `MinPeakHeight = mean + 0.6·std` → nur deutlich herausragende Spitzen (R-Zacken).
- `MinPeakDistance = 0.3·fs` → physiologische Grenze: kein zweiter Herzschlag
  innerhalb von 300 ms (max. ~200 bpm).

### `calculate_rr_intervals(r_locs, fs)`
- RR-Intervall = zeitlicher Abstand zweier R-Zacken in Millisekunden (`diff`).
- **Plausibilitätsprüfung** mit zwei einfachen Bedingungen: (a) feste Grenzen
  300–2000 ms (≈ 30–200 bpm), (b) höchstens 20 % Abweichung vom **lokalen
  Median** (`movmedian` über 5 Intervalle). So werden Ausreißer/Extrasystolen
  erkannt.
- Erkannte Artefakte werden per `interp1(..., 'pchip')` aus den gültigen
  Nachbarn ersetzt und die Anzahl zurückgegeben.

### `interpolate_rr_signal(t_rr, rr, fs_i)`
RR-Werte liegen zu **unregelmäßigen** Zeitpunkten vor, die FFT braucht aber
**gleiche** Abstände. Deshalb Interpolation auf ein festes Raster (Standard
`fs_i = 4 Hz`, ausreichend für das HRV-Band bis 0,4 Hz → Nyquist 2 Hz).
`interp1(..., 'spline')` erzeugt einen glatten Verlauf, `detrend` entfernt einen
linearen Trend (verhindert einen künstlichen Peak bei 0 Hz).

### `apply_window_function(x, window_type)`
Legt das gewählte Fenster über ein Segment. Alle Fenster sind eingebaute
Funktionen: `ones` (Rechteck), `hann`, `hamming`, `blackman`, `kaiser(N,8)`,
`flattopwin`. Rückgabe: gefenstertes Signal **und** Fenstervektor (Letzterer
wird zur Leistungsnormierung gebraucht).

### `calculate_fft(xw, win, fs)`
Bildet das **einseitige, leistungsnormierte** Spektrum (PSD):
`PSD = |FFT|² / (fs · Σwin²)`. Die Normierung auf die Fensterleistung `Σwin²`
sorgt dafür, dass die **integrierte** Bandleistung weitgehend fensterunabhängig
bleibt (fairer Vergleich). Innere Bins werden verdoppelt (einseitiges Spektrum
eines reellen Signals).

### `calculate_hrv_bands(f, psd)`
Integriert die PSD (`trapz` = Fläche unter der Kurve) über die drei Bänder
(Task Force 1996): VLF 0,0033–0,04 Hz, LF 0,04–0,15 Hz, HF 0,15–0,40 Hz.
Berechnet zusätzlich **relative** Leistung (% der Gesamtleistung), normalisierte
Einheiten (n.u.) und das **LF/HF-Verhältnis**.

### `compare_window_functions(...)`
Kern der Untersuchung. Zerlegt das Signal in **5-Minuten-Segmente mit 50 %
Überlappung** (Welch-Verfahren). Für jedes Fenster wird jedes Segment gefenstert,
per FFT in eine PSD überführt, und die PSDs werden **gemittelt** (senkt die
Varianz). Ergebnis pro Fenster: gemittelte PSD, alle Segment-Spektren (für den
Wasserfall), HRV-Kennwerte und Fensterkenngrößen. Ruft am Ende die
Vergleichsplots und die Konsolentabelle auf.

### `window_metrics(win)`
Kenngrößen jedes Fensters, alle aus dem (durch Zero-Padding fein aufgelösten)
Frequenzgang des Fensters mit eingebauten Funktionen bestimmt:
- **PSL** (Peak Side Lobe, dB): höchste Nebenkeule → Maß für Leakage-Unterdrückung.
  Ermittelt als größtes lokales Maximum: `max(findpeaks(WdB))` (der Hauptpeak bei
  Bin 0 ist ein Randpunkt und wird von `findpeaks` nicht mitgezählt).
- **MLW** (Main Lobe Width, Bins): Breite der Hauptkeule → Maß für
  Frequenzauflösung. = 2 × erste Nullstelle, gefunden mit `islocalmin`.
- **ENBW** (`enbw`): äquivalente Rauschbandbreite.

### Plot-Funktionen
- `plot_spectrum` — Einzelspektrum mit farblich hinterlegten HRV-Bändern und
  Kennwert-Textbox.
- `plot_window_characteristics` — Fensterform (Zeit) + Frequenzgang (dB).
- `plot_overlaid_spectra` — alle Fenster überlagert (linear + logarithmisch),
  zeigt direkt Amplituden- und Leakage-Unterschiede.
- `plot_waterfall` — 3D-Wasserfall; nutzt ein **gemeinsames z-Maximum**
  (`zmax_all`), damit alle Fenster auf **identischer Skala** vergleichbar sind.
- `print_comparison_table` — Tabelle der Kenngrößen im Command Window.

---

## 5. Technische Parameter (Abschnitt 3 der Aufgabe)

| Parameter | Wert | Ort im Code |
|---|---|---|
| Interpolations-/Abtastrate | 4 Hz | `FS_INTERP` |
| Segment-/Fensterbreite | 300 s (5 min) | `SEG_LEN_S` |
| Überlappung | 50 % | `OVERLAP` |
| Untersuchte Fenster | rect, hann, hamming, blackman, kaiser, flattop | `WINDOWS` |

---

## 6. Fachliche Bewertung der Fenster (Abschnitt 2.8) + Leitfragen

Grundsätzlicher **Zielkonflikt**: Ein schmaler Hauptpeak (gute
Frequenzauflösung) und stark gedämpfte Nebenkeulen (wenig Leakage) sind nicht
gleichzeitig maximal erreichbar. Jedes Fenster ist ein Kompromiss.

| Fenster | Frequenzauflösung (MLW) | Leakage-Unterdrückung (PSL) | Amplitudengenauigkeit | Rechenaufwand |
|---|---|---|---|---|
| Rechteck | **am besten** (schmalste Hauptkeule) | **am schlechtesten** (~ −13 dB) | mittel | am geringsten |
| Hamming | gut | mittel (~ −43 dB) | gut | gering |
| Hann | mittel | gut (~ −31 dB, fällt schnell ab) | gut | gering |
| Blackman | schlechter (breite Hauptkeule) | **sehr gut** (~ −58 dB) | gut | gering |
| Kaiser (β=8) | einstellbar | sehr gut | gut | gering |
| Flat-Top | **am schlechtesten** (sehr breit) | gut | **am besten** (Amplitude) | gering |

**Leitfrage 1 – Einfluss auf die HRV-Spektralanalyse:**
Das Rechteckfenster zeigt die schärfsten Peaks, „verschmiert" durch starkes
Leakage aber Energie in benachbarte Bänder — das kann besonders die schwache
HF-Leistung verfälschen und das LF/HF-Verhältnis verzerren. Fenster mit guter
Nebenkeulenunterdrückung (Hann, Blackman, Kaiser) liefern **stabilere und
sauberer getrennte** VLF/LF/HF-Werte, dafür etwas breitere Peaks.

**Leitfrage 2 – Bester Kompromiss:**
Für die HRV-Frequenzanalyse ist das **Hann-Fenster** der übliche und
empfehlenswerte Kompromiss: ordentliche Auflösung, gute Leakage-Unterdrückung
und stabile Bandkennwerte. Wenn Leakage besonders kritisch ist, bietet
**Kaiser (β≈8)** einen einstellbaren, sehr guten Kompromiss. Rechteck ist wegen
des Leakage ungeeignet; Flat-Top ist nur sinnvoll, wenn es rein auf die genaue
**Amplitude** eines Peaks ankommt (bei HRV selten).

**Leitfrage 3 – Auswirkungen im 3D-Wasserfall:**
Beim Rechteckfenster erscheinen die Spektren „zackig" und unruhig, mit breiten
Ausläufern zwischen den Bändern (Leakage über die Zeit sichtbar). Hann/Blackman/
Kaiser erzeugen **glattere, ruhigere** Wasserfälle mit klar getrennten Bändern —
gut für die Interpretierbarkeit. Flat-Top verwischt durch die breite Hauptkeule
benachbarte Frequenzen. Dank **einheitlicher z-Skalierung** aller Diagramme sind
diese Unterschiede direkt miteinander vergleichbar.

---

## 7. Erzeugte Ergebnisse (Abschnitt 5 der Aufgabe)

Beim Ausführen entstehen:
- Roh- vs. gefiltertes EKG (Plot),
- EKG mit markierten R-Zacken (Plot),
- RR-Tachogramm (Plot),
- Fenster-Kenngrößen im Zeit- und Frequenzbereich (Plot),
- überlagerte HRV-Spektren linear + logarithmisch (Vergleichsgrafiken),
- ein Einzelspektrum je Fenster (Plot mit HRV-Bändern),
- ein 3D-Wasserfalldiagramm je Fenster (einheitlich skaliert),
- Vergleichstabelle der HRV-Kennwerte + Fensterkenngrößen im Command Window.
