# Vortragstext — FFT-Fensterfunktionen in der HRV-Analyse

**Präsentation:** HRV_Praesentation.pptx · DHBW · MATLAB · Sommersemester 2026
**Vortragende:** Jannik Oeschger · Arvid Haase
**Gesamtdauer:** ca. 8 Minuten · Schwerpunkt: Live-Demo von Code & Plots in MATLAB

**Rollen:** **Jannik** = Person 1 · **Arvid** = Person 2
Regieanweisungen stehen *(kursiv in Klammern)*. In der Live-Demo bedient jeweils die sprechende Person MATLAB.

Grober Redeanteil: Jannik ≈ 3:50 · Arvid ≈ 4:10 — mit Übergabepunkt mitten in der Demo.

---

## 🟦 Folie 1 — Titel · Jannik · (0:00–0:30)

> Hallo zusammen, schön, dass ihr da seid. Ich bin Jannik, das ist Arvid — und in unserem Projekt geht es um **FFT-Fensterfunktionen in der HRV-Analyse**. Die Kernfrage: Wie verändert allein die *Wahl der Fensterfunktion* die Kennwerte der Herzfrequenzvariabilität? Wir halten die Folien bewusst kurz — der Schwerpunkt liegt heute auf einer **Live-Demo** von Code und Ergebnissen direkt in MATLAB.

## 🟦 Folie 2 — Worum geht es? · Arvid · (0:30–1:45)

> Kurz zur Einordnung. Unser Herz schlägt nie exakt gleichmäßig — die Abstände zwischen zwei Schlägen schwanken ständig ein bisschen. Das ist völlig gesund, und genau in diesen Schwankungen stecken **Rhythmen**: die Atmung zum Beispiel, oder die Blutdruckregelung. Diese Herzfrequenzvariabilität — kurz **HRV** — verrät etwas über unser vegetatives Nervensystem.
>
> Um diese Rhythmen sichtbar zu machen, zerlegen wir das Signal mit der **FFT** in seine Frequenzen — im Prinzip wie ein Prisma, das weißes Licht in einen Regenbogen aufteilt, oder wie ein geübtes Ohr, das aus einem Akkord die Einzeltöne heraushört.
>
> Der Haken: Wir haben immer nur einen *endlichen Ausschnitt*. Und der erzeugt einen Effekt namens **Spectral Leakage**, der die Werte verfälscht. **Fensterfunktionen** dämpfen genau dieses Leakage. Unsere Leitfrage — und der rote Faden für die nächsten Minuten — lautet deshalb: **Welches Fenster ist der beste Kompromiss?**

## 🟦 Folie 3 — Die Verarbeitungskette · Jannik · (1:45–3:15)

> Bevor wir überhaupt Fenster vergleichen können, muss aus dem rohen EKG erst ein sauberes Signal werden. Das ist unsere Verarbeitungskette — sechs Schritte:
>
> **Erstens** laden wir das EKG aus einer EDF-Datei. **Zweitens** die Vorverarbeitung: Gleichanteil raus, Baseline-Wanderung durch die Atmung, das 50-Hz-Netzbrummen und hochfrequentes Rauschen. **Drittens** erkennen wir die R-Zacken — die markanten Ausschläge jedes Herzschlags — und berechnen daraus die **RR-Intervalle**, also die Abstände von Schlag zu Schlag. **Viertens** interpolieren wir diese Abstände auf ein gleichmäßiges Raster mit **4 Hertz** — wichtig, weil die FFT gleichmäßige Abstände braucht, die Herzschläge aber unregelmäßig kommen. **Fünftens** kommt das eigentliche Experiment: die FFT, **sechsmal** gerechnet, einmal mit jedem Fenster. Und **sechstens** werten wir die HRV-Bänder aus und erzeugen die Wasserfalldiagramme.
>
> Wichtig ist der **rote Block** in der Mitte: Alles davor ist bei jedem Durchlauf identisch — nur das Fenster in Schritt fünf ändert sich. Das *ist* unser Experiment. Was so ein Fenster überhaupt ist, zeigt euch Arvid.

## 🟦 Folie 4 — Fensterfunktion = ein Kompromiss · Arvid · (3:15–4:15)

> Ein Fenster ist im Grunde eine **Glockenkurve**, mit der wir das Signalstück multiplizieren, damit es an den Rändern sanft auf null geht — dadurch verschwindet der Sprung, der das Leakage verursacht. Aber: Ein Fenster gibt's nicht geschenkt. Jedes Fenster ist ein **Kompromiss** zwischen zwei Eigenschaften.
>
> Links die **Hauptkeule** — der eigentliche Peak einer Frequenz. Ist sie schmal, können wir zwei dicht benachbarte Frequenzen noch trennen: gute **Auflösung**. Rechts die **Nebenkeulen** — das unerwünschte Auslaufen der Energie, also das **Leakage**. Sind sie niedrig, ist das Spektrum sauber.
>
> Und jetzt das Naturgesetz dahinter: **Man kann nicht beides gleichzeitig haben.** Schmalere Hauptkeule → höhere Nebenkeulen, und umgekehrt. Jedes Fenster wählt nur einen anderen Punkt auf dieser Waage. Das **Rechteck** ist der Extremfall „kein Fenster" — super scharf, aber es leckt stark. Das **Flat-Top** ist das Gegenteil — extrem sauber, aber sehr breit. Wir untersuchen sechs Stück: Rechteck, Hann, Hamming, Blackman, Kaiser und Flat-Top. Und wie sich das auswirkt, schauen wir uns am besten direkt in MATLAB an.

## 🟩 Folie 5 — LIVE-DEMO · beide · (4:15–6:45)

**Jannik** *(bedient MATLAB, Schritte 1–4):*

> So — das ist unser Code. Ihr seht: **modular** aufgebaut, jeder Schritt der Kette ist eine eigene Funktion, genau wie in der Aufgabe gefordert. Ich starte ihn einmal.
> *(Skript starten)*
> Während er läuft, ein Blick auf die Konsole: Das EKG ist mit **1000 Hertz** abgetastet, **knapp 40 Minuten** lang. Er findet rund **4.500 R-Zacken** und korrigiert dabei **50 Artefakte** automatisch — einzelne Ausreißer und Fehldetektionen.
> *(auf EKG-Plot)* Hier oben das Signal: grau das Rohsignal, rot nach der Vorverarbeitung — die Baseline ist weg, alles liegt sauber um null, und die R-Zacken sitzen exakt auf den Spitzen.
> *(auf Tachogramm)* Und das ist das **Tachogramm** — die RR-Intervalle über die Zeit. Man sieht richtig, wie die Herzrate lebt, mal schneller, mal langsamer. Das ist die HRV, um die's geht. — Arvid, übernimmst du?

**Arvid** *(bedient MATLAB, Schritte 5–8):*

> Gerne. *(auf Fensterformen)* Hier links die sechs **Fensterformen** im Zeitbereich, rechts ihr Frequenzgang. Schaut aufs **Rechteck**, die blaue Kurve: hohe Nebenkeulen — es leckt am stärksten. Das **Flat-Top**, hellblau: Nebenkeulen ganz unten — sehr sauber, aber die Hauptkeule ist deutlich breiter.
> *(auf dB-Spektrum)* Und jetzt die **wichtigste Grafik** — die überlagerten Spektren in Dezibel. Achtet auf den Boden: Die blaue **Rechteck**-Kurve liegt überall **am höchsten** — das ist der Leakage-„Bodennebel". Die hellblaue **Flat-Top**-Kurve liegt **am tiefsten** — am saubersten.
> *(durch Wasserfälle klicken)* Dasselbe in den **Wasserfällen**: mit dem Rechteck ist die Landschaft zackig und unruhig, mit Hann oder Flat-Top wird sie glatt und die Bänder sind sauber getrennt. Alle nutzen **dieselbe Skala**, damit der Vergleich fair ist.
> *(auf Konsolen-Tabelle)* Und die Tabelle bringt's auf Zahlen — schaut mal besonders auf die **HF-Spalte** und auf **LF/HF**.

*(Fallback: Falls MATLAB zickt — `all_figures.pdf` offen halten und daraus zeigen.)*

## 🟦 Folie 6 — Was wir sehen · Arvid · (6:45–7:30)

> Und genau da steckt die Kernaussage. Unser HRV-Signal hat sehr **starke langsame** Anteile und nur **schwache schnelle**. Das Rechteck lässt durch das Leakage die starke langsame Energie ins schwache **HF-Band** hinüberlecken. Das Ergebnis: Beim Rechteck kommt ein HF von **rund 31** raus, beim Flat-Top nur **rund 8** — also **viermal so viel**, obwohl es dasselbe EKG, dieselbe Person, dieselbe Rechnung ist. Nur das Fenster ist anders.
>
> Und weil HF im Nenner steht, kippt damit auch das **LF/HF-Verhältnis** — von 3,4 auf 5,5. Und das ist eine Zahl, die tatsächlich *medizinisch interpretiert* wird. Die Botschaft ist also: **Die Fensterwahl ist keine Kosmetik — ein falsches Fenster liefert einen falschen Kennwert.**

## 🟦 Folie 7 — Fazit · Jannik · (7:30–7:50)

> Damit zurück zu unserer Leitfrage: Welches Fenster nehmen? Unsere Empfehlung ist das **Hann-Fenster** — der beste Allround-Kompromiss für die HRV: ordentliche Auflösung, wenig Leakage, stabile Bandwerte. Wer feinjustieren will, nimmt **Kaiser mit β ≈ 8**. Das **Rechteck** taugt nur als Negativbeispiel, und das **Flat-Top** nur, wenn es rein auf die exakte Höhe eines Peaks ankommt — bei der HRV also selten.

## 🟦 Folie 8 — Vielen Dank! · beide · (7:50–8:00)

> **Jannik:** Das war's von uns — vielen Dank fürs Zuhören!
> **Arvid:** Wir freuen uns auf eure Fragen.

---

## 🎤 Backup für die Fragerunde

- **Aufnahmedauer:** Datensatz ~40 min (die Aufgabe fordert ≥ 1 h — bewusst so gewählt / verfügbar).
- **Segmentierung:** Fensterbreite 5 min, 50 % Überlappung, **Welch-Mittelung** senkt die Varianz der Schätzung.
- **Bandgrenzen** nach **Task Force 1996**: VLF 0,0033–0,04 Hz · LF 0,04–0,15 Hz · HF 0,15–0,40 Hz.
- **Interpolation:** Spline auf 4 Hz, linearer Trend entfernt (`detrend`), damit langsame Drift die VLF/LF nicht verfälscht.
- **R-Zacken:** `findpeaks` mit Mindesthöhe (Mittelwert + 0,6·Std) und Mindestabstand 0,3 s (≈ 200 bpm-Grenze).
- **Artefaktkorrektur:** Plausibilität 300–2000 ms + gleitender Median (`isoutlier`, movmedian 21), Ersatz per `pchip`.
