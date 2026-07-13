# FFT-Fensterfunktionen einfach erklärt

> Eine Einführung für Leute, die davon **noch nie** gehört haben.
> Ziel: verstehen, worum es in diesem Projekt geht – ohne Vorwissen in Signalverarbeitung.

Wir bauen das Thema Schritt für Schritt auf. Wenn du am Ende die Begriffe
**Spektrum, Leakage, Hauptkeule, Nebenkeule** und **Fensterfunktion** in eigenen
Worten erklären kannst, hat die Datei ihren Zweck erfüllt.

---

## 1. Das große Ganze in drei Sätzen

Unser Herz schlägt nicht wie ein Metronom – die Abstände zwischen den Schlägen
schwanken ständig ein bisschen. In diesen Schwankungen stecken **verborgene
Rhythmen** (z. B. von Atmung und Blutdruckregelung), und die verraten etwas über
unser Nervensystem. Um diese Rhythmen sichtbar zu machen, zerlegen wir das Signal
mit der **FFT** in seine Frequenzen – und genau dabei brauchen wir
**Fensterfunktionen**, damit das Ergebnis nicht verfälscht wird.

Dieses Projekt untersucht: **Welche Fensterfunktion liefert das ehrlichste Ergebnis?**

---

## 2. Was ist überhaupt eine „Frequenz"?

Eine **Frequenz** sagt, *wie oft sich etwas pro Sekunde wiederholt*. Einheit:
**Hertz (Hz)**. 1 Hz = einmal pro Sekunde, 2 Hz = zweimal pro Sekunde, usw.

Ein einfaches, gleichmäßiges Schwingen nennt man eine **Sinuswelle**:

```
  +1 |   .-.        .-.        .-.
     |  /   \      /   \      /   \
   0 | /     \    /     \    /     \
     |/       \  /       \  /       \
  -1 |         '-'        '-'         '-'
     +------------------------------------> Zeit
      <-- eine Wiederholung -->
```

Der Trick der Mathematik: **Jedes** beliebige Signal – egal wie kompliziert –
lässt sich als **Summe vieler solcher Sinuswellen** unterschiedlicher Frequenz
zusammenbauen. Manche schnell, manche langsam, manche stark, manche schwach.

**Analogie – der Musik-Akkord:** Wenn auf dem Klavier ein Akkord erklingt, hörst
du *einen* Klang. Tatsächlich sind das aber mehrere Töne (Frequenzen)
gleichzeitig. Ein geübtes Ohr kann sie heraushören und sagen: „Das ist ein C, ein
E und ein G." Genau das macht die FFT für uns automatisch.

---

## 3. Die FFT – das Signal in seine Frequenzen zerlegen

**FFT** steht für *Fast Fourier Transform* (schnelle Fourier-Transformation). Man
gibt ihr ein Signal (den zeitlichen Verlauf) und bekommt zurück, **welche
Frequenzen mit welcher Stärke** enthalten sind. Dieses Ergebnis heißt
**Spektrum**.

```
   ZEIT-Darstellung                    FREQUENZ-Darstellung (Spektrum)
   (was wir messen)          FFT       (was drinsteckt)
                            ------>
   Signal                             Stärke
     |  /\    /\  /\                     |        |
     | /  \  /  \/  \                    |        |            |
     |/    \/        \                   |    |   |            |
     +----------> Zeit                   +----+---+------------+---> Frequenz
                                              langsam        schnell
```

**Analogie – das Prisma:** Weißes Licht sieht einfarbig aus. Schickt man es durch
ein Prisma, zerfällt es in einen Regenbogen und man sieht, aus welchen Farben
(= Frequenzen des Lichts) es besteht. Die FFT ist das Prisma für unser Signal.

**Bei der HRV** interessieren uns drei Frequenzbereiche („Bänder"):

| Band | Frequenz | Wovon kommt der Rhythmus? |
|------|----------|---------------------------|
| **VLF** | 0,0033 – 0,04 Hz | sehr langsame Prozesse (Hormone, Temperatur) |
| **LF**  | 0,04 – 0,15 Hz  | Blutdruckregelung (Baroreflex), ~0,1 Hz |
| **HF**  | 0,15 – 0,40 Hz  | Atmung (~0,25 Hz ≈ 15 Atemzüge/min) |

Das Verhältnis **LF/HF** wird als grober Indikator für das Gleichgewicht des
vegetativen Nervensystems benutzt. Wenn unsere Analyse LF und HF falsch berechnet,
ist auch dieses Verhältnis falsch – deshalb ist Genauigkeit hier so wichtig.

---

## 4. Das Problem: „Spectral Leakage" (Spektralleckeffekt)

Jetzt kommt der Haken – und der ist der **Kern des ganzen Projekts**.

Die FFT hat eine versteckte Annahme: Sie tut so, als würde sich das kurze Stück
Signal, das man ihr gibt, **unendlich oft nahtlos wiederholen**. Sie klebt in
Gedanken viele Kopien hintereinander.

Bei einem echten Messausschnitt passt das **Ende aber nicht zum Anfang**. An der
Klebestelle entsteht ein **Sprung**:

```
  Ein Ausschnitt ...        ... in Gedanken aneinandergeklebt:

   /\    /\    /\            /\    /\    /\ | /\    /\    /\
  /  \  /  \  /  \          /  \  /  \  /  \|/  \  /  \  /  \
 /    \/    \/    \        /    \/    \/    X    \/    \/    \
                                           ^
                                     SPRUNG! (Ende ≠ Anfang)
```

Dieser künstliche Sprung gehört gar nicht zum echten Signal – aber die FFT „sieht"
ihn und muss ihn mit **zusätzlichen, falschen Frequenzen** erklären. Das Ergebnis:
Die Energie einer einzelnen, sauberen Frequenz wird **breit über die Nachbar­
frequenzen verschmiert**. Dieses „Auslaufen" der Energie heißt **Spectral
Leakage** (Leck­effekt).

**Analogie – der harte Schnitt im Song:** Wenn du ein Lied mittendrin hart
abschneidest und in einer Schleife abspielst, hörst du bei jedem Durchlauf ein
störendes **„Knacken"**. Dieses Knacken ist ein Geräusch, das im Original gar nicht
vorkommt – es entsteht nur durch den Sprung an der Schnittstelle. Leakage ist das­
selbe Phänomen, nur im Spektrum statt im Ton.

**Warum ist das für die HRV schlimm?** Unser HRV-Signal hat sehr **starke langsame**
Anteile (VLF/LF) und nur **schwache schnelle** Anteile (HF). Durch Leakage
„leckt" die starke langsame Energie in das schwache HF-Band hinüber – und dann
sieht HF viel größer aus, als es wirklich ist. (Die konkreten Zahlen dazu kommen
in Abschnitt 9.)

---

## 5. Die Lösung: Fensterfunktionen

Idee: Wenn der **Sprung an den Rändern** das Problem ist, dann sorgen wir dafür,
dass das Signal **an beiden Rändern sanft auf null** geht. Dann passt das Ende zum
Anfang (beide sind ja null), und es gibt keinen Sprung mehr.

Dazu multiplizieren wir das Signalstück mit einer **glockenförmigen Kurve**, die
in der Mitte 1 ist und zu den Rändern hin auf 0 abfällt. Diese Kurve ist die
**Fensterfunktion** (englisch *window*).

```
  Signal (Ausschnitt)   ×   Fenster (Glocke)   =   gefenstertes Signal

   /\  /\  /\  /\             ___                    _
  /  \/  \/  \/  \    ×      /   \          =       / \  /\  _
 /              \          /       \               /   \/  \/ \
 (Ränder "springen")     (0 an den Rändern)     (Ränder sanft auf 0)
```

**Analogie – der DJ-Übergang:** Ein guter DJ schneidet ein Lied nicht hart ab,
sondern **blendet es sanft ein und aus** (Fade-in / Fade-out). Genau das macht ein
Fenster mit unserem Signal – und schon ist das „Knacken" (das Leakage) weg oder
zumindest stark reduziert.

Das **Rechteckfenster** ist übrigens der Spezialfall „**kein** Fenster": man nimmt
den Ausschnitt so, wie er ist (überall 1). Es dient als Vergleichsmaßstab und
zeigt, wie schlimm Leakage *ohne* Gegenmaßnahme ist.

---

## 6. Der Haken an der Lösung: der Kompromiss

Leider ist ein Fenster **kein Gratis-Geschenk**. Um das zu verstehen, schauen wir
uns an, wie eine **einzelne, perfekt saubere Frequenz** im Spektrum aussieht,
nachdem wir ein Fenster benutzt haben. Sie erscheint **nicht** als unendlich dünner
Strich, sondern als eine Form mit einem großen Buckel und kleinen Nebenbuckeln:

```
  Betrag [dB]
    0 dB ->        _____
                  /     \        <-- HAUPTKEULE (engl. main lobe)
                  |     |            = der "echte" Peak.
                  |     |            Ihre BREITE = Frequenzauflösung.
                  |     |
                  |     |
  -13 dB ->  __   |     |   __
            /  \  |     |  /  \     <-- NEBENKEULEN (engl. side lobes)
   ________/    \_|     |_/    \______   = das unerwünschte "Auslaufen".
                                          Ihre HÖHE = wie stark Leakage ist.
            Frequenz  ------>
```

Es gibt also **zwei** Eigenschaften, die gegeneinander arbeiten:

- **Hauptkeulenbreite** → bestimmt die **Frequenzauflösung**
  (Kann ich zwei *dicht benachbarte* Frequenzen noch als zwei getrennte Peaks
  erkennen? Schmale Hauptkeule = ja, breite = sie verschmelzen zu einem.)

- **Nebenkeulenhöhe** → bestimmt die **Leakage-Unterdrückung**
  (Wie viel Energie läuft in fremde Frequenzbereiche aus? Niedrige Nebenkeulen =
  wenig, hohe = viel.)

**Das Naturgesetz dahinter:** Man kann **nicht beides gleichzeitig** optimal haben.
Macht man die Hauptkeule schmaler (bessere Auflösung), werden die Nebenkeulen höher
(mehr Leakage) – und umgekehrt. **Jede Fensterfunktion ist ein Kompromiss** an
einer anderen Stelle dieser Waage.

**Analogie – die Zoom-Linse:** Weitwinkel zeigt dir viel Übersicht, aber wenig
Detail; starker Zoom zeigt viel Detail, aber wenig Übersicht. Du musst dich
entscheiden. Ein Fenster ist die Wahl der „Linse" für dein Spektrum.

### Die drei Kennzahlen im Code

In `window_metrics` misst das Programm diesen Kompromiss mit drei Zahlen:

- **MLW** (*Main Lobe Width*, in **Bins**): Breite der Hauptkeule. **Kleiner = bessere Auflösung.**
- **PSL** (*Peak Side Lobe*, in **dB**): Höhe der größten Nebenkeule. **Negativer = besser** (weniger Leakage).
- **ENBW** (*Equivalent Noise Bandwidth*, in Bins): fasst zusammen, wie sehr ein Fenster die Leistung „verschmiert". **Kleiner = konzentrierter.**

**Kurz erklärt – „Bin":** Die FFT liefert das Spektrum nicht stufenlos, sondern auf
einem festen **Frequenzraster**. Ein **Bin** ist ein Schritt dieses Rasters – die
kleinste Frequenzstufe, die man unterscheiden kann. „Hauptkeule = 2 Bins breit"
heißt also: sie belegt zwei Rasterschritte.

**Kurz erklärt – „Dezibel (dB)":** Eine handliche Schreibweise für Größenverhältnisse,
die riesige Bereiche abdeckt. Faustregel für Höhen: **jede −20 dB bedeuten
„zehnmal kleiner".**

| dB | Nebenkeule ist ... vom Hauptpeak | Bewertung |
|-----|-----------------------------------|-----------|
| −13 dB | ≈ ¼ so hoch | schlecht (viel Leakage) |
| −40 dB | ≈ 1/100 so hoch | gut |
| −60 dB | ≈ 1/1000 so hoch | sehr gut |
| −90 dB | ≈ 1/30 000 so hoch | exzellent |

---

## 7. Die einzelnen Fenster im Überblick

> Sieh dir dazu **Seite 4** in `all_figures.pdf` an: links die Fensterformen (Zeit),
> rechts ihre Frequenzgänge (die Haupt- und Nebenkeulen von oben).

- **Rechteck** – „kein Fenster". **Schmalste** Hauptkeule (beste Auflösung),
  aber **höchste** Nebenkeulen (−13 dB → viel Leakage). Der schlechte Schüler beim
  Leakage, aber der Maßstab.

- **Hann** – sanfte Glocke, geht an den Rändern schön auf 0. Nebenkeulen fallen
  **schnell** ab. Der solide **Allrounder**.

- **Hamming** – ähnlich wie Hann, aber die **erste** Nebenkeule ist noch niedriger
  (−43 dB). Dafür sinkt der „Boden" weiter draußen nicht so schnell.

- **Blackman** – **breitere** Hauptkeule (schlechtere Auflösung), dafür **sehr
  niedrige** Nebenkeulen (−58 dB → wenig Leakage).

- **Kaiser** – das **einstellbare** Fenster: über einen Parameter **β** („beta")
  dreht man selbst am Kompromiss. Mehr β = niedrigere Nebenkeulen, aber breitere
  Hauptkeule. Wir nutzen β = 8 (ähnlich gut wie Blackman).

- **Flat-Top** – **breiteste** Hauptkeule (schlechteste Auflösung) und extrem
  niedrige Nebenkeulen. Seine Spezialität: es misst die **Höhe (Amplitude)** eines
  Peaks besonders **genau**. Deshalb in der Messtechnik/Kalibrierung beliebt – für
  die HRV aber meist zu „breit".

**Eure gemessenen Werte** (aus der Konsolentabelle des Programms) bestätigen die
Theorie perfekt:

| Fenster   | PSL (Leakage) | MLW (Auflösung) | Charakter |
|-----------|---------------|-----------------|-----------|
| Rechteck  | −13 dB 😟      | 2 Bins 😃        | scharf, aber leckt stark |
| Hann      | −32 dB        | 4 Bins          | guter Allrounder |
| Hamming   | −43 dB        | 4 Bins          | erste Nebenkeule sehr niedrig |
| Blackman  | −58 dB        | 6 Bins          | sehr leakage-arm |
| Kaiser β8 | −59 dB        | ~5,5 Bins       | einstellbar, sehr gut |
| Flat-Top  | −93 dB 😃      | 10 Bins 😟       | top bei Amplitude, breit |

(😃 = gut, 😟 = schlecht — man sieht: Wer bei Leakage gewinnt, verliert bei der
Auflösung. Genau der Kompromiss aus Abschnitt 6.)

---

## 8. Das „HRV-Leistungsdichtespektrum" lesen (das Einzel-Diagramm)

Bevor wir Fenster **vergleichen**, schauen wir uns *ein* fertiges Spektrum in Ruhe
an. Das ist das Diagramm mit dem Titel **„HRV-Leistungsdichtespektrum"** (im Code:
`plot_spectrum`, in `all_figures.pdf` z. B. die Einzelseiten pro Fenster).

**Was heißt „Leistungsdichte" (PSD)?**
Die FFT sagt uns, *welche* Frequenzen enthalten sind – die **PSD** (englisch
*Power Spectral Density*) sagt zusätzlich, **wie viel „Leistung" (Stärke) pro
Frequenz** steckt. Das Wort **Dichte** ist wörtlich gemeint, wie bei „Einwohner pro
km²": nicht die Leistung an *einem* exakten Punkt, sondern **pro Frequenz-Abschnitt**
(Einheit ms²/Hz). Der Vorteil:

> Die **Fläche unter der Kurve** in einem Frequenzbereich = die **gesamte Leistung**
> in diesem Bereich. Genau so berechnet der Code die Bandwerte VLF/LF/HF
> (in `calculate_hrv_bands` per Flächenberechnung, `trapz`).

**Wie liest man das Bild?**

```
  Leistungsdichte [ms²/Hz]
     ^
     |####|############|##########|
     |####|############|##########|   <- die schwarze PSD-Kurve verläuft oben drüber
     |VLF |     LF     |    HF    |
     +----+------------+----------+------> Frequenz [Hz]
    0  0.04         0.15        0.40
     (farbig hinterlegt: die drei Bänder)
```

- **x-Achse** = Frequenz, **y-Achse** = Leistungsdichte.
- Die **schwarze Kurve** ist das Spektrum selbst.
- Die **farbigen Flächen** markieren die drei Bänder VLF/LF/HF (im Code:
  `shade_hrv_bands`) – so sieht man sofort, in welchem Band die Kurve „Berge" hat.
- Ein **Berg im HF-Band** (~0,25 Hz) kommt typischerweise von der **Atmung**, ein
  Berg im **LF-Band** (~0,1 Hz) von der **Blutdruckregelung**.
- Die kleine **Textbox** fasst die Kennzahlen zusammen: **LF/HF-Verhältnis** sowie
  die LF- und HF-Leistung (absolut in ms² und als Prozent-Anteil).

**Wofür ist dieses Diagramm gut?** Es ist die **Grundansicht** der ganzen Analyse:
Für *eine* Aufnahme (mit *einem* Fenster) zeigt es auf einen Blick, wo die HRV ihre
Energie hat – eher im schnellen HF-Bereich (Zeichen von Erholung/Ruhe) oder im
langsameren LF-Bereich. Alle folgenden Diagramme bauen darauf auf: Abschnitt 9
**vergleicht** solche Spektren für alle Fenster, Abschnitt 10 zeigt sie **über die
Zeit** (Wasserfall).

---

## 9. Was bedeutet das konkret für die HRV?

Jetzt wird es greifbar. Erinnerung: Das HRV-Signal hat **starke langsame** und nur
**schwache schnelle** Anteile. Schauen wir, wie stark die Fensterwahl das
Endergebnis verändert (ebenfalls eure echten Zahlen):

| Fenster  | HF-Leistung [ms²] | LF/HF-Verhältnis |
|----------|-------------------|------------------|
| Rechteck | **30,9**          | 3,37             |
| Hann     | 16,5              | 4,93             |
| Flat-Top | **7,7**           | 5,45             |

Lies das langsam: Dasselbe EKG, dieselbe Person, dieselbe Rechnung – nur ein
**anderes Fenster** – und die HF-Leistung ist beim Rechteck **viermal so groß** wie
beim Flat-Top!

**Warum?** Genau wegen des Leakage aus Abschnitt 4: Das Rechteck lässt die starke
Energie aus dem langsamen Bereich in das schwache HF-Band **hinüberlecken**. Das
HF-Band wird dadurch **künstlich aufgeblasen** (30,9 statt der ehrlicheren ~7,7).
Und weil HF im Nenner steht, wird auch das **LF/HF-Verhältnis** verzerrt
(3,37 statt 5,45).

> Das ist die Kernbotschaft des Projekts: **Die Wahl der Fensterfunktion verändert
> direkt medizinisch relevante Kennwerte.** Ein leaky Fenster kann eine
> Fehldiagnose-Zahl produzieren.

Am schönsten sieht man das auf **Seite 6** von `all_figures.pdf` (überlagerte
Spektren in dB): Die **blaue** Rechteck-Kurve liegt über den ganzen Bereich
**am höchsten** – das ist der Leakage-„Bodennebel". Die **hellblaue Flat-Top**-Kurve
liegt **am tiefsten** – sie ist am saubersten.

---

## 10. Das 3D-Wasserfalldiagramm

Bis jetzt haben wir *ein* Spektrum über die ganze Messung betrachtet. Aber die
Rhythmen **ändern sich über die Zeit** (mal atmet man schneller, mal ruht man).

Deshalb schneiden wir die Messung in **überlappende 5-Minuten-Stücke**, berechnen
für **jedes** Stück ein eigenes Spektrum und stellen diese Spektren
**hintereinander** dar – wie eine Gebirgslandschaft:

```
   Leistung (Höhe)
       ^          /\
       |     /\  /  \       <- spätere Zeit
       |    /  \/    \
       |   /\        /\     <- frühere Zeit
       |  /  \      /  \
       +-------------------> Frequenz
      (nach hinten: die Zeit)
```

- **x-Achse** = Frequenz, **y-Achse** = Zeit (welches 5-min-Stück),
  **z-Achse (Höhe)** = Stärke.
- So sieht man auf einen Blick, ob ein Rhythmus die ganze Zeit da ist oder kommt
  und geht.

**Der Fenster-Effekt hier:** Mit dem **Rechteck** wird die Landschaft **zackig und
unruhig** (Leakage-„Rauschen", siehe Seite 8). Mit **Hann/Kaiser/Flat-Top** wird
sie **glatter und ruhiger** und die Bänder sind sauber getrennt (Seiten 16, 18) –
viel besser ablesbar. Damit man die Fenster fair vergleichen kann, benutzen alle
Wasserfälle **dieselbe Höhen- und Farbskala**.

---

## 11. Fazit: Welches Fenster nehmen?

Es gibt kein „perfektes" Fenster (Naturgesetz, Abschnitt 6). Aber es gibt einen
**guten Kompromiss** für die HRV:

- **Hann-Fenster = die Standard-Empfehlung.** Ordentliche Auflösung, gute
  Leakage-Unterdrückung, stabile und ehrliche Bandwerte, glatte Wasserfälle.
- **Kaiser (β≈8)**, wenn man den Kompromiss selbst feinjustieren will.
- **Rechteck** nur, um zu **zeigen, wie schlecht** Leakage ohne Fenster ist –
  nicht für echte Ergebnisse.
- **Flat-Top** nur, wenn es rein auf die **exakte Höhe** eines einzelnen Peaks
  ankommt (bei HRV selten).

---

## 12. Mini-Glossar (zum schnellen Nachschlagen)

| Begriff | In einem Satz |
|---------|----------------|
| **Frequenz / Hz** | Wie oft sich etwas pro Sekunde wiederholt. |
| **Sinuswelle** | Die einfachste, gleichmäßige Schwingung – der Grundbaustein. |
| **FFT** | Rechenverfahren, das ein Signal in seine Frequenzen zerlegt. |
| **Spektrum** | Das Ergebnis der FFT: welche Frequenz wie stark enthalten ist. |
| **PSD / Leistungsdichte** | Spektrum, das zeigt, wie viel *Leistung* pro Frequenz steckt. |
| **Spectral Leakage** | Ungewolltes „Auslaufen" von Energie in Nachbarfrequenzen. |
| **Fensterfunktion** | Glockenkurve, die die Ränder sanft auf 0 bringt → weniger Leakage. |
| **Hauptkeule** | Der zentrale Peak einer Frequenz; Breite = Frequenzauflösung. |
| **Nebenkeule** | Kleine Nebenbuckel daneben; Höhe = wie stark man leckt. |
| **Bin** | Ein Schritt im Frequenzraster der FFT. |
| **dB (Dezibel)** | Verhältnis-Schreibweise; −20 dB = „zehnmal kleiner". |
| **HRV** | Herzfrequenzvariabilität – die Schwankung der Herzschlag-Abstände. |
| **VLF / LF / HF** | Die drei HRV-Frequenzbänder (sehr langsam / langsam / schnell). |
| **LF/HF** | Verhältnis der beiden Bänder; grober Indikator fürs Nervensystem. |
| **Wasserfall** | 3D-Bild vieler Spektren über die Zeit hintereinander. |

---

## 13. Und im Code?

Wer das Ganze im Programm nachlesen will:

- Die technische Umsetzung Schritt für Schritt steht in [`DOKUMENTATION.md`](DOKUMENTATION.md).
- Der Code selbst ist [`hrv_fft_analyse.m`](hrv_fft_analyse.m) – jede Fensterfunktion
  wird in `apply_window_function` erzeugt, die Kennzahlen aus Abschnitt 6 in
  `window_metrics`, die HRV-Bänder aus Abschnitt 3/9 in `calculate_hrv_bands`.
- Alle hier erwähnten Diagramme liegen gesammelt in `all_figures.pdf`.
