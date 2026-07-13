# Grundlagen: HRV-Analyse und der Einfluss von Fensterfunktionen

Dieses Dokument erklärt die theoretischen Grundlagen des Projekts – von der
Physiologie des Herzschlags über die Frequenzanalyse bis hin zu den
Fensterfunktionen. Es geht bewusst **nicht** auf den konkreten Code ein,
sondern auf das Thema an sich, so dass man versteht, *warum* die einzelnen
Schritte gemacht werden.

---

## 1. Worum geht es überhaupt? (Die Kernidee)

Unser Herz schlägt nicht wie ein Metronom in exakt gleichmäßigen Abständen.
Selbst bei völliger Ruhe schwankt der zeitliche Abstand zwischen zwei
Herzschlägen ständig – mal ein paar Millisekunden länger, mal kürzer. Diese
natürliche Schwankung nennt man **Herzratenvariabilität (HRV)**, im Englischen
*Heart Rate Variability*.

Das Entscheidende: Diese Schwankung ist **kein Zufall und kein Fehler**, sondern
ein Fenster in die Steuerung des Körpers. *Wie stark* und *in welchem Rhythmus*
der Herzschlag schwankt, verrät, wie das vegetative (autonome) Nervensystem
gerade arbeitet – also grob gesagt, wie das Verhältnis von **Anspannung** und
**Entspannung** aussieht.

- **Hohe HRV** → meist ein Zeichen für gute Erholung, Entspannung, Fitness.
- **Niedrige HRV** → oft Zeichen für Stress, Belastung oder Erschöpfung.

Ziel der Analyse ist es, diese Schwankungen messbar zu machen und in ihre
verschiedenen **Rhythmen (Frequenzen)** zu zerlegen.

---

## 2. Physiologische Grundlagen

### 2.1 Wie entsteht der Herzschlag?

Das Herz besitzt einen eigenen „Taktgeber", den **Sinusknoten**. Er erzeugt von
selbst elektrische Impulse, die den Herzmuskel zum Zusammenziehen bringen. Der
Sinusknoten allein würde relativ gleichmäßig takten – die eigentliche
Feinregelung kommt von außen, über das **vegetative Nervensystem**.

### 2.2 Das vegetative Nervensystem – Gas und Bremse

Das vegetative (autonome) Nervensystem arbeitet unbewusst und hat zwei
Gegenspieler, die permanent auf den Sinusknoten einwirken:

| System | Wirkung auf das Herz | Bildlich | Reaktionsgeschwindigkeit |
|---|---|---|---|
| **Sympathikus** | beschleunigt den Herzschlag | „Gaspedal" | langsam (Sekunden) |
| **Parasympathikus** (Vagusnerv) | verlangsamt den Herzschlag | „Bremse" | sehr schnell (unter 1 Sekunde) |

Der **Sympathikus** wird bei Belastung, Stress oder Aufregung aktiv
(„Kampf oder Flucht"). Der **Parasympathikus** dominiert bei Ruhe, Erholung und
Verdauung („Rest and Digest").

Weil beide gleichzeitig und ständig gegeneinander regeln, entsteht das
charakteristische „Zittern" im Herzrhythmus – die HRV. Genau dieses
Zusammenspiel wollen wir sichtbar machen.

### 2.3 Die Atmung als sichtbarster Einfluss

Der auffälligste Effekt ist die **respiratorische Sinusarrhythmie (RSA)**:

- Beim **Einatmen** wird der Herzschlag **schneller**.
- Beim **Ausatmen** wird er **langsamer**.

Diese Schwankung im Atemtakt wird fast vollständig vom Parasympathikus
(Vagusnerv) gesteuert. Sie ist der Grund, warum sich die Atmung später als
deutlicher Ausschlag im hochfrequenten Bereich (HF-Band) des Spektrums
wiederfindet.

---

## 3. Vom EKG zum HRV-Signal

Bevor man die HRV analysieren kann, braucht man erst einmal die genauen
Zeitpunkte der Herzschläge. Der Weg dahin sieht so aus:

```
EKG-Signal  →  R-Zacken finden  →  Abstände messen  →  Tachogramm
```

### 3.1 Das EKG und die R-Zacke

Ein **EKG (Elektrokardiogramm)** zeichnet die elektrische Aktivität des Herzens
auf. Jeder Herzschlag hinterlässt ein typisches Muster (den QRS-Komplex). Der
mit Abstand höchste, schmalste Ausschlag darin ist die **R-Zacke**. Sie ist gut
erkennbar und markiert sehr präzise den Zeitpunkt eines Herzschlags.

### 3.2 RR-Intervalle und das Tachogramm

Der zeitliche Abstand zwischen zwei aufeinanderfolgenden R-Zacken heißt
**RR-Intervall** (auch NN-Intervall für „normal-to-normal") und wird in
**Millisekunden** angegeben.

- Ein RR-Intervall von 1000 ms entspricht einem Puls von 60 Schlägen/Minute.
- Ein RR-Intervall von 800 ms entspricht 75 Schlägen/Minute.

Trägt man alle RR-Intervalle nacheinander auf, erhält man das **Tachogramm** –
die Kurve, die die Schwankung des Herzrhythmus über die Zeit zeigt. **Dieses
Tachogramm ist das eigentliche Signal, das wir analysieren.**

### 3.3 Artefakte

Nicht jeder erkannte Schlag ist echt: Bewegungen, Muskelzucken oder eine
verpasste R-Zacke erzeugen unrealistische Werte (z. B. ein RR-Intervall von
200 ms = 300 Schläge/Minute, physiologisch unmöglich). Solche **Artefakte**
müssen erkannt und korrigiert werden, sonst verfälschen sie das ganze Spektrum.

---

## 4. Zeitbereich vs. Frequenzbereich

Es gibt zwei grundsätzliche Sichtweisen auf die HRV.

### 4.1 Zeitbereich (einfache Statistik)

Man kann direkt aus den RR-Intervallen statistische Kennwerte berechnen, z. B.:

- **SDNN** – die Standardabweichung aller RR-Intervalle (Maß für die
  Gesamtvariabilität).
- **RMSSD** – Maß für die schnellen Schwankungen von Schlag zu Schlag (stark
  vom Parasympathikus geprägt).

Diese Werte sagen aber nur, *wie stark* der Rhythmus schwankt – **nicht**, in
welchen *Rhythmen* er das tut.

### 4.2 Frequenzbereich (die spannende Sichtweise)

Interessanter ist die Frage: Aus welchen **Rhythmen** setzt sich die Schwankung
zusammen? Gibt es eine langsame Welle alle 10 Sekunden? Eine schnelle Welle im
Atemtakt? Genau das beantwortet die **Frequenzanalyse**.

> **Analogie:** Ein Musikakkord klingt für das Ohr wie ein einziger Ton. In
> Wahrheit besteht er aus mehreren einzelnen Tönen (Frequenzen). Die
> Frequenzanalyse ist wie ein „musikalisches Ohr", das den Akkord in seine
> Einzeltöne zerlegt.

---

## 5. Die Fourier-Transformation (FFT)

### 5.1 Die Grundidee

Die **Fourier-Transformation** beruht auf einer mächtigen Erkenntnis: *Jedes*
Signal lässt sich als Summe vieler einfacher Schwingungen (Sinus- und
Kosinuswellen) unterschiedlicher Frequenz darstellen. Die Transformation dreht
den Blickwinkel:

- **Vorher (Zeitbereich):** „Wie ändert sich das Signal über die Zeit?"
- **Nachher (Frequenzbereich):** „Aus welchen Frequenzen besteht das Signal, und
  wie stark ist jede davon vertreten?"

Die **FFT (Fast Fourier Transform)** ist einfach ein schnelles Rechenverfahren,
um das am Computer zu berechnen.

### 5.2 Das Ergebnis: das Leistungsdichtespektrum (PSD)

Das Ergebnis ist ein **Spektrum**: eine Kurve, die für jede Frequenz zeigt, wie
viel „Energie" (Leistung) in diesem Rhythmus steckt. Genauer spricht man vom
**Leistungsdichtespektrum** (engl. *Power Spectral Density, PSD*). Ein hoher
Ausschlag bei einer bestimmten Frequenz bedeutet: Dieser Rhythmus ist im
Herzsignal stark vertreten.

### 5.3 Eine wichtige Voraussetzung: gleichmäßige Abtastung

Die FFT setzt voraus, dass die Messwerte in **exakt gleichen Zeitabständen**
vorliegen. Die RR-Intervalle kommen aber genau *nicht* gleichmäßig – sie
entstehen ja bei jedem (unregelmäßigen) Herzschlag. Deshalb muss das
Tachogramm vorher auf ein gleichmäßiges Zeitraster umgerechnet werden
(**Interpolation / Resampling**, z. B. auf 4 Werte pro Sekunde).

> **Nyquist-Regel:** Man muss mindestens doppelt so schnell abtasten wie die
> höchste interessierende Frequenz. Bei 4 Hz Abtastung sind Frequenzen bis 2 Hz
> darstellbar – das reicht für die HRV (max. 0,4 Hz) locker aus.

---

## 6. Die HRV-Frequenzbänder

Das Spektrum wird in drei physiologisch definierte **Frequenzbänder**
aufgeteilt (nach dem internationalen Standard der Task Force von 1996). Die
Leistung in jedem Band entspricht der **Fläche unter der Spektrumkurve** in
diesem Bereich.

| Band | Frequenzbereich | Entsprechender Rhythmus | Physiologische Bedeutung |
|---|---|---|---|
| **VLF** (very low frequency) | 0,0033 – 0,04 Hz | > 25 Sekunden pro Welle | Langsame Prozesse: Temperaturregelung, Hormone, Gefäßtonus. Nicht vollständig geklärt. |
| **LF** (low frequency) | 0,04 – 0,15 Hz | ca. 7–25 Sekunden pro Welle | Blutdruckregelung (Baroreflex, „Mayer-Wellen" bei ~0,1 Hz). Mischung aus Sympathikus **und** Parasympathikus. |
| **HF** (high frequency) | 0,15 – 0,40 Hz | ca. 2,5–7 Sekunden pro Welle | Atmung (RSA). Fast reiner Parasympathikus (Vagus). |

### 6.1 Wie man von Frequenz auf Zeit umrechnet

Frequenz und Periodendauer sind Kehrwerte voneinander: `Zeit = 1 / Frequenz`.

- 0,25 Hz → eine Welle alle 4 Sekunden → passt zu einer ruhigen Atmung von
  ca. 15 Atemzügen pro Minute (deshalb HF = Atmung).
- 0,1 Hz → eine Welle alle 10 Sekunden → typische Blutdruck-„Mayer-Welle"
  (deshalb im LF-Band).

### 6.2 Das LF/HF-Verhältnis

Aus den beiden Bändern LF und HF bildet man oft das **LF/HF-Verhältnis**. Es
wird klassisch als Maß für die **„sympatho-vagale Balance"** interpretiert:

- **hohes LF/HF** → eher Richtung Anspannung/Sympathikus.
- **niedriges LF/HF** → eher Richtung Entspannung/Parasympathikus.

> **Hinweis:** Diese Interpretation ist verbreitet, aber in der Fachwelt auch
> umstritten, weil das LF-Band nicht rein sympathisch ist. Für dieses Projekt
> reicht das LF/HF-Verhältnis als anschaulicher, gängiger Kennwert.

---

## 7. Das Fenster-Problem (der Kern des Projekts)

Jetzt kommt der eigentliche Untersuchungsgegenstand des Projekts: die
**Fensterfunktionen**.

### 7.1 Warum entsteht überhaupt ein Problem?

Die FFT geht insgeheim davon aus, dass sich das analysierte Signalstück
**endlos periodisch wiederholt** – als würde man das Stück aneinanderreihen. In
der Realität schneiden wir aber einfach einen endlichen Ausschnitt heraus.
Dadurch passen Anfang und Ende meist nicht zusammen, und es entstehen an den
**Rändern harte Sprünge (Kanten)**.

Diese künstlichen Kanten sind im Originalsignal gar nicht vorhanden. Die FFT
„sieht" sie aber und deutet sie als zusätzliche Frequenzen. Das Ergebnis:
Energie „verschmiert" in Frequenzbereiche, in die sie nicht gehört. Dieser
Effekt heißt **spektrale Leckage (Leakage)** und verfälscht das Spektrum –
gerade schwache Signale können von der Leckage starker Nachbarn überdeckt
werden.

### 7.2 Die Lösung: Fensterfunktionen

Ein **Fenster** ist eine Gewichtungskurve, mit der man das Signalstück
multipliziert. Sie ist in der Mitte ≈ 1 und läuft zu den Rändern hin sanft auf
0 zu. Dadurch werden die störenden Randsprünge **weich ausgeblendet**, und die
Leckage sinkt deutlich.

> **Analogie:** Ein Musikstück, das nicht abrupt abbricht, sondern am Ende sanft
> ausgeblendet wird (Fade-out) – ohne störendes „Knacken".

### 7.3 Der grundlegende Zielkonflikt

Es gibt **kein perfektes Fenster**. Jedes ist ein Kompromiss zwischen zwei
gegensätzlichen Eigenschaften:

- **Frequenzauflösung** (Schärfe): Wie gut lassen sich zwei dicht benachbarte
  Frequenzen noch trennen? → beschrieben durch die **Hauptkeulenbreite**.
- **Leckage-Unterdrückung**: Wie gut werden falsche Nachbarfrequenzen
  unterdrückt? → beschrieben durch die **Nebenkeulendämpfung**.

Die Regel lautet:

> **Wer die Ränder stark abdämpft, unterdrückt die Leckage gut – verliert aber
> an Schärfe. Wer wenig abdämpft, bleibt scharf – lässt aber mehr Leckage zu.**

### 7.4 Die untersuchten Fenster

| Fenster | Charakter | Stärke | Schwäche |
|---|---|---|---|
| **Rechteck (rect)** | gar keine Abdämpfung | schärfste Auflösung | stärkste Leckage |
| **Hann** | glockenförmig | guter Allrounder | mittel |
| **Hamming** | ähnlich Hann | gute Nebenkeulendämpfung | leichte Rest-Leckage |
| **Blackman** | breiter, stark gedämpft | sehr geringe Leckage | breitere Hauptkeule |
| **Kaiser** | einstellbar (Parameter β) | flexibel anpassbar | je nach Einstellung |
| **Flat-Top** | sehr breite, flache Spitze | beste Amplitudentreue | schlechteste Auflösung |

### 7.5 Kennzahlen zur Beschreibung eines Fensters

Um Fenster objektiv zu vergleichen, nutzt man drei Maße:

- **MLW (Main Lobe Width, Hauptkeulenbreite):** Wie breit ist der zentrale
  Ausschlag? → **klein = scharfe Auflösung.**
- **PSL (Peak Side Lobe, höchste Nebenkeule):** Wie hoch sind die störenden
  Nebenausschläge (in dB unter dem Hauptausschlag)? → **stark negativ = wenig
  Leckage.**
- **ENBW (Equivalent Noise Bandwidth):** Effektive Bandbreite des Fensters; ein
  Maß dafür, wie das Fenster Rauschen „einsammelt".

---

## 8. Das Welch-Verfahren (stabilere Spektren)

Ein einzelnes FFT-Spektrum ist oft stark **verrauscht** und zappelig. Um ein
ruhigeres, verlässlicheres Ergebnis zu bekommen, verwendet man das
**Welch-Verfahren**:

1. Das Gesamtsignal wird in mehrere **überlappende Abschnitte** zerlegt
   (z. B. jeweils 5 Minuten lang, mit 50 % Überlappung).
2. Jeder Abschnitt wird **gefenstert** und einzeln per FFT in sein Spektrum
   umgerechnet.
3. Alle Einzelspektren werden **gemittelt**.

> **Analogie:** Ein Foto bei wenig Licht ist verrauscht. Macht man zehn Fotos
> und mittelt sie, wird das Bild ruhig und klar. Genau das bewirkt das Mitteln
> der Spektren – das zufällige Rauschen mittelt sich weg, die echten Anteile
> bleiben.

**Kompromiss:** Kürzere Abschnitte → mehr Mittelungen → ruhigeres Spektrum,
aber **gröbere** Frequenzauflösung. Es gilt also wieder abzuwägen.

---

## 9. Die eigentliche Forschungsfrage des Projekts

Alle bisherigen Bausteine führen auf eine zentrale Frage zu:

> **Wie stark hängen die berechneten HRV-Werte (VLF, LF, HF, LF/HF) davon ab,
> welche Fensterfunktion man wählt?**

Das ist praktisch relevant, weil in Studien und Geräten unterschiedliche Fenster
verwendet werden. Wenn allein die Wahl des Fensters die Ergebnisse spürbar
verändert, sind Messwerte **nur dann vergleichbar**, wenn dieselbe Methode
verwendet wurde. Das Projekt macht diesen Einfluss sichtbar, indem es dieselben
Herzdaten mit sechs verschiedenen Fenstern auswertet und die Ergebnisse
gegenüberstellt.

---

## 10. Zusammenfassung in einem Absatz

Das Herz schlägt durch das Zusammenspiel von Sympathikus (Gas) und
Parasympathikus (Bremse) leicht unregelmäßig – diese **Herzratenvariabilität**
spiegelt den Zustand des vegetativen Nervensystems wider. Aus dem EKG bestimmt
man die genauen Herzschlag-Zeitpunkte, daraus die **RR-Intervalle** (das
Tachogramm) und zerlegt dieses Signal mit der **FFT** in seine Frequenzanteile.
Das Spektrum wird in die Bänder **VLF, LF und HF** aufgeteilt, deren Verhältnis
etwas über Anspannung und Entspannung aussagt. Weil man das Signal für die FFT
in endliche Stücke schneiden muss, entsteht **spektrale Leckage**, die man mit
**Fensterfunktionen** dämpft – jede davon ein Kompromiss zwischen Schärfe und
Leckage. Das Projekt untersucht, **wie stark die Wahl der Fensterfunktion die
HRV-Ergebnisse beeinflusst.**

---

## Glossar (Kurzüberblick)

| Begriff | Bedeutung |
|---|---|
| **HRV** | Herzratenvariabilität – die Schwankung der Herzschlag-Abstände |
| **EKG** | Elektrokardiogramm – Aufzeichnung der elektrischen Herzaktivität |
| **R-Zacke** | höchster Ausschlag pro Herzschlag im EKG, markiert den Schlagzeitpunkt |
| **RR-Intervall** | zeitlicher Abstand zwischen zwei R-Zacken (in ms) |
| **Tachogramm** | Kurve aller RR-Intervalle über die Zeit |
| **FFT** | Fast Fourier Transform – zerlegt ein Signal in seine Frequenzen |
| **PSD** | Leistungsdichtespektrum – Energie je Frequenz |
| **Leakage** | spektrale Leckage – Verschmieren von Energie in falsche Frequenzen |
| **Fensterfunktion** | Gewichtungskurve, die die Signalränder weich ausblendet |
| **Welch-Verfahren** | Mittelung mehrerer Segment-Spektren für ein ruhigeres Ergebnis |
| **VLF / LF / HF** | die drei HRV-Frequenzbänder |
| **Sympathikus / Parasympathikus** | die zwei Gegenspieler des vegetativen Nervensystems |
