%% ========================================================================
%  HRV-ANALYSE: EINFLUSS UNTERSCHIEDLICHER FFT-FENSTERFUNKTIONEN
%  ========================================================================
%
%  Diese Datei enthaelt die VOLLSTAENDIGE, in sich geschlossene
%  Implementierung des Projekts "Einfluss unterschiedlicher
%  FFT-Fensterfunktionen auf die HRV-Analyse".
%
%  Alle Teilfunktionen sind als LOKALE FUNKTIONEN am Ende dieser Datei
%  enthalten (moeglich ab MATLAB R2016b). Es ist also nur EINE Datei
%  noetig -- einfach in MATLAB oeffnen und ausfuehren (F5 / "Run").
%
%  ------------------------------------------------------------------------
%  WAS MACHT DAS PROGRAMM?
%  ------------------------------------------------------------------------
%  Die Herzfrequenzvariabilitaet (HRV) beschreibt die natuerliche
%  Schwankung der zeitlichen Abstaende aufeinanderfolgender Herzschlaege
%  (der sogenannten RR-Intervalle). Diese Schwankung enthaelt
%  diagnostisch wertvolle Information ueber das vegetative Nervensystem.
%  Man analysiert sie ueblicherweise im FREQUENZBEREICH und unterscheidet
%  drei Frequenzbaender (VLF, LF, HF).
%
%  Um vom Frequenzbereich zu sprechen, braucht man eine FFT. Da ein
%  endliches Signalstueck analysiert wird, entstehen bei der FFT
%  spektrale Leckeffekte ("Leakage"). FENSTERFUNKTIONEN daempfen diese
%  Effekte -- allerdings auf Kosten der Frequenzaufloesung. Dieses
%  Programm vergleicht sechs gaengige Fenster (Rechteck, Hann, Hamming,
%  Blackman, Kaiser, Flat-Top) und zeigt ihren Einfluss auf die
%  HRV-Spektralanalyse.
%
%  ------------------------------------------------------------------------
%  WICHTIG - ZWEI VERSCHIEDENE BEDEUTUNGEN VON "FENSTER"
%  ------------------------------------------------------------------------
%    1) ANALYSEFENSTER (Segment): ein ZEITABSCHNITT des Signals, hier
%       5 Minuten lang. Die lange Aufzeichnung wird in solche Abschnitte
%       zerlegt (mit 50 % Ueberlappung), die einzeln ausgewertet und dann
%       gemittelt werden (Welch-Methode). -> Parameter SEG_LEN_S, OVERLAP.
%    2) FENSTERFUNKTION (Hann, Hamming, ...): eine Gewichtungskurve, die
%       INNERHALB jedes Abschnitts auf das Signal multipliziert wird, um
%       das spektrale Leakage zu reduzieren. -> Parameter WINDOWS.
%    Kurz: Erst in 5-min-Abschnitte schneiden (1), dann auf jeden Abschnitt
%    eine Fensterfunktion legen (2) und die FFT rechnen.
%
%  ------------------------------------------------------------------------
%  VERARBEITUNGSKETTE (jeweils eine lokale Funktion)
%  ------------------------------------------------------------------------
%    load_ecg_data          EKG laden (EDF-Datei ODER synthetische Daten)
%        v
%    preprocess_ecg         Stoerungen filtern (DC, Baseline, 50 Hz, Rauschen)
%        v
%    detect_r_peaks         R-Zacken finden (vereinfachter Pan-Tompkins)
%        v
%    calculate_rr_intervals RR-Intervalle bilden + Artefakte korrigieren
%        v
%    interpolate_rr_signal  ungleichmaessige RR-Reihe -> gleichmaessiges Raster
%        v
%    compare_window_functions  je Fenster ein Welch-Spektrum berechnen
%      ( nutzt: apply_window_function, calculate_fft, calculate_hrv_bands )
%        v
%    plot_spectrum / plot_waterfall   Ergebnisse je Fenster darstellen
%
%  ------------------------------------------------------------------------
%  BEDIENUNG
%  ------------------------------------------------------------------------
%    * Ohne eigene Daten: einfach ausfuehren. Es wird automatisch ein
%      synthetisches 1-stuendiges Test-EKG mit BEKANNTER HRV erzeugt
%      (LF-Anteil bei 0.10 Hz, HF-Anteil bei 0.25 Hz). An den Spektren
%      laesst sich so die Korrektheit der Kette pruefen.
%    * Mit eigener EDF-Datei: unten die Variable EDF_FILE auf den Pfad
%      setzen (EDF-Import benoetigt MATLAB R2020b oder neuer).
%
%  ------------------------------------------------------------------------
%  BENOETIGTE TOOLBOX
%  ------------------------------------------------------------------------
%    Signal Processing Toolbox (butter, filtfilt, findpeaks, hann,
%    hamming, blackman, kaiser, flattopwin). Der 50-Hz-Notch ist bewusst
%    OHNE iirnotch implementiert (manueller Biquad), um die Abhaengigkeiten
%    gering zu halten.
%  ========================================================================

clear;        % Arbeitsbereich leeren
clc;          % Befehlsfenster leeren
close all;    % alle Abbildungen schliessen

%% ----------------------- EINSTELLBARE PARAMETER -------------------------
% Diese Werte steuern die gesamte Analyse und koennen frei angepasst werden.

EDF_FILE   = '';      % Pfad zur EDF-Datei. '' => synthetische Testdaten.
FS_INTERP  = 4;       % Resampling-Frequenz des RR-Signals [Hz].
                      %   4 Hz ist Standard in der HRV-Literatur: die
                      %   hoechste interessierende Frequenz (HF-Band) endet
                      %   bei 0.40 Hz, 4 Hz liegt also weit ueber dem
                      %   Nyquist-Kriterium (>0.80 Hz) -> kein Aliasing.
SEG_LEN_S  = 300;     % Laenge eines ANALYSEFENSTERS (Segments) [s] = 5 min.
                      %   ACHTUNG - nicht mit der FensterFUNKTION verwechseln!
                      %   Hier ist ein ZEITABSCHNITT gemeint: die lange
                      %   Aufzeichnung wird in 5-min-Stuecke zerlegt, und JEDES
                      %   Stueck wird einzeln spektral ausgewertet. 5 min ist
                      %   der Standard der HRV-Kurzzeitanalyse (Task Force 1996):
                      %   lang genug, um die tiefen LF-Frequenzen aufzuloesen,
                      %   und kurz genug, dass die Herzrate im Abschnitt
                      %   naeherungsweise stationaer bleibt (Voraussetzung der FFT).
OVERLAP    = 0.5;     % Ueberlappung aufeinanderfolgender Segmente (0.5 = 50 %).
                      %   Beispiel: Segment 1 = 0-5 min, Segment 2 = 2.5-7.5 min,
                      %   Segment 3 = 5-10 min, ...  Durch die Ueberlappung
                      %   entstehen aus denselben Daten mehr Segmente; ihre
                      %   Spektren werden gemittelt (Welch-Methode), was das
                      %   Ergebnis glaettet und die Schwankung verringert.
WINDOWS    = {'rect','hann','hamming','blackman','kaiser','flattop'};
                      % Die zu vergleichenden FENSTERFUNKTIONEN. Diese werden
                      %   INNERHALB jedes Segments angewandt (siehe oben) -- das
                      %   ist die zweite, voellig andere Bedeutung von "Fenster".

%% --------------------------- 1. EKG LADEN -------------------------------
fprintf('=== 1. EKG-Daten laden ===\n');
[ecg_raw, fs] = load_ecg_data(EDF_FILE);
t = (0:numel(ecg_raw)-1).' / fs;     % Zeitvektor [s] (Spaltenvektor)

%% ------------------------ 2. VORVERARBEITUNG ----------------------------
fprintf('\n=== 2. Vorverarbeitung ===\n');
ecg = preprocess_ecg(ecg_raw, fs);

% --- Roh- gegen gefiltertes Signal zeigen (erste 10 s zur Anschauung)
figure('Color','w','Position',[100 100 950 500]);
sel = t <= 10;                       % Maske: nur die ersten 10 Sekunden
subplot(2,1,1);
plot(t(sel), ecg_raw(sel), 'Color', [0.6 0.6 0.6]);
title('EKG - Rohsignal (erste 10 s)'); ylabel('Amplitude');
grid on; box on; xlim([0 10]);
subplot(2,1,2);
plot(t(sel), ecg(sel), 'Color', [0.85 0.2 0.2]);
title('EKG - nach Vorverarbeitung'); xlabel('Zeit [s]'); ylabel('Amplitude');
grid on; box on; xlim([0 10]);

%% --------------------------- 3. R-ZACKEN --------------------------------
fprintf('\n=== 3. R-Zacken-Erkennung ===\n');
r_locs = detect_r_peaks(ecg, fs);

% --- erkannte R-Zacken im Signal markieren (erste 10 s)
figure('Color','w','Position',[100 100 950 350]);
plot(t(sel), ecg(sel), 'k'); hold on;
r_sel = r_locs(r_locs <= 10*fs);     % nur R-Zacken in den ersten 10 s
plot(r_sel/fs, ecg(r_sel), 'ro', 'MarkerFaceColor','r', 'MarkerSize',5);
title('Erkannte R-Zacken (erste 10 s)');
xlabel('Zeit [s]'); ylabel('Amplitude');
legend({'EKG','R-Zacken'}, 'Location','northeast');
grid on; box on; xlim([0 10]); hold off;

%% ------------------------- 4. RR-INTERVALLE -----------------------------
fprintf('\n=== 4. RR-Intervalle + Artefaktkorrektur ===\n');
[t_rr, rr, n_art] = calculate_rr_intervals(r_locs, fs);

% --- Tachogramm: RR-Intervall ueber der Zeit
figure('Color','w','Position',[100 100 950 350]);
plot(t_rr, rr, '-', 'Color', [0.1 0.4 0.8]); hold on;
plot(t_rr, rr, '.', 'Color', [0.1 0.4 0.8], 'MarkerSize', 6);
title(sprintf('Tachogramm (RR-Intervalle) - %d Artefakte korrigiert', n_art));
xlabel('Zeit [s]'); ylabel('RR-Intervall [ms]');
grid on; box on; hold off;

%% ------------------------- 5. INTERPOLATION -----------------------------
fprintf('\n=== 5. Interpolation (Resampling) ===\n');
[~, sig, fs_i] = interpolate_rr_signal(t_rr, rr, FS_INTERP);

%% ----------------------- 6. FENSTERVERGLEICH ----------------------------
fprintf('\n=== 6. Vergleich der Fensterfunktionen ===\n');
results = compare_window_functions(sig, fs_i, WINDOWS, SEG_LEN_S, OVERLAP);

%% --------------------- 7. EINZELDARSTELLUNGEN ---------------------------
fprintf('\n=== 7. Einzelspektren und Wasserfalldiagramme ===\n');
for iw = 1:numel(results)
    R = results(iw);
    plot_spectrum(R.f, R.psd, R.name, R.hrv);            % 2D-Spektrum
    plot_waterfall(R.f, R.seg_psd, R.t_seg, R.name);     % 3D-Wasserfall
end

fprintf('\n=== Analyse abgeschlossen ===\n');


%% ========================================================================
%  AB HIER: LOKALE FUNKTIONEN
%  (in MATLAB ab R2016b in derselben Skriptdatei erlaubt)
%  ========================================================================


% =========================================================================
function [ecg, fs] = load_ecg_data(filename)
%LOAD_ECG_DATA  Liest ein EKG-Signal aus einer EDF-Datei.
%   [ecg, fs] = LOAD_ECG_DATA(filename) liest den ersten Signalkanal der
%   EDF-Datei FILENAME und gibt das Signal ECG (Spaltenvektor) sowie die
%   Abtastrate FS (Hz) zurueck.
%
%   Wird kein Dateiname uebergeben oder existiert die Datei nicht, so wird
%   automatisch ein synthetisches Test-EKG erzeugt (GENERATE_SYNTHETIC_ECG),
%   damit das Programm in jedem Fall lauffaehig ist.
%
%   EDF (European Data Format) ist das Standardformat fuer biomedizinische
%   Zeitreihen (EEG, EKG, ...). Der EDF-Import nutzt EDFREAD/EDFINFO und
%   benoetigt MATLAB R2020b oder neuer.

    % Fallback auf synthetische Daten, falls keine gueltige Datei vorliegt
    if nargin < 1 || isempty(filename) || exist(filename, 'file') ~= 2
        if nargin >= 1 && ~isempty(filename)
            warning('Datei "%s" nicht gefunden -> synthetische Testdaten.', filename);
        else
            fprintf('Kein EDF-Dateiname uebergeben -> synthetische Testdaten.\n');
        end
        [ecg, fs] = generate_synthetic_ecg(3600, 250);   % 1 h @ 250 Hz
        return;
    end

    % --- Abtastrate aus dem EDF-Header bestimmen ------------------------
    % Im EDF-Header stehen Dauer eines Datensatzes und Anzahl der Samples
    % pro Datensatz; daraus folgt die Abtastrate.
    info   = edfinfo(filename);
    recDur = seconds(info.DataRecordDuration);     % Dauer eines Datensatzes [s]
    fs     = double(info.NumSamples(1)) / recDur;  % Samples/s des 1. Kanals

    % --- Signaldaten lesen ----------------------------------------------
    % edfread liefert eine Timetable mit einer Zeile pro Datensatz; jede
    % Zelle enthaelt die Samples dieses Datensatzes. Diese werden zu einem
    % durchgehenden Vektor aneinandergehaengt.
    tt    = edfread(filename);
    vname = tt.Properties.VariableNames{1};        % erster Signalkanal
    col   = tt.(vname);
    if iscell(col)
        ecg = cell2mat(col);
    else
        ecg = col(:);
    end
    ecg = double(ecg(:));

    fprintf('EDF geladen: "%s" | Kanal "%s" | fs = %.1f Hz | Dauer = %.1f min\n', ...
            filename, vname, fs, numel(ecg)/fs/60);
end


% =========================================================================
function [ecg, fs, t] = generate_synthetic_ecg(duration, fs)
%GENERATE_SYNTHETIC_ECG  Erzeugt ein synthetisches EKG mit bekannter HRV.
%   [ecg, fs, t] = GENERATE_SYNTHETIC_ECG(duration, fs) liefert ein
%   kuenstliches EKG der Laenge DURATION (s) bei Abtastrate FS (Hz).
%
%   Der Sinn dieser Funktion ist die VALIDIERUNG: In das Signal werden
%   gezielt zwei HRV-Komponenten "eingebaut" -- eine LF-Schwingung bei
%   0.10 Hz und eine HF-Schwingung bei 0.25 Hz. Funktioniert die gesamte
%   Analysekette korrekt, MUESSEN im Endergebnis genau bei diesen
%   Frequenzen Spektralpeaks erscheinen.
%
%   Zusaetzlich werden typische Stoerungen ueberlagert, damit die
%   Vorverarbeitung etwas zu tun hat (DC-Offset, Baseline-Wanderung durch
%   Atmung, langsame Drift, 50-Hz-Netzbrummen, hochfrequentes Rauschen).

    if nargin < 1 || isempty(duration), duration = 3600; end
    if nargin < 2 || isempty(fs),       fs       = 250;  end

    rng(42);   % feste Zufallssaat -> reproduzierbare Ergebnisse

    % ----------------------------------------------------------------
    % 1) RR-Intervallreihe mit eingebetteter LF/HF-Modulation aufbauen
    % ----------------------------------------------------------------
    % Idee: Das Grund-RR-Intervall (mean_rr) wird langsam moduliert. Die
    % Modulation besteht aus einer LF- und einer HF-Sinusschwingung. So
    % entsteht eine Herzschlagfolge mit definierter HRV.
    mean_rr = 0.857;            % mittleres RR-Intervall [s]  (~70 bpm)
    f_lf    = 0.10;  a_lf = 0.040;   % LF: Frequenz [Hz] / Modulationstiefe [s]
    f_hf    = 0.25;  a_hf = 0.025;   % HF: Frequenz [Hz] / Modulationstiefe [s]

    % Schlagzeitpunkte iterativ erzeugen: jeder neue Schlag liegt um das
    % (modulierte) aktuelle RR-Intervall nach dem vorherigen.
    t_beats = zeros(1, ceil(duration/mean_rr) + 100);
    t_beats(1) = 0;
    k = 1;
    while t_beats(k) < duration
        tcur = t_beats(k);
        rr   = mean_rr ...
             + a_lf * sin(2*pi*f_lf*tcur) ...
             + a_hf * sin(2*pi*f_hf*tcur) ...
             + 0.004 * randn();          % winzige zufaellige Schwankung
        k = k + 1;
        t_beats(k) = tcur + rr;
    end
    t_beats = t_beats(1:k-1);            % letzten (ueber duration) verwerfen

    % ----------------------------------------------------------------
    % 2) EKG-Kurve aus einer QRS-T-Vorlage an jedem Schlag aufbauen
    % ----------------------------------------------------------------
    % Ein realer EKG-Schlag besteht aus P-Welle, QRS-Komplex und T-Welle.
    % Hier wird er vereinfacht aus fuenf Gauss-Glocken zusammengesetzt.
    t   = (0:1/fs:duration-1/fs).';      % Zeitvektor (Spalte)
    ecg = zeros(numel(t), 1);

    % Vorlagen-Loben:  [Offset (s), Amplitude (mV), Breite sigma (s)]
    lobes = [ -0.160,  0.10, 0.020;      % P-Welle
              -0.025, -0.12, 0.008;      % Q-Zacke
               0.000,  1.00, 0.009;      % R-Zacke (groesste Auslenkung)
               0.025, -0.18, 0.008;      % S-Zacke
               0.180,  0.28, 0.040 ];    % T-Welle

    for ib = 1:numel(t_beats)
        tb = t_beats(ib);
        % Aus Geschwindigkeitsgruenden nur ein lokales Zeitfenster um den
        % Schlag herum bearbeiten (nicht das ganze Signal).
        i0 = max(1,        floor((tb - 0.35) * fs) + 1);
        i1 = min(numel(t), floor((tb + 0.40) * fs) + 1);
        if i1 < i0, continue; end
        loc = t(i0:i1);
        for il = 1:size(lobes, 1)
            off = lobes(il, 1);  amp = lobes(il, 2);  sig = lobes(il, 3);
            ecg(i0:i1) = ecg(i0:i1) + amp * exp(-((loc - (tb + off)).^2) ./ (2*sig^2));
        end
    end

    % ----------------------------------------------------------------
    % 3) Stoerungen ueberlagern (zum Test der Vorverarbeitung)
    % ----------------------------------------------------------------
    ecg = ecg + 0.30 * sin(2*pi*0.20*t);          % Baseline-Wanderung (Atmung)
    ecg = ecg + 0.15 * sin(2*pi*0.05*t + 1.0);    % langsame Drift
    ecg = ecg + 0.08 * sin(2*pi*50.0*t);          % 50-Hz-Netzbrummen
    ecg = ecg + 0.02 * randn(numel(t), 1);        % hochfrequentes Rauschen
    ecg = ecg + 0.5;                              % Gleichanteil (DC-Offset)

    fprintf(['Synthetisches EKG erzeugt: %.0f min @ %.0f Hz | %d Schlaege | ', ...
             'LF=%.2f Hz, HF=%.2f Hz\n'], duration/60, fs, numel(t_beats), f_lf, f_hf);
end


% =========================================================================
function ecg_filt = preprocess_ecg(ecg, fs, mains_freq)
%PREPROCESS_ECG  Filtert ein EKG-Rohsignal.
%   ecg_filt = PREPROCESS_ECG(ecg, fs) entfernt der Reihe nach:
%     1. den Gleichanteil  (Mittelwertabzug),
%     2. die Baseline-Wanderung  (Butterworth-Hochpass, 0.5 Hz),
%     3. das Netzbrummen          (Notch-Filter, Standard 50 Hz),
%     4. hochfrequentes Rauschen  (Butterworth-Tiefpass, 40 Hz).
%
%   Alle Filter werden NULLPHASIG mit FILTFILT angewendet. filtfilt
%   filtert vorwaerts und rueckwaerts; dadurch entsteht KEINE
%   Phasenverschiebung. Das ist hier entscheidend, denn eine
%   Phasenverschiebung wuerde die R-Zacken zeitlich verschieben und damit
%   die spaeter berechneten RR-Intervalle verfaelschen.

    if nargin < 3 || isempty(mains_freq), mains_freq = 50; end

    ecg = double(ecg(:));
    nyq = fs / 2;             % Nyquist-Frequenz (halbe Abtastrate)

    % --- 1) Gleichanteil (DC) entfernen ---------------------------------
    % Der konstante Offset traegt keine Information und stoert nachfolgende
    % Filter; einfacher Mittelwertabzug genuegt.
    ecg_filt = ecg - mean(ecg);

    % --- 2) Baseline-Wanderung: Hochpass 0.5 Hz -------------------------
    % Langsame Schwankungen (Atmung, Bewegung, Elektroden-Drift) liegen
    % unterhalb ~0.5 Hz und werden mit einem Hochpass entfernt.
    [b, a]   = butter(2, 0.5/nyq, 'high');
    ecg_filt = filtfilt(b, a, ecg_filt);

    % --- 3) Netzbrummen: schmaler Notch bei mains_freq ------------------
    % Das 50-Hz-Brummen des Stromnetzes wird mit einem schmalen Kerbfilter
    % (Guete Q=30) gezielt herausgefiltert, ohne benachbarte Frequenzen
    % stark zu beeinflussen.
    [bn, an] = notch_coeffs(mains_freq, fs, 30);
    ecg_filt = filtfilt(bn, an, ecg_filt);

    % --- 4) Hochfrequentes Rauschen: Tiefpass 40 Hz ---------------------
    % Die diagnostisch relevante EKG-Energie liegt unter ~40 Hz; daher
    % wird hoeherfrequentes Rauschen (z.B. Muskelartefakte) gedaempft.
    [b, a]   = butter(4, 40/nyq, 'low');
    ecg_filt = filtfilt(b, a, ecg_filt);
end


% =========================================================================
function [b, a] = notch_coeffs(f0, fs, Q)
%NOTCH_COEFFS  Koeffizienten eines RBJ-Bandstop-Biquads (Kerbfilter).
%   Entwurf nach der bekannten "Audio EQ Cookbook"-Formel
%   (Robert Bristow-Johnson). Diese manuelle Variante wird bewusst statt
%   IIRNOTCH verwendet, damit der Code auch ohne diese Toolbox-Funktion
%   (und unter GNU Octave) laeuft.
%
%   f0 = Kerbfrequenz [Hz], fs = Abtastrate [Hz], Q = Guetefaktor
%   (groesseres Q => schmalere, gezieltere Kerbe).
    w0    = 2*pi*f0/fs;          % normierte Kreisfrequenz
    alpha = sin(w0) / (2*Q);     % Bandbreitenparameter
    b = [ 1,          -2*cos(w0),  1         ];   % Zaehler
    a = [ 1 + alpha,  -2*cos(w0),  1 - alpha ];   % Nenner
    b = b / a(1);                % auf a(1) normieren
    a = a / a(1);
end


% =========================================================================
function r_locs = detect_r_peaks(ecg, fs)
%DETECT_R_PEAKS  Automatische Detektion der R-Zacken im EKG.
%   r_locs = DETECT_R_PEAKS(ecg, fs) liefert die Abtastindizes der
%   erkannten R-Zacken.
%
%   Das Verfahren ist eine vereinfachte Variante des klassischen
%   PAN-TOMPKINS-Algorithmus. Die Grundidee: den QRS-Komplex durch eine
%   Kette von Transformationen so hervorheben, dass er als deutlicher,
%   gut detektierbarer Energieimpuls erscheint.

    ecg = double(ecg(:));
    nyq = fs / 2;

    % --- 1) Bandpass 5-15 Hz --------------------------------------------
    % Die Hauptenergie des QRS-Komplexes liegt etwa in diesem Band. Der
    % Bandpass unterdrueckt P-/T-Wellen und Rest-Rauschen.
    [b, a] = butter(2, [5/nyq, 15/nyq], 'bandpass');
    f = filtfilt(b, a, ecg);

    % --- 2) Ableitung ---------------------------------------------------
    % Die zeitliche Ableitung betont die steilen Flanken des QRS-Komplexes.
    d = [f(1); diff(f)];

    % --- 3) Quadrierung -------------------------------------------------
    % Quadrieren macht alle Werte positiv und ueberhoeht grosse (steile)
    % Werte zusaetzlich -> der QRS-Komplex sticht noch deutlicher hervor.
    sq = d .^ 2;

    % --- 4) gleitende Mittelung ueber 150 ms ----------------------------
    % Bildet eine glatte Energie-Einhuellende, in der jeder QRS-Komplex zu
    % genau einem breiten Maximum verschmilzt.
    w     = max(1, round(0.150 * fs));
    integ = conv(sq, ones(w, 1) / w, 'same');

    % --- 5) Peaksuche ueber der Schwelle --------------------------------
    % Schwelle relativ zum Mittelwert der Einhuellenden. Der Mindestabstand
    % von 0.25 s (=> max. 240 bpm) verhindert Doppelerkennungen eines QRS.
    thr     = 0.4 * mean(integ);
    minDist = round(0.25 * fs);
    [~, pk] = findpeaks(integ, 'MinPeakHeight', thr, ...
                               'MinPeakDistance', minDist);

    % --- 6) Feinjustierung ----------------------------------------------
    % Die Einhuellende ist gegenueber der echten R-Zacke leicht verschoben.
    % Daher wird im EKG im Fenster +/- 50 ms um jeden gefundenen Punkt das
    % tatsaechliche lokale Maximum (die R-Spitze) gesucht.
    wref  = round(0.05 * fs);
    r_locs = zeros(numel(pk), 1);
    for i = 1:numel(pk)
        lo = max(1,          pk(i) - wref);
        hi = min(numel(ecg), pk(i) + wref);
        [~, rel] = max(ecg(lo:hi));
        r_locs(i) = lo + rel - 1;
    end

    r_locs = unique(r_locs);    % Duplikate entfernen, aufsteigend sortiert
    fprintf('R-Zacken erkannt: %d\n', numel(r_locs));
end


% =========================================================================
function [t_rr, rr, n_artifacts] = calculate_rr_intervals(r_locs, fs)
%CALCULATE_RR_INTERVALS  RR-Intervalle aus R-Zacken mit Artefaktkorrektur.
%   [t_rr, rr, n_artifacts] = CALCULATE_RR_INTERVALS(r_locs, fs)
%
%   Das RR-Intervall ist der zeitliche Abstand zweier aufeinanderfolgender
%   R-Zacken -- die zentrale Groesse der HRV-Analyse.
%
%   Ausgaben:
%     t_rr        : Zeitpunkte der RR-Intervalle (s)
%     rr          : RR-Intervalle in Millisekunden (artefaktkorrigiert)
%     n_artifacts : Anzahl ersetzter unplausibler Intervalle
%
%   Eine einzelne falsch erkannte oder uebersehene R-Zacke wuerde ein
%   stark abweichendes RR-Intervall erzeugen und das Spektrum verfaelschen.
%   Deshalb werden unplausible Werte erkannt und ersetzt.

    r_locs = r_locs(:);
    t_r    = r_locs / fs;          % R-Zacken-Zeitpunkte [s]
    rr     = diff(t_r) * 1000;     % RR-Intervalle [ms]  (diff = Differenzen)
    t_rr   = t_r(2:end);           % Zeitstempel je Intervall (2. R-Zacke)

    % --- Plausibilitaetspruefung ----------------------------------------
    % (a) physiologischer Bereich: 300 ms (200 bpm) .. 2000 ms (30 bpm)
    ok  = (rr > 300) & (rr < 2000);
    % (b) zusaetzlich: nicht mehr als 50 % Abweichung vom Median
    %     (faengt Ausreisser ab, die zwar im Bereich liegen, aber zum
    %      restlichen Verlauf nicht passen).
    med = median(rr(ok));
    ok  = ok & (abs(rr - med) < 0.5 * med);
    bad = ~ok;
    n_artifacts = sum(bad);

    % --- Artefakte ersetzen ---------------------------------------------
    % Unplausible Intervalle werden aus den gueltigen Nachbarwerten per
    % "pchip" (formerhaltende kubische Interpolation) ersetzt, damit die
    % Zeitreihe luekenlos bleibt.
    if n_artifacts > 0 && sum(ok) >= 2
        rr(bad) = interp1(t_rr(ok), rr(ok), t_rr(bad), 'pchip', 'extrap');
    end

    fprintf(['RR-Intervalle: %d | Artefakte korrigiert: %d | ', ...
             'mittl. RR = %.1f ms (%.1f bpm) | SDNN = %.1f ms\n'], ...
             numel(rr), n_artifacts, mean(rr), 60000/mean(rr), std(rr));
end


% =========================================================================
function [t_u, sig, fs_i] = interpolate_rr_signal(t_rr, rr, fs_i)
%INTERPOLATE_RR_SIGNAL  Gleichmaessig abgetastetes HRV-Signal erzeugen.
%   [t_u, sig, fs_i] = INTERPOLATE_RR_SIGNAL(t_rr, rr, fs_i)
%
%   PROBLEM: Die RR-Intervalle liegen zu den (unregelmaessigen)
%   Herzschlag-Zeitpunkten vor. Die FFT setzt aber ein GLEICHMAESSIG
%   abgetastetes Signal voraus. Daher muss die RR-Reihe zunaechst auf ein
%   gleichmaessiges Zeitraster interpoliert ("resampled") werden.
%
%   Es wird eine kubische SPLINE-Interpolation verwendet: sie liefert
%   einen glatten, stetig differenzierbaren Verlauf und erzeugt damit
%   weniger kuenstliche Spektralartefakte als eine lineare Interpolation.
%   Danach wird mit DETREND ein linearer Trend (inkl. Mittelwert) entfernt
%   -- ein Gleich- oder Trendanteil wuerde sonst das niederfrequente Ende
%   des Spektrums dominieren.

    if nargin < 3 || isempty(fs_i), fs_i = 4; end

    t_rr = t_rr(:);
    rr   = rr(:);

    % gleichmaessiges Zeitraster vom ersten bis zum letzten RR-Zeitpunkt
    t_u = (t_rr(1) : 1/fs_i : t_rr(end)).';

    % kubische Spline-Interpolation auf dieses Raster
    sig = interp1(t_rr, rr, t_u, 'spline');

    % linearen Trend + Mittelwert entfernen
    sig = detrend(sig, 1);

    fprintf('Interpoliert: %d Werte @ %.0f Hz | Dauer = %.0f s\n', ...
            numel(sig), fs_i, t_u(end) - t_u(1));
end


% =========================================================================
function [xw, win] = apply_window_function(x, window_type)
%APPLY_WINDOW_FUNCTION  Wendet eine FFT-Fensterfunktion auf ein Segment an.
%   [xw, win] = APPLY_WINDOW_FUNCTION(x, window_type) multipliziert das
%   Signalsegment X elementweise mit der Fensterfunktion und liefert das
%   gefensterte Signal XW SOWIE den Fenstervektor WIN selbst zurueck.
%   WIN wird spaeter zur Leistungs-Normierung der FFT gebraucht
%   (siehe CALCULATE_FFT).
%
%   HINTERGRUND: Die FFT nimmt implizit an, dass sich das endliche Segment
%   periodisch fortsetzt. Passen Anfang und Ende nicht zusammen, entsteht
%   ein "Sprung", der sich im Spektrum als LEAKAGE (Verschmieren von
%   Energie auf Nachbarfrequenzen) zeigt. Fensterfunktionen daempfen die
%   Segmentraender weich gegen Null und reduzieren so das Leakage -- zum
%   Preis einer breiteren Hauptkeule (= geringere Frequenzaufloesung).
%
%   Den reinen Fenstervektor erhaelt man durch Uebergabe eines
%   Einsen-Signals:  [~, w] = apply_window_function(ones(N,1), 'hann');
%
%   Unterstuetzte Fenster: 'rect','hann','hamming','blackman','kaiser','flattop'.

    x = double(x(:));
    N = numel(x);

    switch lower(window_type)
        case 'rect'                       % Rechteck = gar keine Fensterung
            win = ones(N, 1);
        case 'hann'
            win = hann(N);
        case 'hamming'
            win = hamming(N);
        case 'blackman'
            win = blackman(N);
        case 'kaiser'
            win = kaiser(N, 8);           % beta = 8: starke Nebenkeulendaempfung
        case 'flattop'
            win = flattopwin(N);
        otherwise
            error('apply_window_function:unknownWindow', ...
                  'Unbekannter Fenstertyp "%s".', window_type);
    end

    win = win(:);
    xw  = x .* win;                       % elementweise Fensterung
end


% =========================================================================
function [f, psd] = calculate_fft(xw, win, fs)
%CALCULATE_FFT  Einseitige Leistungsdichte (PSD) eines gefensterten Signals.
%   [f, psd] = CALCULATE_FFT(xw, win, fs) berechnet aus dem bereits
%   gefensterten Segment XW (zugehoeriger Fenstervektor WIN, Abtastrate FS)
%   die einseitige spektrale Leistungsdichte PSD und den Frequenzvektor F.
%
%   Berechnung:
%       PSD(f) = |FFT(xw)|^2 / (fs * U),     U = sum(win.^2)
%
%   WARUM die Division durch U = sum(win.^2)?
%   Die Fensterung verringert die Signalenergie. Die Normierung auf die
%   Fensterleistung U gleicht das aus, sodass die INTEGRIERTE Leistung
%   (und damit die HRV-Bandleistung) weitgehend UNABHAENGIG vom gewaehlten
%   Fenster ist. Genau das ist physikalisch korrekt: Die in einem Band
%   enthaltene Leistung ist eine Eigenschaft des Signals, nicht des
%   Analysewerkzeugs. Die Fensterunterschiede zeigen sich folglich nicht
%   in der Bandleistung, sondern in Aufloesung, Leakage und Peak-Form.
%
%   WARUM "einseitig" und der Faktor 2?
%   Das Spektrum eines reellen Signals ist symmetrisch; die negativen
%   Frequenzen tragen dieselbe Leistung wie die positiven. Wir behalten nur
%   die positive Haelfte (0..Nyquist) und verdoppeln daher die inneren Bins,
%   damit die Gesamtleistung erhalten bleibt (DC und Nyquist nicht).

    xw  = double(xw(:));
    win = double(win(:));
    N   = numel(xw);

    X = fft(xw);
    X = X(1:floor(N/2)+1);              % nur 0 .. Nyquist behalten

    U   = sum(win .^ 2);                % Fensterleistung
    psd = (abs(X) .^ 2) / (fs * U);

    if N > 2
        psd(2:end-1) = 2 * psd(2:end-1);   % innere Bins verdoppeln
    end

    f = (0:floor(N/2)).' * (fs / N);    % zugehoerige Frequenzen [Hz]
end


% =========================================================================
function hrv = calculate_hrv_bands(f, psd)
%CALCULATE_HRV_BANDS  Spektrale HRV-Parameter in den Frequenzbaendern.
%   hrv = CALCULATE_HRV_BANDS(f, psd) integriert die Leistungsdichte PSD
%   ueber den drei klassischen HRV-Baendern und gibt eine Struktur zurueck.
%
%   Frequenzbaender (Standard nach Task Force 1996):
%     VLF : 0.0033 - 0.04 Hz   (very low frequency)
%     LF  : 0.04   - 0.15 Hz   (low frequency,  u.a. Barorezeptor-Aktivitaet)
%     HF  : 0.15   - 0.40 Hz   (high frequency, atmungsgekoppelt, vagal)
%
%   Felder der Ausgabestruktur:
%     .VLF/.LF/.HF        absolute Leistung je Band            [ms^2]
%     .Total              Gesamtleistung                       [ms^2]
%     .VLF_rel/.LF_rel/.HF_rel   relative Leistung             [%]
%     .LF_nu/.HF_nu       normalisierte Einheiten (nur LF/HF)  [-]
%     .LF_HF              Verhaeltnis LF/HF (sympatho-vagale Balance)
%
%   Die Leistung in einem Band ist die FLAECHE unter der PSD-Kurve in
%   diesem Frequenzbereich -- numerisch per Trapezregel (TRAPZ) bestimmt.

    f   = f(:);
    psd = psd(:);

    % Bandgrenzen [untere, obere] in Hz
    bands = struct('VLF', [0.0033 0.04], ...
                   'LF',  [0.04   0.15], ...
                   'HF',  [0.15   0.40]);

    hrv   = struct();
    names = fieldnames(bands);
    for i = 1:numel(names)
        rng_b = bands.(names{i});
        m     = (f >= rng_b(1)) & (f < rng_b(2));   % Maske fuer dieses Band
        if nnz(m) >= 2
            hrv.(names{i}) = trapz(f(m), psd(m));    % Flaeche unter der PSD
        else
            hrv.(names{i}) = 0;
        end
    end

    hrv.Total = hrv.VLF + hrv.LF + hrv.HF;

    % relative Leistung (Anteil an der Gesamtleistung, in Prozent)
    if hrv.Total > 0
        hrv.VLF_rel = 100 * hrv.VLF / hrv.Total;
        hrv.LF_rel  = 100 * hrv.LF  / hrv.Total;
        hrv.HF_rel  = 100 * hrv.HF  / hrv.Total;
    else
        hrv.VLF_rel = 0; hrv.LF_rel = 0; hrv.HF_rel = 0;
    end

    % normalisierte Einheiten (n.u.): Anteil von LF bzw. HF an (LF+HF)
    lf_hf = hrv.LF + hrv.HF;
    if lf_hf > 0
        hrv.LF_nu = 100 * hrv.LF / lf_hf;
        hrv.HF_nu = 100 * hrv.HF / lf_hf;
    else
        hrv.LF_nu = 0; hrv.HF_nu = 0;
    end

    % sympatho-vagale Balance
    if hrv.HF > 0
        hrv.LF_HF = hrv.LF / hrv.HF;
    else
        hrv.LF_HF = NaN;
    end
end


% =========================================================================
function results = compare_window_functions(sig, fs_i, windows, seg_len_s, overlap)
%COMPARE_WINDOW_FUNCTIONS  Vergleicht FFT-Fensterfunktionen fuer die HRV.
%   results = COMPARE_WINDOW_FUNCTIONS(sig, fs_i, windows, seg_len_s, overlap)
%   berechnet fuer jedes Fenster aus WINDOWS das gemittelte HRV-Spektrum
%   und stellt die Ergebnisse gegenueber.
%
%   VERFAHREN (WELCH-METHODE):
%   Das Signal wird in Segmente der Laenge SEG_LEN_S mit OVERLAP
%   Ueberlappung zerlegt. Jedes Segment wird gefenstert, in eine PSD
%   ueberfuehrt und die PSDs werden GEMITTELT. Die Mittelung ueber mehrere
%   Segmente reduziert die Varianz (das "Rauschen") der Spektralschaetzung
%   deutlich -- man erhaelt ein stabileres, glatteres Spektrum, allerdings
%   auf Kosten der Frequenzaufloesung (kuerzere Segmente => groebere
%   Aufloesung). Das ist der klassische Bias-Varianz-Kompromiss.
%
%   Ausgabe: Struktur-Array (ein Eintrag je Fenster) mit Feldern
%     .name .f .psd .seg_psd .t_seg .hrv .metrics
%
%   Zusaetzlich werden drei Vergleichsabbildungen erzeugt und eine
%   Vergleichstabelle in der Konsole ausgegeben.

    if nargin < 3 || isempty(windows)
        windows = {'rect','hann','hamming','blackman','kaiser','flattop'};
    end
    if nargin < 4 || isempty(seg_len_s), seg_len_s = 300; end
    if nargin < 5 || isempty(overlap),   overlap   = 0.5; end

    sig = double(sig(:));

    % --- Segmentlaenge in Samples -----------------------------------------
    % seg_len_s ist in Sekunden gegeben (z.B. 300 s = 5 min). Multipliziert
    % mit der Abtastrate fs_i (z.B. 4 Hz) ergibt das die Segmentlaenge in
    % Abtastwerten:  300 s * 4 Hz = 1200 Samples pro 5-min-Segment.
    N   = round(seg_len_s * fs_i);
    if N > numel(sig)
        warning(['Signal kuerzer als ein Segment (%d < %d Samples) -> ', ...
                 'es wird das gesamte Signal als ein Segment verwendet.'], ...
                 numel(sig), N);
        N = numel(sig);
    end

    % --- Schrittweite zwischen den Segment-Startpunkten -------------------
    % overlap ist der Ueberlappungsanteil (0.5 = 50 %). Die Startpunkte
    % ruecken jeweils um (1-overlap)*N weiter. Bei 50 % Ueberlappung also
    % um N/2 -> Segment 1 ab 0 s, Segment 2 ab 2.5 min, Segment 3 ab 5 min ...
    step = max(1, round(N * (1 - overlap)));

    % --- Startindizes aller Segmente --------------------------------------
    % Vom ersten Sample bis dorthin, wo gerade noch ein volles Segment der
    % Laenge N hineinpasst (numel(sig)-N+1). So entstehen mehrere, sich
    % ueberlappende 5-min-Abschnitte, ueber die spaeter gemittelt wird.
    starts = 1:step:(numel(sig) - N + 1);
    if isempty(starts), starts = 1; end
    t_seg = (starts - 1) / fs_i;            % Segment-Startzeiten [s]
    nSeg  = numel(starts);

    fprintf(['\nFenstervergleich: %d Segmente a %.0f s (%.0f%% Ueberlappung), ', ...
             '%d Fenster\n'], nSeg, seg_len_s, overlap*100, numel(windows));

    % Ergebnis-Struktur vorbereiten
    results = struct('name', {}, 'f', {}, 'psd', {}, 'seg_psd', {}, ...
                     't_seg', {}, 'hrv', {}, 'metrics', {});

    % --- je Fenster: Welch-Mittelung ueber alle Segmente ----------------
    for iw = 1:numel(windows)
        wname = windows{iw};
        % Fenstervektor einmal vorab erzeugen (gilt fuer alle Segmente)
        [~, win] = apply_window_function(ones(N, 1), wname);

        seg_psd = [];
        for is = 1:nSeg
            s0  = starts(is);
            seg = sig(s0 : s0 + N - 1);
            seg = seg - mean(seg);          % Mittelwert je Segment abziehen
                                            % (entfernt lokalen DC-Anteil)
            xw  = seg .* win;               % fenstern
            [f, psd] = calculate_fft(xw, win, fs_i);
            if isempty(seg_psd)
                seg_psd = zeros(numel(f), nSeg);
            end
            seg_psd(:, is) = psd;           % PSD dieses Segments speichern
        end

        psd_avg = mean(seg_psd, 2);         % Welch-Mittelung ueber Segmente
        hrv     = calculate_hrv_bands(f, psd_avg);
        metrics = window_metrics(win);      % Kennzahlen des Fensters

        results(iw).name    = wname;
        results(iw).f       = f;
        results(iw).psd     = psd_avg;
        results(iw).seg_psd = seg_psd;
        results(iw).t_seg   = t_seg;
        results(iw).hrv     = hrv;
        results(iw).metrics = metrics;
    end

    % --- Vergleichsdarstellungen + Konsolentabelle ----------------------
    plot_window_characteristics(windows, N);   % Form + Frequenzgang
    plot_overlaid_spectra(results, false);      % ueberlagerte Spektren (linear)
    plot_overlaid_spectra(results, true);       % dito in dB (zeigt Leakage)
    print_comparison_table(results);            % Tabelle ins Befehlsfenster
end


% =========================================================================
function m = window_metrics(win)
%WINDOW_METRICS  Charakteristische Kenngroessen einer Fensterfunktion.
%   Liefert eine Struktur mit
%     .PSL   hoechster Nebenkeulenpegel (peak side lobe level) [dB]
%            -> je negativer, desto besser die Leakage-Unterdrueckung
%     .MLW   Hauptkeulenbreite (main lobe width)               [Bins]
%            -> je kleiner, desto besser die Frequenzaufloesung
%     .ENBW  aequivalente Rauschbandbreite                     [Bins]
%
%   Die Werte werden aus dem stark "gezeropaddeten" (kuenstlich hoeher
%   aufgeloesten) Betragsspektrum des Fensters bestimmt.
    win = win(:);
    N   = numel(win);

    % ENBW: equivalent noise bandwidth in Frequenzbins
    m.ENBW = N * sum(win.^2) / (sum(win)^2);

    % hochaufgeloestes Betragsspektrum des Fensters (starkes Zero-Padding)
    Nfft = 65536;
    W    = abs(fft(win, Nfft));
    W    = W(1:floor(Nfft/2)+1);
    W    = W / W(1);                    % auf den DC-Wert (= Hauptkeulen-Max.) normieren
    WdB  = 20*log10(W + eps);           % in Dezibel
    fbin = (0:numel(W)-1).' * N / Nfft; % Frequenzachse in DFT-Bins der Laenge N

    % Erste echte Nullstelle der Hauptkeule = erstes LOKALES MINIMUM, das
    % mindestens 0.2 Bins von DC entfernt liegt. Die Bedingung "lokales
    % Minimum" (statt "erster Anstieg") macht die Suche robust -- sonst
    % wuerde z.B. der wellige, flache Scheitel des Flat-Top-Fensters
    % faelschlich als Nullstelle erkannt.
    nullidx = [];
    for i = 2:numel(W)-1
        if W(i) < W(i-1) && W(i) < W(i+1) && fbin(i) > 0.2
            nullidx = i; break;
        end
    end

    if ~isempty(nullidx)
        m.MLW = 2 * fbin(nullidx);          % volle Hauptkeulenbreite [Bins]
        m.PSL = max(WdB(nullidx:end));      % hoechste Nebenkeule [dB]
    else
        m.MLW = NaN; m.PSL = NaN;
    end
end


% =========================================================================
function plot_window_characteristics(windows, N)
%PLOT_WINDOW_CHARACTERISTICS  Zeitform und Frequenzgang der Fenster.
%   Zeigt links die Fensterformen im Zeitbereich und rechts ihren
%   Frequenzgang in dB. Im Frequenzgang erkennt man direkt den
%   Kompromiss: schmale Hauptkeule (gute Aufloesung) gegen tiefe
%   Nebenkeulen (wenig Leakage).
    colors = lines(numel(windows));
    figure('Color','w','Position',[100 100 1000 420]);

    % --- (a) Zeitbereich
    subplot(1,2,1); hold on;
    for iw = 1:numel(windows)
        [~, win] = apply_window_function(ones(N,1), windows{iw});
        plot((0:N-1)/(N-1), win, 'LineWidth', 1.4, 'Color', colors(iw,:));
    end
    xlabel('normierte Zeit'); ylabel('Amplitude');
    title('Fensterfunktionen - Zeitbereich');
    legend(windows, 'Location','south', 'Interpreter','none');
    grid on; box on; ylim([0 1.05]); hold off;

    % --- (b) Frequenzgang in dB
    subplot(1,2,2); hold on;
    Nfft = 4096;
    for iw = 1:numel(windows)
        [~, win] = apply_window_function(ones(N,1), windows{iw});
        W   = abs(fft(win, Nfft)); W = W(1:Nfft/2+1);
        WdB = 20*log10(W/max(W) + eps);
        fb  = (0:Nfft/2) / Nfft * N;         % Frequenz in FFT-Bins
        plot(fb, WdB, 'LineWidth', 1.3, 'Color', colors(iw,:));
    end
    xlabel('Frequenz [Bins]'); ylabel('Betrag [dB]');
    title('Fensterfunktionen - Frequenzgang');
    legend(windows, 'Location','northeast', 'Interpreter','none');
    grid on; box on; xlim([0 10]); ylim([-120 5]); hold off;
end


% =========================================================================
function plot_overlaid_spectra(results, logScale)
%PLOT_OVERLAID_SPECTRA  Ueberlagerte HRV-Spektren aller Fenster.
%   logScale = false : lineare Darstellung (Peak-Hoehen vergleichen)
%   logScale = true  : dB-Darstellung (macht den Leakage-Boden sichtbar)
    colors = lines(numel(results));
    figure('Color','w'); hold on;

    fmax = 0.5;                          % dargestellter Frequenzbereich [Hz]
    for iw = 1:numel(results)
        f   = results(iw).f;
        psd = results(iw).psd;
        if logScale
            y = 10*log10(psd + eps);
        else
            y = psd;
        end
        plot(f, y, 'LineWidth', 1.3, 'Color', colors(iw,:));
    end

    % HRV-Bandgrenzen als senkrechte Hilfslinien (plot-basiert, damit es
    % sowohl in MATLAB als auch in Octave funktioniert)
    yl = ylim;
    for fb = [0.04 0.15 0.40]
        plot([fb fb], yl, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
    end
    ylim(yl);

    xlim([0 fmax]); xlabel('Frequenz [Hz]');
    if logScale
        ylabel('Leistungsdichte [dB]');
        title('Ueberlagerte HRV-Spektren (logarithmisch) - Leakage-Vergleich');
    else
        ylabel('Leistungsdichte [ms^2/Hz]');
        title('Ueberlagerte HRV-Spektren (linear)');
    end
    legend({results.name}, 'Location','northeast', 'Interpreter','none');
    grid on; box on; hold off;
end


% =========================================================================
function print_comparison_table(results)
%PRINT_COMPARISON_TABLE  Tabellarischer Vergleich im Befehlsfenster.
    fprintf('\n=================== Vergleich der Fensterfunktionen ===================\n');
    fprintf('%-9s | %7s %6s %6s | %8s %8s %7s | %6s\n', ...
            'Fenster', 'PSL', 'MLW', 'ENBW', 'LF', 'HF', 'LF/HF', 'VLF');
    fprintf('%-9s | %7s %6s %6s | %8s %8s %7s | %6s\n', ...
            '', '[dB]', '[Bin]', '[Bin]', '[ms^2]', '[ms^2]', '[-]', '[ms^2]');
    fprintf('----------------------------------------------------------------------\n');
    for iw = 1:numel(results)
        m = results(iw).metrics;
        h = results(iw).hrv;
        fprintf('%-9s | %7.1f %6.2f %6.2f | %8.1f %8.1f %7.2f | %6.1f\n', ...
                results(iw).name, m.PSL, m.MLW, m.ENBW, ...
                h.LF, h.HF, h.LF_HF, h.VLF);
    end
    fprintf('======================================================================\n');
    fprintf(['Hinweis: PSL = hoechste Nebenkeule (je negativer, desto weniger ', ...
             'Leakage),\n         MLW = Hauptkeulenbreite (je kleiner, desto ', ...
             'bessere Aufloesung),\n         ENBW = aequivalente Rauschbandbreite.\n']);
end


% =========================================================================
function plot_spectrum(f, psd, window_name, hrv)
%PLOT_SPECTRUM  Stellt ein HRV-Leistungsdichtespektrum dar.
%   PLOT_SPECTRUM(f, psd, window_name) zeichnet die PSD im HRV-relevanten
%   Bereich 0..0.5 Hz und hinterlegt die drei Baender (VLF/LF/HF) farblich.
%   Mit dem optionalen Argument HRV wird zusaetzlich ein Kennwerte-Textfeld
%   (LF/HF, Bandleistungen) eingeblendet.

    if nargin < 3 || isempty(window_name), window_name = ''; end

    figure('Color','w'); hold on;

    fmax = 0.5;
    ymax = 1.05 * max(psd(f <= fmax));
    if ~isfinite(ymax) || ymax <= 0, ymax = 1; end

    shade_hrv_bands(ymax);               % HRV-Baender farblich hinterlegen
    plot(f, psd, 'k', 'LineWidth', 1.3); % eigentliches Spektrum

    xlim([0 fmax]); ylim([0 ymax]);
    xlabel('Frequenz [Hz]'); ylabel('Leistungsdichte [ms^2/Hz]');
    title(sprintf('HRV-Leistungsdichtespektrum - Fenster: %s', window_name));
    grid on; box on;
    legend({'VLF (0.0033-0.04 Hz)', 'LF (0.04-0.15 Hz)', ...
            'HF (0.15-0.40 Hz)', 'PSD'}, 'Location','northeast');

    % optionales Kennwerte-Textfeld
    if nargin >= 4 && ~isempty(hrv)
        txt = sprintf(['LF/HF = %.2f\nLF = %.1f ms^2 (%.0f%%)\n', ...
                       'HF = %.1f ms^2 (%.0f%%)'], ...
                       hrv.LF_HF, hrv.LF, hrv.LF_rel, hrv.HF, hrv.HF_rel);
        text(0.97*fmax, 0.95*ymax, txt, ...
             'HorizontalAlignment','right', 'VerticalAlignment','top', ...
             'BackgroundColor','w', 'EdgeColor',[0.6 0.6 0.6], 'FontSize',9);
    end
    hold off;
end


% =========================================================================
function shade_hrv_bands(ymax)
%SHADE_HRV_BANDS  Hinterlegt die drei HRV-Frequenzbaender farblich.
    bands = [0.0033 0.04;      % VLF
             0.04   0.15;      % LF
             0.15   0.40];     % HF
    cols  = [0.85 0.90 0.98;   % hellblau   (VLF)
             0.80 0.95 0.85;   % hellgruen  (LF)
             0.99 0.88 0.82];  % hellorange (HF)
    for i = 1:size(bands,1)
        x = bands(i,:);
        patch([x(1) x(2) x(2) x(1)], [0 0 ymax ymax], cols(i,:), ...
              'EdgeColor','none', 'FaceAlpha',0.7);
    end
end


% =========================================================================
function plot_waterfall(f, seg_psd, t_seg, window_name)
%PLOT_WATERFALL  3D-Wasserfalldarstellung der zeitlichen Spektrenfolge.
%   PLOT_WATERFALL(f, seg_psd, t_seg, window_name) zeichnet die zeitliche
%   Entwicklung des HRV-Spektrums: jede "Linie" des Wasserfalls ist das
%   Spektrum eines Analysesegments. So wird sichtbar, ob die LF-/HF-Anteile
%   ueber die Aufzeichnung stabil bleiben oder schwanken.
%
%   Eingaben:
%     f           : Frequenzvektor [Hz]
%     seg_psd     : Matrix der Segment-PSDs [numel(f) x nSeg]
%     t_seg       : Startzeiten der Segmente [s]
%     window_name : Name des Fensters (fuer den Titel)

    if nargin < 4 || isempty(window_name), window_name = ''; end

    f   = f(:);
    sel = f <= 0.5;                      % nur HRV-relevanter Bereich
    fz  = f(sel);
    Z   = seg_psd(sel, :).';             % Zeilen = Segmente, Spalten = Frequenz

    nSeg = size(Z, 1);
    if nargin < 3 || isempty(t_seg), t_seg = 1:nSeg; end
    t_seg = t_seg(:);

    figure('Color','w');
    waterfall(fz, t_seg, Z);             % 3D-Wasserfall
    colormap(turbo);
    xlabel('Frequenz [Hz]'); ylabel('Segment-Startzeit [s]');
    zlabel('Leistungsdichte [ms^2/Hz]');
    title(sprintf('HRV-Wasserfalldiagramm - Fenster: %s', window_name));
    grid on; box on; view(40, 30);       % Blickwinkel

    % --- HRV-Bandgrenzen als senkrechte Ebenen markieren ---------------
    % Jede Ebene wird in zwei Dreiecke zerlegt (Dreiecks-Patches sind in
    % MATLAB und Octave gleichermassen robust).
    hold on;
    zmax = max(Z(:)); if ~isfinite(zmax) || zmax <= 0, zmax = 1; end
    ymin = min(t_seg); ymax = max(t_seg); if ymax == ymin, ymax = ymin + 1; end
    for fb = [0.04, 0.15, 0.40]          % LF-Beginn, LF/HF-Grenze, HF-Ende
        patch([fb fb fb], [ymin ymax ymax], [0 0 zmax], [0.3 0.3 0.3], ...
              'FaceAlpha',0.12, 'EdgeColor',[0.3 0.3 0.3], 'LineStyle','--');
        patch([fb fb fb], [ymin ymax ymin], [0 zmax zmax], [0.3 0.3 0.3], ...
              'FaceAlpha',0.12, 'EdgeColor','none');
    end
    hold off;
end
