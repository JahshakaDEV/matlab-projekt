%% ========================================================================
%  HRV ANALYSIS: INFLUENCE OF DIFFERENT FFT WINDOW FUNCTIONS
%  ------------------------------------------------------------------------
%  Modular layout: every processing step is implemented as its own (local)
%  function at the end of this file.
%
%  Processing chain:
%    load_ecg_data -> preprocess_ecg -> detect_r_peaks ->
%    calculate_rr_intervals -> interpolate_rr_signal ->
%    compare_window_functions ( apply_window_function, calculate_fft,
%    calculate_hrv_bands ) -> plot_spectrum / plot_waterfall
%
%  Requires: Signal Processing Toolbox.
%  ========================================================================

clear; clc; close all;

%% ----------------------- Parameters ------------------------------------
EDF_FILE   = '16-49-27.EDF';
FS_INTERP  = 4;       % resampling frequency of the RR signal [Hz]
SEG_LEN_S  = 300;     % analysis window length [s] = 5 min
OVERLAP    = 0.5;     % segment overlap
WINDOWS    = {'rect','hann','hamming','blackman','kaiser','flattop'};

%% --------------------------- 1. Load ECG -------------------------------
% Aufgabe 2.1

fprintf('=== 1. EKG-Daten laden ===\n');
[ecg_raw, fs] = load_ecg_data(EDF_FILE);
t = (0:numel(ecg_raw)-1).' / fs; % Zeilenvektor der Sample Indizes -1 (0,1,2,3,..,N-1)
% . transporniert den Zeilenvektor zum Spaltenvektor / fs teilt den Index durch die Abtastrate,
% wandelt Also Sample Nr. N in Zeitpunkt N/fs Sekunden um

%% ------------------------ 2. Preprocessing -----------------------------
fprintf('\n=== 2. Vorverarbeitung ===\n');
ecg = preprocess_ecg(ecg_raw, fs);

% Raw vs. filtered signal (excerpt 10-20 s)
figure('Color','w','Position',[100 100 950 500]);
sel = t >= 0;                                    % select the whole signal (xlim below only zooms the view)
subplot(2,1,1);
plot(t(sel), ecg_raw(sel), 'Color', [0.6 0.6 0.6]);
title('EKG - Rohsignal (Ausschnitt 10-20 s)'); ylabel('Amplitude');
grid on; box on; xlim([10 20]);
subplot(2,1,2);
plot(t(sel), ecg(sel), 'Color', [0.85 0.2 0.2]);
title('EKG - nach Vorverarbeitung (Ausschnitt 10-20 s)'); xlabel('Zeit [s]'); ylabel('Amplitude');
grid on; box on; xlim([10 20]);

%% --------------------------- 3. R peaks --------------------------------
% Aufgabe 2.2

fprintf('\n=== 3. R-Zacken-Erkennung ===\n');
r_locs = detect_r_peaks(ecg, fs);

figure('Color','w','Position',[100 100 950 350]);
plot(t(sel), ecg(sel), 'k'); hold on;
r_sel = r_locs; % all detected R peaks (xlim below only limits the visible view)
plot(r_sel/fs, ecg(r_sel), 'ro', 'MarkerFaceColor','r', 'MarkerSize',5);
title('Erkannte R-Zacken');
xlabel('Zeit [s]'); ylabel('Amplitude');
legend({'EKG','R-Zacken'}, 'Location','northeast');
grid on; box on; xlim([0 10]); hold off;

%% ------------------------- 4. RR intervals -----------------------------
fprintf('\n=== 4. RR-Intervalle + Artefaktkorrektur ===\n');
[t_rr, rr, n_art] = calculate_rr_intervals(r_locs, fs);

figure('Color','w','Position',[100 100 950 350]);
plot(t_rr, rr, '-', 'Color', [0.1 0.4 0.8]); hold on;
plot(t_rr, rr, '.', 'Color', [0.1 0.4 0.8], 'MarkerSize', 6);
title(sprintf('Tachogramm (RR-Intervalle) - %d Artefakte korrigiert', n_art));
xlabel('Zeit [s]'); ylabel('RR-Intervall [ms]');
grid on; box on; hold off;

%% ------------------------- 5. Interpolation ----------------------------
% Aufgabe 2.3

fprintf('\n=== 5. Interpolation (Resampling) ===\n');
[~, sig, fs_i] = interpolate_rr_signal(t_rr, rr, FS_INTERP);

%% ----------------------- 6. Window comparison --------------------------
% Aufgabe 2.4
fprintf('\n=== 6. Vergleich der Fensterfunktionen ===\n');
results = compare_window_functions(sig, fs_i, WINDOWS, SEG_LEN_S, OVERLAP);

%% --------------------- 7. Individual plots -----------------------------
fprintf('\n=== 7. Einzelspektren und Wasserfalldiagramme ===\n');

% Common z maximum across all windows .
zmax_all = 0;
for iw = 1:numel(results)
    sel_f  = results(iw).f <= 0.5;
    Zi     = results(iw).seg_psd(sel_f, :);
    zmax_all = max(zmax_all, max(Zi(:)));
end

for iw = 1:numel(results)
    R = results(iw);
    plot_spectrum(R.f, R.psd, R.name, R.hrv);
    plot_waterfall(R.f, R.seg_psd, R.t_seg, R.name, zmax_all);
end

% Export all open figures to a PDF (one figure per page).
pdfName = 'all_figures.pdf';
if isfile(pdfName), delete(pdfName); end          % remove the old export (otherwise pages are appended)

figHandles = findall(0, 'Type', 'figure');
[~, order] = sort([figHandles.Number]);           % logical order (ECG first)
figHandles = figHandles(order);
for i = 1:numel(figHandles)
    exportgraphics(figHandles(i), pdfName, 'Append', true);
end

fprintf('\n=== Analyse abgeschlossen ===\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% =========================================================================
function [ecg, fs] = load_ecg_data(filename)
%LOAD_ECG_DATA  Read an ECG signal and its sampling rate from an EDF file.
%   [ECG, FS] = LOAD_ECG_DATA(FILENAME) reads the sampling rate from the
%   EDF header and joins the record blocks of the first signal channel into
%   one continuous vector.
%
%   Input:
%     FILENAME - path/name of the EDF file, e.g. '16-49-27.EDF'.
%
%   Output:
%     ECG - raw ECG signal as a double column vector.
%     FS  - sampling rate in Hz, derived from the header.

info   = edfinfo(filename); % Creates edfinfo object, contain information such as file size, number of data records, number of signals, and number of samples.
recDur = seconds(info.DataRecordDuration); % Milli in Sekunden
fs     = double(info.NumSamples(1)) / recDur; % Abtastrate auslesen

tt    = edfread(filename); % eig EKG Signal auslesen --> timetable
vname = tt.Properties.VariableNames{1}; % Name der ersten Spalte / ersten Kanals.
col   = tt.(vname); % extrahiert die Spalte aus dem <Timetable --> Rohdaten erster Kanal
if iscell(col), ecg = cell2mat(col); else, ecg = col(:); end
% Ist col ein Cell Array --> cell2mat fügt alle Blöcke zu einem durchgehenden Vektor zusammen
% Ist es bereits numerisch --> col(:) erzwingt nur einen Spaltenvektor
ecg = double(ecg(:)); % Stellt sicher das EDC ein Spaltenvektor mit double ist,
% double wird für FFT gebraucht

fprintf('EDF geladen: "%s" | Kanal "%s" | fs = %.1f Hz | Dauer = %.1f min\n', ...
    filename, vname, fs, numel(ecg)/fs/60);
end

% =========================================================================
function ecg_filt = preprocess_ecg(ecg, fs)
%PREPROCESS_ECG  Remove DC, baseline drift, 50 Hz hum and high-frequency noise.
%   ECG_FILT = PREPROCESS_ECG(ECG, FS) applies, in this order, DC removal, a
%   0.5 Hz high-pass, a 50 Hz notch and a 40 Hz low-pass. All filters are
%   zero-phase (filtfilt), so the R peaks are not shifted in time.
%
%   Input:
%     ECG - raw ECG signal (vector).
%     FS  - sampling rate of ECG in Hz.
%
%   Output:
%     ECG_FILT - filtered ECG, same length as ECG.

ecg_filt = ecg - mean(ecg); % Entfernt den DC-Offset, Zentrierung um 0
ecg_filt = highpass(ecg_filt,0.5, fs); % Hochpassfilter, Entfernt sehr langsame Schwankungen, Baseline Wander
% bspw durch Atmung oder Elektrodenbewegung entstehen und sonst die R-Zacken Erkennung Stören. Alles unter 0.5 wird gedämpft

wo = 50/(fs/2); % fs/2 ist Nyquist Frequenz, Matlab Filterfunktionen erwarten normierte Frequenzen im Bereich 0-1
% wobei 1 = Nyquist entsprcucht. wo ist also 50Hz normiert auf diesen Bereicht
bw = wo/40;v% Q = 40, Bandbreite, ein hohes Q bedeutet eine sehr schamle scharfe Kerbe, es soll wirklich nur exakt 50Hz entfernen werden
[b,a] = iirnotch(wo, bw); % Entwirft die Filterkoeffizienten (b, a) für dieses schmale Kerbfilter (IIR = Infinite Impulse Response)
% bei genau der Frequenz wo mit Bandbreite bw.
ecg_filt = filtfilt(b, a, ecg_filt); % wendet den Filter an, 1. vorwärts 1. rückwärts, wegen Phasenverschiebung, R-Zacken werden nicht verschoben

ecg_filt = lowpass(ecg_filt, 40, fs); % Tiefpass mit Grenzfrequewnz 40Hz, Entfernt hochfrequentes Rauschen (z.B. Muskelzittern/EMG-Störungen,
% hochfrequente Messartefakte) oberhalb 40 Hz, da die relevante EKG-Information (inkl. QRS-Komplex) darunter liegt.
end

% =========================================================================
function locs = detect_r_peaks(ecg, fs)
%DETECT_R_PEAKS  Locate the R peaks of an ECG with findpeaks.
%   LOCS = DETECT_R_PEAKS(ECG, FS) accepts a peak only if it satisfies two
%   criteria: a minimum height (so only clearly protruding peaks count as an
%   R peak) and a minimum distance (physiologically there is no second beat
%   within 300 ms).
%
%   Input:
%     ECG - preprocessed ECG signal (vector).
%     FS  - sampling rate of ECG in Hz.
%
%   Output:
%     LOCS - sample indices of the detected R peaks.

schwelle       = mean(ecg) + 0.6*std(ecg); % Mittelwert des Signals + 0.6*der Standardabweichung
RR_min_samples = round(0.3*fs); % 0.3 s -> upper limit of ~200 bpm
[~, locs] = findpeaks(ecg, ...
    'MinPeakHeight', schwelle, 'MinPeakDistance', RR_min_samples); % alle lokale Maxima im ecg, die über schwelle liegen
% und mindesetens RR_min_samples Samples voneinander entferne sidn. Peak Werte werden verworfen ~, nur die Positionen locs
fprintf('R-Zacken erkannt: %d\n', numel(locs));
end

% =========================================================================
function [t_rr, rr, n_artifacts] = calculate_rr_intervals(r_locs, fs)
%CALCULATE_RR_INTERVALS  Derive RR intervals from R peaks and repair artefacts.
%   [T_RR, RR, N_ARTIFACTS] = CALCULATE_RR_INTERVALS(R_LOCS, FS) computes the
%   beat-to-beat distances and flags a value as an artefact if it lies outside
%   300-2000 ms OR is a statistical outlier with respect to the moving median
%   (isoutlier). Flagged values are replaced by pchip interpolation from their
%   valid neighbours.
%
%   Input:
%     R_LOCS - sample indices of the R peaks (see DETECT_R_PEAKS).
%     FS     - sampling rate the indices refer to, in Hz.
%
%   Output:
%     T_RR        - time stamp of each interval in s (that of the later beat).
%     RR          - RR intervals in ms, artefacts already corrected.
%     N_ARTIFACTS - number of corrected values.
r_locs = r_locs(:); % erzwingt Spaltenvektor
t_r  = r_locs / fs; % Wandelt die Sample-Indizes der R-Zacken in Zeitpunkte (sek) um.
rr   = diff(t_r) * 1000; % [ms] different zwischen Aufeinanderfolgenden Zeitpunkte, also benachbarte Herzschläge
t_rr = t_r(2:end); % Ordnet jedem RR-Intervall einem Zeitstempel zu, der spätere Zeitpunkt der beiden Schläge

% 300-2000 ms corresponds to ~30-200 bpm. The MAD-based movmedian check on top
% of it adapts to the local spread, so it also catches single missed/false
% beats in passages where the rate itself fluctuates.
implausible = (rr < 300) | (rr > 2000); % Alle werte außerhalb 30-200bpm raus
outlier     = isoutlier(rr, 'movmedian', 21); % Erkennt statistische Ausreißer relativ zu einem gleitenden Median,
% über ein Fenster von 21 Werte.
bad = implausible | outlier; % kobiniert beide Kriterien als bad
n_artifacts = sum(bad); % zählt wie viele RR-Werte als Artefakt markiert wurden

if n_artifacts > 0 && sum(~bad) >= 2
    rr(bad) = interp1(t_rr(~bad), rr(~bad), t_rr(bad), 'pchip', 'extrap');
    %  ersetzt werte mit interp1 mit der Methode 'pchip' (stückweise kubische Hermite-Interpolation,
    % glatt und ohne Überschwingen) schätzt die Werte an den Zeitpunkten der Artefakte (t_rr(bad))
    % anhand der gültigen Nachbarwerte (t_rr(~bad), rr(~bad)). 'extrap' erlaubt Extrapolation,
    % falls ein Artefakt am Rand liegt, wo keine Nachbarwerte auf beiden Seiten existieren.
end

fprintf(['RR-Intervalle: %d | Artefakte: %d | mittl. RR = %.1f ms ', ...
    '(%.1f bpm) | SDNN = %.1f ms\n'], numel(rr), n_artifacts, ...
    mean(rr), 60000/mean(rr), std(rr));
end

% =========================================================================
function [t_u, sig, fs_i] = interpolate_rr_signal(t_rr, rr, fs_i)
%INTERPOLATE_RR_SIGNAL  Resample the RR series onto a uniform grid and detrend it.
%   [T_U, SIG, FS_I] = INTERPOLATE_RR_SIGNAL(T_RR, RR, FS_I) spline-interpolates
%   the RR values onto an evenly spaced time grid and removes the linear trend.
%   The FFT assumes uniform sampling, but the RR values arrive at the irregular
%   times of the heartbeats.
%
%   Input:
%     T_RR - time stamps of the RR intervals in s.
%     RR   - RR intervals in ms.
%     FS_I - target sampling rate in Hz (optional, default 4).
%
%   Output:
%     T_U  - uniform time grid in s.
%     SIG  - interpolated, linearly detrended RR signal in ms.
%     FS_I - the sampling rate actually used, in Hz.
% Diese Funktion resampelt das unregelmäßig abgetastete RR-Intervall-Signal
% auf ein gleichmäßiges Zeitraster (Grundvoraussetzung für die FFT) und entfernt den linearen Trend.
if nargin < 3 || isempty(fs_i), fs_i = 4; end % falls kein 3. Arg --> 4Hz als Standardwert
t_rr = t_rr(:);  rr = rr(:); % Erzwingt Spaltenvektoren

t_u = (t_rr(1) : 1/fs_i : t_rr(end)).'; % Erzeugt das neue gleichmäßige Zeitraster: beginnend beim ersten RR-Zeitstempel
% endend beim letzten, mit Schrittweite 1/fs_i Sekunden. .' transponiert das Ergebnis zu einem Spaltenvektor
sig = interp1(t_rr, rr, t_u, 'spline'); % Imterpoliert die RR-Werte, die ursprünglich zu unregelmäßigen Zeitpunkten vorliegen,
% auf die neuen gleichmäßigen Zeitpunkte t_u. spline sorgt für einen glatten Kurvenveerlauf zwischen den Stützstellen
sig = detrend(sig, 1); % Entfernt einen linearen Trend aus dem Signal (1 --> linear, 0 wäre nur DC Offset), Ein langsamer
% linearer Anstieg/Abfall der HF über die Aufnahme hinweg hat keine Bedeutung für die interessierenden HRV-Frequenzenbänder
% würde aber als starke, tieffrequente Komponente die FFT verfäschen, deshalb entfernen

fprintf('Interpoliert: %d Werte @ %.0f Hz | Dauer = %.0f s\n', ...
    numel(sig), fs_i, t_u(end)-t_u(1));
end

% =========================================================================
function [xw, win] = apply_window_function(x, window_type)
%APPLY_WINDOW_FUNCTION  Multiply a segment by a named window function.
%   [XW, WIN] = APPLY_WINDOW_FUNCTION(X, WINDOW_TYPE) returns the windowed
%   signal XW and the window vector WIN itself (WIN is needed for the power
%   normalisation in CALCULATE_FFT). Windows taper the segment edges and thus
%   reduce spectral leakage, but they widen the main lobe (lower resolution).
%
%   Input:
%     X           - segment to be windowed (vector of length N).
%     WINDOW_TYPE - 'rect', 'hann', 'hamming', 'blackman', 'kaiser' or
%                   'flattop'; anything else raises an error.
%
%   Output:
%     XW  - windowed segment, column vector of length N.
%     WIN - window vector used, column vector of length N.
%
%   To obtain the bare window vector:
%     [~, w] = apply_window_function(ones(N,1), 'hann');

x = double(x(:));
N = numel(x);
switch lower(window_type)
    case 'rect',     win = ones(N,1);     % no windowing
    case 'hann',     win = hann(N);
    case 'hamming',  win = hamming(N);
    case 'blackman', win = blackman(N);
    case 'kaiser',   win = kaiser(N, 8);  % beta = 8
    case 'flattop',  win = flattopwin(N);
    otherwise
        error('apply_window_function:unknownWindow', ...
            'Unbekannter Fenstertyp "%s".', window_type);
end
win = win(:);
xw  = x .* win;
end

% =========================================================================
function [f, psd] = calculate_fft(xw, win, fs)
%CALCULATE_FFT  Compute the one-sided, power-normalised PSD of a windowed segment.
%   [F, PSD] = CALCULATE_FFT(XW, WIN, FS) evaluates
%   PSD = |FFT(XW)|^2 / (FS*U) with U = sum(WIN.^2). Normalising by the window
%   power U makes the integrated (band) power largely independent of the window.
%   The inner bins are doubled (one-sided spectrum of a real signal).
%
%   Input:
%     XW  - windowed segment (vector of length N, see APPLY_WINDOW_FUNCTION).
%     WIN - the window vector that was applied, same length as XW.
%     FS  - sampling rate of the segment in Hz.
%
%   Output:
%     F   - frequency axis in Hz, floor(N/2)+1 points from 0 to FS/2.
%     PSD - power spectral density in signal units^2/Hz, same size as F.

xw  = double(xw(:));  win = double(win(:));
N   = numel(xw);

X = fft(xw);
X = X(1:floor(N/2)+1);
U   = sum(win .^ 2);
psd = (abs(X).^2) / (fs * U);
if N > 2, psd(2:end-1) = 2*psd(2:end-1); end

f = (0:floor(N/2)).' * (fs / N);
end

% =========================================================================
function hrv = calculate_hrv_bands(f, psd)
%CALCULATE_HRV_BANDS  Integrate a PSD into the VLF/LF/HF bands and derive HRV measures.
%   HRV = CALCULATE_HRV_BANDS(F, PSD) computes the band power as the area under
%   the PSD (Trapez-Formel):
%   VLF 0.0033-0.04 Hz, LF 0.04-0.15 Hz, HF 0.15-0.40 Hz.
%
%   Input:
%     F   - frequency axis in Hz (see CALCULATE_FFT).
%     PSD - power spectral density belonging to F, in ms^2/Hz.
%
%   Output:
%     HRV - struct with the absolute band powers VLF, LF, HF and Total in ms^2,
%           the relative shares VLF_rel/LF_rel/HF_rel in % of Total, the
%           normalised units LF_nu/HF_nu in % of (LF+HF), and the ratio LF_HF
%           (NaN if HF is zero).

f = f(:);  psd = psd(:);
bands = struct('VLF',[0.0033 0.04], 'LF',[0.04 0.15], 'HF',[0.15 0.40]);

hrv = struct();
names = fieldnames(bands);
for i = 1:numel(names)
    r = bands.(names{i});
    m = (f >= r(1)) & (f < r(2));
    if nnz(m) >= 2, hrv.(names{i}) = trapz(f(m), psd(m));
    else,           hrv.(names{i}) = 0; end
end

hrv.Total = hrv.VLF + hrv.LF + hrv.HF;
if hrv.Total > 0
    hrv.VLF_rel = 100*hrv.VLF/hrv.Total;
    hrv.LF_rel  = 100*hrv.LF /hrv.Total;
    hrv.HF_rel  = 100*hrv.HF /hrv.Total;
else
    hrv.VLF_rel = 0; hrv.LF_rel = 0; hrv.HF_rel = 0;
end

lf_hf = hrv.LF + hrv.HF;                  % normalised units
if lf_hf > 0
    hrv.LF_nu = 100*hrv.LF/lf_hf;  hrv.HF_nu = 100*hrv.HF/lf_hf;
else
    hrv.LF_nu = 0; hrv.HF_nu = 0;
end

if hrv.HF > 0, hrv.LF_HF = hrv.LF/hrv.HF; else, hrv.LF_HF = NaN; end
end

% =========================================================================
function results = compare_window_functions(sig, fs_i, windows, seg_len_s, overlap)
%COMPARE_WINDOW_FUNCTIONS  Estimate a Welch spectrum per window function and compare them.
%   RESULTS = COMPARE_WINDOW_FUNCTIONS(SIG, FS_I, WINDOWS, SEG_LEN_S, OVERLAP)
%   splits the signal into segments, windows each segment and forms its PSD,
%   then averages the PSDs (Welch), which lowers the variance of the estimate.
%   It also draws the comparison plots and prints the console table.
%
%   Input:
%     SIG       - uniformly sampled RR signal (see INTERPOLATE_RR_SIGNAL).
%     FS_I      - sampling rate of SIG in Hz.
%     WINDOWS   - cell array of window names (optional, default all six).
%     SEG_LEN_S - segment length in s (optional, default 300).
%     OVERLAP   - segment overlap as a fraction 0..1 (optional, default 0.5).
%
%   Output:
%     RESULTS - 1-by-numel(WINDOWS) struct array with the fields
%               .name    window name
%               .f       frequency axis in Hz
%               .psd     Welch-averaged PSD in ms^2/Hz
%               .seg_psd PSD of every segment, one column per segment
%               .t_seg   start time of every segment in s
%               .hrv     band measures (see CALCULATE_HRV_BANDS)
%               .metrics window measures (see WINDOW_METRICS)

if nargin < 3 || isempty(windows) % falls kein 3. Arg --> alle 6 Fenster als Standard
    windows = {'rect','hann','hamming','blackman','kaiser','flattop'};
end
if nargin < 4 || isempty(seg_len_s), seg_len_s = 300; end % Standard Segmentlaenge 300s = 5min
if nargin < 5 || isempty(overlap),   overlap   = 0.5; end % Standard Ueberlappung 50%

sig = double(sig(:)); % erzwingt double Spaltenvektor
N   = round(seg_len_s * fs_i); % Segmentlaenge in Samples = Dauer[s] * Abtastrate[Hz]
if N > numel(sig) % Signal kuerzer als ein Segment abfangen
    warning('Signal kuerzer als ein Segment -> ganzes Signal als ein Segment.');
    N = numel(sig); % dann ganzes Signal als ein Segment
end
step   = max(1, round(N*(1-overlap))); % Versatz zwischen Segmentanfaengen, bei 50% = halbe Segmentlaenge
starts = 1:step:(numel(sig)-N+1); % Startindizes aller Segmente, so dass jedes noch komplett ins Signal passt
if isempty(starts), starts = 1; end % Sicherheit: mind. ein Segment ab Index 1
t_seg  = (starts-1) / fs_i; % Startzeit jedes Segments in Sekunden
nSeg   = numel(starts); % Anzahl der Segmente

fprintf('\nFenstervergleich: %d Segmente a %.0f s (%.0f%% Ueberlappung), %d Fenster\n', ...
    nSeg, seg_len_s, overlap*100, numel(windows));

results = struct('name',{}, 'f',{}, 'psd',{}, 'seg_psd',{}, ... % leeres Ergebnis-Struct-Array vorbereiten
    't_seg',{}, 'hrv',{}, 'metrics',{});
for iw = 1:numel(windows) % ueber jede Fensterfunktion iterieren
    wname = windows{iw}; % Name des aktuellen Fensters
    [~, win] = apply_window_function(ones(N,1), wname); % nur den Fenstervektor holen (Eingabe = Einsen)

    seg_psd = zeros(floor(N/2)+1, nSeg); % Speicher: eine PSD-Spalte pro Segment
    for is = 1:nSeg % ueber jedes Segment iterieren
        s0  = starts(is); % Startindex dieses Segments
        seg = sig(s0:s0+N-1); % Segment aus dem Signal ausschneiden
        seg = seg - mean(seg);              % local DC would leak into the VLF band
        [f, psd] = calculate_fft(seg.*win, win, fs_i); % Fenster anwenden und PSD berechnen
        seg_psd(:, is) = psd; % PSD dieses Segments abspeichern
    end

    psd_avg = mean(seg_psd, 2);             % Welch averaging
    results(iw).name    = wname; % Fenstername
    results(iw).f       = f; % Frequenzachse in Hz
    results(iw).psd     = psd_avg; % gemittelte PSD
    results(iw).seg_psd = seg_psd; % alle Einzelsegment-PSDs (fuer Wasserfall)
    results(iw).t_seg   = t_seg; % Startzeiten der Segmente
    results(iw).hrv     = calculate_hrv_bands(f, psd_avg); % HRV-Bandleistungen (VLF/LF/HF)
    results(iw).metrics = window_metrics(win); % Fenster-Guetemasse (PSL/MLW/ENBW)
end

plot_window_characteristics(windows, N); % Fensterformen Zeit- und Frequenzbereich plotten
plot_overlaid_spectra(results, false); % ueberlagerte Spektren linear
plot_overlaid_spectra(results, true); % ueberlagerte Spektren logarithmisch (dB)
print_comparison_table(results); % Vergleichstabelle in die Konsole drucken
end

% =========================================================================
function m = window_metrics(win)
%WINDOW_METRICS  Measure the leakage, resolution and noise bandwidth of a window.
%   M = WINDOW_METRICS(WIN) derives all three figures from the frequency
%   response of the window using built-in functions.
%
%   Input:
%     WIN - window vector of length N (see APPLY_WINDOW_FUNCTION).
%
%   Output:
%     M - struct with the fields
%         .PSL  peak side lobe in dB relative to the main lobe (leakage measure)
%         .MLW  main lobe width in FFT bins (resolution measure)
%         .ENBW equivalent noise bandwidth in FFT bins

win  = win(:);  N = numel(win); % Spaltenvektor erzwingen, Fensterlaenge N
m.ENBW = enbw(win);                          % equivalent noise bandwidth

% Frequency response of the window, finely resolved by zero padding
Nfft = 8192; % starkes Zero-Padding fuer feine Frequenzaufloesung
W    = abs(fft(win, Nfft)); % Betrag des Frequenzgangs des Fensters
W    = W(1:Nfft/2+1); % nur die einseitige (positive) Haelfte behalten
fbin = (0:Nfft/2).' * N / Nfft;              % frequency axis in FFT bins
WdB  = 20*log10(W/max(W) + eps);             % normalised, main peak = 0 dB

% First null = first local minimum beyond 0.5 bins (the main lobe is always
% wider, so ripple in the flat top of the flat-top window does not interfere).
nulls     = fbin(islocalmin(W)); % Frequenzen aller lokalen Minima (Nullstellen)
nulls     = nulls(nulls > 0.5); % nur die jenseits 0.5 Bins (Hauptkeule ausschliessen)
firstNull = nulls(1); % erste Nullstelle = Rand der Hauptkeule

m.MLW = 2 * firstNull;                       % main lobe width = 2 x first null
m.PSL = max(WdB(fbin > firstNull));          % highest side lobe (everything beyond the main lobe)
end

% =========================================================================
function plot_window_characteristics(windows, N)
%PLOT_WINDOW_CHARACTERISTICS  Plot the window shapes (time) next to their frequency responses (dB).
%   PLOT_WINDOW_CHARACTERISTICS(WINDOWS, N) opens one figure with two subplots
%   that show why the windows differ in leakage and resolution.
%
%   Input:
%     WINDOWS - cell array of window names.
%     N       - window length in samples the windows are built with.

colors = lines(numel(windows)); % eine unterscheidbare Farbe pro Fenster
figure('Color','w','Position',[100 100 1000 420]);

subplot(1,2,1); hold on; % linker Plot: Zeitbereich (Fensterform)
for iw = 1:numel(windows) % ueber jedes Fenster
    [~, win] = apply_window_function(ones(N,1), windows{iw}); % Fenstervektor holen
    plot((0:N-1)/(N-1), win, 'LineWidth', 1.4, 'Color', colors(iw,:)); % Form ueber normierte Zeit 0..1 plotten
end
xlabel('normierte Zeit'); ylabel('Amplitude');
title('Fensterfunktionen - Zeitbereich');
legend(windows, 'Location','south', 'Interpreter','none');
grid on; box on; ylim([0 1.05]); hold off;

subplot(1,2,2); hold on; % rechter Plot: Frequenzgang in dB
Nfft = 4096; % Zero-Padding fuer glatten Frequenzgang
for iw = 1:numel(windows) % ueber jedes Fenster
    [~, win] = apply_window_function(ones(N,1), windows{iw}); % Fenstervektor holen
    W   = abs(fft(win, Nfft)); W = W(1:Nfft/2+1); % Betragsspektrum, einseitig
    WdB = 20*log10(W/max(W) + eps); % auf Hauptkeule normiert, in dB (eps gegen log(0))
    fb  = (0:Nfft/2) / Nfft * N; % Frequenzachse in FFT-Bins
    plot(fb, WdB, 'LineWidth', 1.3, 'Color', colors(iw,:)); % Hauptkeule + Nebenkeulen plotten
end
xlabel('Frequenz [Bins]'); ylabel('Betrag [dB]');
title('Fensterfunktionen - Frequenzgang');
legend(windows, 'Location','northeast', 'Interpreter','none');
grid on; box on; xlim([0 10]); ylim([-120 5]); hold off;
end

% =========================================================================
function plot_overlaid_spectra(results, logScale)
%PLOT_OVERLAID_SPECTRA  Plot the HRV spectra of all windows on top of each other.
%   PLOT_OVERLAID_SPECTRA(RESULTS, LOGSCALE) opens one figure with the HRV band
%   limits marked, so the windows can be compared directly.
%
%   Input:
%     RESULTS  - struct array from COMPARE_WINDOW_FUNCTIONS.
%     LOGSCALE - false for a linear PSD axis, true for a dB axis; the dB view
%                exposes the leakage floor, the linear view the band powers.

colors = lines(numel(results)); % eine Farbe pro Fenster
figure('Color','w'); hold on;
% Zoom the linear view to just past the HF band (0.40 Hz) -> the empty range
% 0.45-0.5 Hz is dropped. The dB view shows the whole leakage floor.
if logScale, fmax = 0.5; else, fmax = 0.45; end % obere Frequenzgrenze je nach Ansicht
for iw = 1:numel(results) % ueber jedes Fenster
    f = results(iw).f;  psd = results(iw).psd; % Frequenzachse und gemittelte PSD holen
    if logScale, y = 10*log10(psd + eps); else, y = psd; end % dB-Ansicht oder lineare PSD
    plot(f, y, 'LineWidth', 1.3, 'Color', colors(iw,:)); % Spektrum ueberlagert plotten
end
yl = ylim;                               % freeze the limits before drawing the band markers
for fb = [0.04 0.15 0.40] % senkrechte gestrichelte Linien an den HRV-Bandgrenzen
    plot([fb fb], yl, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off'); % HandleVisibility off = nicht in Legende
end
ylim(yl); xlim([0 fmax]); xlabel('Frequenz [Hz]');
if logScale % Achsenbeschriftung/Titel je nach Ansicht
    ylabel('Leistungsdichte [dB]');
    title('Ueberlagerte HRV-Spektren (logarithmisch) - Leakage-Vergleich');
else
    ylabel('Leistungsdichte [ms^2/Hz]');
    title('Ueberlagerte HRV-Spektren (linear)');
end
legend({results.name}, 'Location','northeast', 'Interpreter','none'); % Fensternamen als Legende
grid on; box on; hold off;
end

% =========================================================================
function print_comparison_table(results)
%PRINT_COMPARISON_TABLE  Print the window and HRV measures as a table in the command window.
%   PRINT_COMPARISON_TABLE(RESULTS) lists PSL, MLW and ENBW next to the band
%   powers of every window, so the numbers behind the plots can be read off.
%
%   Input:
%     RESULTS - struct array from COMPARE_WINDOW_FUNCTIONS.

fprintf('\n=================== Vergleich der Fensterfunktionen ===================\n');
fprintf('%-9s | %7s %6s %6s | %8s %8s %7s | %6s\n', ... % Kopfzeile mit Spaltennamen
    'Fenster','PSL','MLW','ENBW','LF','HF','LF/HF','VLF');
fprintf('%-9s | %7s %6s %6s | %8s %8s %7s | %6s\n', ... % zweite Kopfzeile mit Einheiten
    '','[dB]','[Bin]','[Bin]','[ms^2]','[ms^2]','[-]','[ms^2]');
fprintf('----------------------------------------------------------------------\n');
for iw = 1:numel(results) % ueber jedes Fenster eine Zeile drucken
    m = results(iw).metrics;  h = results(iw).hrv; % Guetemasse und HRV-Werte holen
    fprintf('%-9s | %7.1f %6.2f %6.2f | %8.1f %8.1f %7.2f | %6.1f\n', ... % Zahlen formatiert ausgeben
        results(iw).name, m.PSL, m.MLW, m.ENBW, h.LF, h.HF, h.LF_HF, h.VLF);
end
fprintf('======================================================================\n');
end

% =========================================================================
function plot_spectrum(f, psd, window_name, hrv)
%PLOT_SPECTRUM  Plot one HRV power spectral density with the VLF/LF/HF bands shaded.
%   PLOT_SPECTRUM(F, PSD, WINDOW_NAME, HRV) opens one figure for a single
%   window function.
%
%   Input:
%     F           - frequency axis in Hz.
%     PSD         - power spectral density belonging to F, in ms^2/Hz.
%     WINDOW_NAME - window name for the title (optional).
%     HRV         - band measures (optional); if given, the LF/HF values are
%                   shown in a text box inside the plot.

if nargin < 3 || isempty(window_name), window_name = ''; end % Standard: leerer Fenstername
figure('Color','w'); hold on;

fmax = 0.45;    % up to just past the HF band (0.40 Hz); the empty range 0.45-0.5 Hz is dropped
ymax = 1.05 * max(psd(f <= fmax)); % obere y-Grenze = 5% ueber dem groessten sichtbaren PSD-Wert
if ~isfinite(ymax) || ymax <= 0, ymax = 1; end % Sicherheit falls PSD leer/0

shade_hrv_bands(ymax); % VLF/LF/HF-Baender farbig hinterlegen
plot(f, psd, 'k', 'LineWidth', 1.3); % eigentliche PSD-Kurve in schwarz

xlim([0 fmax]); ylim([0 ymax]);
xlabel('Frequenz [Hz]'); ylabel('Leistungsdichte [ms^2/Hz]');
title(sprintf('HRV-Leistungsdichtespektrum - Fenster: %s', window_name));
grid on; box on;
legend({'VLF (0.0033-0.04 Hz)','LF (0.04-0.15 Hz)','HF (0.15-0.40 Hz)','PSD'}, ...
    'Location','northeast');

if nargin >= 4 && ~isempty(hrv) % nur wenn HRV-Werte uebergeben wurden
    txt = sprintf('LF/HF = %.2f\nLF = %.1f ms^2 (%.0f%%)\nHF = %.1f ms^2 (%.0f%%)', ... % Textbox-Inhalt zusammenbauen
        hrv.LF_HF, hrv.LF, hrv.LF_rel, hrv.HF, hrv.HF_rel);
    % Place the text box at half height: the legend sits at the top right
    % (northeast), and above the PSD the right half is empty -> no overlap.
    text(0.97*fmax, 0.55*ymax, txt, 'HorizontalAlignment','right', ... % Textbox auf halber Hoehe rechts platzieren
        'VerticalAlignment','top', 'BackgroundColor','w', ...
        'EdgeColor',[0.6 0.6 0.6], 'FontSize',9);
end
hold off;
end

% =========================================================================
function shade_hrv_bands(ymax)
%SHADE_HRV_BANDS  Shade the VLF/LF/HF bands of the current axes in colour.
%   SHADE_HRV_BANDS(YMAX) draws one patch per band across the full plot height.
%
%   Input:
%     YMAX - upper edge of the patches, i.e. the y limit of the axes.

bands = [0.0033 0.04; 0.04 0.15; 0.15 0.40]; % Frequenzgrenzen der drei HRV-Baender [Hz]
cols  = [0.85 0.90 0.98; 0.80 0.95 0.85; 0.99 0.88 0.82]; % je eine RGB-Farbe pro Band
for i = 1:size(bands,1) % ueber jedes Band
    x = bands(i,:); % untere/obere Frequenzgrenze dieses Bandes
    patch([x(1) x(2) x(2) x(1)], [0 0 ymax ymax], cols(i,:), ... % farbiges Rechteck ueber volle Hoehe
        'EdgeColor','none', 'FaceAlpha',0.7); % keine Umrandung, halbtransparent
end
end

% =========================================================================
function plot_waterfall(f, seg_psd, t_seg, window_name, zmax_common)
%PLOT_WATERFALL  Plot the segment spectra as a 3D waterfall with the HRV band limits.
%   PLOT_WATERFALL(F, SEG_PSD, T_SEG, WINDOW_NAME, ZMAX_COMMON) draws one line
%   per segment, which makes the stability of the spectrum over time visible.
%
%   Input:
%     F           - frequency axis in Hz.
%     SEG_PSD     - PSD of every segment, one column per segment.
%     T_SEG       - start time of every segment in s (optional, default 1:nSeg).
%     WINDOW_NAME - window name for the title (optional).
%     ZMAX_COMMON - upper limit of the z and colour axis (optional); pass the
%                   same value for every window to make the plots comparable.

if nargin < 4 || isempty(window_name), window_name = ''; end % Standard: leerer Fenstername
f   = f(:); % Spaltenvektor erzwingen
sel = f <= 0.5; % nur Frequenzen bis 0.5 Hz zeigen (darueber ist nichts Relevantes)
fz  = f(sel); % gefilterte Frequenzachse
Z   = seg_psd(sel, :).';                 % rows = segments
nSeg = size(Z, 1); % Anzahl der Segmente
if nargin < 3 || isempty(t_seg), t_seg = 1:nSeg; end % Standard: Segmentnummer statt Zeit
t_seg = t_seg(:); % Spaltenvektor erzwingen

figure('Color','w');
waterfall(fz, t_seg, Z); % 3D-Wasserfall: eine Linie pro Segment
colormap(turbo); % Farbverlauf fuer die Hoehe
xlabel('Frequenz [Hz]'); ylabel('Segment-Startzeit [s]');
zlabel('Leistungsdichte [ms^2/Hz]');
title(sprintf('HRV-Wasserfalldiagramm - Fenster: %s', window_name));
grid on; box on; view(40, 30); % Blickwinkel (Azimut 40, Elevation 30)

% Common z maximum: either passed in or taken from this diagram.
if nargin >= 5 && ~isempty(zmax_common) && isfinite(zmax_common) && zmax_common > 0 % gemeinsames z-Maximum uebergeben?
    zmax = zmax_common; % dann dieses verwenden --> alle Wasserfaelle vergleichbar
    zlim([0 zmax]); clim([0 zmax]);      % uniform z and colour scale (clim = current replacement for caxis)
else
    zmax = max(Z(:)); % sonst Maximum aus diesem Diagramm
end

% Mark each HRV band limit as an upright, semi-transparent rectangle
% (rectangle in the plane x = fb, spanned over the time and z axis).
hold on;
if ~isfinite(zmax) || zmax <= 0, zmax = 1; end % Sicherheit gegen ungueltiges zmax
ymin = min(t_seg); ymax = max(t_seg); if ymax == ymin, ymax = ymin+1; end % y-Bereich (Zeitachse) der Flaechen
for fb = [0.04 0.15 0.40] % je Bandgrenze eine aufrechte Flaeche
    patch([fb fb fb fb], [ymin ymin ymax ymax], [0 zmax zmax 0], ... % Rechteck in der Ebene x = fb
        [0.4 0.4 0.4], 'FaceAlpha',0.12, 'EdgeColor','none'); % grau, sehr transparent, ohne Rand
end
hold off;
end
