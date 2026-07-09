%% ========================================================================
%  HRV-ANALYSE: EINFLUSS UNTERSCHIEDLICHER FFT-FENSTERFUNKTIONEN
%  ------------------------------------------------------------------------
%  Modularer Aufbau: alle Teilschritte sind als eigene (lokale) Funktionen
%  am Ende der Datei umgesetzt.
%
%  Verarbeitungskette:
%    load_ecg_data -> preprocess_ecg -> detect_r_peaks ->
%    calculate_rr_intervals -> interpolate_rr_signal ->
%    compare_window_functions ( apply_window_function, calculate_fft,
%    calculate_hrv_bands ) -> plot_spectrum / plot_waterfall
%
%
%  Benoetigt: Signal Processing Toolbox.
%  ========================================================================

clear; clc; close all;

%% ----------------------- Parameter -------------------------------------
EDF_FILE   = '16-49-27.EDF';
FS_INTERP  = 4;       % Resampling-Frequenz des RR-Signals [Hz]
SEG_LEN_S  = 300;     % Analysefensterlaenge [s] = 5 min
OVERLAP    = 0.5;     % Segment-Ueberlappung
WINDOWS    = {'rect','hann','hamming','blackman','kaiser','flattop'};

%% --------------------------- 1. EKG laden ------------------------------
% Aufgabe 2.1

fprintf('=== 1. EKG-Daten laden ===\n');
[ecg_raw, fs] = load_ecg_data(EDF_FILE);
t = (0:numel(ecg_raw)-1).' / fs;

%% ------------------------ 2. Vorverarbeitung ---------------------------
fprintf('\n=== 2. Vorverarbeitung ===\n');
ecg = preprocess_ecg(ecg_raw, fs);

% Roh- vs. gefiltertes Signal (Ausschnitt 10-20 s)
figure('Color','w','Position',[100 100 950 500]);
sel = t >= 0;                                    % gesamtes Signal auswaehlen (xlim unten zoomt nur die Ansicht)
subplot(2,1,1);
plot(t(sel), ecg_raw(sel), 'Color', [0.6 0.6 0.6]);
title('EKG - Rohsignal (Ausschnitt 10-20 s)'); ylabel('Amplitude');
grid on; box on; xlim([10 20]);
subplot(2,1,2);
plot(t(sel), ecg(sel), 'Color', [0.85 0.2 0.2]);
title('EKG - nach Vorverarbeitung (Ausschnitt 10-20 s)'); xlabel('Zeit [s]'); ylabel('Amplitude');
grid on; box on; xlim([10 20]);

%% --------------------------- 3. R-Zacken -------------------------------
% Aufgabe 2.2

fprintf('\n=== 3. R-Zacken-Erkennung ===\n');
r_locs = detect_r_peaks(ecg, fs);

figure('Color','w','Position',[100 100 950 350]);
plot(t(sel), ecg(sel), 'k'); hold on;
r_sel = r_locs;                                  % alle erkannten R-Zacken (xlim unten begrenzt nur die sichtbare Ansicht)
plot(r_sel/fs, ecg(r_sel), 'ro', 'MarkerFaceColor','r', 'MarkerSize',5);
title('Erkannte R-Zacken');
xlabel('Zeit [s]'); ylabel('Amplitude');
legend({'EKG','R-Zacken'}, 'Location','northeast');
grid on; box on; xlim([0 10]); hold off;

%% ------------------------- 4. RR-Intervalle ----------------------------
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

%% ----------------------- 6. Fenstervergleich ---------------------------
% Aufgabe 2.4
fprintf('\n=== 6. Vergleich der Fensterfunktionen ===\n');
results = compare_window_functions(sig, fs_i, WINDOWS, SEG_LEN_S, OVERLAP);

%% --------------------- 7. Einzeldarstellungen --------------------------
fprintf('\n=== 7. Einzelspektren und Wasserfalldiagramme ===\n');

% Gemeinsames z-Maximum ueber alle Fenster (Anforderung 2.7: einheitliche
% Skalierung -> alle Wasserfalldiagramme werden vergleichbar).
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

% Alle offenen Figures als PDF exportieren (eine Figure pro Seite).
pdfName = 'all_figures.pdf';
if isfile(pdfName), delete(pdfName); end          % alten Export entfernen (sonst wird angehaengt)

figHandles = findall(0, 'Type', 'figure');
[~, order] = sort([figHandles.Number]);           % logische Reihenfolge (EKG zuerst)
figHandles = figHandles(order);
for i = 1:numel(figHandles)
    exportgraphics(figHandles(i), pdfName, 'Append', true);
end

fprintf('\n=== Analyse abgeschlossen ===\n');


%% ========================================================================
%  LOKALE FUNKTIONEN
%  ========================================================================

% =========================================================================
function [ecg, fs] = load_ecg_data(filename)
%LOAD_ECG_DATA  Liest ein EKG aus einer EDF-Datei.

% Abtastrate aus dem EDF-Header
info   = edfinfo(filename);
recDur = seconds(info.DataRecordDuration);
fs     = double(info.NumSamples(1)) / recDur;

% Signaldaten (erster Kanal) zu einem durchgehenden Vektor zusammenfuegen
tt    = edfread(filename);
vname = tt.Properties.VariableNames{1};
col   = tt.(vname);
if iscell(col), ecg = cell2mat(col); else, ecg = col(:); end
ecg = double(ecg(:));

fprintf('EDF geladen: "%s" | Kanal "%s" | fs = %.1f Hz | Dauer = %.1f min\n', ...
    filename, vname, fs, numel(ecg)/fs/60);
end

% =========================================================================
function ecg_filt = preprocess_ecg(ecg, fs)
%PREPROCESS_ECG  Filterung: DC, Baseline-Hochpass, 50-Hz-Notch, Tiefpass.
%   Alle Filter nullphasig (filtfilt) -> keine Phasenverschiebung, die
%   R-Zacken werden zeitlich nicht verschoben.

ecg_filt = ecg - mean(ecg);                       % Gleichanteil entfernen

ecg_filt = highpass(ecg_filt,0.5, fs);            % High-pass für Entfernen der Baseline-Wanderung

wo = 50/(fs/2);                                   % 50 ist die Frequenz
bw = wo/40;                                       % 40 ist die Güte des Filters
[b,a] = iirnotch(wo, bw);
ecg_filt = filtfilt(b, a, ecg_filt);              % 50 Hz Netzbrummen entfernen

ecg_filt = lowpass(ecg_filt, 40, fs);             % Tiefpass Filter ab 40Hz
end

% =========================================================================
function locs = detect_r_peaks(ecg, fs)
%DETECT_R_PEAKS  Automatische R-Zacken-Erkennung mit findpeaks.
%   Zwei Kriterien: Mindesthoehe (nur deutlich herausragende Spitzen zaehlen
%   als R-Zacke) und Mindestabstand (physiologisch kein zweiter Schlag < 300 ms).

schwelle       = mean(ecg) + 0.6*std(ecg);        % Schwelle: nur Peaks deutlich ueber dem Mittelwert
RR_min_samples = round(0.3*fs);                   % 0.3 s Mindestabstand -> schliesst Doppeldetektion aus (max. ~200 bpm)
[~, locs] = findpeaks(ecg, ...                    % Rueckgabe: Sample-Indizes der Peaks (Amplituden ignoriert)
    'MinPeakHeight', schwelle, 'MinPeakDistance', RR_min_samples);
fprintf('R-Zacken erkannt: %d\n', numel(locs));
end

% =========================================================================
function [t_rr, rr, n_artifacts] = calculate_rr_intervals(r_locs, fs)
%CALCULATE_RR_INTERVALS  RR-Intervalle [ms] mit Plausibilitaetspruefung.
%   Als Artefakt gilt: ausserhalb 300-2000 ms ODER statistischer Ausreisser
%   gegenueber dem gleitenden Median (isoutlier). Solche Werte werden per
%   pchip aus den gueltigen Nachbarn ersetzt.

r_locs = r_locs(:);
t_r  = r_locs / fs;
rr   = diff(t_r) * 1000;       % RR [ms]
t_rr = t_r(2:end);

% Plausibilitaetspruefung (beide Schritte mit nativen Funktionen):
%  1) physiologisch unmoegliche Werte ausserhalb 300-2000 ms (~30-200 bpm),
%  2) statistische Ausreisser gegenueber dem gleitenden Median
%     (isoutlier, Methode 'movmedian' ueber 21 Intervalle). Die MAD-basierte
%     Pruefung passt sich der lokalen Streuung an -> faengt einzelne
%     verpasste/falsche Schlaege auch in schwankenden Abschnitten.
implausible = (rr < 300) | (rr > 2000);
outlier     = isoutlier(rr, 'movmedian', 21);
bad = implausible | outlier;
n_artifacts = sum(bad);

if n_artifacts > 0 && sum(~bad) >= 2
    rr(bad) = interp1(t_rr(~bad), rr(~bad), t_rr(bad), 'pchip', 'extrap');
end

fprintf(['RR-Intervalle: %d | Artefakte: %d | mittl. RR = %.1f ms ', ...
    '(%.1f bpm) | SDNN = %.1f ms\n'], numel(rr), n_artifacts, ...
    mean(rr), 60000/mean(rr), std(rr));
end

% =========================================================================
function [t_u, sig, fs_i] = interpolate_rr_signal(t_rr, rr, fs_i)
%INTERPOLATE_RR_SIGNAL  RR-Reihe auf gleichmaessiges Raster (Spline).
%   Die FFT setzt gleichmaessige Abtastung voraus; die RR-Werte liegen aber
%   zu unregelmaessigen Zeitpunkten vor. Danach lineare Trendbereinigung.

if nargin < 3 || isempty(fs_i), fs_i = 4; end
t_rr = t_rr(:);  rr = rr(:);

t_u = (t_rr(1) : 1/fs_i : t_rr(end)).';
sig = interp1(t_rr, rr, t_u, 'spline');
sig = detrend(sig, 1);

fprintf('Interpoliert: %d Werte @ %.0f Hz | Dauer = %.0f s\n', ...
    numel(sig), fs_i, t_u(end)-t_u(1));
end

% =========================================================================
function [xw, win] = apply_window_function(x, window_type)
%APPLY_WINDOW_FUNCTION  Fensterung eines Segments.
%   Gibt das gefensterte Signal XW UND den Fenstervektor WIN zurueck (WIN
%   wird zur Leistungsnormierung in calculate_fft gebraucht). Fenster
%   daempfen die Segmentraender und reduzieren so das spektrale Leakage,
%   verbreitern aber die Hauptkeule (geringere Aufloesung).
%   Reinen Fenstervektor: [~,w] = apply_window_function(ones(N,1),'hann');

x = double(x(:));
N = numel(x);
switch lower(window_type)
    case 'rect',     win = ones(N,1);     % keine Fensterung
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
%CALCULATE_FFT  Einseitige, leistungsnormierte PSD.
%   PSD = |FFT(xw)|^2 / (fs*U) mit U = sum(win.^2). Die Normierung auf die
%   Fensterleistung U macht die integrierte (Band-)Leistung weitgehend
%   fensterunabhaengig. Innere Bins werden verdoppelt (einseitiges Spektrum
%   eines reellen Signals).

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
%CALCULATE_HRV_BANDS  Leistung in VLF/LF/HF + LF/HF, relativ, n.u.
%   Baender (Task Force 1996): VLF 0.0033-0.04, LF 0.04-0.15, HF 0.15-0.40 Hz.
%   Bandleistung = Flaeche unter der PSD (Trapezregel).

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

lf_hf = hrv.LF + hrv.HF;                  % normalisierte Einheiten
if lf_hf > 0
    hrv.LF_nu = 100*hrv.LF/lf_hf;  hrv.HF_nu = 100*hrv.HF/lf_hf;
else
    hrv.LF_nu = 0; hrv.HF_nu = 0;
end

if hrv.HF > 0, hrv.LF_HF = hrv.LF/hrv.HF; else, hrv.LF_HF = NaN; end
end

% =========================================================================
function results = compare_window_functions(sig, fs_i, windows, seg_len_s, overlap)
%COMPARE_WINDOW_FUNCTIONS  Welch-Spektren je Fenster + Vergleich.
%   Das Signal wird in Segmente (Laenge seg_len_s, Ueberlappung overlap)
%   zerlegt; je Segment wird gefenstert und eine PSD gebildet, die PSDs
%   werden gemittelt (Welch). Das senkt die Varianz der Schaetzung.
%   Rueckgabe: Struktur-Array je Fenster (.name .f .psd .seg_psd .t_seg
%   .hrv .metrics). Zusaetzlich Vergleichsplots + Konsolentabelle.

if nargin < 3 || isempty(windows)
    windows = {'rect','hann','hamming','blackman','kaiser','flattop'};
end
if nargin < 4 || isempty(seg_len_s), seg_len_s = 300; end
if nargin < 5 || isempty(overlap),   overlap   = 0.5; end

sig = double(sig(:));
N   = round(seg_len_s * fs_i);
if N > numel(sig)
    warning('Signal kuerzer als ein Segment -> ganzes Signal als ein Segment.');
    N = numel(sig);
end
step   = max(1, round(N*(1-overlap)));
starts = 1:step:(numel(sig)-N+1);
if isempty(starts), starts = 1; end
t_seg  = (starts-1) / fs_i;
nSeg   = numel(starts);

fprintf('\nFenstervergleich: %d Segmente a %.0f s (%.0f%% Ueberlappung), %d Fenster\n', ...
    nSeg, seg_len_s, overlap*100, numel(windows));

results = struct('name',{}, 'f',{}, 'psd',{}, 'seg_psd',{}, ...
    't_seg',{}, 'hrv',{}, 'metrics',{});
for iw = 1:numel(windows)
    wname = windows{iw};
    [~, win] = apply_window_function(ones(N,1), wname);

    seg_psd = zeros(floor(N/2)+1, nSeg);    % Spektren aller Segmente vorab reservieren (FFT-Laenge = floor(N/2)+1)
    for is = 1:nSeg
        s0  = starts(is);                   % Startindex des Segments
        seg = sig(s0:s0+N-1);               % Segment der Laenge N ausschneiden
        seg = seg - mean(seg);              % lokalen DC abziehen
        [f, psd] = calculate_fft(seg.*win, win, fs_i);
        seg_psd(:, is) = psd;               % PSD dieses Segments in Spalte is ablegen
    end

    psd_avg = mean(seg_psd, 2);             % Welch-Mittelung
    results(iw).name    = wname;
    results(iw).f       = f;
    results(iw).psd     = psd_avg;
    results(iw).seg_psd = seg_psd;
    results(iw).t_seg   = t_seg;
    results(iw).hrv     = calculate_hrv_bands(f, psd_avg);
    results(iw).metrics = window_metrics(win);
end

plot_window_characteristics(windows, N);
plot_overlaid_spectra(results, false);
plot_overlaid_spectra(results, true);
print_comparison_table(results);
end

% =========================================================================
function m = window_metrics(win)
%WINDOW_METRICS  Kenngroessen: PSL [dB], MLW [Bins], ENBW [Bins].
%   PSL = hoechste Nebenkeule (Leakage-Mass), MLW = Hauptkeulenbreite
%   (Aufloesungs-Mass), ENBW = aequivalente Rauschbandbreite. Alle drei
%   werden mit eingebauten Funktionen aus dem Frequenzgang bestimmt.

win  = win(:);  N = numel(win);
m.ENBW = enbw(win);                          % aequivalente Rauschbandbreite

% Frequenzgang des Fensters, fein aufgeloest durch Zero-Padding
Nfft = 8192;
W    = abs(fft(win, Nfft));
W    = W(1:Nfft/2+1);
fbin = (0:Nfft/2).' * N / Nfft;              % Frequenzachse in FFT-Bins
WdB  = 20*log10(W/max(W) + eps);             % normiert, Hauptpeak = 0 dB

% Erste Nullstelle = erstes lokales Minimum jenseits von 0.5 Bins (die
% Hauptkeule ist stets breiter; so stoert eine Welligkeit im flachen
% Flat-Top-Scheitel nicht).
nulls     = fbin(islocalmin(W));
nulls     = nulls(nulls > 0.5);
firstNull = nulls(1);

m.MLW = 2 * firstNull;                       % Hauptkeulenbreite = 2 x erste Nullstelle
m.PSL = max(WdB(fbin > firstNull));          % hoechste Nebenkeule (alles nach der Hauptkeule)
end

% =========================================================================
function plot_window_characteristics(windows, N)
%PLOT_WINDOW_CHARACTERISTICS  Fensterform (Zeit) und Frequenzgang (dB).
colors = lines(numel(windows));
figure('Color','w','Position',[100 100 1000 420]);

subplot(1,2,1); hold on;
for iw = 1:numel(windows)
    [~, win] = apply_window_function(ones(N,1), windows{iw});
    plot((0:N-1)/(N-1), win, 'LineWidth', 1.4, 'Color', colors(iw,:));
end
xlabel('normierte Zeit'); ylabel('Amplitude');
title('Fensterfunktionen - Zeitbereich');
legend(windows, 'Location','south', 'Interpreter','none');
grid on; box on; ylim([0 1.05]); hold off;

subplot(1,2,2); hold on;
Nfft = 4096;
for iw = 1:numel(windows)
    [~, win] = apply_window_function(ones(N,1), windows{iw});
    W   = abs(fft(win, Nfft)); W = W(1:Nfft/2+1);
    WdB = 20*log10(W/max(W) + eps);
    fb  = (0:Nfft/2) / Nfft * N;
    plot(fb, WdB, 'LineWidth', 1.3, 'Color', colors(iw,:));
end
xlabel('Frequenz [Bins]'); ylabel('Betrag [dB]');
title('Fensterfunktionen - Frequenzgang');
legend(windows, 'Location','northeast', 'Interpreter','none');
grid on; box on; xlim([0 10]); ylim([-120 5]); hold off;
end

% =========================================================================
function plot_overlaid_spectra(results, logScale)
%PLOT_OVERLAID_SPECTRA  Ueberlagerte HRV-Spektren (linear oder dB).
colors = lines(numel(results));
figure('Color','w'); hold on;
fmax = 0.5;
for iw = 1:numel(results)
    f = results(iw).f;  psd = results(iw).psd;
    if logScale, y = 10*log10(psd + eps); else, y = psd; end
    plot(f, y, 'LineWidth', 1.3, 'Color', colors(iw,:));
end
yl = ylim;                               % Bandgrenzen (plot-basiert)
for fb = [0.04 0.15 0.40]
    plot([fb fb], yl, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
end
ylim(yl); xlim([0 fmax]); xlabel('Frequenz [Hz]');
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
%PRINT_COMPARISON_TABLE  Vergleichstabelle im Befehlsfenster.
fprintf('\n=================== Vergleich der Fensterfunktionen ===================\n');
fprintf('%-9s | %7s %6s %6s | %8s %8s %7s | %6s\n', ...
    'Fenster','PSL','MLW','ENBW','LF','HF','LF/HF','VLF');
fprintf('%-9s | %7s %6s %6s | %8s %8s %7s | %6s\n', ...
    '','[dB]','[Bin]','[Bin]','[ms^2]','[ms^2]','[-]','[ms^2]');
fprintf('----------------------------------------------------------------------\n');
for iw = 1:numel(results)
    m = results(iw).metrics;  h = results(iw).hrv;
    fprintf('%-9s | %7.1f %6.2f %6.2f | %8.1f %8.1f %7.2f | %6.1f\n', ...
        results(iw).name, m.PSL, m.MLW, m.ENBW, h.LF, h.HF, h.LF_HF, h.VLF);
end
fprintf('======================================================================\n');
end

% =========================================================================
function plot_spectrum(f, psd, window_name, hrv)
%PLOT_SPECTRUM  HRV-Leistungsdichtespektrum mit schattierten Baendern.
%   Optionales hrv-Argument blendet LF/HF-Kennwerte als Textfeld ein.

if nargin < 3 || isempty(window_name), window_name = ''; end
figure('Color','w'); hold on;

fmax = 0.5;
ymax = 1.05 * max(psd(f <= fmax));
if ~isfinite(ymax) || ymax <= 0, ymax = 1; end

shade_hrv_bands(ymax);
plot(f, psd, 'k', 'LineWidth', 1.3);

xlim([0 fmax]); ylim([0 ymax]);
xlabel('Frequenz [Hz]'); ylabel('Leistungsdichte [ms^2/Hz]');
title(sprintf('HRV-Leistungsdichtespektrum - Fenster: %s', window_name));
grid on; box on;
legend({'VLF (0.0033-0.04 Hz)','LF (0.04-0.15 Hz)','HF (0.15-0.40 Hz)','PSD'}, ...
    'Location','northeast');

if nargin >= 4 && ~isempty(hrv)
    txt = sprintf('LF/HF = %.2f\nLF = %.1f ms^2 (%.0f%%)\nHF = %.1f ms^2 (%.0f%%)', ...
        hrv.LF_HF, hrv.LF, hrv.LF_rel, hrv.HF, hrv.HF_rel);
    % Textbox auf halbe Hoehe setzen: die Legende sitzt oben rechts (northeast),
    % oberhalb der PSD ist die rechte Haelfte leer -> keine Ueberlappung mehr.
    text(0.97*fmax, 0.55*ymax, txt, 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'BackgroundColor','w', ...
        'EdgeColor',[0.6 0.6 0.6], 'FontSize',9);
end
hold off;
end

% =========================================================================
function shade_hrv_bands(ymax)
%SHADE_HRV_BANDS  Hinterlegt VLF/LF/HF farblich.
bands = [0.0033 0.04; 0.04 0.15; 0.15 0.40];
cols  = [0.85 0.90 0.98; 0.80 0.95 0.85; 0.99 0.88 0.82];
for i = 1:size(bands,1)
    x = bands(i,:);
    patch([x(1) x(2) x(2) x(1)], [0 0 ymax ymax], cols(i,:), ...
        'EdgeColor','none', 'FaceAlpha',0.7);
end
end

% =========================================================================
function plot_waterfall(f, seg_psd, t_seg, window_name, zmax_common)
%PLOT_WATERFALL  3D-Wasserfall der Segmentspektren mit HRV-Bandgrenzen.
%   Jede Linie = Spektrum eines Segments -> zeitliche Stabilitaet sichtbar.
%   zmax_common (optional) setzt eine fuer alle Fenster einheitliche
%   z- und Farbskala (Anforderung 2.7 "einheitliche Skalierung").

if nargin < 4 || isempty(window_name), window_name = ''; end
f   = f(:);
sel = f <= 0.5;
fz  = f(sel);
Z   = seg_psd(sel, :).';                 % Zeilen = Segmente
nSeg = size(Z, 1);
if nargin < 3 || isempty(t_seg), t_seg = 1:nSeg; end
t_seg = t_seg(:);

figure('Color','w');
waterfall(fz, t_seg, Z);
colormap(turbo);
xlabel('Frequenz [Hz]'); ylabel('Segment-Startzeit [s]');
zlabel('Leistungsdichte [ms^2/Hz]');
title(sprintf('HRV-Wasserfalldiagramm - Fenster: %s', window_name));
grid on; box on; view(40, 30);

% Gemeinsames z-Maximum: entweder uebergeben oder aus diesem Diagramm.
if nargin >= 5 && ~isempty(zmax_common) && isfinite(zmax_common) && zmax_common > 0
    zmax = zmax_common;
    zlim([0 zmax]); clim([0 zmax]);      % einheitliche z- und Farbskala (clim = aktueller Ersatz fuer caxis)
else
    zmax = max(Z(:));
end

% HRV-Bandgrenzen je als senkrechtes, halbtransparentes Rechteck markieren
% (Rechteck in der Ebene x = fb, aufgespannt ueber Zeit- und z-Achse).
hold on;
if ~isfinite(zmax) || zmax <= 0, zmax = 1; end
ymin = min(t_seg); ymax = max(t_seg); if ymax == ymin, ymax = ymin+1; end
for fb = [0.04 0.15 0.40]
    patch([fb fb fb fb], [ymin ymin ymax ymax], [0 zmax zmax 0], ...
        [0.4 0.4 0.4], 'FaceAlpha',0.12, 'EdgeColor','none');
end
hold off;
end
