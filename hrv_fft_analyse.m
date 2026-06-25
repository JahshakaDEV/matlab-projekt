%% ========================================================================
%  HRV-ANALYSE: EINFLUSS UNTERSCHIEDLICHER FFT-FENSTERFUNKTIONEN
%  ------------------------------------------------------------------------
%  Modularer Aufbau: alle Teilschritte sind als eigene (lokale) Funktionen
%  am Ende der Datei umgesetzt (moeglich ab MATLAB R2016b). Einfach diese
%  Datei oeffnen und ausfuehren.
%
%  Verarbeitungskette:
%    load_ecg_data -> preprocess_ecg -> detect_r_peaks ->
%    calculate_rr_intervals -> interpolate_rr_signal ->
%    compare_window_functions ( apply_window_function, calculate_fft,
%    calculate_hrv_bands ) -> plot_spectrum / plot_waterfall
%
%  Bedienung:
%    * EDF_FILE = '' : es wird ein synthetisches 1-h-Test-EKG mit bekannter
%      HRV erzeugt (LF 0.10 Hz, HF 0.25 Hz) -> Pipeline pruefbar.
%    * EDF_FILE = 'pfad.edf' : eigene Aufzeichnung (EDF-Import ab R2020b).
%
%  Benoetigt: Signal Processing Toolbox. Der 50-Hz-Notch ist bewusst ohne
%  iirnotch implementiert (manueller Biquad).
%  ========================================================================

clear; clc; close all;

%% ----------------------- Parameter -------------------------------------
EDF_FILE   = '11-18-37.EDF';      % '' => synthetische Testdaten
FS_INTERP  = 4;       % Resampling-Frequenz des RR-Signals [Hz]
SEG_LEN_S  = 300;     % Analysefensterlaenge [s] = 5 min
OVERLAP    = 0.5;     % Segment-Ueberlappung
WINDOWS    = {'rect','hann','hamming','blackman','kaiser','flattop'};

%% --------------------------- 1. EKG laden ------------------------------
fprintf('=== 1. EKG-Daten laden ===\n');
[ecg_raw, fs] = load_ecg_data(EDF_FILE);
t = (0:numel(ecg_raw)-1).' / fs;

%% ------------------------ 2. Vorverarbeitung ---------------------------
fprintf('\n=== 2. Vorverarbeitung ===\n');
ecg = preprocess_ecg(ecg_raw, fs);

% Roh- vs. gefiltertes Signal (erste 10 s)
figure('Color','w','Position',[100 100 950 500]);
sel = t <= 360;
subplot(2,1,1);
plot(t(sel), ecg_raw(sel), 'Color', [0.6 0.6 0.6]);
title('EKG - Rohsignal (erste 10 s)'); ylabel('Amplitude');
grid on; box on; xlim([190 200]);
subplot(2,1,2);
plot(t(sel), ecg(sel), 'Color', [0.85 0.2 0.2]);
title('EKG - nach Vorverarbeitung'); xlabel('Zeit [s]'); ylabel('Amplitude');
grid on; box on; xlim([190 200]);

%% --------------------------- 3. R-Zacken -------------------------------
fprintf('\n=== 3. R-Zacken-Erkennung ===\n');
r_locs = detect_r_peaks(ecg, fs);

figure('Color','w','Position',[100 100 950 350]);
plot(t(sel), ecg(sel), 'k'); hold on;
r_sel = r_locs(r_locs <= 10*fs);
plot(r_sel/fs, ecg(r_sel), 'ro', 'MarkerFaceColor','r', 'MarkerSize',5);
title('Erkannte R-Zacken (erste 10 s)');
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
fprintf('\n=== 5. Interpolation (Resampling) ===\n');
[~, sig, fs_i] = interpolate_rr_signal(t_rr, rr, FS_INTERP);

%% ----------------------- 6. Fenstervergleich ---------------------------
fprintf('\n=== 6. Vergleich der Fensterfunktionen ===\n');
results = compare_window_functions(sig, fs_i, WINDOWS, SEG_LEN_S, OVERLAP);

%% --------------------- 7. Einzeldarstellungen --------------------------
fprintf('\n=== 7. Einzelspektren und Wasserfalldiagramme ===\n');
for iw = 1:numel(results)
    R = results(iw);
    plot_spectrum(R.f, R.psd, R.name, R.hrv);
    plot_waterfall(R.f, R.seg_psd, R.t_seg, R.name);
end

fprintf('\n=== Analyse abgeschlossen ===\n');


%% ========================================================================
%  LOKALE FUNKTIONEN
%  ========================================================================

% =========================================================================
function [ecg, fs] = load_ecg_data(filename)
%LOAD_ECG_DATA  Liest ein EKG aus einer EDF-Datei oder erzeugt Testdaten.
%   Ohne (gueltige) Datei wird auf ein synthetisches Test-EKG
%   zurueckgegriffen, damit das Programm immer lauffaehig ist.
%   EDF-Import benoetigt MATLAB R2020b+.

    if nargin < 1 || isempty(filename) || exist(filename, 'file') ~= 2
        if nargin >= 1 && ~isempty(filename)
            warning('Datei "%s" nicht gefunden -> synthetische Testdaten.', filename);
        else
            fprintf('Kein EDF-Dateiname uebergeben -> synthetische Testdaten.\n');
        end
        [ecg, fs] = generate_synthetic_ecg(3600, 250);   % 1 h @ 250 Hz
        return;
    end

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
function ecg_filt = preprocess_ecg(ecg, fs, mains_freq)
%PREPROCESS_ECG  Filterung: DC, Baseline-Hochpass, 50-Hz-Notch, Tiefpass.
%   Alle Filter nullphasig (filtfilt) -> keine Phasenverschiebung, die
%   R-Zacken werden zeitlich nicht verschoben.

    if nargin < 3 || isempty(mains_freq), mains_freq = 50; end
    ecg = double(ecg(:));
    nyq = fs / 2;

    ecg_filt = ecg - mean(ecg);                       % 1) DC entfernen

    [b,a]    = butter(2, 0.5/nyq, 'high');            % 2) Baseline-Hochpass
    ecg_filt = filtfilt(b, a, ecg_filt);

    [bn,an]  = notch_coeffs(mains_freq, fs, 30);      % 3) 50-Hz-Notch
    ecg_filt = filtfilt(bn, an, ecg_filt);

    [b,a]    = butter(4, 40/nyq, 'low');              % 4) HF-Rauschen-Tiefpass
    ecg_filt = filtfilt(b, a, ecg_filt);
end

% =========================================================================
function [b, a] = notch_coeffs(f0, fs, Q)
%NOTCH_COEFFS  RBJ-Bandstop-Biquad (Kerbfilter), ohne iirnotch.
    w0    = 2*pi*f0/fs;
    alpha = sin(w0) / (2*Q);
    b = [ 1,         -2*cos(w0),  1         ];
    a = [ 1 + alpha, -2*cos(w0),  1 - alpha ];
    b = b / a(1);  a = a / a(1);
end

% =========================================================================
function r_locs = detect_r_peaks(ecg, fs)
%DETECT_R_PEAKS  R-Zacken-Detektion (vereinfachter Pan-Tompkins).
%   Bandpass 5-15 Hz -> Ableitung -> Quadrierung -> gleitende Mittelung
%   -> Schwellwert-Peaksuche -> Feinjustierung auf das EKG-Maximum.

    ecg = double(ecg(:));
    nyq = fs / 2;

    [b,a] = butter(2, [5/nyq, 15/nyq], 'bandpass');   % QRS-Band betonen
    f  = filtfilt(b, a, ecg);
    d  = [f(1); diff(f)];                             % Ableitung
    sq = d .^ 2;                                      % Quadrierung
    w     = max(1, round(0.150*fs));
    integ = conv(sq, ones(w,1)/w, 'same');           % Energie-Einhuellende

    thr     = 0.4 * mean(integ);
    minDist = round(0.25 * fs);                       % max. ~240 bpm
    [~, pk] = findpeaks(integ, 'MinPeakHeight', thr, 'MinPeakDistance', minDist);

    % Feinjustierung: echtes R-Maximum im Fenster +/- 50 ms suchen
    wref   = round(0.05 * fs);
    r_locs = zeros(numel(pk), 1);
    for i = 1:numel(pk)
        lo = max(1, pk(i)-wref);  hi = min(numel(ecg), pk(i)+wref);
        [~, rel] = max(ecg(lo:hi));
        r_locs(i) = lo + rel - 1;
    end
    r_locs = unique(r_locs);
    fprintf('R-Zacken erkannt: %d\n', numel(r_locs));
end

% =========================================================================
function [t_rr, rr, n_artifacts] = calculate_rr_intervals(r_locs, fs)
%CALCULATE_RR_INTERVALS  RR-Intervalle [ms] mit Plausibilitaetspruefung.
%   Plausibel: 300-2000 ms und |RR-Median| < 50 %. Unplausible Werte
%   werden per pchip aus den gueltigen Nachbarn ersetzt.

    r_locs = r_locs(:);
    t_r  = r_locs / fs;
    rr   = diff(t_r) * 1000;       % RR [ms]
    t_rr = t_r(2:end);

    ok  = (rr > 300) & (rr < 2000);
    med = median(rr(ok));
    ok  = ok & (abs(rr - med) < 0.5*med);
    bad = ~ok;
    n_artifacts = sum(bad);

    if n_artifacts > 0 && sum(ok) >= 2
        rr(bad) = interp1(t_rr(ok), rr(ok), t_rr(bad), 'pchip', 'extrap');
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

        seg_psd = [];
        for is = 1:nSeg
            s0  = starts(is);
            seg = sig(s0:s0+N-1);
            seg = seg - mean(seg);              % lokalen DC abziehen
            [f, psd] = calculate_fft(seg.*win, win, fs_i);
            if isempty(seg_psd), seg_psd = zeros(numel(f), nSeg); end
            seg_psd(:, is) = psd;
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
%   PSL = hoechste Nebenkeule (Leakage), MLW = Hauptkeulenbreite
%   (Aufloesung), ENBW = aequivalente Rauschbandbreite.

    win = win(:);  N = numel(win);
    m.ENBW = enbw(win);

    Nfft = 65536;
    W    = abs(fft(win, Nfft));
    W    = W(1:floor(Nfft/2)+1);
    W    = W / W(1);
    WdB  = 20*log10(W + eps);
    fbin = (0:numel(W)-1).' * N / Nfft;

    % erste echte Nullstelle = erstes lokales Minimum (>0.2 Bins von DC);
    % robust auch beim Flat-Top (welliger flacher Scheitel)
    nullidx = [];
    for i = 2:numel(W)-1
        if W(i) < W(i-1) && W(i) < W(i+1) && fbin(i) > 0.2
            nullidx = i; break;
        end
    end
    if ~isempty(nullidx)
        m.MLW = 2 * fbin(nullidx);
        m.PSL = max(WdB(nullidx:end));
    else
        m.MLW = NaN;  m.PSL = NaN;
    end
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
        text(0.97*fmax, 0.95*ymax, txt, 'HorizontalAlignment','right', ...
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
function plot_waterfall(f, seg_psd, t_seg, window_name)
%PLOT_WATERFALL  3D-Wasserfall der Segmentspektren mit HRV-Bandgrenzen.
%   Jede Linie = Spektrum eines Segments -> zeitliche Stabilitaet sichtbar.

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

    % Bandgrenzen als senkrechte Ebenen (je 2 Dreiecke -> MATLAB+Octave)
    hold on;
    zmax = max(Z(:)); if ~isfinite(zmax) || zmax <= 0, zmax = 1; end
    ymin = min(t_seg); ymax = max(t_seg); if ymax == ymin, ymax = ymin+1; end
    for fb = [0.04 0.15 0.40]
        patch([fb fb fb], [ymin ymax ymax], [0 0 zmax], [0.3 0.3 0.3], ...
              'FaceAlpha',0.12, 'EdgeColor',[0.3 0.3 0.3], 'LineStyle','--');
        patch([fb fb fb], [ymin ymax ymin], [0 zmax zmax], [0.3 0.3 0.3], ...
              'FaceAlpha',0.12, 'EdgeColor','none');
    end
    hold off;
end
