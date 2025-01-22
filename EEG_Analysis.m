% Clear Workspace and Load Data
clear; clc; close all;

%% Load EEG data
data = importdata('eeg_signal.txt');
time = data(:,1);
signal = data(:,2);
fs = 128;  % Sampling frequency

% Remove mean, centre signal
signal = signal - mean(signal);

% Plot Raw EEG Signal
figure('Name', 'Raw EEG Signal - Time Domain', 'NumberTitle', 'off');
plot(time, signal);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('EEG Signal');

%% Segment the data (10-second segments)
segment_length = fs * 10;  
num_segments = floor(length(signal) / segment_length);

% Define the segment indices to extract
selected_segments = [1, 2, 3, floor(num_segments/2)-1, floor(num_segments/2), floor(num_segments/2)+1, num_segments-2, num_segments-1, num_segments];

% Figure for spectrograms
figure('Name', 'EEG Signal Spectrograms', 'NumberTitle', 'off');
for i = 1:length(selected_segments)
    idx = selected_segments(i);
    segment = signal((idx-1)*segment_length+1 : idx*segment_length);

    % Compute spectrogram
    [s, f, t] = spectrogram(segment, fs);

    % Plot spectrogram as a subplot
    subplot(3, 3, i);  % 3x3 grid for selected segments
    imagesc(t, f, log(abs(s)));
    axis xy;
    title(['Segment ', num2str(idx)]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;
end
sgtitle('Spectrograms of Selected EEG Segments');

%% Split the segments into 3 separate groups (first, middle, last)
segment_groups = {selected_segments(1:3), selected_segments(4:6), selected_segments(7:9)};

% Create a figures for subplots
window = hamming(256);
noverlap = 128;
nfft = 512;
for fig_num = 1:3
    figure;
    sgtitle(['EEG Segments Group ', num2str(fig_num), ': Raw Signal and PSD']);

    for i = 1:length(segment_groups{fig_num})
        idx = segment_groups{fig_num}(i);
        segment_start = (idx-1)*segment_length + 1;
        segment_end = idx*segment_length;
        segment_time = time(segment_start:segment_end);
        segment = signal(segment_start:segment_end);

        % Compute Welch's power spectral density (PSD)
        window = hamming(256);
        noverlap = 128;
        nfft = 512;
        [pxx, f] = pwelch(segment, window, noverlap, nfft, fs);

        % Plot raw signal (left column)
        subplot(3, 2, (i-1)*2 + 1);
        plot(segment_time, segment);
        title(['Segment ', num2str(idx), ' - Time Domain']);
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;

        % Plot PSD (right column)
        subplot(3, 2, (i-1)*2 + 2);
        plot(f, 10*log10(pxx));
        title(['Segment ', num2str(idx), ' - PSD']);
        xlabel('Frequency (Hz)');
        ylabel('Power/Frequency (dB/Hz)');
        grid on;
    end
end
sgtitle('Power Spectral Density (PSD) of Selected EEG Segments')

%% Define frequency bands of interest (in Hz)
freq_bands = [0.5 4.0; 4.0 8.0; 8.0 13.0];
band_labels = {'0.5-4 Hz (Delta)', '4-8 Hz (Theta)', '8-13 Hz (Alpha)'};

% Selected segments indices
selected_segments = [1, 2, 3, num_segments/2-1, num_segments/2, num_segments/2+1, num_segments-2, num_segments-1, num_segments];

% Initialize storage for absolute and relative power per segment
num_segments_selected = length(selected_segments);
absolute_power = zeros(num_segments_selected, size(freq_bands, 1));
relative_power = zeros(num_segments_selected, size(freq_bands, 1));

% Compute absolute and relative power for each segment
for seg = 1:num_segments_selected
    idx = selected_segments(seg);
    segment_start = (idx-1) * segment_length + 1;
    segment_end = idx * segment_length;
    segment = signal(segment_start:segment_end);

    % Compute Welch's power spectral density (PSD) for the segment
    window = hamming(256);
    noverlap = 128;
    nfft = 512;
    [pxx, f] = pwelch(segment, window, noverlap, nfft, fs);

    % Frequency resolution for integration
    freq_res = f(2) - f(1);

    % Compute absolute band power using Simpson's rule (approximated with trapz)
    for i = 1:size(freq_bands, 1)
        band_idx = (f >= freq_bands(i, 1)) & (f <= freq_bands(i, 2));
        absolute_power(seg, i) = trapz(f(band_idx), pxx(band_idx));
    end

    % Compute total power for the segment
    total_power = trapz(f, pxx);

    % Calculate relative band power as a percentage
    relative_power(seg, :) = (absolute_power(seg, :) / total_power) * 100;
end

% Display the results for each segment in table format
for seg = 1:num_segments_selected
    disp(['Results for Segment ', num2str(selected_segments(seg))]);
    results_table = table(band_labels', absolute_power(seg, :)', relative_power(seg, :)', ...
        'VariableNames', {'Frequency_Range_Hz', 'Absolute_Band_Power', 'Relative_Band_Power_Percent'});
    disp(results_table);
    fprintf('\n');
end