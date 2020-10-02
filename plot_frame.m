function [jth, pks, i1, i2] = plot_frame(jth, pks, figpath, fpath, fname, hann_window, noverlap, NFFT, Fs, k, tt, xx, th, nn_a, manual_check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_frame.m
%
% Plots the spectogram and time series including detection numbers and
% saves the figure as JPG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole segment

% Create output folder if not existent
if ~isfolder(figpath)
    mkdir(figpath);
end
% Hide plots if no manual correction is selected to make it faster
if manual_check
    figure()
else
    figure('visible','off')
end
clf
splpath = split(fpath,filesep);

% Spectrogram plot
subplot(2,1,1)
title(strcat(splpath{end} ,filesep, fname), 'Interpreter', 'none')
spectrogram(xx, hann_window, noverlap, NFFT*2, Fs, 'yaxis')
set(gca, 'YLim', [90 180], 'XTick', [], 'XLabel', [])
caxis ([-30 30])
colorbar('off')
title(strcat(splpath{end} ,filesep, fname), 'Interpreter', 'none')

% Time series plot
subplot(2,1,2)
cnums = arrayfun(@num2str, k:k+length(jth)-1, 'UniformOutput', 0); % cell array of click numbers
plot(tt,xx,'b', ...
     tt([1 end]),th*nn_a*[1 1],'r--', ...
     tt([1 end]),nn_a*[1 1],'k--', ...
     tt(jth),pks,'Or') % Comment if plot gets crowded    
text(tt(jth), pks, cnums, 'FontSize', 8)

% Delete points manually
if manual_check
    [jth_rm, pks_rm] = remove_detections(tt, jth, pks, k, xx, th, nn_a);
    jth = jth_rm;
    pks = pks_rm;
    i1 = jth - 256;
    i2 = jth + 255;
end
if isempty(jth)
    close(gcf)
    return
end    
xlabel('Time [s]')
ylabel('Amplitude [counts]')
xticks = get(gca, 'Xtick');
xtlabel = cellstr(datestr(seconds(xticks), 'MM:SS.FFF'));
set(gca, 'XTick', xticks, 'XTickLabel',xtlabel)
ylim([-max(pks)-10 max(pks)+10]) 
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, [figpath filesep 'Clicks' num2str(k) ...
       '-' num2str(k+length(jth)-1) '-plt.jpg']);
close(gcf)