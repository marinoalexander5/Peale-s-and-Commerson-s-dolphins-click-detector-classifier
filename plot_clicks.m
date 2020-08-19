function hp2 = plot_clicks(jth, i1, i2, ii, click, k, tt, tukey_window, pks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_clicks.m
%
% Create subplots of the time series of each detected click including detection
% numbers and saves the figure as JPG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create minimum size square matrix to locate subplots
subp_shape = NaN(ceil(sqrt(length(jth)))) ; % ver donde ubicar para evitar repetir en cada loop
hp2 = figure(2);
subplot(size(subp_shape,1),size(subp_shape,2), ii-k+1)
% trying to plot only 256 samples for closer time zoom resolution
plot(tt(i1(ii-k+1)+128:i2(ii-k+1)-128),click(129:512-128),'b', ...
    tt(i1(ii-k+1)+128:i2(ii-k+1)-128), tukey_window(129:512-128)*(pks(ii-k+1)+1), '--k');        
axis tight
xticks = get(gca, 'Xtick');
xtlabel = cellstr(datestr(seconds(xticks(ceil(length(xticks)/2))), 'MM:SS.FFF'));
set(gca, 'XTick', xticks(ceil(length(xticks)/2)), 'XTickLabel',xtlabel)
%                 set(gca,'fontsize',10)
%                 xlabel('Time [s]')
%                 ylabel('Amplitude [counts]')
set(gca, 'YTick', [])
title(['Click ', num2str(ii)])