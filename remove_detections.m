function [jth_rm, pks_rm] = remove_detections(tt, jth, pks, k, xx, th, nn_a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% remove_detections.m
%
% Remove detections by clicking on plot.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nClick on detection marker to remove, press ENTER to finish')
confirm = 'N';
while (confirm == 'N' || confirm == 'n') 
    subplot(2,1,2)
    cnums = arrayfun(@num2str, k:k+length(jth)-1, 'UniformOutput', 0); % cell array of click numbers
    plot(tt,xx,'b', ...
             tt([1 end]),th*nn_a*[1 1],'r--', ...
             tt([1 end]),nn_a*[1 1],'k--', ...
             tt(jth),pks,'Or', 'MarkerSize', 5) % Comment if plot gets crowded  
        text(tt(jth), pks, cnums, 'FontSize', 10)
    xlabel('Time [s]')
    ylabel('Amplitude [counts]')
    xticks = get(gca, 'Xtick');
    xtlabel = cellstr(datestr(seconds(xticks), 'MM:SS.FFF'));
    set(gca, 'XTick', xticks, 'XTickLabel',xtlabel)
    ylim([-max(pks)-10 max(pks)+10]) 
    set(gcf, 'Position', get(0, 'Screensize'));
    [ptsx, ~] = getpts;
    if any(ptsx)
%         rm_dets = ~ismembertol(tt(jth), ptsx, 5e-3);
        rm_idx = knnsearch(tt(jth),ptsx,'k',1);
        jth_rm = jth;
        pks_rm = pks;
        jth_rm(rm_idx) = [];
        pks_rm(rm_idx) = [];
        subplot(2,1,2)
        cnums = arrayfun(@num2str, k:k+length(jth_rm)-1, 'UniformOutput', 0); % cell array of click numbers
        plot(tt,xx,'b', ...
             tt([1 end]),th*nn_a*[1 1],'r--', ...
             tt([1 end]),nn_a*[1 1],'k--', ...
             tt(jth_rm),pks_rm,'Or', 'MarkerSize', 5) % Comment if plot gets crowded  
        text(tt(jth_rm), pks_rm, cnums, 'FontSize', 10)
        xlabel('Time [s]')
        ylabel('Amplitude [counts]')
        xticks = get(gca, 'Xtick');
        xtlabel = cellstr(datestr(seconds(xticks), 'MM:SS.FFF'));
        set(gca, 'XTick', xticks, 'XTickLabel',xtlabel)
        ylim([-max(pks)-10 max(pks)+10]) 
        set(gcf, 'Position', get(0, 'Screensize'));
    else
        jth_rm = jth;
        pks_rm = pks;
    end
    % Confirm detections
    confirm = input('\nConfirm Detections (Y/N): ','s');
    if isempty(confirm)
        confirm = 'y';
    end
end