function th = input_threshold(files, index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input_threshold.m
%
% Asks user to input detection SNR threshold. Plots signal and input value
% to confirm value or enter a new one.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import entire file (should be shorter than 10 s)
[ff,fs]= audioread(files{index},'native');

% Setup filter
fc1 = 90000/(fs/2); % Frecuencia de corte inferior en Hz % 50 k para incluir BBPs
fc2 = 160000/(fs/2); % Frecuencia de corte superior en Hz
[N,~] = buttord(fc1,fc2,0.1,60);
[b,a] = butter(N,[fc1 fc2]);

% Setup signal
ff = ff - mean(ff);
xx = filtfilt(b,a,double(ff));
tt = (0:length(xx)-1)'/fs;
nn_a = rms(xx);
       
confirm = 'N';
while confirm == 'N'
    th = input('\nEnter linear SNR ratio or press enter for defautl (default=10):  ');
    if isempty(th)
        th = 10;    
    end
    % Check plot
    plot(tt,xx,'b', ...
    tt([1 end]),th*nn_a*[1 1],'r--', ...
    tt([1 end]),nn_a*[1 1],'k--') 
    % Confirm choice from plot
    confirm = input('Confirm Threshold (Y/N): ','s');
end
clf(gcf)

