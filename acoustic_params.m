function [pf, cf , dur10dB, bw3dB, bw10dB] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB)
% function [pf,cf,dur10dB, bw3dB, bw10dB] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB )
%
% click: Time domain click signal  
% Fs: Sample rate
% NFFT: number of DFT points
% magnitude: click spectrum
% magnitudedB: click spectrum in dB
%
% Calculates   pf: peak frequency 
%              cf: centroid frequency
%              dur: 10dB duration  
%              3dbbw: 3dB bandwidth  
%              10dbbw: 10dB bandwidth

persistent freq_kHz 
freq_kHz = Fs/2*linspace(0,1,NFFT/2)*1e-3;
%
% Peak Frequency
[max_f,pos_f] = max(magnitudedB(2:end)); % 2:end para evitar continua
pos_f = pos_f + 1;
pf = freq_kHz(pos_f);
% Centroid Frequency
cf = sum(magnitude.*freq_kHz')/sum(magnitude);
% 10dB Duration
clickenv = abs(hilbert(click));   % Envolvente de señal
[max_t,pos_t] = max(clickenv);
clickenvflip = flipud(clickenv);
dur10dB_first = pos_t - (find(clickenvflip(end-pos_t+1:end)<=(max_t/(10^(10/20))), 1, 'first')) + 1; % 10^(10/20) = 10dB 'linear'
dur10dB_last = pos_t + (find(clickenv(pos_t+1:end)<=(max_t/(10^(10/20))), 1, 'first'));
if isempty(dur10dB_last) == 1; 
    dur10dB_last = length(click);
end
if isempty(dur10dB_first) == 1; 
    dur10dB_first = 1;
end
dur10dB = (dur10dB_last - dur10dB_first)/Fs*1e6; 
% 3dB Bandwidth
magnitudedBflip = flipud(magnitudedB); 
bw3dB_first = pos_f - (find(magnitudedBflip(end-pos_f+1:end)<=(max_f-3), 1, 'first')) + 1;
bw3dB_last = pos_f + (find(magnitudedB(pos_f+1:end)<=(max_f-3), 1, 'first'));
bw3dB = freq_kHz(bw3dB_last)-freq_kHz(bw3dB_first);
% 10dB Bandwidth
bw10dB_first = pos_f - (find(magnitudedBflip(end-pos_f+1:end)<=(max_f-10), 1, 'first')) + 1;
bw10dB_last = pos_f + (find(magnitudedB(pos_f+1:end)<=(max_f-10), 1, 'first'));
bw10dB = freq_kHz(bw10dB_last)-freq_kHz(bw10dB_first);
end

% Alexander Marino 11/2017
