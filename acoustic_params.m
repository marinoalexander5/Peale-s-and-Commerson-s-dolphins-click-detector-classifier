function [pf, cf , dur10db, bw3db, bw10db, bwrms, Q3db, Qrms] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedb)
% function [pf,cf,dur10db, bw3db, bw10db] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedb )
%
% click: Time domain click signal  
% Fs: Sample rate
% NFFT: number of DFT points
% magnitude: click spectrum
% magnitudedb: click spectrum in db
%
% Calculates   pf: peak frequency 
%              cf: centroid frequency
%              dur: 10db duration  
%              3dbbw: 3db bandwidth  
%              10dbbw: 10db bandwidth

persistent freq_kHz
freq_kHz = Fs*(0:NFFT/2)/NFFT*1e-3;
%

% Peak Frequency
[max_f,pos_f] = max(magnitudedb(2:end)); % 2:end para evitar continua
pos_f = pos_f + 1;
pf = freq_kHz(pos_f);

% Centroid Frequency
cf = sum(abs(magnitude.^2).*freq_kHz')/sum(abs(magnitude.^2)); % (new)
% cf = sum(magnitude.*freq_kHz')/sum(magnitude);

% 10db Duration
clickenv = abs(hilbert(click));   % Envolvente de se�al
[max_t,pos_t] = max(clickenv);
clickenvflip = flipud(clickenv);
dur10db_first = pos_t - (find(clickenvflip(end-pos_t+1:end)<=(max_t/(10^(10/20))), 1, 'first')) + 1; % 10^(10/20) = 10db 'linear'
dur10db_last = pos_t + (find(clickenv(pos_t+1:end)<=(max_t/(10^(10/20))), 1, 'first'));
if isempty(dur10db_last) == 1
    dur10db_last = length(click);
end
if isempty(dur10db_first) == 1 
    dur10db_first = 1;
end
dur10db = (dur10db_last - dur10db_first)/Fs*1e6; 

% 3db Bandwidth
magnitudedbflip = flipud(magnitudedb); 
bw3db_first = pos_f - (find(magnitudedbflip(end-pos_f+1:end)<=(max_f-3), 1, 'first')) + 1;
bw3db_last = pos_f + (find(magnitudedb(pos_f+1:end)<=(max_f-3), 1, 'first'));
bw3db = freq_kHz(bw3db_last)-freq_kHz(bw3db_first);

% 10db Bandwidth
bw10db_first = pos_f - (find(magnitudedbflip(end-pos_f+1:end)<=(max_f-10), 1, 'first')) + 1;
bw10db_last = pos_f + (find(magnitudedb(pos_f+1:end)<=(max_f-10), 1, 'first'));
bw10db = freq_kHz(bw10db_last)-freq_kHz(bw10db_first);

% RMS Bandwidth
bwrms = sqrt(sum(((freq_kHz'-cf).^2).*(abs(magnitude).^2))/sum(abs(magnitude).^2));

% Q3db
Q3db = cf/bw3db;

% Qrms
Qrms = cf/bwrms;

% Alexander Marino 11/2017
