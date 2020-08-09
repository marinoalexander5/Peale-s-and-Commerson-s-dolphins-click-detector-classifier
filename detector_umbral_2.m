function [jth, i1, i2, pks] = detector_umbral_2(yy_a, nn_a, th, dt_max, Fs, tt)
%%%
%%%   [jth,i1,i2] = detector_umbral_2(snr_a,tt,th,dt_max,yy_a,nn_a,to)
%%%   
%%%   Peak detection function for click extraction
%%%   
%%%   snr_a: instant signal to noise ratio vector for whole signal (sample
%%%         by sample)
%%%   tt: time vector of the lenght of the signal
%%%   th: detector threshold (in times, linear ratio)
%%%   dt_max: minimum time window between detections
%%%   yy_a: envelope of the signal
%%%   nn_a: estimated noise rms value
%%%   to: initial time indicator
%%%
%%%   Retunrs:
%%%   jth: vector of max peaks above threshold
%%%   i1: vector of detection starting indices
%%%   i2: vector of detection ending indices

[pks, locs] = findpeaks(abs(yy_a), Fs, 'MinPeakDistance', dt_max, 'MinPeakHeight', th*nn_a); % locs in time(tt) scale
jth = round(locs*Fs); % turning to samples for easier indexing
if any(jth)
    i1 = jth - 256;
    i2 = jth + 255;
else
    i1 = [];
    i2 = [];
end

