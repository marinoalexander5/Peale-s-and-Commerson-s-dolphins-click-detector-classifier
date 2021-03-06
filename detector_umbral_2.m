function [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to, Fs)
%%%
%%%   [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to)
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
%%%   jth: 
%%%   i1: vector of detection starting indices
%%%   i2: vector of detection ending indices

[pks, locs] = finpeakes(snr_a, Fs, 'MinPeakDistance', dt_max, 'MinPeakHeight', th, 'Threshold', 1);
