function [tfnum, FT, hydrophone] = load_tf(Fs, NFFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% load_tf.m
%
% Load transfer function according to specified hydrophone
% Adapted from triton (cetus.ucsd.edu, SIO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fax = Fs*(0:NFFT/2)/NFFT*1e-3;
hydro_list = {'Reson','Antarctic_array','SoundTrap HF300'};
tfnum = input('\nSelect hydrophone used for recordings\n1 = Reson dip hydrophone\n2 = Antarctic Array\n3 = SoundTrap ST300\nNumber: '); 
if isempty(tfnum)
    tfnum = 0;
end
while ~any((1:length(hydro_list)) == tfnum)
   fprintf('Invalid number, try again\n\n') 
   tfnum = input('\nSelect hydrophone used for recordings\n1 = Reson dip hydrophone\n2 = Antarctic Array\n3 = SoundTrap ST300\nNumber: ');
   if isempty(tfnum)
    tfnum = 0;
   end
end
switch tfnum
    case 1 % Hidrofono Cethus (HF573)
        tfid = fopen('HF573_101213_invSensit.tf','r');
        [A,~] = fscanf(tfid,'%f %f',[2,inf]);
        fclose(tfid);
    case 2  % Arreglo hidrofonos SIO Antartida (HF629-ch6) %%%% PEDIR A BRUCE LA QUE CORRESPONDE
        tfid = fopen('HF631_140122_sig1_invSensit.tf','r');
        [A,~] = fscanf(tfid,'%f %f',[2,inf]);
%         FT = 0;
        fclose(tfid);
    case 3
        hydrophone = hydro_list{tfnum};
        FT = 0; % zeros(1, NFFT/2+1)
        return
        % No transfer function for SoundTrap, flat over all BW
end
freq = A(1,:);
uppc = A(2,:);    % [dB re uPa(rms)^2/counts^2]
[~,ia,ic] = unique(freq);
if length(ia) ~= length(ic)
    freq = freq(ia);
    uppc = uppc(ia);
end
Ptf = interp1(freq,uppc,fax,'linear','extrap');
FT = 10*log10(Ptf);
hydrophone = hydro_list{tfnum};