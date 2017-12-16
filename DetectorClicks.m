% Antes de correr seleccionar: Umbral(th), Función transferencia según hidrófono 
clf
close all
clear
clc
addpath(genpath(pwd))
% files = dir('*.wav');
File = struct([]); % Preallocation salida de detector
ai = audioinfo(files(1).name);
Fs = ai.SampleRate; % Frecuencia de muestreo
nbits = ai.BitsPerSample; % Número de bits por muestra
% 
% Filtro pasa banda
fc1=90000/(Fs/2); % Frecuencia de corte inferior en Hz
fc2=240000/(Fs/2); % Frecuencia de corte superior en Hz
[N,Wn]=buttord(fc1,fc2,0.1,60);
[b,a]=butter(N,[fc1 fc2]);
%
clip = (2^(nbits))/2-1; % Determinar nivel de saturación
th = ;  % Determinar umbral SNR
dt_max = 2/1000; % Ventana de 2 ms entre cada detección
NFFT = 512; % puntos de FFT
window = hanning(150); % ventana hanning fft
SVMmat = zeros(0,6); % SVM matrix input preallocation
j = 0; % indice de detecciones totales
%
for index = 1 : length(files);
    fprintf('File %d/%d\n Importando %s\n', index, length(files), files(index).name);
    to = 1; % Tiempo inicial del segmento a analizar (en muestras)
    k = 1; % Indice iteración detecciones
    Det = struct('filename', (char),'clickn', 0,'signal', [],'itime', 0,'spectrum', zeros(256,1),'snr',0,'pfreq',0,'cfreq',0,'dur10dB',0,'bw3dB',0,'bw10dB',0);% Preallocation click data structure
    Ts = audioinfo(files(index).name);
    while to <= Ts.TotalSamples;
        if to <= Ts.TotalSamples-Fs;
            [ff,Fs]= audioread(files(index).name,[to to+Fs-1],'native');  % Fragmento de 1 segundo de archivo WAV
        else
            [ff,Fs]= audioread(files(index).name,[to Ts.TotalSamples],'native');
        end
        ff = ff - mean(ff); % Remover offset DC
        xx = filtfilt(b,a,double(ff));
        tt = (to:to+length(xx)-1)'/Fs; % Vector temporal
        yy_a = hilbert(xx);   % Envolvente de señal filtrada
        nn_a = sqrt(mean(abs(yy_a).^2)); % Estimación de ruido
        snr_a = abs(yy_a)/nn_a; % Relación señal / ruido
        %
        [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to);
        % Redefinir segmento temporal si click queda en el límite de la división o no hay detecciones
        if isempty(jth);
            to = to+Fs;
            continue
        elseif ~any(abs(yy_a(1:jth(1)))<nn_a) && to>1;
            to = to-300;
            continue           
        elseif i2(end) > length(xx);
            jth(end)=[]; i1(end)=[]; i2(end)=[];
            if to > Ts.TotalSamples-Fs;
                to = to-(length(xx)-i1(end));
            end
        end
        %%  Almacenamiento de información de cada click          
%         Calibración función de transferencia
%         [FT,TXT,RAW]=xlsread('función de transferencia.xls',2); % Hidrófono Cethus
%         FT = FT(:,6);
        [FT,TXT,RAW]=xlsread('función de transferencia.xls',4);  % Arreglo hidrófonos SIO Antártida
        FT = FT(:,8); FT = FT(1:256);
        for ii = k : k + length(jth) - 1;
            % Click Señal Temporal 
            click = xx(i1(ii-k+1):i2(ii-k+1));
            % Eliminar clicks saturados
            if isempty(find(click > clip, 1)) == 0;
                Det(ii).signal = 'clipping';
                continue
            end            
            Det(ii).signal = click;
            % Filename
            Det(ii).filename = files(index).name; 
            % Número de Click
            Det(ii).clickn = ii;
            % Tiempo Inicial de Click
            Det(ii).itime = datestr(tt(i1(ii-k+1))/(24*60*60), 'DD:HH:MM:SS.FFF');
            % Almacenar índices de número de archivo y número de detección
            indices(j+1,:) = [index ii];
            % Relación Señal a Ruido del click
            clickrms = sqrt(mean(abs(click).^2));
            Det(ii).snr = clickrms/nn_a;
            % Espectro de Click
            magnitude = 2*abs(fft(click.*window,NFFT)./NFFT);
            magnitude = magnitude(1:NFFT/2);
            magnitudedB = (20*log10(magnitude)) + FT;
            Det(ii).spectrum = magnitudedB;
            % Parámetros Acústicos
            [pfreq, cfreq, dur10dB, bw3dB, bw10dB ] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB);
            Det(ii).pfreq = pfreq;
            Det(ii).cfreq = cfreq;
            Det(ii).dur10dB = dur10dB;
            Det(ii).bw3dB = bw3dB;
            Det(ii).bw10dB = bw10dB;
            %
            j = j+1;
        end
        k = k + length(jth);
        to = to+Fs;
    end
      SVMmat = vertcat(SVMmat,[[Det.snr]' [Det.pfreq]' [Det.cfreq]' [Det.dur10dB]' [Det.bw3dB]' [Det.bw10dB]']); % Matriz de datos a clasificar
    File(index).Det = Det;
end
SVMmat( ~any(SVMmat,2), : ) = [];
SVMmat_indexada = [indices SVMmat];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alexander Marino 11/2017
