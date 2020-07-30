% Antes de correr seleccionar: Umbral(th), Función transferencia según hidrófono 
clf
close all
clear
clc
addpath(genpath(pwd))
files = dir('*.wav');
File = struct([]); % Preallocation salida de detector
ai = audioinfo(files(1).name);
Fs = ai.SampleRate; % Frecuencia de muestreo
nbits = ai.BitsPerSample; % Número de bits por muestra
% 
% Filtro pasa banda
fc1=100000/(Fs/2); % Frecuencia de corte inferior en Hz
fc2=160000/(Fs/2); % Frecuencia de corte superior en Hz
[N,Wn]=buttord(fc1,fc2,0.1,60);
[b,a]=butter(N,[fc1 fc2]);
%
clip = (2^(nbits))/2-1; % Determinar nivel de saturación
th = 5;  % Determinar umbral SNR (linear ratio, not dB)
dt_max = 5/1000; % Ventana de 5 ms entre cada detección
NFFT = 512; % puntos de FFT
han_window = hanning(150); % ventana hanning fft
SVMmat = zeros(0,6); % SVM matrix input preallocation
j = 0; % indice de detecciones totales
%
% Definir funcion de transferencia
tfnum = input('Select hydrophone used for recordings\n 1 = Reson dip hydrophone            2 = Antarctic Array\nNumber: ');
while ~any([1, 2] == tfnum);
   fprintf('Invalid number, try again\n\n') 
   tfnum = input('Select hydrophone used for recordings\n 1 = Reson dip hydrophone            2 = Antarctic Array\nNumber: ');
end
%%
counter = 0;
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
%         fprintf('DEBUGG: SNR from singal: %d/ from envelope %d\n', rms(xx), nn_a);
        snr_a = abs(yy_a)/nn_a; % Relación señal / ruido
        %
%%      %%%%%%%%%%%%%%%%% DEBUGGING START
        figure (1)
        plot(tt, yy_a, [tt(1) tt(end)], [th th], '--r')
        [pks, locs] = findpeaks(abs(yy_a), Fs, 'MinPeakDistance', dt_max, 'MinPeakHeight', th*nn_a); % locs in time(tt) scale
        % plot test
        if any(pks);
            hold on
            plot(counter+locs , pks, 'Or');
            plot ([0 10],[th th], '--r')
        end
        % crear branck para trabaajar en una version alternativa
        % arreglar esta parte para ver como ajustar correctamente ventana tukey centrada en pks y locs(unificar locs(tiempo) con nro muestras (+256 y -256))
        for i = (1:length(locs));
            plot((counter+locs(i)-256/Fs:1/Fs:counter+locs(i)+255/Fs),tukeywin(512)*(pks(i)+1), '--k') % Centred 512 samples Tuckey window in click peak 
            hold on
        end
        counter = counter + 1;
        to = to+Fs;
    end
end
    
%%        %%%%%%%%%%%% DEBUGING END
        [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to);
%         % Redefinir segmento temporal si click queda en el límite de la división o no hay detecciones
%         if isempty(jth);
%             to = to+Fs;
%             continue
%         elseif ~any(abs(yy_a(1:jth(1)))<nn_a) && to>1;
%             to = to-300;
%             continue           
%         elseif i2(end) > length(xx);
%             jth(end)=[]; i1(end)=[]; i2(end)=[];
%             if to > Ts.TotalSamples-Fs;
%                 to = to-(length(xx)-i1(end));
%             end
%         end
%         %%  Almacenamiento de información de cada click
% %         Calibración función de transferencia
%         if tfnum == 1;
%             [FT,TXT,RAW]=xlsread('función de transferencia.xls',2); % Hidrófono Cethus
%             FT = FT(:,6);
%         elseif tnum == 2;
%             [FT,TXT,RAW]=xlsread('función de transferencia.xls',4);  % Arreglo hidrófonos SIO Antártida
%             FT = FT(:,8); FT = FT(1:256);
%         end
%         for ii = k : k + length(jth) - 1;
%             % Click Señal Temporal 
%             click = xx(i1(ii-k+1):i2(ii-k+1));
%             % Eliminar clicks saturados
%             if isempty(find(click > clip, 1)) == 0;
%                 Det(ii).signal = 'clipping';
%                 continue
%             end            
%             Det(ii).signal = click;
%             % Filename
%             Det(ii).filename = files(index).name; 
%             % Número de Click
%             Det(ii).clickn = ii;
%             % Tiempo Inicial de Click
%             Det(ii).itime = datestr(tt(i1(ii-k+1))/(24*60*60), 'DD:HH:MM:SS.FFF');
%             % Almacenar índices de número de archivo y número de detección
%             indices(j+1,:) = [index ii];
%             % Relación Señal a Ruido del click
%             clickrms = sqrt(mean(abs(click).^2));
%             Det(ii).snr = clickrms/nn_a;
%             % Espectro de Click
%             magnitude = 2*abs(fft(click.*hann_window,NFFT)./NFFT);
%             magnitude = magnitude(1:NFFT/2);
%             magnitudedB = (20*log10(magnitude)) + FT;
%             Det(ii).spectrum = magnitudedB;
%             % Parámetros Acústicos
%             [pfreq, cfreq, dur10dB, bw3dB, bw10dB ] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB);
%             Det(ii).pfreq = pfreq;
%             Det(ii).cfreq = cfreq;
%             Det(ii).dur10dB = dur10dB;
%             Det(ii).bw3dB = bw3dB;
%             Det(ii).bw10dB = bw10dB;
%             %
%             j = j+1;
%         end
%         k = k + length(jth);
%         to = to+Fs;
%     end
%       SVMmat = vertcat(SVMmat,[[Det.snr]' [Det.pfreq]' [Det.cfreq]' [Det.dur10dB]' [Det.bw3dB]' [Det.bw10dB]']); % Matriz de datos a clasificar
%     File(index).Det = Det;
%     %% Si hay detecciones agregar filename a .txt
%     files_clicks = fopen('files-with-clicks.txt', 'a');
%     f_content = fscanf(files_clicks, 's');
%     if isempty(strfind(files(index).name, f_content)); % ver si file ya está incluido en la lista
%         fprintf(files_clicks, '%s \n', files(index).name);
%     end
%     fclose(files_clicks);
% end
% SVMmat( ~any(SVMmat,2), : ) = [];
% SVMmat_indexada = [indices SVMmat];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alexander Marino 11/2017
