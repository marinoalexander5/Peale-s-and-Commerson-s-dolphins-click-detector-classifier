% Antes de correr seleccionar: Umbral(th), Función transferencia según hidrófono 
clf
close all
clear
clc
addpath(genpath(pwd))
warning('off', 'signal:findpeaks:largeMinPeakHeight'); % for findpeaks not detecting and OpenGL plotting
%% Ask if plots need to be checked for debugging
answer = input('Check detection plots? [y/n]:  ', 's');
while ~any(strcmp({'y', 'n'},answer))
  fprintf('Please answer y or n\n') 
  answer = input('Check detection plots? [y/n]:  ', 's');
end
switch answer
    case 'y'
        plot_check = 1;
    case 'n'
        plot_check = 0;
end
%%
files = dir('*.wav');
%
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
hann_window = hanning(NFFT);
tukey_window = tukeywin(NFFT); % ventana tukey fft
SVMmat = zeros(0,6); % SVM matrix input preallocation
j = 0; % indice de detecciones totales
%
% Definir funcion de transferencia
hydro_list = {'Reson','Antarctic_array'};
tfnum = input('Select hydrophone used for recordings\n 1 = Reson dip hydrophone            2 = Antarctic Array\nNumber: ');
while ~any([1, 2] == tfnum) || isempty(tfnum);
   fprintf('Invalid number, try again\n\n') 
   tfnum = input('Select hydrophone used for recordings\n 1 = Reson dip hydrophone            2 = Antarctic Array\nNumber: ');
end
% % Calibración función de transferencia
% if tfnum == 1;
%     [FT,TXT,RAW]=xlsread('función de transferencia.xls',2); % Hidrófono Cethus
%     FT = FT(:,6);
% elseif tfnum == 2;
%     [FT,TXT,RAW]=xlsread('función de transferencia.xls',4);  % Arreglo hidrófonos SIO Antártida
%     FT = FT(:,8); FT = FT(1:256);
% end
%%
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
        % Uncomment función seleccionada (conservando ambas para comparación de desempeño)
%         [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to); % i1, i2 in samples
        [jth,i1,i2,pks] = detector_umbral_2(yy_a, nn_a, th, dt_max, Fs, tt); % jth, i1, i2 in samples
%       % DEBUGGING PLOTS
        if plot_check;
            % Whole segment
            figure (1)
            plot(tt, yy_a, [tt(1) tt(end)], [th*nn_a th*nn_a], '--r')
            if any(jth);
                hold on
                plot(to/Fs+jth/Fs , pks, 'Or');
            end
            for i = (1:length(jth));
                plot((to/Fs+(jth(i)-256)/Fs:1/Fs:to/Fs+(jth(i)+255)/Fs),tukey_window*(pks(i)+1), '--k') % Centred 512 samples Tuckey window in click peak 
                hold on
            end
        end%
        if isempty(jth);
            to = to+Fs;
            continue
        % Redefinir segmento temporal si click queda en el límite de la división o no hay detecciones
%         elseif ~any(abs(yy_a(1:jth(1)))<nn_a) && to>1;
%             to = to-300;
%             continue           
%         elseif i2(end) > length(xx);
%             jth(end)=[]; i1(end)=[]; i2(end)=[];
%             if to > Ts.TotalSamples-Fs;
%                 to = to-(length(xx)-i1(end));
%             end
        end
%         %%  Almacenamiento de información de cada click
        FT = 0; %%%%%%%%%%%%%% temporal para DEBBUGGING
        for ii = k : k + length(jth) - 1;
            % Click Señal Temporal 
            click = xx(i1(ii-k+1):i2(ii-k+1)).*tukey_window;
            % Eliminar clicks saturados (guardar tiempo para ICI)
            if any(click > clip);
%                 Det(ii).signal = 'clipping';
                Det(ii).filename = files(index).name; 
                Det(ii).itime = datestr(tt(floor(i1(ii-k+1)))/(24*60*60), 'DD:HH:MM:SS.FFF');
                Det(ii).hydrophone = hydro_list{tfnum};
                % replace the rest with Nan
                fn = fieldnames(Det);
                for f = 1:numel(fn)
                    if (isempty(Det(ii).(fn{f})) )
                        Det(ii).(fn{f}) = nan;
                    end
                end
                continue
            end            
            Det(ii).signal = click;
            % Filename
            Det(ii).filename = files(index).name; 
            % Número de Click
            Det(ii).clickn = ii;
            % Tiempo Inicial de Click
            Det(ii).itime = datestr(tt(floor(i1(ii-k+1)))/(24*60*60), 'DD:HH:MM:SS.FFF');
            % Almacenar índices de número de archivo y número de detección
            indices(j+1,:) = [index ii];
            % Relación Señal a Ruido del click
            csum = cumsum(click.^2); % Ventana del 95% de energia
            [~,rmswin_start] = min(abs(csum-0.025*csum(end)));
            [~,rmswin_end] = min(abs(csum-0.975*csum(end)));
            clickrms = rms(click(rmswin_start:rmswin_end));
            Det(ii).snr = clickrms/nn_a;
            % Espectro de Click
            magnitude = 2*abs(fft(click,NFFT)./NFFT);
            magnitude = magnitude(1:NFFT/2);
            magnitudedB = (20*log10(magnitude)) + FT;
            Det(ii).spectrum = magnitudedB;
            % Parámetros Acústicos
            [pfreq, cfreq, dur10dB, bw3dB, bw10dB ] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB);
            % note: probably smarter to use arrfun to round
            Det(ii).pfreq = round(pfreq, 2);
            Det(ii).cfreq = round(cfreq, 2);
            Det(ii).dur10dB = round(dur10dB, 2);
            Det(ii).bw3dB = round(bw3dB, 2);
            Det(ii).bw10dB = round(bw10dB, 2);
            % Hidrófono usado
            Det(ii).hydrophone = hydro_list{tfnum};
            j = j+1;
        end
        k = k + length(jth);
        to = to+Fs;
    end
    File(index).Det = Det;
    %% Write ouput files
    % Create output folder if not existent
    if ~isdir('../detector-output/delfin-austral');
        mkdir('../detector-output/delfin-austral');
    end
    % Si hay detecciones agregar filename a .txt
    files_clicks = fopen('../detector-output/delfin-austral/files-with-clicks.txt', 'a');
    f_content = fscanf(files_clicks, 's');
    if isempty(strfind(files(index).name, f_content)); % ver si file ya está incluido en la lista
        fprintf(files_clicks, '%s \n', files(index).name);
    end
    fclose(files_clicks);
    % Crear .xls para exportar data
    % Separate table for char parameters
    T1 = cell2table({Det.filename}', 'VariableNames', {'filename'});
    T2 = cell2table({Det.itime}', 'VariableNames', {'time_in_file'});
    T3 = cell2table({Det.hydrophone}', 'VariableNames', {'hydrophone'});
    T4 = table([Det.clickn]', [Det.snr]', [Det.pfreq]', [Det.cfreq]', [Det.dur10dB]', [Det.bw3dB]', [Det.bw10dB]');
    T4.Properties.VariableNames = {'click_num' 'snr' 'pfreq' 'cfreq' 'dur10db' 'bw3db' 'bw10db'};
    % Concatenate tables and write to file
    T = [T1, T2, T3, T4];
    writetable( T, '../delfin-austral/detector-output.xls');
%     'date' 'time_in_day' ;
    %
%     SVMmat = vertcat(SVMmat,[[Det.clickn]' [Det.itime]' [Det.snr]' [Det.pfreq]' [Det.cfreq]' [Det.dur10dB]' [Det.bw3dB]' [Det.bw10dB]' ]); % Matriz de datos a clasificar

end
% SVMmat( ~any(SVMmat,2), : ) = [];
% SVMmat_indexada = [indices SVMmat];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alexander Marino 11/2017
