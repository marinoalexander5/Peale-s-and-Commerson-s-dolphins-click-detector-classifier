% fix plotting

% Antes de correr seleccionar: Umbral(th)
% Assuming all files in directory have same sample rate for speed, if not
%   the case it needs to be inscluded in the loop
% Recheck detection indexin to make script more readable
clf
close all
clear
clc
addpath(genpath(pwd))
warning('off', 'signal:findpeaks:largeMinPeakHeight'); % for findpeaks not detecting and OpenGL plotting
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
warning('off', 'MATLAB:xlswrite:AddSheet');
% Ask for folder containing files
inpath = uigetdir(pwd,'Select Location of WAV files');
addpath(inpath);
% files = dir(fullfile(inpath,'*.wav'));
files = getAllFiles (inpath);
% [pathstr, name, ext] = fileparts(inpath);

% Ask for output folder 
outpath = uigetdir(pwd,'Select Location of Detector output');

% Ask for dataset identifier
dataset_name = input('Enter dataset name identifier (for output excel sheet name):  ','s');

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
File = struct([]); % Preallocation salida de detector
ai = audioinfo(files{1}); 
Fs = ai.SampleRate; % Frecuencia de muestreo
nbits = ai.BitsPerSample; % Número de bits por muestra

% Filtro pasa banda
fc1=100000/(Fs/2); % Frecuencia de corte inferior en Hz
fc2=160000/(Fs/2); % Frecuencia de corte superior en Hz
[N,Wn]=buttord(fc1,fc2,0.1,60);
[b,a]=butter(N,[fc1 fc2]);

clip = (2^(nbits))/2-1; % Determinar nivel de saturación
th = 16; % 24 dB  % Determinar umbral SNR (linear ratio, not dB)
dt_max = 5/1000; % Ventana de 5 ms entre cada detección
NFFT = 512; % puntos de FFT
hann_window = hanning(NFFT);
tukey_window = tukeywin(NFFT); % ventana tukey fft
SVMmat = zeros(0,6); % SVM matrix input preallocation
j = 0; % indice de detecciones totales


%% Start detection
%% Definir funcion de transferencia
hydro_list = {'Reson','Antarctic_array','SoundTrap HF300'};
tfnum = input('Select hydrophone used for recordings\n1 = Reson dip hydrophone\n2 = Antarctic Array\n3 = SoundTrap ST300\nNumber: '); 
if isempty(tfnum);
    tfnum = 0;
end
while ~any((1:length(hydro_list)) == tfnum);
   fprintf('Invalid number, try again\n\n') 
   tfnum = input('Select hydrophone used for recordings\n1 = Reson dip hydrophone\n2 = Antarctic Array\n3 = SoundTrap ST300\nNumber: ');
   if isempty(tfnum);
    tfnum = 0;
   end
end
% % Calibración función de transferencia
% if tfnum == 1;
%     [FT,TXT,RAW]=xlsread('función de transferencia.xls',2); % Hidrófono Cethus
%     FT = FT(:,6);
% elseif tfnum == 2;
%     [FT,TXT,RAW]=xlsread('función de transferencia.xls',4);  % Arreglo hidrófonos SIO Antártida
%     FT = FT(:,8); FT = FT(1:256);

% elseif tnum == 3;
%     AGREGAR SOUNDTRAP TRANSFER FUNCTION
% end

%%
for index = 1 : length(files);
    fprintf('File %d/%d\n Importando %s\n', index, length(files), files{index});
    to = 1; % Tiempo inicial del segmento a analizar (en muestras)
    k = 1; % Indice iteración detecciones
    Det = struct('filename', (char),'clickn', 0,'signal', [],'itime', 0,'spectrum', zeros(256,1),'snr',0,'pfreq',0,'cfreq',0,'dur10dB',0,'bw3dB',0,'bw10dB',0);% Preallocation click data structure
    Ts = audioinfo(files{index});
    
    while to <= Ts.TotalSamples;
        if to <= Ts.TotalSamples-Fs;
            [ff,Fs]= audioread(files{index},[to to+Fs-1],'native');  % Fragmento de 1 segundo de archivo WAV
        else
            [ff,Fs]= audioread(files{index},[to Ts.TotalSamples],'native');
        end
        
        ff = ff - mean(ff); % Remover offset DC
        xx = filtfilt(b,a,double(ff));
        tt = (to:to+length(xx)-1)'/Fs; % Vector temporal
        yy_a = hilbert(xx);   % Envolvente de señal filtrada
        nn_a = sqrt(mean(abs(yy_a).^2)); % Estimación de ruido
%         fprintf('DEBUGG: SNR from signal: %d/ from envelope %d\n', rms(xx), nn_a);
        snr_a = abs(yy_a)/nn_a; % Relación señal / ruido
        %
        % Uncomment función seleccionada (conservando ambas para comparación de desempeño)
%         [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to); % i1, i2 in samples
        [jth,i1,i2,pks] = detector_umbral_2(yy_a, nn_a, th, dt_max, Fs, tt); % jth, i1, i2 in samples
%%      Performance check plots
        if plot_check;
            % Whole segment
            if ~isempty(jth);
                figpath = strcat(outpath,'\', dataset_name, '\File', ...
                    num2str(index));
                % Create output folder if not existent
                if ~isdir(figpath);
                    mkdir(figpath);
                end
                figure
                clf
                hp1 = plot(tt,xx,'b',tt([1 end]),th*nn_a*[1 1],'r--',tt(jth),pks,'Or',tt([1 end]),nn_a*[1 1],'k--');        
%                 set(hp(2:4),'linewidth',2)
                set(gca,'fontsize',20)
                xlabel('Time [s]')
                ylabel('Amplitude [counts]')
                ylim([-max(pks)-10 max(pks)+10]) 
                title(files{index});
                saveas(gcf, [figpath '-Frame' num2str(j+1) '.jpg']);
                close(hp1)
            end
        end
        
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
        
%%      Almacenamiento de información de cada click
        FT = 0; %%%%%%%%%%%%%% temporal para DEBBUGGING
        for ii = k : k + length(jth) - 1;
            % Hidrófono usado
            Det(ii).hydrophone = hydro_list{tfnum};
            
            % Click Señal Temporal 
            click = xx(i1(ii-k+1):i2(ii-k+1)).*tukey_window;
            Det(ii).signal = click;
            
            % Filename
            Det(ii).filename = files{index}; 
            
            % Detection date and time UTC
            Det(ii).itime = seconds(tt(i1(ii-k+1)));
            Det(ii).date_time = wavname2date (Det(ii).filename); % Returns datenum
            Det(ii).date_time = datetime(Det(ii).date_time, 'ConvertFrom', 'datenum');
            Det(ii).date_time = Det(ii).date_time + Det(ii).itime; % Apply time position to date ref
            Det(ii).date_time.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            
            % Tiempo Inicial de Click relative to file
            Det(ii).itime.Format = 'hh:mm:ss.SSS';
            
            % Número de Click
            Det(ii).clickn = ii;
%             % Almacenar índices de número de archivo y número de detección
%             indices(j+1,:) = [index ii];
            
            % Eliminar clicks saturados (guardar tiempo para ICI)
            if any(click > clip);
%                 Det(ii).signal = 'clipping';
                % replace (snr,spectrum,acoustic params) with NaN
                fn = fieldnames(Det);
                for f = 1:numel(fn)
                    if (isempty(Det(ii).(fn{f})) )
                        Det(ii).(fn{f}) = nan;
                    end
                end
                continue
            end
            
            % Subplot with all detections in 1s frame 
            if plot_check;
                % Create minimum size square matrix to locate subplots
                subp_shape = NaN(ceil(sqrt(length(jth)))) ;
                hp2 = figure;
                subplot(size(subp_shape,1),size(subp_shape,2), ii-k+1)
                plot(tt(i1(ii-k+1):i2(ii-k+1)),click,'b', ...
                    tt(i1(ii-k+1):i2(ii-k+1)), tukey_window*(pks(ii-k+1)+1), '--k');        
%                 set(hp(2:4),'linewidth',2)
                set(gca,'XTickLabel',{num2str(tt(i1(ii-k+1)),'%.2f'), ...
                    num2str(tt(jth(ii-k+1))),num2str(tt(i2(ii-k+1)),'%.2f')}, ...
                    'XTick',[tt(i1(ii-k+1)) tt(jth(ii-k+1)) tt(i2(ii-k+1))])
                set(gca,'fontsize',12)
                xlabel('Time [s]')
                ylabel('Amplitude [counts]')
                ylim([-max(pks)-10 max(pks)+10]) 
                title(files{index});
            end

            
            % Relación Señal a Ruido del click
            csum = cumsum(click.^2); % Ventana del 95% de energia
            [~,rmswin_start] = min(abs(csum-0.025*csum(end)));
            [~,rmswin_end] = min(abs(csum-0.975*csum(end)));
            clickrms = rms(click(rmswin_start:rmswin_end));
            Det(ii).snr = round(clickrms/nn_a, 2);
            
            % Espectro de Click
            magnitude = 2*abs(fft(click,NFFT)./NFFT);
            magnitude = magnitude(1:NFFT/2);
            magnitudedB = (20*log10(magnitude)) + FT;
            Det(ii).spectrum = magnitudedB;
            
            % Parámetros Acústicos
            [pfreq, cfreq, dur10dB, bw3dB, bw10dB ] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB);
            % note: probably smarter way to round multiple fields
            Det(ii).pfreq = round(pfreq, 2);
            Det(ii).cfreq = round(cfreq, 2);
            Det(ii).dur10dB = round(dur10dB, 2);
            Det(ii).bw3dB = round(bw3dB, 2);
            Det(ii).bw10dB = round(bw10dB, 2);
            
            j = j+1;
        end
        saveas(gcf, strcat(outpath,'\', dataset_name, '\File', ...
            num2str(index),'\Frame', num2str(j+1) ,'-clicks', ...
            num2str(k), '-', num2str(k+length(jth)),'.jpg'));
        close(hp2)
        k = k + length(jth);
        to = to+Fs;
    end
    File(index).Det = Det;
    %% Write ouput files
    % Create output folder if not existent
%     if ~isdir(strcat(outpath,'\', dataset_name));
%         mkdir(strcat(outpath,'\', dataset_name));
%     end
    % Si hay detecciones agregar filename a .txt
    files_clicks = fopen(fullfile(strcat(outpath,'\', dataset_name), 'files-with-clicks.txt'), 'a');
    f_content = fscanf(files_clicks, 's');
    if isempty(strfind(files{index}, f_content)); % ver si file ya está incluido en la lista
        fprintf(files_clicks, '%s \n', files{index});
    end
    fclose(files_clicks);
    % Crear .xls para exportar data
    % Separate table for char parameters
    T0 = cell2table({Det.date_time}', 'VariableNames', {'date_time'});
    T1 = cell2table({Det.filename}', 'VariableNames', {'filename'});
    T2 = cell2table({Det.itime}', 'VariableNames', {'time_in_file'});
    T3 = cell2table({Det.hydrophone}', 'VariableNames', {'hydrophone'});
    T4 = table([Det.clickn]', [Det.snr]', [Det.pfreq]', [Det.cfreq]', [Det.dur10dB]', [Det.bw3dB]', [Det.bw10dB]');
    T4.Properties.VariableNames = {'click_num' 'snr' 'pfreq' 'cfreq' 'dur10db' 'bw3db' 'bw10db'};
    % Concatenate tables and write to file
    T = [T0, T1, T2, T3, T4];
    writetable( T, fullfile(strcat(outpath,'\', dataset_name), 'clicks-features.xls'), 'Sheet', dataset_name);
    %
%     SVMmat = vertcat(SVMmat,[[Det.clickn]' [Det.itime]' [Det.snr]' [Det.pfreq]' [Det.cfreq]' [Det.dur10dB]' [Det.bw3dB]' [Det.bw10dB]' ]); % Matriz de datos a clasificar

end
% SVMmat( ~any(SVMmat,2), : ) = [];
% SVMmat_indexada = [indices SVMmat];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alexander Marino 11/2017
