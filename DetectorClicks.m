% TODO:
    % added corrected cf to acoustic_params and bwrms. Check array
    % dimensios to confirm and added bwrms when calling from main script.    
        
% Commit message:
%       Fixed excel output to fit in same file - different pages
%       changed dB to lowercase for easier handling 
%       prompt for threshold input
%       Added bwrms, Qrms and Q3db
%       Changed centroid frequency calculation

% Assuming all files in directory have same sample rate for speed, if not
%   the case it needs to be included in the loop
% Recheck detection indexing to make script more readable
clf
close all
clear
clc
addpath(genpath(pwd))
warning('off', 'signal:findpeaks:largeMinPeakHeight'); % for findpeaks not detecting and OpenGL plotting
warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
warning('off', 'MATLAB:xlswrite:AddSheet');
% Temporary disable figures to avoid popping up and improve speed(not working)
% figure(1, 'Visible', 'off')
% figure(2, 'Visible', 'off')
% Ask for folder containing files
inpath = uigetdir(pwd,'Select Location of WAV files');
addpath(genpath(inpath));
files = getAllFiles (inpath);

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
nbits = ai.BitsPerSample; % N�mero de bits por muestra

% Filtro pasa banda
fc1=90000/(Fs/2); % Frecuencia de corte inferior en Hz % 50 k para incluir BBPs
fc2=160000/(Fs/2); % Frecuencia de corte superior en Hz
[N,Wn]=buttord(fc1,fc2,0.1,60);
[b,a]=butter(N,[fc1 fc2]);


%% Ask for SNR threshold (linear ratio, not dB) - default 10 
th = input('Enter linear SNR ratio or press enter for defautl (default=10):  ');
if isempty(th)
    th = 10;    
end

clip = (2^(nbits))/2-1; % Determinar nivel de saturaci�n
dt_max = 0.5/1000; % Ventana de 5 ms entre cada detecci�n
NFFT = 512; % puntos de FFT
hann_window = hanning(NFFT);
tukey_window = tukeywin(NFFT); % ventana tukey fft
noverlap = 0.5 * NFFT;
SVMmat = zeros(0,6); % SVM matrix input preallocation
j = 0; % indice de detecciones totales

%% Start detection
%% Definir funcion de transferencia
[tfnum, FT, hydrophone] = load_tf(Fs, NFFT);
%%
for index = 1 : length(files)
    fprintf('File %d/%d\n Importando %s\n', index, length(files), files{index});
    to = 1; % Tiempo inicial del segmento a analizar (en muestras)
    k = 1; % Indice iteraci�n detecciones
    Det = struct('filepath', (char),'filename', (char),'clickn', 0,'signal', [],'itime', 0,'spectrum', zeros(256,1),'snr',0,'pfreq',0,'cfreq',0,'dur10dB',0,'bw3dB',0,'bw10dB',0);% Preallocation click data structure
    Ts = audioinfo(files{index});
    [fpath,fname,~] = fileparts(files{index});
    figpath = strcat(outpath,filesep, dataset_name, filesep, fname);
    while to <= Ts.TotalSamples
        if to <= Ts.TotalSamples-Fs
            [ff,Fs]= audioread(files{index},[to to+Fs-1],'native');  % Fragmento de 1 segundo de archivo WAV
        else
            [ff,Fs]= audioread(files{index},[to Ts.TotalSamples],'native');
        end
        
        ff = ff - mean(ff); % Remover offset DC
        xx = filtfilt(b,a,double(ff));
        tt = (to:to+length(xx)-1)'/Fs; % Vector temporal
        yy_a = hilbert(xx);   % Envolvente de se�al filtrada
        nn_a = sqrt(mean(abs(yy_a).^2)); % Estimaci�n de ruido
        snr_a = abs(yy_a)/nn_a; % Relaci�n se�al / ruido
        %
        % Uncomment funci�n seleccionada (conservando ambas para comparaci�n de desempe�o)
%         [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to); % i1, i2 in samples
        [jth,i1,i2,pks] = detector_umbral_2(yy_a, nn_a, th, dt_max, Fs); % jth, i1, i2 in samples
        if isempty(jth)
            to = to+Fs;
            continue
        end
%%      Performance check plots (probably better way to organize parameters)
        if plot_check
            plot_frame(jth, pks, figpath, fpath, fname, hann_window, noverlap, NFFT, Fs, k, tt, xx, th, nn_a);
        end        
%%      Almacenamiento de informaci�n de cada click
        FT = 0; %%%%%%%%%%%%%% temporal para DEBBUGGING
        for ii = k : k + length(jth) - 1
            % Hidr�fono usado
            Det(ii).hydrophone = hydrophone;
            
            % Click Se�al Temporal 
            try
                click = xx(i1(ii-k+1):i2(ii-k+1)).*tukey_window;
            catch
                warning('Click cut off by segment. Saving just for ICI.');
                click = [xx(i1(ii-k+1):i2(ii-k+1)); clip+1];
            end
            Det(ii).signal = click;
            
            % Filepath
            Det(ii).filepath = fpath;
            
            % Filename
            Det(ii).filename = fname; 
            
            % Detection date and time UTC
            Det(ii).itime = seconds(tt(i1(ii-k+1)));
            Det(ii).date_time = wavname2date (Det(ii).filename); % Returns datenum
            Det(ii).date_time = datetime(Det(ii).date_time, 'ConvertFrom', 'datenum');
            Det(ii).date_time = Det(ii).date_time + Det(ii).itime; % Apply time position to date ref
            Det(ii).date_time.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            
            % Tiempo Inicial de Click relative to file
            Det(ii).itime.Format = 'hh:mm:ss.SSS';
            
            % N�mero de Click
            Det(ii).clickn = ii;
%             % Almacenar �ndices de n�mero de archivo y n�mero de detecci�n
%             indices(j+1,:) = [index ii];
            
            % Eliminar clicks saturados (guardar tiempo para ICI)
            if any(click > clip)
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
            if plot_check
                hp2 = plot_clicks(jth, i1, i2, ii, click, k, tt, tukey_window, pks);
            end

            
            % Relaci�n Se�al a Ruido del click
            csum = cumsum(click.^2); % Ventana del 95% de energia
            [~,rmswin_start] = min(abs(csum-0.025*csum(end)));
            [~,rmswin_end] = min(abs(csum-0.975*csum(end)));
            click_rms = rms(click(rmswin_start:rmswin_end));
            try
                noise_rms = rms(xx(i1(ii-k+1)-500:i1(ii-k+1)));
            catch
                warning('Click cut off by segment. Careful with estimated SNR.')
                noise_rms = rms(xx(1:i1(ii-k+1)));
            end
            click_rms = click_rms - noise_rms;
            Det(ii).snr = round(click_rms/noise_rms, 2);
            
            % Espectro de Click
            P2 = abs(fft(click,NFFT)/NFFT);
            magnitude = P2(1:NFFT/2+1);
            magnitude(2:end-1) = 2*magnitude(2:end-1);
            magnitudedB = (10*log10(magnitude)) + FT; % 20log10 power spectrum
            Det(ii).spectrum = magnitudedB;
            
            % Par�metros Ac�sticos
            [pfreq, cfreq , dur10dB, bw3db, bw10db, bwrms, Q3db, Qrms] = acoustic_params(click, Fs, NFFT, magnitude, magnitudedB);
            % note: probably smarter way to round multiple fields
            Det(ii).pfreq = round(pfreq, 2);
            Det(ii).cfreq = round(cfreq, 2);
            Det(ii).dur10dB = round(dur10db, 2);
            Det(ii).bw3dB = round(bw3db, 2);
            Det(ii).bw10dB = round(bw10db, 2);
            Det(ii).bwrms = round(bwrms, 2);
            Det(ii).Q3db = round(Q3db, 2);
            Det(ii).Qrms = round(Qrms, 2);
            
            j = j+1;
        end
        saveas(hp2, [figpath filesep 'Clicks' num2str(k) ...
                       '-' num2str(k+length(jth)-1) '-subplt.jpg']);
        clf(hp2)
        k = k + length(jth);
        to = to+Fs;
    end
    if k > 1  
        File(index).Det = Det;
        %% Write ouput files
        % Si hay detecciones agregar filename a .txt
        files_clicks = fopen(fullfile(strcat(outpath, filesep, dataset_name), 'files-with-clicks.txt'), 'a');
        f_content = fscanf(files_clicks, 's');
        if isempty(strfind(files{index}, f_content)) % ver si file ya est� incluido en la lista
            fprintf(files_clicks, '%s \n', files{index});
        end
        fclose(files_clicks);
        % Crear .xls para exportar data
        % Separate table for char parameters
        %%%% smarter to tramsform and filter columns once last file is finished
            % T = struct2table(struct2array(File))
            % T = T(:,{'date_time' 'filepath' 'filename' 'itime' 'hydrophone' 'clickn' 'snr' 'pfreq' 'cfreq' 'dur10dB' 'bw3dB' 'bw10dB'})
        T0 = cell2table({Det.date_time}', 'VariableNames', {'date_time'});
        T1 = cell2table({Det.filepath}', 'VariableNames', {'filepath'});
        T2 = cell2table({Det.filename}', 'VariableNames', {'filename'});
        T3 = cell2table({Det.itime}', 'VariableNames', {'time_in_file'});
        T4 = cell2table({Det.hydrophone}', 'VariableNames', {'hydrophone'});
        T5 = table([Det.clickn]', [Det.snr]', [Det.pfreq]', [Det.cfreq]', [Det.dur10dB]', [Det.bw3dB]', [Det.bw10dB]', [Det.bwrms]', [Det.Q3db]', [Det.Qrms]');
        T5.Properties.VariableNames = {'click_num' 'snr' 'pfreq' 'cfreq' 'dur10db' 'bw3db' 'bw10db' 'bwrms' 'Q3db' 'Qrms'};
        % Concatenate tables and write to file
        T = [T0, T1, T2, T3, T4, T5];
        writetable( T, fullfile(outpath, 'clicks-features.xls'), 'Sheet', dataset_name);
        % Save output struct to .mat file
        save(fullfile(outpath, dataset_name, 'Detections.mat'),'File')        
    %     SVMmat = vertcat(SVMmat,[[Det.clickn]' [Det.itime]' [Det.snr]' [Det.pfreq]' [Det.cfreq]' [Det.dur10dB]' [Det.bw3dB]' [Det.bw10dB]' ]); % Matriz de datos a clasificar
    end
end
% SVMmat( ~any(SVMmat,2), : ) = [];
% SVMmat_indexada = [indices SVMmat];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alexander Marino 08/2020
