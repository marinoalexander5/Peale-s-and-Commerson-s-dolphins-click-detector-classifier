function [jth,i1,i2] = detector_umbral(snr_a,tt,th,dt_max,yy_a,nn_a,to)
ith = find (snr_a > th); 
jth=ith;
p = 1; 
for a = 2:length(ith);
    if tt(jth(a))-tt(jth(p))> dt_max; % a:muestra actual  / p:muestra previa
    p = a;
    else
        jth(a)= 0;
    end
end
jth(jth==0)=[]; % Conservar solamente primer cruce de umbral para cada detección
%
% Ventana de duración fija para cada click (150 muestras) = 300 µs
i1 = zeros(length(jth),1); % Preallocation
if isempty(jth) || (~any(abs(yy_a(1:jth(1)))<nn_a) && to>1); % Redefinir segmento temporal si hay detección incompleta al inicio
    i2=[];
    return
else if ~any(abs(yy_a(1:jth(1)))<nn_a);
        jth(1)=[]; i1(1)=[];
    end
    yy_aflip = flipud(abs(yy_a)); 
    for n = 1:length(jth); % Calcular tiempos iniciales y finales de detecciones
        cs_nn = find(yy_aflip(end-jth(n)+1:end)<= nn_a, 1,'first'); % Primer cruce por debajo de nivel de ruido previo al primer cruce de umbral (click start)
        i1(n) = jth(n) - cs_nn + 1;  
    end
    correct = jth-i1 > 100;
    i1(correct) = jth(correct) - 20;
    i2 = i1+149; % Fin del click 
end


% Alexander Marino 11/2017