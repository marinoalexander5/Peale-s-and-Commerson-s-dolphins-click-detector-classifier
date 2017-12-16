% % % % Clasificador: preparado de datos y busqueda de mejores parametros % % % %
% s:scaled
% r:random
% l:labeled
%% Union matrices
All_data=[Aus_data;To_data];
%% Escalado de set de datos entre [0,1]
Min_Max_scale = [min(All_data); max(All_data)];
for i=1:size(Min_Max_scale,2);
All_data_scale(:,i) = (All_data(:,i) - Min_Max_scale(1,i))/(Min_Max_scale(2,i)-Min_Max_scale(1,i));
end
Aus_data_scale = All_data_scale(1:size(Aus_data,1),:);
To_data_scale = All_data_scale(size(Aus_data,1)+1:end,:);
%% Randomizar
randmat_Aus = [randi(size(Aus_data,1),[size(Aus_data,1) 1]) Aus_data_scale];
randmat_To = [randi(size(To_data,1),[size(To_data,1) 1]) To_data_scale];
Aus_data_s=sortrows(randmat_Aus);
To_data_s=sortrows(randmat_To);
Aus_data_sr=Aus_data_s(:,2:end);
To_data_sr=To_data_s(:,2:end);
%% Labels
Aus_label=ones(size(Aus_data,1),1);
To_label=(-1)*ones(size(To_data,1),1);
%% ploteo SVM 1 vs 1
x=[Aus_data_sr; To_data_sr];
group=[repmat('Delfín Austral', [size(Aus_data,1),1]);repmat('Tonina Overa  ', [size(To_data,1),1])]; %[Aus_label; To_label]; 
feature = (1:5)';
feature_label = {'pf' 'cf' 'dur10dB' 'bw3dB' 'bw10dB'};
for k=1:5
for i=2:5      
    figure ()
    %figure (feature(1))
    %subplot(2,2,i-1)
    gscatter(x(:,feature(1)),x(:,feature(i)),group);
    title([feature_label{feature(1)} ' vs ' feature_label{feature(i)}]);
    h = get(gca,'Children');
    set(gca,'Children',[h(2) h(1)]);
if i==5;
    feature = circshift(feature,-1);
end
end
end
%% Dividir train y test
splitnum_Aus = ceil(size(Aus_data,1)*0.8); % Ver porcentaje de training + validation y de testing
splitnum_To = ceil(size(To_data,1)*0.8);
train_data = sparse([Aus_data_sr(1:splitnum_Aus,:); To_data_sr(1:splitnum_To,:)]);
train_label = [Aus_label(1:splitnum_Aus,:); To_label(1:splitnum_To,:)];
test_data = sparse([Aus_data_sr(splitnum_Aus+1:end,:); To_data_sr(splitnum_To+1:end,:)]);
test_label = [Aus_label(splitnum_Aus+1:end,:); To_label(splitnum_To+1:end,:)];
%% libsvm
libsvmwrite('.\train_data', train_label, train_data)
libsvmwrite('.\test_data', test_label, test_data)
[train_label, train_data] = libsvmread('.\train_data');
%% Mejor parametro
bestcv = 0;
n = 0:14;
m = 0:12;
accuracy_mat = zeros(numel(n),numel(m));
for i = 1:10;
for c = (n*9)+2; 
  for g = (m*5)+4; 
    cmd = ['-v 5 -c ', num2str(c), ' -g ', num2str(g)];
    cv = svmtrain2(train_label, train_data, cmd);
    if (cv >= bestcv),
      bestcv = cv; bestc = c; bestg = g;
    end
    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', c, g, cv, bestc, bestg, bestcv);
    accuracy_mat ((c-2)/9+1,(g-4)/5+1,i) = cv;
  end
end
end
mean_acc = mean(accuracy_mat,3); % matriz de medias de porcentaje correcto
sd_acc = std(accuracy_mat,[],3); % matriz de desvio estandar de procentaje correcto
%% Testear Modelo
[test_label, test_data] = libsvmread('.\test_data');
model = svmtrain2(train_label, train_data, ['-t 2 -c ', num2str(bestc), ' -g ', num2str(bestg)]);
[predict_label_L, accuracy_L, values_L] = svmpredict(test_label, test_data, model);
if (model.Label(1) == -1)
    values_L = -values_L;
end
accuracy_L

% Alexander Marino 11/2017

