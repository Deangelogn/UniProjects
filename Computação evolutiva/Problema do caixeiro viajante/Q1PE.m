%PE

clc;
clear;
close all;

filename = 'exePE6.mat';

fileID = fopen('cidades.txt','r');
C = textscan(fileID,'%f %f');
fclose(fileID);

cityCoords = [C{1} C{2}]; % Coordenadas das cidades
numcity = size(cityCoords,1); % Número de cidades
pop_size = 90; % Tamanho da população
pop = zeros(pop_size, numcity); 
maxIt = 100000;
numEval = 0; % Número de avaliações
mFit = 5e10; %melhor fitness;
q = 2; % Número de indivíduos no torneio
evalLimit = 500000;

% Inicializando população
for i=1:pop_size
    pop(i,:) = randperm(numcity);
end

for i=1:maxIt
    i
    
    % Mutação -----------------------------------------------------------
    for j=1:pop_size
        k = randi(numcity,1,2);
        k = sort(k);
        pop(end+1,:) = MutacaoCadeiaReversa(k,pop(j,:));
        %pop(end,k(1)) = pop(j,k(2));
        %pop(end,k(2)) = pop(j,k(1));
    end
    %--------------------------------------------------------------------
    
    %Fitness da população -----------------------------------------------
    fp = zeros(size(pop,1),1); % fp = fitness da população
    numEval = numEval + length(fp);
    for j=1:size(pop,1);
        fp(j) = TotalDistance(cityCoords(pop(j,:),:));
    end
   %--------------------------------------------------------------------
   
    %Seleção por torneio ------------------------------------------------
    [np,fnp]=TorneioPE( q,pop,fp,pop_size );
    % -------------------------------------------------------------------
   
   
    % média, máximo e mínimo do fitness da nova população----------------
     meanP(i) = mean(fnp);
     maxP(i) = max(fnp);
     minP(i) = min(fnp);
    % -------------------------------------------------------------------
    
     if(i>1)
     if(minP(i-1)<minP(i))
        disp('bugg');
        break;
     end
    end
    
    if(numEval>evalLimit)
        break;
    end
    
     % A nova população se torna a população atual------------------------
     pop = np;
    %---------------------------------------------------------------------
    
end

[v, minIdx]=min(fnp);
best = cityCoords(pop(minIdx,:),:);
best(end+1,:) = best(1,:);

figure (1);
plot(C{1}, C{2},'*');
grid on
hold on
plot(best(:,1), best(:,2));
start = scatter(best(1,1),best(1,2),'g');
leng = [start]
legend(leng,'inicio');

figure(2);
plot(maxP,'r-');

figure(3);
plot(meanP,'g-');

figure(4)
plot(minP,'b-');

figure(5)
hold on
%plot(maxP,'r-');
plot(meanP,'r-');
plot(minP,'b-');
xlim([0 i]);
%ylim([0 meanP(1)]);
legend('Fitness médio', 'Melhor fitness');

figure(5)
hold on
plot(meanP(1:500),'r-');
plot(minP(1:500),'b-');
legend('Fitness médio', 'Melhor fitness');


save(filename);

