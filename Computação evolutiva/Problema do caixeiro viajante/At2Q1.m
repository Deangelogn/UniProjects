    %At2Q1

clc;
clear;
close all;

filename = 'exe1.mat';

fileID = fopen('cidades.txt','r');
C = textscan(fileID,'%f %f');
fclose(fileID);

cityCoords = [C{1} C{2}]; % Coordenadas das cidades
numcity = size(cityCoords,1); % N�mero de cidades
pop_size = 100; % Tamanho da popula��o
pop = zeros(pop_size, numcity); 
crossProb = 0.8; % probabilidade de crossover
mutProb = 0.2; % probabilidade de muta��o
maxIt = 500000;
numEval = 0; % N�mero de avalia��es
mFit = 5e10; %melhor fitness;
q = 2; % N�mero de indiv�duos no torneio
evalLimit = 500000;


% Inicializando popula��o
for i=1:pop_size
    pop(i,:) = randperm(numcity);
end

for i=1:maxIt
    % Crossover -----------------------------------------------------------
    popProb = rand(size(pop,1),1); % probabilidade da popula��o;
    crossIdx = find(popProb>(1-crossProb)); % vetor com os indices dos indiv�duos que v�o para o crossover
    
    while(length(crossIdx)>0)
        if(length(crossIdx)==1)
            crossIdx(end+1) = randi(size(pop,1),1); 
        end
        k = randi(numcity,1,2);
        a = randperm(length(crossIdx),2);
        k = sort(k);
        [s1 s2] = CrossoverPosition(pop(a(1),:), pop(a(2),:),k);
        pop(end+1,:) = s1;
        pop(end+1,:) = s2;
        crossIdx([a(1),a(2)])=[];
    end
    %---------------------------------------------------------------------
    
    % Muta��o -----------------------------------------------------------
    popProb = rand(size(pop,1),1);
    mutacaoIdx = find(popProb>(1-mutProb));
    while(length(mutacaoIdx)>0)
        k = randi(numcity,1,2);
        k = sort(k);
        a = randperm(length(mutacaoIdx),1);
        pop(a,:) = MutacaoCadeiaReversa(k,pop(a,:));
        mutacaoIdx(a)=[];
    end
    %--------------------------------------------------------------------
    
    %Fitness da popula��o -----------------------------------------------
    fp = zeros(size(pop,1),1); % fp = fitness da popula��o
    numEval = numEval + length(fp);
    for j=1:size(pop,1);
        fp(j) = TotalDistance(cityCoords(pop(j,:),:));
    end
    %--------------------------------------------------------------------
    
    % Melhor individuo --------------------------------------------------
    if(min(fp) < mFit)
        [v, minIdx]=min(fp);
        mFit = v;
        mIndv = pop(minIdx,:);
        %mIndv
        %pause;
    end
    % -------------------------------------------------------------------
    
    %Sele��o por torneio ------------------------------------------------
    [np,fnp]=SelecaoTorneio( q,pop,fp,pop_size );
    % -------------------------------------------------------------------
    
    %Verificar se o melhor indiv�duo est� na nova popula��o---------------
    if (v<min(fnp))
        np(end+1,:)=mIndv;
        fnp(end+1)=v;
    end
    %---------------------------------------------------------------------
    
    % m�dia, m�ximo e m�nimo do fitness da nova popula��o----------------
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
    
     % A nova popula��o se torna a popula��o atual------------------------
     pop = np;
    %---------------------------------------------------------------------
end

best = cityCoords(mIndv,:);
best(end+1,:) = best(1,:);
figure (1);
plot(C{1}, C{2},'*');
grid on
hold on
plot(best(:,1), best(:,2));
start = scatter(best(1,1),best(1,2),'g');
leng = [start];
legend(leng,'origem');

figure(2);
plot(maxP,'r-');

figure(3);
plot(meanP,'g-');

figure(4)
plot(minP,'b-');

figure(5)
hold on
%plot(maxP,'r-');
plot(meanP,'g-');
plot(minP,'b-');
xlim([0 i]);
%ylim([0 meanP(1)]);
legend('Fitness m�dio', 'Melhor fitness');

save(filename);



