%PE

clc;
clear;
close all;

filename = 'exePE6.mat';

fileID = fopen('cidades.txt','r');
C = textscan(fileID,'%f %f');
fclose(fileID);

cityCoords = [C{1} C{2}]; % Coordenadas das cidades
numcity = size(cityCoords,1); % N�mero de cidades
pop_size = 90; % Tamanho da popula��o
pop = zeros(pop_size, numcity); 
maxIt = 100000;
numEval = 0; % N�mero de avalia��es
mFit = 5e10; %melhor fitness;
q = 2; % N�mero de indiv�duos no torneio
evalLimit = 500000;

% Inicializando popula��o
for i=1:pop_size
    pop(i,:) = randperm(numcity);
end

for i=1:maxIt
    i
    
    % Muta��o -----------------------------------------------------------
    for j=1:pop_size
        k = randi(numcity,1,2);
        k = sort(k);
        pop(end+1,:) = MutacaoCadeiaReversa(k,pop(j,:));
        %pop(end,k(1)) = pop(j,k(2));
        %pop(end,k(2)) = pop(j,k(1));
    end
    %--------------------------------------------------------------------
    
    %Fitness da popula��o -----------------------------------------------
    fp = zeros(size(pop,1),1); % fp = fitness da popula��o
    numEval = numEval + length(fp);
    for j=1:size(pop,1);
        fp(j) = TotalDistance(cityCoords(pop(j,:),:));
    end
   %--------------------------------------------------------------------
   
    %Sele��o por torneio ------------------------------------------------
    [np,fnp]=TorneioPE( q,pop,fp,pop_size );
    % -------------------------------------------------------------------
   
   
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
legend('Fitness m�dio', 'Melhor fitness');

figure(5)
hold on
plot(meanP(1:500),'r-');
plot(minP(1:500),'b-');
legend('Fitness m�dio', 'Melhor fitness');


save(filename);

