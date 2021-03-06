clear;
clc;
close all;

% fit = valor de fitness
% indv = indiv�duo

% Os indiv�duos s�o representados por um vetor coluna de bin�rios, que � do
% mesmo tamnho do vetor vals, ou seja, cada indiv�duo possui dimens�es de
% 1x300. Cada bit do vetor indica se aquele n�mero faz parte do grupo ou
% nao. Exemplo, se no indece i do indiv�duo o bit for igual a 1, ent�o o
% valor corresponde no vetor vals faz parte daquele indiv�duo, caso o bit
% do indece i for igual a zero, ent�o o valor correspondete em vals nao
% pertence ao indiv�duo.

% A fun��o de fitness escolhida foi f(x)=1/(1+(total-2*sum(eval(x)))), onde
% total � o somat�rio de todos os valores da parti��o e x representa um
% indiv�duo. Para resolver o problema das parti��es, n�o �, extritamente,
% necess�rio dividir os n�meros em dois grupos (A e B), basta saber se um
% grupo qualquer � capaz de se cancelar e de cancelar o outro grupo. Caso
% essa condi��o seja satisfeita o total = 0.

% total = A + B
% A - B = 0 => A = B
% total = A + A => total = 2A =>  total - 2A = 0

% Como estamos tentando busca a minima diferen�a entre os grupos, ou seja,
% zero, ent�o para que a equa��o de fitness n�o divirja � necess�rio
% adicionar 1 ao divisor. Com a equa��o de fitness podemos inferir que,
% quanto maior a diferen�a entre os dois grupos, menor ser� o valor de
% fitness, enquanto que para diferen�as pequenas o valor de fitness
% aumenta. Caso o diferen�a entre os grupos seja 0, o valor de fitness
% atinge o valor m�ximo, que � de 1.

% Dois crit�rios de parada foram estabelecidos, que s�o: caso o fitness
% atinjao valor 1, que ja indica que a melhor solu��o de todas foi
% encontrada, ou caso o n�mero de itara��es seja maior que 10.000.
% Originalmente, o espa�o de solu��o consiste de (2^300)-1 poss�veis
% combina�oes de cadeia, ent�o 10.000 itera��es pode ser visto como uma
% valor bem relevante para a busca de solu��es.


filename = 'exe21.mat';

%lendo valores
fileID = fopen('numeros_particao.txt','r');
formatSpec = '%f';
vals = fscanf(fileID,formatSpec);
fclose(fileID);

%vals=1:55;
%vals = vals(1:10);
l = length(vals);

% Par�metros 
n = 500; % N�mero de indiv�duos  
maxIt = 5000; % N�mero m�ximo de itera��es
total = sum(vals); % Soma de todos os valores
m = 0.2; % Probabilidade de uma muta��o acontecer
c = 0.8; % Probabilidade de Crossover
q = 2; % N�mero de indiv�duos para o torneio
d = 0.5; % fator de diversidade

%p = randi([0 1], n,length(vals)); % Inicializa��o rand�mica da popula��o

p = rand(n,l);
p(p>0.5)=1;
p(p<0.5)=0;

p = logical(p); 


T = 10;
contM =0; % Contador de muta��es

mFit =0; % Melhor valor de fitness
bFit = zeros(maxIt,1);
mIndv = logical(1:l); % Melhor indiv�duo;


numEval = 0;
meanP = zeros(1,maxIt);
maxP = zeros(1,maxIt);
minP = zeros(1,maxIt);





for i=1:maxIt
    i
    mp = rand(n,1); % Probabilidades de sofrer muta��es e crossOver da popula��o

    % Crossover------------------------------
    crossIdx = find(mp>(1-c)); %crossIdx = vetor de indices dos invid�vuos selecionados
    
    while(length(crossIdx)>0)
        k = randi([1,l-1],2);
        if(length(crossIdx)==1)
            crossIdx(end+1) = randi(n,1); 
        end
        a = randperm(length(crossIdx),2);
        k = sort(k);
        f = crossoverOX(k(1),k(2),p([a(1),a(2)],:));
        p(end+1,:) = f(1,:);
        p(end+1,:) = f(2,:);
        crossIdx([a(1),a(2)])=[];
    end
    %-------------------------------------
    
    
    
    %Muta��o------------------------------------ 
    mp = rand(size(p,1),1);
    mutacaoIdx = find(mp>(1-m));
    
    while(length(mutacaoIdx)>0)
        k = randi(l,1,2);
        k = sort(k);
        a = randperm(length(mutacaoIdx),1);
        %p(a,:) = MutacaoCadeiaReversa(k,p(a,:));
        p(a,:) = MutacaoPontual(k,p(a,:));
        mutacaoIdx(a)=[];
    end
    %-------------------------------------
    
    %Fitness da popula��o ------------------------------------
    fp = zeros(size(p,1),1); % fp = fitness da popula��o
    numEval = numEval + length(fp);
    for j=1:size(p,1);
        fp(j) = 1/(1+abs(total-(2*sum(vals(p(j,:) ) ) ) ) );
    end
    
    %fp
    %pause;
    %---------------------------------------
    
    % Melhor individuo ---------------------
    if(max(fp) > mFit)
        [v, maxIdx]=max(fp);
        mFit = v;
        mIndv = p(maxIdx,:);
        %mIndv
        %pause;
    end
   
    % ---------------------------------------
    
    %Sele��o por torneio -------------------------
    [np,fnp]=SelecaoTorneio( q,p,fp,n );
    % ------------------------------------------
    
    %Sele��o por rank ----------------------
    %[np,fnp]=SelecaoRank(fp,n,p);
    %----------------------------------------
   % size(unique(np,'rows'),1)
   % pause;
    %Sele��o por rank -------------------------
%     np = logical(zeros(n,length(vals))); % Nova gera��o de indiv�duos
%     [P,I]=sort(fp,'descend');
%     fnp = fp(1:n);
%     for j=1:n;
%         np(j,:) = p(I(j),:);
%     end
    %-----------------------------------------
%     if(size(unique(np,'rows'),1)/n < d)
%         mutacaoIdx = randperm(n,n-size(unique(np,'rows'),1));
%    %     size(unique(np,'rows'),1)
%         while(length(mutacaoIdx)>0)
%               a = randperm(length(mutacaoIdx),1);
%               k = randperm(l,1);
%               np(a,:) = MutacaoInvBit(k,np(a,:));
%             %np(a,:) = MutacaoPontual(k,np(a,:));
%             mutacaoIdx(a)=[];
%         end
% %         size(unique(np,'rows'),1)
% %         pause;
%     end
    
    %Verificar se o melhor indiv�duo est� na novapopula��o
    if (v>max(fnp))
        np(end+1,:)=mIndv;
        fnp(end+1)=v;
    end
    
    %Euristica
    if(mod(i,T)==0)
        for j=1:l
           aux = mIndv;
            if(aux(j)==1)
                aux(j)=0;
            else
                aux(j)=1;
            end
            nv = 1/(1+abs(total-(2*sum(vals(aux) ) ) ) );
            if (nv > v)
                v=nv;
            end
        end
        
      
        
    end
    
    
    %nnz(ismember(p,mIndv,'rows'))
    %nnz(ismember(p,mIndv,'rows'))
    % m�dia, m�ximo e m�nimo do fitness da nova popula��o 
     meanP(i) = mean(fnp);
     maxP(i) = max(fnp);
     minP(i) = min(fnp);
    if(i>1)
     if(maxP(i)<maxP(i-1))
        disp('bugg');
        break;
     end
    end
    % A nova popula��o se torna a popula��o atual
     p = np;
end

diff = (1-mFit)/mFit;

figure(1);
%plot(meanP,'g-');
    %hold on
plot(maxP,'r-');
%plot(minP,'b-');

figure(2);
plot(meanP,'g-');

figure(3)
plot(minP,'b-');

    figure(4)
    hold on
    plot(maxP,'r-');
    plot(meanP,'g-');
    plot(minP,'b-');

save(filename);


% Fun��o de  Fitness
%fit = 1/(1+abs(total-2*sum(vals(indv))))