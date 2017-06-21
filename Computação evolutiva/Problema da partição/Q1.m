clear;
clc;
close all;

% fit = valor de fitness
% indv = indivíduo

% Os indivíduos são representados por um vetor coluna de binários, que é do
% mesmo tamnho do vetor vals, ou seja, cada indivíduo possui dimensões de
% 1x300. Cada bit do vetor indica se aquele número faz parte do grupo ou
% nao. Exemplo, se no indece i do indivíduo o bit for igual a 1, então o
% valor corresponde no vetor vals faz parte daquele indivíduo, caso o bit
% do indece i for igual a zero, então o valor correspondete em vals nao
% pertence ao indivíduo.

% A função de fitness escolhida foi f(x)=1/(1+(total-2*sum(eval(x)))), onde
% total é o somatório de todos os valores da partição e x representa um
% indivíduo. Para resolver o problema das partições, não é, extritamente,
% necessário dividir os números em dois grupos (A e B), basta saber se um
% grupo qualquer é capaz de se cancelar e de cancelar o outro grupo. Caso
% essa condição seja satisfeita o total = 0.

% total = A + B
% A - B = 0 => A = B
% total = A + A => total = 2A =>  total - 2A = 0

% Como estamos tentando busca a minima diferença entre os grupos, ou seja,
% zero, então para que a equação de fitness não divirja é necessário
% adicionar 1 ao divisor. Com a equação de fitness podemos inferir que,
% quanto maior a diferença entre os dois grupos, menor será o valor de
% fitness, enquanto que para diferenças pequenas o valor de fitness
% aumenta. Caso o diferença entre os grupos seja 0, o valor de fitness
% atinge o valor máximo, que é de 1.

% Dois critérios de parada foram estabelecidos, que são: caso o fitness
% atinjao valor 1, que ja indica que a melhor solução de todas foi
% encontrada, ou caso o número de itarações seja maior que 10.000.
% Originalmente, o espaço de solução consiste de (2^300)-1 possíveis
% combinaçoes de cadeia, então 10.000 iterações pode ser visto como uma
% valor bem relevante para a busca de soluções.


filename = 'exe21.mat';

%lendo valores
fileID = fopen('numeros_particao.txt','r');
formatSpec = '%f';
vals = fscanf(fileID,formatSpec);
fclose(fileID);

%vals=1:55;
%vals = vals(1:10);
l = length(vals);

% Parâmetros 
n = 500; % Número de indivíduos  
maxIt = 5000; % Número máximo de iterações
total = sum(vals); % Soma de todos os valores
m = 0.2; % Probabilidade de uma mutação acontecer
c = 0.8; % Probabilidade de Crossover
q = 2; % Número de indivíduos para o torneio
d = 0.5; % fator de diversidade

%p = randi([0 1], n,length(vals)); % Inicialização randômica da população

p = rand(n,l);
p(p>0.5)=1;
p(p<0.5)=0;

p = logical(p); 


T = 10;
contM =0; % Contador de mutações

mFit =0; % Melhor valor de fitness
bFit = zeros(maxIt,1);
mIndv = logical(1:l); % Melhor indivíduo;


numEval = 0;
meanP = zeros(1,maxIt);
maxP = zeros(1,maxIt);
minP = zeros(1,maxIt);





for i=1:maxIt
    i
    mp = rand(n,1); % Probabilidades de sofrer mutações e crossOver da população

    % Crossover------------------------------
    crossIdx = find(mp>(1-c)); %crossIdx = vetor de indices dos invidívuos selecionados
    
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
    
    
    
    %Mutação------------------------------------ 
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
    
    %Fitness da população ------------------------------------
    fp = zeros(size(p,1),1); % fp = fitness da população
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
    
    %Seleção por torneio -------------------------
    [np,fnp]=SelecaoTorneio( q,p,fp,n );
    % ------------------------------------------
    
    %Seleção por rank ----------------------
    %[np,fnp]=SelecaoRank(fp,n,p);
    %----------------------------------------
   % size(unique(np,'rows'),1)
   % pause;
    %Seleção por rank -------------------------
%     np = logical(zeros(n,length(vals))); % Nova geração de indivíduos
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
    
    %Verificar se o melhor indivíduo está na novapopulação
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
    % média, máximo e mínimo do fitness da nova população 
     meanP(i) = mean(fnp);
     maxP(i) = max(fnp);
     minP(i) = min(fnp);
    if(i>1)
     if(maxP(i)<maxP(i-1))
        disp('bugg');
        break;
     end
    end
    % A nova população se torna a população atual
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


% Função de  Fitness
%fit = 1/(1+abs(total-2*sum(vals(indv))))