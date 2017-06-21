close all;
clc;
clear;

filename = 'exe1.mat';

%lendo valores
file = load('knapsack.mat');
p = file.p;
v = file.v;
c = file.c;
opt = file.opt;

%representação das soluções - cadeias binárias
n = 100; % número de indivíduos
cr = 0.8; % probabilidade de crossover 
m = 0.2; % probabilidade de mutação
maxIt = 1000; % número máximo de iterações
opt = logical(opt);

l = length(v);

mFit =0; % Melhor valor de fitness
mIndv = logical(1:l); % Melhor indivíduo;

po = randi([0 1], n,l);
po = logical(po);

q = 2;
numEval = 0;


for i=1:maxIt
    i
    %crossOver ----------------------
    pc = rand(n,1);
    crossIdx = find(pc>(1-cr));
    
    while(length(crossIdx)>0)
        k = randi([1,l-1],2);
        if(length(crossIdx)==1)
            crossIdx(end+1) = randi(n,1); 
        end
        a = randperm(length(crossIdx),2);
        k = sort(k);
        f = crossoverOX(k(1),k(2),po([a(1),a(2)],:));
        po(end+1,:) = f(1,:);
        po(end+1,:) = f(2,:);
        crossIdx([a(1),a(2)])=[];
    end
    %-----------------------------------
    
    %Mutação------------------------------------ 
    pm = rand(n,1);
    mutacaoIdx = find(pm>(1-m));
    
    while(length(mutacaoIdx)>0)
        k = randi(l,1,2);
        k = sort(k);
        a = randperm(length(mutacaoIdx),1);
        %po(a,:) = MutacaoCadeiaReversa(k,po(a,:));
        po(a,:) = MutacaoPontual(k,po(a,:));
        mutacaoIdx(a)=[];
    end
    %-------------------------------------
    
     %Fitness da população ------------------------------------
     fp = zeros(size(po,1),1); % fp = fitness da população
     numEval = numEval + length(fp);
     for j=1:size(po,1);
        if(sum(p(po(j,:)))>c);
            fp(j)=0;
        else
            fp(j) = sum(v(po(j,:)));
        end
     end
     %----------------------------------------------------------
    % Melhor individuo ---------------------
    if(max(fp) > mFit)
        [mv, maxIdx]=max(fp);
        mFit = mv;
        mIndv = po(maxIdx,:);
        
    end
    % ---------------------------------------
    
    %Seleção por torneio -------------------------
    [np,fnp]=SelecaoTorneio( q,po,fp,n );
    % ------------------------------------------
    
    % Salvaciosnismo -----------------------------
    if (mv>max(fnp))
        disp('entrou');
        np(j+1,:)=mIndv;
        fnp(end+1)=mv;
    end
    % --------------------------------------------
    
    % médio, máximo e mínimo fitness ---------------
    meanP(i) = mean(fnp);
    maxP(i) = max(fnp);
    minP(i) = min(fnp);
    optval(i) = sum(v(opt));
    if(i>1)
     if(maxP(i)<maxP(i-1));
        disp('bugg');
        break;
     end
    end
    %----------------------------------------------
    
    po = np;
    
    
    if(mFit >= sum(v(logical(opt))))
    disp('melhor resposta');
    break;
    end
    
    if(numEval>=100000)
        break;
    end
    
end

numEval = 100000;
for j=1:20
    mr = 0;
for i=1:numEval
   fit = logical(randi([0 1],1,24)); 
   if(sum(p(fit))>c);
        fp =0;
    else
        fp = sum(v(fit));
   end
   
   if(fp>mr)
       mr = fp;
   end
   
end
    obt(j)=mr;
end

% figure(1);
% plot(maxP,'r-');
% 
% figure(2);
% plot(meanP,'g-');
% 
% figure(3)
% plot(minP,'b-');
% 
% figure(4)
% hold on
% plot(optval,'k-');
% plot(maxP,'r-');
% plot(meanP,'g-');
% plot(minP,'b-');

save(filename);

