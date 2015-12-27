% Inputs : 
% M : weight matrix (the matrix is symetrized with
% the sum of weights in both directions)
% s : 1 = Recursive computation
%   : 0 = Just one level computation
% self : 1 = Use self weights
%        0 = Do not use self weights
% debug   : 1 = outputs some debug messages
% verbose : 1 = outputs some messages
%
% Output :
% COMTY, structure with the following information
% for each level i :
%   COMTY.COM{i} : vector of community IDs (sorted by community sizes)
%   COMTY.SIZE{i} : vector of community sizes
%   COMTY.MOD(i) : modularity of clustering
%   COMTY.Niter(i) : Number of iteration before convergence
%
function [COMTY ending] = my_cluster_algorithm(M,s,self,debug,verbose)

if nargin < 1
  error('not enough argument');
end

if nargin < 2
  s = 1;
end

if nargin < 3
  self = 1;
end

if nargin < 4
  debug = 0;
end

if nargin < 5
  verbose = 0;
end

S = size(M);
N = S(1);

ddebug = 0;
ending = 0;

M = M + M';
if (self == 0)
  M((N+1).*[0:N-1]+1) = 0;
end
M2 = M;
M2((N+1).*[0:N-1]+1) = 0;

m = sum(sum(M));
Niter = 1;

if m==0 | N == 1
  fprintf('No more possible decomposition\n');
  ending = 1;
  COMTY = 0;
  return;
end

% Main loop
K = sum(M); % Sum of wieght incident to node i
SumTot = sum(M);
SumIn = diag(M); % Sum of weight inside community i
COM = 1:S(1); % Community of node i
for k=1:N
  Neighbor{k} = find(M(k,:));
end

sCost = 10;
gain = 1;
while (gain == 1)
  Cost = zeros(1,N);
  gain = 0;
  for i=1:N
    Ci = COM(i);
    NB = Neighbor{i};
    G = zeros(1,N); % Gain vector
    best_increase = -1;
    Cnew = Ci;
    COM(i) = -1;
    SumTot(Ci) = SumTot(Ci) - K(i);
    CNj1 = find(COM==Ci);
    SumIn(Ci) = SumIn(Ci) - 2*sum(M(i,CNj1)) - M(i,i);
    for j=1:length(NB)
      Cj = COM(NB(j));
      if (G(Cj) == 0)
        CNj = find(COM==Cj);
        Ki_in = 2*sum(M(i,CNj));
        G(Cj) = Ki_in/m - 2*K(i)*SumTot(Cj)/(m*m);
        if (ddebug)
          fprintf('Gaim for comm %d => %g\n',Cj-1,G(Cj));
        end
        if G(Cj) > best_increase;
          best_increase = G(Cj);
          Cnew_t = Cj;
        end
      end
    end
    if best_increase > 0
      Cnew = Cnew_t;
      if (debug)
        fprintf('Move %d => %d\n',i-1,Cnew-1);
      end
      Cost(i) = best_increase;
    end
    Ck = find(COM==Cnew);
    SumIn(Cnew) = SumIn(Cnew) + 2*sum(M(i,Ck));
    SumTot(Cnew) = SumTot(Cnew) + K(i);
    COM(i) = Cnew;
    if (Cnew ~= Ci)
      gain = 1;
    end
    
  end
  Niter = Niter + 1;
end

Niter = Niter - 1;


COMu = unique(COM);

ClusterAdj=M;
for i=1:length(COMu)
    icluster =setdiff(find(COM==COMu(i)),COMu(i));
    for j =1:length(icluster)
        ClusterAdj(COMu(i),:)=ClusterAdj(COMu(i),:)+ClusterAdj(icluster(j),:);
        ClusterAdj(:,COMu(i))=ClusterAdj(:,COMu(i))+ClusterAdj(:,icluster(j));
     
    end
end


ClusterAdj = subgraph(ClusterAdj,COMu);


ClusterAdj = symmetrize(ClusterAdj);

Degree = degrees(ClusterAdj);

AverageDegree = averageDegree(ClusterAdj);

[Sorteddegree,Index] = sort(Degree,'descend');
p = std(Degree);
S = size(ClusterAdj);
ClusterN = S(1);
i =1;
Seeds=[];
while Degree(Index(i))>AverageDegree+p
    Seeds =[Seeds,Index(i)];
    i=i+1;
end

for i=1:ClusterN
    ClusterNeighbor{i} = setdiff(find(ClusterAdj(i,:)),i);
end

Seeds_Length = length(Seeds);

if(Seeds_Length==0)
    Seeds =[Index(1)];
    Seeds_Length = 1;
end
%Seeds
ClusterK = sum(ClusterAdj); % Sum of wieght incident to node i
ClusterSumTot = sum(ClusterAdj); %Sum of the weights of the links incident to nodes in community i
ClusterSumIn = diag(ClusterAdj); % Sum of weight inside community i

% init the nodes community
ClusterCOM = [];
i=1;
while (i<=ClusterN)
   ClusterCOM=[ClusterCOM,-1];
   i=i+1;
end

for i=1:Seeds_Length
    ClusterCOM(Seeds(i)) = Seeds(i);
end

ClusterDoneN=0;
Clustergain=1;
loop=0;
while (Clustergain==1)
    ClusterDoneN=0;
    Clustergain = 0;
    loop=loop+1;
    for i = 1 :ClusterN
       if(debug==1)
            fprintf('loop %d\n',loop);
       end
       if(length(find(Seeds==i))>0)
            continue;
       end
       Ci = ClusterCOM(i);
       NB = ClusterNeighbor{i};
       best_increase = -1;
       Cnew = Ci;
       G = zeros(1,ClusterN); % Gain vector
             
       if(Ci~=-1)
            ClusterCOM(i) = -1;
            ClusterSumTot(Ci) = ClusterSumTot(Ci) - ClusterK(i);
            CNj1 = find(ClusterCOM==Ci);
            ClusterSumIn(Ci) = ClusterSumIn(Ci) - 2*sum(ClusterAdj(i,CNj1)) - ClusterAdj(i,i);
            if(debug==1)
                fprintf('node index of %d has already in cluster %d\n',i,Ci);
             end
        end
             
        for j=1:length(NB)
             Cj = ClusterCOM(NB(j));
             if(Cj==-1)
                 if(debug==1)
                     fprintf('this neighbor %d did not belongs to a cluster\n',NB(j));
                 end
                 continue;
             end
             if(G(Cj)==0)
                 CNj = find(ClusterCOM==Cj);
                 Ki_in = 2*sum(ClusterAdj(i,CNj));
                 G(Cj) = Ki_in/m -2* ClusterK(i)*ClusterSumTot(Cj)/(m*m);
                 if (G(Cj) > best_increase)
                     best_increase = G(Cj);
                     Cnew_t = Cj;
                 end
             end
          end
          if best_increase >0
               Cnew = Cnew_t;
          end
             
           if(Cnew~=-1)
                 Ck = find(ClusterCOM==Cnew);
                 ClusterSumIn(Cnew) = ClusterSumIn(Cnew) + 2*sum(ClusterAdj(i,Ck));
                 ClusterSumTot(Cnew) = ClusterSumTot(Cnew) + ClusterK(i);
                 ClusterCOM(i) = Cnew;
            end
            if (Cnew ~= Ci)
               Clustergain = 1;
            else
                 ClusterDoneN=ClusterDoneN+1;
            end
      end
         
      if(debug ==1)
           fprintf('node %d is done\n',i);
      end
end

ClusterCOMu=unique(ClusterCOM);

for i =1:length(ClusterCOMu)
    if(ClusterCOMu(i)~=-1)
        Combin=find(ClusterCOM==ClusterCOMu(i));
        for j=1:length(Combin)
            COM(find(COM==COMu(Combin(j))))=COMu(Combin(1));
        end
    end
end



[COM COMSIZE] = reindex_com(COM);
COMTY.COM{1} = COM;
COMTY.SIZE{1} = COMSIZE;
COMTY.MOD(1) = compute_modularity(COM,M);
COMTY.AVGSIZE = mean(COMSIZE);
COMTY.Niter(1) = Niter;


end

 



% Re-index community IDs
function [C Ss] = reindex_com(COMold)

C = zeros(1,length(COMold));
COMu = unique(COMold);
S = zeros(1,length(COMu));
for l=1:length(COMu)
    S(l) = length(COMold(COMold==COMu(l)));
end
[Ss INDs] = sort(S,'descend');

for l=1:length(COMu)
    C(COMold==COMu(INDs(l))) = l;
end
end

%Compute modulartiy
function MOD = compute_modularity(C,Mat)

m = sum(sum(Mat));
MOD = 0;
COMu = unique(C);
for j=1:length(COMu)
    Cj = find(C==COMu(j));
    Ec = sum(sum(Mat(Cj,Cj)));
    Et = sum(sum(Mat(Cj,:)));
    if Et>0
        MOD = MOD + Ec/m-(Et/m)^2;
    end
end

end


