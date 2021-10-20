clear all
clc

%% Examples
adj.G = [0 1 0 0 0 0 0 0 0; %1
               0 0 1 0 0 0 0 0 0; %2
               0 1 0 0 0 0 0 0 0; %3
               1 0 0 0 1 0 0 1 0; %4
               0 1 0 0 0 1 0 1 0; %5
               0 0 1 0 0 0 0 0 1; %6
               0 0 0 0 0 1 0 0 0; %7 
               0 0 0 0 0 1 0 0 1; %8
               0 0 0 0 1 0 0 0 0;];%9
adj.H=[0 0; %1
               0 0; %2
               0 0; %3
               0 0; %4
               0 0; %5
               1 1; %6
               1 1; %7 
               0 0; %8
               0 0;];%9];  
           
 adj.R =[];
 Gmatrix = adj.G;
varargin{1} =  adj;
  L = length(adj.G);

targetmodule=[]; % Synthesis for full network identifiability, thus all modules are target modules
for i=1:L
     inputi=find(Gmatrix(i,:));
     for j=1:length(inputi)
        targetmodule = vertcat(targetmodule,[inputi(j) i]);
     end
 end
 

%% Signal allocation
 [vararginOut] =  modidfsynthesis_fullm(varargin,targetmodule,1);

[TestOld,NonIdentiOld] = modidfanalysis_fullm(varargin,targetmodule); % Verify identifiability of the old set
[Test,NonIdenti] = modidfanalysis_fullm(vararginOut,targetmodule);    % Verify identifiability of the new set



%% Plots
% Plot the original model set
figure
ExtendGraph = adj2adjm(adj);
Graph=digraph((ExtendGraph)',{'w1','w2','w3','w4','w5','w6','w7','w8','w9','e1','e2'});
H=plot(Graph,'Layout','force');
% Plot the new model set after signal allocation
ExtendGraphNew = adj2adjm(vararginOut{1});
GraphNew=digraph((ExtendGraphNew)',{'w1','w2','w3','w4','w5','w6','w7','w8','w9','e1','e2','r1','r2'});
figure
Hnew=plot(GraphNew,'Layout','force');
