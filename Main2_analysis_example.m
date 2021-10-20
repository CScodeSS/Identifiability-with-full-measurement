clear all
clc
% Example of testing generic identifiability of a set of target modules

% Basic paramters
adj.G = [0 0 0 0; %1
         1 0 0 0; %2
         0 1 0 0; %3
         1 1 1 0;]; %4
adj.H = [1;0;0;0];
adj.R = [0;2;0;0];
varargin{1}=adj;

%% Fixed module case with G_{42} known 

% adjfix.G = [0 0 0 0; %1
%          0 0 0 0; %2
%          0 0 0 0; %3
%          0 1 0 0;]; %4
% adjfix.H = [0;0;0;0];
% adjfix.R = [0;0;0;0];
% varargin{2}=adjfix;

%%  Conduct verification
targetmodule=[1 4;3 4;1 2];
[Test,NonIdenti] = modidfanalysis_fullm(varargin,targetmodule);

% Plot the results
ExtendGraph = adj2adjm(adj);
Graph=digraph((ExtendGraph)',{'w1','w2','w3','w4','e','r'});
H=plot(Graph,'Layout','force');
%layout(H)
for i=1:size(NonIdenti,1)
    index = NonIdenti(i,:);
    highlight(H,[index(1) index(2)],'EdgeColor','r')    % Non-identifiable target modules are marked in red
end
if ~isempty(NonIdenti)
    identifiablemodule = setdiff(targetmodule,NonIdenti,'rows');
else
    identifiablemodule = targetmodule;
end
    for j=1:size(identifiablemodule,1)
        index2= identifiablemodule(j,:);
      highlight(H,[index2(1) index2(2)],'EdgeColor','g')    % identifiable target modules are marked in green
    end