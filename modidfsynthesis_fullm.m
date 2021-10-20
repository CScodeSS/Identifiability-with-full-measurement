function [vararginOut] = modidfsynthesis_fullm(varargin,moduleindex,fixr)
% [vararginOut] = modidfsynthesis_fullm(varargin,moduleindex,fixr)
%
% Allocate additional actuators to achieve generic identifiability of a
% subset of modules, where the modules can be from different MISO
% subsystems, for the sitution of full measurement. 
%
% Inputs: varargin: 1-dimensional or 2-dimensional cell array, where 
%                   1. varargin{1} contains the adjacency matrix object of a network and must be specified;        
%                   2. varargin{2} contains the adjacency matrix object of the fixed modules:
%                           a. If not specified, then the defaul is that no
%                              fixed modules exist
%         moduleindex: m x 2 matrix that contains the indices of target modules, where 
%                      1. m is the total number of target modules
%                      2. Each row vector [i j] denotes a directed edge from w_i to w_j and thus the module G_{ji}  
%         fixr : Binary variable:
%                a. fixr=1 if all the transfer functions in R matrix corresponding to the allocated r
%                   signals are known
%                b. fixr = 0 (default if not specified): The above transfer
%                   functions are unknown, but only known if one r is directly
%                   allocated at the considered output.
% Outputs: vararginOut: The outputed network with the allocated r signals;
%                       The same format of the input varargin    
%
% Reference:
% (a) S. Shi, X. Cheng and P.M.J. Van den Hof, "Generic identifiability of subnetworks in a linear dynamic network: the full measurement case",  
%     arXiv preprint arXiv:2008.01495,2020.
%
%   Author:  Shengling Shi
%            Control Systems Group
%            Eindhoven University of Technology.
%   Version: 1.1 
%   Date:    26- Aug-2021; name changed by Paul, 2-9-2021; change counting
%   by Shi, 3-9-2021
%

%% Check inputs
 if ~(size(moduleindex,2) == 2)
  error('the dimension of the input moduleindex is not correct')
end  
adjobj = varargin{1};   
testmatrix = adjobj.G;    
for n=1:size(moduleindex,1)
    testindex = moduleindex(n,:);
    if ~(testmatrix(testindex(2),testindex(1))==1)
          error('the specified target modules do not exist in the network')  
    end
end
if ~exist('fixr','var') || isempty(fixr)
     fixr=0; 
end


%% Basic build up
% Build the extended graph of the network

[A,~] = adj2adjm(adjobj);
L = size(adjobj.G,1); % number of nodes
K = size(adjobj.R,2); % number of external excitations
p = size(adjobj.H,2); % number of white noises

if length(varargin)<2
    adjobjfix.G=zeros(L,L);
    adjobjfix.H=zeros(L,p);
    adjobjfix.R=zeros(L,K);
else
    adjobjfix = varargin{2};
end

% Check the modules belong to what MISO subsystems
output= unique(moduleindex(:,2));


%%


numMISO = length(output);     % Number of different MISO subsystems that contain the target modules
input = cell(numMISO,1);
for n=1:numMISO
%    ind = find(output(i)==moduleindex(:,2)); 
     int = moduleindex(:,1);
%    input{i} = int(ind);       % Input indices for each MISO subsystem
     input{n} = int(output(n)==moduleindex(:,2));
end
GraphExten= digraph(A');


%% Signal allocation for modules in the same MISO subsystem separately
Gmax = adjobj.G;
Gmaxfix = adjobjfix.G;
Rmax = adjobj.R;
Rmaxfix = adjobjfix.R;
Hmax = adjobj.H;
Hmaxfix = adjobjfix.H;
RHmatrix = [Hmax Rmax]+[Hmaxfix Rmaxfix];
TestMISO = zeros(numMISO,1);
NonIdenti=[];
for n = 1:numMISO
    wj=output(n);
    inneighbor = find(Gmax(wj,:));
    inneighborfix = find(Gmaxfix(wj,:));
    Wj = setdiff(inneighbor,inneighborfix );
    Wjbar = setdiff(input{n},inneighborfix ); 
    if ~isempty(Wjbar)
        X=1:1:K+p;
        Xex = find(RHmatrix(wj,:)==1);
        Xj = setdiff( X,Xex);      
        Xj = Xj+L;      % Indices of intially present external signals, which have no unknown edge to wj, in extended graph 
        
        % Step 1 of Algorithm 1 in the paper
        NplusWjbar=[];
       for kk0=1:length(Wjbar) 
           NplusWjbar= union(NplusWjbar,successors(GraphExten,Wjbar(kk0)));
       end
       adjobj.R=Rmax;
       [A,~] = adj2adjm(adjobj); 
        [~,D,~] = mindismaxpath(A',union(Xj,NplusWjbar),setdiff(Wj,Wjbar));
        % Step 2
         PP=union(D,Wjbar);
        [b,~,P]= mindismaxpath(A',Xj,PP);
        % Step 3
        if b< length(PP)
        GraphOnly=digraph(P);
        NodeUnExcited = PP;
            for kk1=1:length(Xj)
                node1= Xj(kk1);
                for kk2=1:length(PP)
                    node2 = PP(kk2);
                    paths=allpaths(GraphOnly,node1,node2);
                    if ~isempty(paths)
                        if length(paths)>1
                           error('Multiple paths')
                        else
                           pathvec=paths{:};
                           endnode=pathvec(end);
                           NodeUnExcited = setdiff(NodeUnExcited,endnode);
                        end
                    end
                end
            end
           RnewColumns=zeros(L,length(NodeUnExcited)); 
           RnewFix =[];
           for kk3=1:length(NodeUnExcited)
               RnewColumns(NodeUnExcited(kk3),kk3)=1;
               if  NodeUnExcited(kk3)== wj || fixr==1
                   RnewFix =[RnewFix RnewColumns(:,kk3)];
               else
                   RnewFix = [RnewFix zeros(L,1)];
               end
           end
            Rmax= [Rmax RnewColumns];
            if ~isempty(RnewFix)
               Rmaxfix=  [Rmaxfix RnewFix];
            end
         RHmatrix = [Hmax Rmax]+[Hmaxfix Rmaxfix];   
        end
    end
    K = size(Rmax,2);
end
adjobj.R=Rmax;
vararginOut{1} = adjobj;
if ~isempty(find(Rmaxfix,1)) % replace by any(Rmaxfix(:));
    adjobjfix.R = Rmaxfix;
    vararginOut{2}=adjobjfix;
end

end

