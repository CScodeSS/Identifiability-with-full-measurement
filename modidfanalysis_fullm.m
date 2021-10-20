function [Test,NonIdenti] = modidfanalysis_fullm(varargin,moduleindex)
% Verify the generic identifiability of a subset of modules in a dynamic
% network, where the modules can be from different MISO
% subsystems
%
% [Test,NonIdenti] = subidentianaly(varargin,moduleindex)
%
% Inputs: varargin: 1-dimensional or 2-dimensional cell array, where 
%                   1. varargin{1} contains the adjacency matrix object of a network and must be specified;        
%                   2. varargin{2} contains the adjacency matrix object of the fixed modules:
%                           a. If not specified, then the defaul is that no
%                              fixed modules exist
%         moduleindex: m x 2 matrix that contains the indices of target modules, where 
%                      1. m is the total number of modules to be verfied
%                      2. Each row vector [i j] denotes a directed edge from w_i to w_j and thus the module G_{ji}  
% Outputs: Test: Logic variable, 1 if the target modules are identifiable, and 0 otherwise 
%          NonIdenti: z x 2 matrix that contains indices of non-identifiable modules; Empty if Test=1. 
%
% Reference:
% (a) S. Shi, X. Cheng and P.M.J. Van den Hof, "Generic identifiability of subnetworks in a linear dynamic network: the full measurement case",  
%     arXiv preprint arXiv:2008.01495,2020.
%
%   Author:  Shengling Shi
%            Control Systems Group
%            Eindhoven University of Technology.
%   Version: 1.1 
%   Date:    05- Aug-2021
%
%   Note: Rely on the functions "mindismaxpath" and "adj2adjm"

%% Check inputs
if ~(nargin == 2)
  error('number of input arguments is not correct')
end  
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


output= unique(moduleindex(:,2));
numMISO = length(output);     % Number of different MISO subsystems that contain the modules
input = cell(numMISO,1);
for n=1:numMISO
   int = moduleindex(:,1);
   input{n} = int(output(n)==moduleindex(:,2));       % Input indices for each MISO subsystem
end


%% verify identifiability for modules in the same MISO subsystem separately
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
        Xj = Xj+L;      % Indices of external signals, which have no unknown edge to wj, in extended graph       
        [b1,~,~]= mindismaxpath(A',Xj,Wjbar);
        [b2,~,~]= mindismaxpath(A',Xj,Wj); 
        [b3,~,~]= mindismaxpath(A',Xj,setdiff(Wj,Wjbar)); 

        if (b1==length(Wjbar)) && (b2 == b1+b3)     % Identifiability conditions of the modules in one MISO subsystem
            TestMISO(n)=1;
        else              % check which module in the MISO is not identifiable
            for k =1: length(Wjbar)
              [b1,~,~]= mindismaxpath(A',Xj,Wjbar(k));
              [b3,~,~]= mindismaxpath(A',Xj,setdiff(Wj,Wjbar(k))); 
               if (b1==1) && (b2 == b1+b3)           % Single module identifiability test
               else                      
                   NonIdenti = [NonIdenti ; Wjbar(k)  wj];          
               end
            end
        end
    else
        TestMISO(n)=1;  
    end
end
Test = all(TestMISO);
end

