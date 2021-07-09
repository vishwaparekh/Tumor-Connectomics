function tcf  =  getTCF(img ,mask, t)
% img: Input image (2D for a single image, 3D for multiparametric MRI)
% mask: mask for the roi
% t = threshold for connectomics analysis

if length(size(img))==2
    
    XY_norm = img(:);
    
elseif length(size(img))==3
    XY_norm = zeros(size(img,1)*size(img,2),size(img,3));
    for i = 1:size(img,3)
        XY_norm(:,i)=normalizeImage(img);
    end
end

indROI = mask(:);




kIso = round(t*100);
X_norm = XY_norm(indROI,:);
% Distance matrix
C1 = squareform(pdist(X_norm(:,:)));
C1 = C1./sqrt(size(X_norm,2));
C1(isnan(C1)) = 0;

% Threshold the image for computing graph metrics
D1 = C1<=t;
D1 = min(D1,D1');

% Create Graph
G = graph(D1);

% Compute Node Strength
nodeStrength = mean(C1);    
Z = X_norm(:,1);
Z(indROI) = normalizeImage(nodeStrength);
Z(~indROI) = 0;
tcf.nodeStrength = reshape(Z,[size(img,1), size(img,2)]);



% Compute Degree Centrality
degCent = sum(D1);
degCent = degCent./length(degCent); %Normalize by number of nodes
Z = X_norm(:,1);
Z(indROI) = normalizeImage(degCent);
Z(~indROI) = 0;
tcf.degreeCentrality = reshape(Z,[size(img,1), size(img,2)]);


% Compute Betweenness Centrality
nNodes = length(degCent);
bC = centrality(G,'betweenness');
bC = bC./((nNodes-1)*(nNodes-2)/2);
Z = X_norm(:,1);
Z(indROI) = normalizeImage(bC);
Z(~indROI) = 0;
tcf.betweennessCentrality = reshape(Z,[size(img,1), size(img,2)]);

% Get connected components in the graph
B = conncomp(G);
tcf.NumConnComponents=B;

% Compute Eigenvector Centrality
eC = centrality(G,'eigenvector');
Z = X_norm(:,1);
Z(indROI) = normalizeImage(eC);
Z(~indROI) = 0;
tcf.eigenvectorCentrality = reshape(Z,[size(img,1), size(img,2)]);

% Compute Clustering Coefficient
clusteringCoefficient = clustercoeffs(double(D1)); 
Z = X_norm(:,1);
Z(indROI) = normalizeImage(clusteringCoefficient);
Z(~indROI) = 0;
tcf.clusteringCoefficient = reshape(Z,[size(img,1), size(img,2)]);

% Compute average path length
spD  = distances(G);
Z = X_norm(:,1);
Z(indROI) = normalizeImage(mean(spD));
Z(~indROI) = 0;
tcf.avgPathLength = reshape(Z,[size(img,1), size(img,2)]);

        