function [cls_idx]= setPatchIdx(X, vec)

% 这段代码重新写了 by jbwei

% f2          =   size(X,1);
% set         =   1:size(X, 2);
% L           =   size(set,2);
% 
% for j = 1 : L
%     dis   =   (vec(1, :) -  X(1, set(j))).^2;
%     for i = 2 : f2
%         dis  =  dis + (vec(i, :)-X(i, j)).^2;
%     end
%     [val ind]      =   min( dis );
%     cls_idx( j )   =   ind;
% end
% cls_idx = cls_idx';

[winSize, nPatches] = size(X);
if (size(vec,1)~=winSize)
	fprintf('Error:: window size of patches (%d) does not match the preset vec (%d)!', winSize, size(vec,1));
end

nClusters = size(vec, 2);
distance = zeros(nClusters, nPatches);
for iCluster=1:nClusters
	subMat = repmat(vec(:,iCluster), 1, nPatches);
	distance(iCluster, :) = sum((X - subMat).^2, 1);
end

[maxValue, cls_idx] = min(distance, [], 1);
cls_idx = cls_idx;
