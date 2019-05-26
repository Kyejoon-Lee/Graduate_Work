function [C] = read_neo_result(G,graph,nparts,alpha,beta,sigma)

fileName = strcat(graph,'_clust_',num2str(nparts),'_alpha_',alpha,'_beta_',beta,'_sigma_',sigma);
fid = fopen(fileName);

n = size(G,1);
noNodes=0;
noElem=0;
tline = fgets(fid);
while ischar(tline)
    noNodes=noNodes+1;
    clustInd = strread(tline,'%d ');
    noElem = noElem + length(clustInd);
    tline = fgets(fid);
end
fclose(fid);

assert(n==noNodes, 'no. of nodes is not equal to the no. of lines.\n');

rowNodes = zeros(noElem,1);
colClusters = zeros(noElem,1);
values = ones(noElem,1);
fid = fopen(fileName);
tline = fgets(fid);
nodeID=0;
sind = 1;
while ischar(tline)
    nodeID = nodeID+1;
    clustInd = strread(tline,'%d ');
    if isempty(clustInd)
        tline = fgets(fid);
        continue;
    end
    clustInd = clustInd + 1; % cluster index starts from 0
    eind = sind+length(clustInd)-1;
    colClusters(sind:eind) = clustInd;
    rowNodes(sind:eind) = repmat(nodeID,length(clustInd),1);
    sind = eind+1;
    tline = fgets(fid);
end

C = sparse(rowNodes,colClusters,values,n,nparts);

assert(size(C,2)==nparts, 'no. of columns of C is not equal to the no. of clusters');

%fprintf('graclus_overlap result... no. of clusters: %d\n',nparts);
%fprintf('no. of assigned nodes: %d out of %d (%6.4f %%)\n',nnz(sum(C,2)),n,nnz(sum(C,2))/n*100);

%graclus_overlap = C;
%save(strcat('graclus_overlap_',graph,'_C.mat'),'graclus_overlap','-v7.3');

end