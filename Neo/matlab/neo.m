function [C] = neo(A,k,alpha,beta,initC)
%% NEO-K-Means Matlab interface

if nargin < 2
    error('Please specify an adjacency matrix and the number of clusters.');
end

if nnz(abs(A-round(A)))~=0
    A = round(1000*A);
end

if nargin < 5
    one_level = 0;
    initC = [];
else
    one_level = 1;
end

if ~issparse(A)
    A = sparse(A);
end

if ~issparse(initC)
    initC = sparse(initC);
end
% binary assignment matrix
initC=double(logical(initC));

[out,out2]=neo_mex(A,nnz(A),k,alpha,beta,one_level,initC',size(initC,1),nnz(initC));
tt = ([out out2]~=0);
out = out(sum(tt,2)==2);
out2 = out2(sum(tt,2)==2);
C=sparse(out,out2,ones(length(out),1),size(A,1),k);

end