%% compile
clear;
neo_compile_script

%% load input
load ../data/karateClub.mat

%% multi-level
% usage: neo(input,k,alpha,beta);
C_multi=neo(sparse(A),2,0.2,0);
save('karate_neo.mat','C_multi');

% %% one-level
% % usage: neo(input,k,alpha,beta,initial_C);
% C_one=neo(A,2,0.2,0,C_multi);
% save('karate_neo.mat','C_one','-append');