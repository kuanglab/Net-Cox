%% Net-Cox with co-expression network
clc
clear all

% load the Bonome data
% Data: a m-by-n matrix, m is the number of genes, n is the number of
% samples. Note: the sample dimesion of Data is sorted by the survival time (Times) 
% Times: survival time (in month), applied to sort the Data
% d: censor information, 1 is uncensored, 0 is censored.
load Bonome_SloanKettering

%Normalize the data by dividing square root of column sum and row sum
Data = Normal_M(Data);

%Calculate the co-expression network. 
S = getS(Data);



% parameter lambda balance the total likelihood and the network constraint, the larger the lambda the
% more weight on the network constraint
lambda = [1e-5,1e-4];

% alpha weighting the network matrix and the identity matrix in the network
% constraint, the smaller the alpha, the more we trust on network
% information. alpha should belong to (0,1].
alpha = [0.01,0.5,0.95];
bta_coexpression = NetCox(Data,lambda,alpha,d, S);

% [k l m] = size(bta_coexpression), k is the number of lambda, l is the number of alpha,
% m is the number of the genes. if setting beta =
% reshape(bta_coexpression(2,1,:),[],1), it return the coefficients learned
% from the algorithm when lambda equal to 1e-4 and alpha equal to 0.01.

