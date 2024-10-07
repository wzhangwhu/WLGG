warning off
clear
clc


addpath('data')
addpath('.\function')
% addpath('.\function\GenLouvain1.2')


% dataset = {'Test_1_mECS', 'Test_2_Kold', 'Test_3_Pollen', 'Test_4_Usoskin','Test_5_Goolam','Test_6_Braindata','Test_7_Treutlin','Test_8_Ting','Test_10_Ginhoux','Test_11_Deng','Test_12_Buettner', 'Test_14_hcancer','Test_15_human', 'Test_16_islet','Test_Zheng','CellLines','yan_2013','Test_GSE52583_LESdata_refine','Test_GSE36552_HEEdata_refine','Test_17_Breast_cancer_GSE75688','Test_Biase3celltypes','Test_biase_GSE57249','Test_camp2_data','Wang-pancreas','Test_Kumar','Test_Engel'} ;
% % dataset = {'Test_16_islet'};
% 
% % 选lambda1 and lambda2
% dataset = { 'Test_7_Treutlin','Test_8_Ting', 'Test_Zheng','Test_1_Zeisel_big','Enge_data','Test_2_Melanoma_micro'};
% % dataset = {'Test_7_Treutlin'}
% dataset = { 'Test_2_Kold','Test_6_Braindata','Test_7_Treutlin','Test_8_Ting','Test_15_human', 'Test_16_islet','Test_Zheng','Test_GSE52583_LESdata_refine','Enge_data','Test_2_Melanoma_micro','Quake_10x_Spleen'};
% dataset = {'Quake_10x_Bladder_data'};
% 
% dataset = {'Test_1_mECS', 'Test_2_Kold', 'Test_3_Pollen', 'Test_4_Usoskin','Test_5_Goolam','Test_6_Braindata','Test_7_Treutlin','Test_8_Ting','Test_10_Ginhoux','Test_11_Deng','Test_12_Buettner', 'Test_14_hcancer','Test_15_human', 'Test_16_islet','Test_Zheng','CellLines','yan_2013','Test_GSE52583_LESdata_refine','Test_GSE36552_HEEdata_refine','Test_17_Breast_cancer_GSE75688','Test_Biase3celltypes','Test_biase_GSE57249','Test_camp2_data','Wang-pancreas','Test_Kumar','Test_Engel','Test_1_Zeisel_big','Test_2_Macosko_big','Test_3_Tasic_big','Quake_10x_Limb_Muscle','Quake_10x_Bladder_data','PBMC_4k_10X','Test_3_PBMC_4k_10X','Enge_data','Melanoma_mico','Quake_10x_Spleen','Test_lake_data','Young_data_big'} ;
% dataset = {'Test_2_Macosko_big','Test_3_Tasic_big','Test_lake_data'};
% dataset = {'reprocess_LES_log'};
% dataset = {'reprocess_LES_subdivide_E18'};
%% Paper used datasets
% dataset = {'Test_2_Kold','Test_7_Treutlin','Test_8_Ting','Test_15_human', 'Test_16_islet','Test_Zheng','Test_GSE52583_LESdata_refine','Test_17_Breast_cancer_GSE75688','Test_1_Zeisel_big','Quake_10x_Bladder_data','Enge_data','Melanoma_mico','Quake_10x_Spleen','GSE115978_preprocessed'} ;
dataset = {'Test_Zheng'} ;

NMI_collect_newcombine=[]; 
ARI_collect_newcombine=[];
ACC_collect_newcombine=[];

%%
%parameters
miu = 0.01;
lambda1 = 3e-2;
lambda2 = 8e-3;
lambda3 = 3e-3;   % random
max_iter = 180;
% lambdap1=[ 5e-4 1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 7e-2 1e-1 5e-1 1 5 10 20 50 100];
%%
for ii = 1:length(dataset)
    ii

    load( dataset{ii});
    [m,n] = size(in_X);

    fea = in_X;                         % cellnum*genenum
    gnd=true_labs(:);                   %转变为列向量
    nn_class =length(unique(gnd));
    
    
    % ---------- Initilization for X  -------- %     初始化X  耗时开始
    fea = double(fea);  %%
    select_sample = [];
    select_gnd    = [];
    for i = 1:nn_class  %% 类的数目
        idx = find(gnd == i);
        idx_sample    = fea(idx,:);
        select_sample = [select_sample;idx_sample];
        select_gnd    = [select_gnd;gnd(idx)];
    end
    fea = select_sample';                      % 转置: genenum*cellnum
    fea = fea./repmat(sqrt(sum(fea.^2)),[size(fea,1) 1]);
    gnd = select_gnd;
     fea = log2(fea+1);
    
    X = fea;
    
    clear fea select_gnd select_sample idx
    %
    tic()
    % ---------- Initilization for Z  -------- %
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 10;
%   options.WeightMode = 'Binary';
    options.WeightMode = 'Cosine';
    Z = constructW(X',options);
    Z = full(Z);  
    Z1 = Z-diag(diag(Z));
      
    clear Z
    
        %---------- Optimization for Z  -------- % 
    [Z] = LGG_impro_original_2(X,Z1,lambda1,lambda2,lambda3,max_iter,miu);
    %     [Z] = LGG(X,Z1,lambda1,lambda2,max_iter,miu);
    Z= Z - diag(diag(Z));
    Z = (Z+Z')/2;
    
%%        ---------- Clustering ------------%
        [result_label, kerNS]= SpectralClustering(Z, nn_class);
    
    %------------ eigengap --------------%    2023.04.10
%{  
    Z(isnan(Z))=0;
    c=[];
    [predict_cluster] = clusteringCells_estimate_clusternum(Z,c,2:15,'showFigure');
    [result_label, kerNS]= SpectralClustering(Z, predict_cluster);
%}
    %------------ Louvain --------------%   
   %{
    k = full(sum(Z));
    twom = sum(k);
    C = @(v) Z(:,v) - k'*k(v)/twom;
    limit = 10000;
    verbose = 1;
    [S,Q] = genlouvain(C,limit,verbose,0);
    result_label = S(:);
    %}
    true_cluster = length(unique(gnd));
    predict_cluster = length(unique(result_label));
    
    time_WLGG=toc();
   
    %%    % ---------- evaluation ------- %
    NMI=Cal_NMI_newused(gnd, result_label);
    ARI=Contingency_ARI_newused(gnd, result_label);
    ACC=ACC_ClusteringMeasure(gnd, result_label);
        fprintf(['NMI_for_ ' dataset{ii} ' is %f\n'],NMI)
     fprintf(['ARI_for_ ' dataset{ii} ' is %f\n'],ARI)
     fprintf(['ACC_for_ ' dataset{ii} ' is %f\n'],ACC)
%     fprintf(['true_cluster_for_ ' dataset{ii} ' is %f\n'],true_cluster)
%     fprintf(['predict_cluster_for_ ' dataset{ii} ' is %f\n'],predict_cluster)
    fprintf(['time_for_' dataset{ii} ' is %f\n'],time_WLGG)


    
    collect_result_NMI(ii)=NMI;
    collect_result_ARI(ii)=ARI;
    collect_result_ACC(ii)=ACC;
    collect_result_predict_clutser(ii)=predict_cluster;
    collect_result_time(ii)=time_WLGG;
end
%
xlswrite(['LGG_SC_parameter_result_NMI_MATLAB2017b' '.xlsx'],collect_result_NMI)
xlswrite(['LGG_SC_parameter_result_ARI_MATLAB2017b' '.xlsx'],collect_result_ARI)
xlswrite(['LGG_SC_parameter_result_ACC_MATLAB2017b' '.xlsx'],collect_result_ACC)
 xlswrite(['GCF_Louvain_new_result_Time_MATLAB2017b' '.xlsx'],collect_result_time)
%}
% save parameter_WLGG_impro2_all_datasets_NMI.mat collect_result_NMI
% save parameter_WLGG_impro2_all_datasets_ARI.mat collect_result_ARI
% save parameter_WLGG_impro2_all_datasets_ACC.mat collect_result_ACC
% save parameter_WLGG_impro2_all_datasets_predict_clutser.mat collect_result_predict_clutser





%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the
% clustering of the nodes using the spectral clustering algorithm of
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
% Modified @ Chong You, 2015
%--------------------------------------------------------------------------

function [groups, kerNS] = SpectralClustering(CKSym,n)

warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans
REPlic = 100; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}

DN = diag(1./sqrt(sum(CKSym)+eps) );  %eps=2.2204e-16
LapN = speye(N) - DN * CKSym * DN;  % speye(N)生成N*N对角元素为1的矩阵  构建拉普拉斯矩阵
[~,~,vN] = svd(LapN);
kerN = vN(:,N-n+1:N);
%kerN = vN(:,N-12:N);
normN = sum(kerN .^2, 2) .^.5;
kerNS = bsxfun(@rdivide, kerN, normN + eps);  %% 矩阵kerN 每一行除以normN+eps

%-------------
%Y = pdist(kerNS);
%Z = linkage(Y);
%groups = cluster(Z,'maxclust',n);
kerNS = real(kerNS);
groups = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
end






