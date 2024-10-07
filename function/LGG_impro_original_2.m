function [Z] = LGG_impro_original_2(X,Z_ini,lambda1,lambda2,lambda3,max_iter,miu)
% L:Lasso constraint  GG:Gaussian Graph
[m,n] = size(X);
max_miu = 1e8;        % µü´ú´ÎÊý
tol  = 1e-6;
tol2 = 1e-2;
rho = 1.1;
C1 = zeros(size(Z_ini));

% tic()
for iter = 1:max_iter
    if iter == 1        
        Z = Z_ini;
        U = Z_ini;
    end
    Z_old = Z;
    U_old = U;
       
    % -------- Update Z --------- %
    % calculate D
    D = L2_distance_1(X,X);

%     % calculate P by KNN with theta 
%     k = 10;
%     theta = 1.2;
%     H = construct_knnGraph(X', k);  % construct KNN graph
%     P = zeros(size(H))+theta;       % initialize P
%     for i = 1:length(H)
%         IDX = find(H(i,:));         % find the neighbor subscripts of ith cell 
%         P(i,IDX) = 1;               % let the neighbors' value to 1
%     end
   % calculate P by KNN with Gaussian kernel function 
   k = 10;
   sigma = 1;
   H = construct_knnGraph(X', k);  % construct KNN graph
%    P = zeros(size(H));             % initialize P
   P = ones(size(H));
   G = exp(-D/(2*sigma^2));
   for i = 1:length(H)
       IDX = find(H(i,:));         % find the neighbor subscripts of ith cell 
       P(i,IDX) = G(i,IDX);        % let the neighbors' value to gaussian
   end

   % calculate Z
    M1 = U-C1/miu;
    G1 = 2*X'*X+2*lambda1*eye(n)+miu*eye(n);
    Z = G1\(2*X'*X+miu*M1-lambda3*(P.*D));

    % optimization
%     Z(logical(eye(size(Z)))) = 0;
%     Z = Z-diag(diag(Z));
%     Z = (Z+Z')/2;


    % -------- U ------------ %
%     tic()
   M2 = Z+C1/miu;
   G2 = max(M2-lambda2/miu,0)+min(M2+lambda2/miu,0);
   U = max(G2,0);
%     timeU=toc()
    % -------- Update C1 miu -------- %
    L1 = Z-U;
    C1 = C1+miu*L1;
    
    LL1 = norm(Z-Z_old,'fro');
    LL2 = norm(U-U_old,'fro');
    SLSL = max(max(LL1,LL2))/norm(X,'fro');
    if miu*SLSL < tol2
        miu = min(rho*miu,max_miu);
    end
    % --------- obj ---------- %
%         leq1 = max(abs(L1(:)));
    stopC = abs(L1(:));
    %     stopC = max(leq1,max(abs(L3(:))));
    if stopC < tol
        iter
        break;
    end
end
end
