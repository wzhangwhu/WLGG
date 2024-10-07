function [H] = construct_knnGraph(X, k,t)
% X : n*d
% t RBF中超参
% k
% 构建RBF相似度矩阵

% addpath('E:\scRNA_seq_test_data')
% load('Test_7_Treutlin.mat')
% X = in_X;
% k = 10;

[n,dim] = size(X);
Dis =EuDist2(X, X);
t = sum(sum(Dis))/n/n;
% A = exp(-Dis/(2*t^2));

[dumb, idx] = sort(Dis, 1);
A_1 = zeros(n);
A_1 = A_1 + diag(ones(n,1));
for i = 1:n
    id = idx(2:k+1,i);
    di = Dis(id, i);
    A_1(id,i) = exp(-di/(2*t^2));
end
% H = exp(-W/(2*t^2));
H = A_1;
H = H-eye(n);
H = H';
% W_e = sum(H);
% % 构建d_v 点的度
% for i = 1:n
%     temp = 0;
%     for j = 1:n
%         if H(i,j) ~= 0
%             temp = temp + H(i,j)*W_e(j);
%         end
%     end
%     d_v(i) = temp; 
% end
% 
% d_e = sum(H);
% 
% D_v = diag(d_v');
% D_e = diag(d_e');
% W = diag(W_e');
% I = diag(ones(n,1));
% L = I - D_v^(-0.5)*H*W*D_e^(-1)*H'*D_v^(-0.5);

end



