function [V, D, mc, com] = mykpca(data, com, ratio, ker, para)
% Kernel Principal Component Analysis. Components are sorted according to
% the principal values. 
if nargin < 4
    ker = 'linear';
    para = [];
end

% Initial kernel matrix.
K = kf(ker, para, data, data);
mc = mean(K);
K = kcenter(K, mc);

% KPCA components.
[V, D] = eig(K);
% Sort eigenvalues in descending order.
D = diag(D);
nzidx = find(D>1e-6);
D = D(nzidx);
V = V(:, nzidx');
[D, D_idx] = sort(D); D = flipud(D);
D_idx = flipud(D_idx);
if ~isempty(ratio)
    acc_var = 0; 
    for i = 1:length(nzidx)
        acc_var = acc_var + D(i);
        if acc_var >= ratio*sum(D)
            break;
        end
    end
    com = i;
else
    if com == -1    % For full PCA. 
        com = length(nzidx);
    elseif com > length(nzidx)
        com = length(nzidx);
    end
end
D = D(1:com); D = diag(D); % For projection. 
D_idx = D_idx(1:com); V = V(:, D_idx);



















