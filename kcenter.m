function K = kcenter(K, mc)
% Center the kernel matrix.
rows = size(K, 1);
cols = size(K, 2);
w = repmat(mc, [rows 1]);
K = K - w;
mr = mean(K, 2);
w = repmat(mr, [1, cols]);
K = K - w;
