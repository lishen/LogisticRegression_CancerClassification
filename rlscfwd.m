function [y, outputs] = rlscfwd(net, test_data)
% Regularized Least-Squares Classification testing.
outputs = test_data*net.w + net.bias;
y = outputs;
y(find(y >= 0)) = 1;
y(find(y < 0)) = -1;
