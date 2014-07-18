function [y, outputs] = rlsc0fwd(net, test_data)
% Regularized Least-Squares Classification testing.
if strcmp(net.ker, 'linear')
    outputs = test_data*net.w + net.bias;
else
    K = kf(net.ker, net.para, test_data, net.train_data);
    outputs = K*net.c + net.bias;
end
y = outputs;
y(find(y >= 0)) = 1;
y(find(y < 0)) = -1;
