function net = rlsc0train(net, train_data, train_label)
% Regularized Least-Squares Classification training.
[net.samples net.vars] = size(train_data);

R = eye(net.samples) - ones(net.samples)/net.samples;
K = kf(net.ker, net.para, train_data, train_data);
L = net.samples*net.lamb*eye(net.samples);
net.c = (R*K + L)\(R*train_label);
net.bias = mean(train_label) - mean(K*net.c);

if strcmp(net.ker, 'linear')
    net.w = train_data'*net.c;
else
    net.train_data = train_data;
end

