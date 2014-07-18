function net = rlsctrain(net, train_data, train_label)
% Regularized Least-Squares Classification training.
% Assuming all variables are mutually orthogonal. 
[net.samples net.vars] = size(train_data);

% Calculate matrix R.
R = zeros(1+net.vars, 1+net.vars);
R(1, 1) = net.samples;
for s = 1:net.vars
    R(1, s+1) = sum(train_data(:, s));
    R(s+1, 1) = R(1, s+1);
    R(s+1, s+1) = 1 + net.samples*net.lamb;
end

% Variable weights and bias.
w1 = R\([ones(net.samples, 1) train_data]'*train_label);
net.w = w1(2:net.vars+1);
net.bias = w1(1);

