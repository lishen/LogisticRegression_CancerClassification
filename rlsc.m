function net = rlsc(lamb)
% Regularized Least-Squares Classification. Construct the structure of net. 
net.lamb = lamb;
net.vars = 0;
net.samples = 0;
net.bias = [];
net.w = [];
