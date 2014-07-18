function net = rlsc0(ker, para, lamb)
% Regularized Least-Squares Classification. Construct the structure of net. 
net.lamb = lamb;
net.vars = 0;
net.samples = 0;
net.bias = [];
net.c = [];
net.w = [];
net.ker = ker;
net.para = para;
net.train_data = [];
