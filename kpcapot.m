function tstV = kpcapot(test_data, train_data, mc, V, D, ker, para)
% Projection of test samples in KPCA space.
if nargin < 6
    ker = 'linear';
    para = [];
end
tstK = kf(ker, para, test_data, train_data);
tstK = kcenter(tstK, mc);
tstV = tstK*V*inv(D);