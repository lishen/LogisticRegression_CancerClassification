function tstT = kplspot(test_data, train_data, K, mc, U, T, ker, para)
% Projection of test samples in KPLS space. 
if nargin < 7
    ker = 'linear';
    para = [];
end
tstK = kf(ker, para, test_data, train_data);
tstK = kcenter(tstK, mc);
tstT = tstK*U*inv(T'*K*U);