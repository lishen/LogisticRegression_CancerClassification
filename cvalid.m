function [cverr, acmop, acmy]= cvalid(classifier, nfold, data, label, fsrt)
% Perform cross-validation and return the CV error, accumulated outputs and class labels. 
% Start training.
samples = size(data, 1);
if nfold == -1 % LOOCV
    nfold = samples;
else
    randidx = randperm(samples);
    data = data(randidx, :);
    label = label(randidx, :);
end

cverr = 0; acmop = []; acmy = [];
for i = 1:nfold
    [train_data train_label test_data test_label] = cvpar(data, label, nfold, i);
    if exist('fsrt') && ~isempty(fsrt)
        if strcmp(fsrt.method, 'PLS')
            [trnT U K mc com] = mykpls(train_data, train_label, fsrt.types, fsrt.com, fsrt.ratio, fsrt.ker, fsrt.para);
            tstT = kplspot(test_data, train_data, K, mc, U, trnT, fsrt.ker, fsrt.para);
            train_data = trnT; test_data = tstT;
        elseif strcmp(fsrt.method, 'PCA')
            [trnT D mc com] = mykpca(train_data, fsrt.com, fsrt.ratio, fsrt.ker, fsrt.para);
            tstT = kpcapot(test_data, train_data, mc, trnT, D, fsrt.ker, fsrt.para);
            train_data = trnT; test_data = tstT;
        end
    end
    [y output] = csroutput(classifier, train_data, train_label, test_data);
    cverr = cverr + nnz(y - test_label);
    acmop = [acmop; output];
    acmy = [acmy; test_label];
end





