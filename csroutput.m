function [y, output] = csroutput(classifier, train_data, train_label, test_data)
% Prediction results of classifiers. 
algo = classifier.algo;
ker = classifier.ker;
para = classifier.para;
C = classifier.C;
switch algo
case 'svm'
    % package by Schwaighofer.
    net = svm0(size(train_data, 2), ker, para, C);
    net = svmtrain(net, train_data, train_label);
    [y output] = svmfwd(net, test_data);
case 'svm1'
    % package by Steve Gunn. modified by Li Shen.
    [alpha, bias] = mysvc(train_data, train_label, ker, para, C);
    [y output] = mysvcfwd(train_data, train_label, test_data, ker, para, alpha, bias);
case 'rlsc'
    net = rlsc(C);
    net = rlsctrain(net, train_data, train_label);
    [y output] = rlscfwd(net, test_data);
case 'rlsc0'
    net = rlsc0(ker, para, C);
    net = rlsc0train(net, train_data, train_label);
    [y output] = rlsc0fwd(net, test_data);
end
