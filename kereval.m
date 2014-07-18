clc
clear all
warning off;
%--------------------- Load datasets-----------------------%
alldata = struct('file',{'lungCancer', 'prostate_TumorVSNormal', 'breastCancer', 'AMLALL'},...
    'path',{'data'},...
    'name',{'Lung Cancer', 'Prostate Cancer', 'Breast Cancer', 'ALL & AML'},...
    'genes',{12533, 12600, 24481, 7129},...
    'types',{2, 2, 2, 2},...
    'label',{{'Mesothelioma', 'ADCA'}, {'Tumor', 'oNrmal'}, {'relapse', 'onn-relapse'}, {'ALL', 'AML'}},...
    'trains',{32, 102, 78, 38}, 'train_label',{1}, 'train_data',{1},...
    'tests',{149, 34, 19, 34}, 'test_label',{1}, 'test_data',{1});

cvdata = struct('file',{'colonTumor', 'ovarian', 'centralNervousSystem_outcome'},...
    'path',{'data'},...
    'name',{'Colon Tumor', 'Ovarian Cancer', 'Central Nervous System'},...
    'genes',{2000, 15154, 7129},...
    'types',{2, 2, 2},...
    'label',{{'engative', 'positive'}, {'oNrmal', 'Cancer'}, {'Class1', 'Class0'}},...
    'trains',{62, 253, 60}, 'train_label',{1}, 'train_data',{1});

seldata = 3; cvtag = 1;
nfold = 3; % -1 for LOOCV
iters = 100; % Iterations for random partition.

fsrt.method = 'PCA';
fsrt.com = []; fsrt.ratio = 0.85;   % 0.99 for PLS label variance. 0.85 for PCA variable variance.
fsrt.ker = 'linear'; fsrt.para = [];

nlker = 'rbf';
P_ary = [1 1e2 1e4 1e6];
% P_ary = [];

classifier.algo = 'svm1'; % svm, svm1, rlsc, rlsc0.
classifier.ker = 'linear';  % Kernel.
classifier.para = [];   % Kernel parameters.
classifier.C = []; 
C_ary = [10000 100 1 0.01 0.0001];
% C_ary = [0.0001 0.01 1 100 10000];

usesig = 0; % Use sigmoid or not.

conlvl = 0.95;  % T-test for gene selection.
selnum = [];

fid = fopen('results/KER/kpcasvm_cns.txt', 'w');
% fid = 1;

if ~cvtag && seldata ~= 2
    datasrt = loadcd(alldata(seldata));
    train_data = datasrt.train_data; train_label = datasrt.train_label;
    test_data = datasrt.test_data; test_label = datasrt.test_label;
    trains = datasrt.trains; tests = datasrt.tests; samples = trains + tests;
    comb_data = [train_data; test_data];
    comb_label = [train_label; test_label];
else
    if ~cvtag && seldata == 2
        datasrt = loadcd(alldata(seldata));
    else
        datasrt = loadcv(cvdata(seldata));
    end
    samples = datasrt.trains;
    comb_data = datasrt.train_data;
    comb_label = datasrt.train_label;
end

comb_data = comb_data - min(min(comb_data)) + 1;
comb_data = log(comb_data);

if ~isempty(conlvl)
    grank = trank(comb_data, comb_label, conlvl);
    comb_data = comb_data(:, grank);
    selnum = length(grank);
elseif ~isempty(selnum)
    grank = trank(comb_data, comb_label);
    comb_data = comb_data(:, grank(1:selnum));
end

name = datasrt.name;
genes = datasrt.genes;
types = datasrt.types;
fsrt.types = types;

fprintf(fid, '\n\n----------------------------------------------------------------\n');
fprintf(fid, 'Classification of %s.\n', name);
fprintf(fid, 'Total samples: %d.\n', samples);
fprintf(fid, 'Classifier: %s, Kernel: %s.\n', classifier.algo, classifier.ker);
if usesig
    fprintf(fid, 'Use sigmoid output.\n');
end
fprintf(fid, 'Feature Selection: %s, Kernel: %s.\n', fsrt.method, fsrt.ker);
if nfold == -1
    fprintf(fid, 'LOOCV is employed\n');
else
    fprintf(fid, '%d-fold cross-validation is employed\n', nfold);
end
if ~isempty(selnum)
    fprintf(fid, 'Preliminary selected genes: %d.\n', selnum);
end

lstart = cputime;

partrains = ceil(samples*2/3);
for i = 1:iters
    [train_data train_label test_data test_label] = mypartition(comb_data, comb_label, partrains);
    [train_data test_data] = datanorm(train_data, test_data);   % Mean = 0, Var = 1. 
    
    minerr = Inf;
    % Grid kernel parameters selection.
    for j = 1:length(P_ary)+1
        % If use dimension reduction, set kernel parameters for it.
        % If do not use dimension reduction, set kernel parameters for SVM.
        if ~strcmp(fsrt.method, 'NO')
            if j == length(P_ary)+1
                fsrt.ker = 'linear';
                fsrt.para = [];
            else
                fsrt.ker = nlker;
                fsrt.para = P_ary(j);
            end
        else
            if j == length(P_ary)+1
                classifier.ker = 'linear';
                classifier.para = [];
            else
                classifier.ker = nlker;
                classifier.para = P_ary(j);
            end
        end
        % Grid regularization parameter selection.
        for k = 1:length(C_ary)
            classifier.C = C_ary(k);
            if usesig
                [A B err] = myplatt(classifier, nfold, train_data, train_label, fsrt);
            else
                err = cvalid(classifier, nfold, train_data, train_label, fsrt);
            end
            % Choose parameters when a smaller CV error is found.
            if err <= minerr
                if strcmp(fsrt.method, 'NO')
                    kerused = classifier.ker; paraused = classifier.para;
                else
                    kerused = fsrt.ker; paraused = fsrt.para;
                end
                if usesig
                    Aused = A; Bused = B;
                end
                Cused = classifier.C;
                minerr = err;
            end
        end
    end
    classifier.C = Cused;
    % Classifier trained on training data and test on testing data, use parameters from grid selection. 
    if ~isempty(fsrt)
        if strcmp(fsrt.method, 'PLS')
            [trnT U K mc com] = mykpls(train_data, train_label, types, fsrt.com, fsrt.ratio, kerused, paraused);
            tstT = kplspot(test_data, train_data, K, mc, U, trnT, kerused, paraused);
            train_data = trnT; test_data = tstT;
        elseif strcmp(fsrt.method, 'PCA')
            [trnT D mc com] = mykpca(train_data, fsrt.com, fsrt.ratio, kerused, paraused);
            tstT = kpcapot(test_data, train_data, mc, trnT, D, kerused, paraused);
            train_data = trnT; test_data = tstT;
        end
    else
        classifier.ker = kerused; classifier.para = paraused; 
    end
    [y output] = csroutput(classifier, train_data, train_label, test_data);
    if usesig
        y = sigmoid(output, Aused, Bused);
    end
    terr = nnz(y - test_label);
    if isempty(paraused)
        paraused = -1;
    end
    if strcmp(fsrt.method, 'NO')
        com = -1;
    end
    fprintf(fid, 'Iter: %d, ker = %s, para = %.4f, C = %.6f, com = %d, cverr = %d, terr = %d.\n', i, kerused, paraused, Cused, com, minerr, terr);
    com_ary(i) = com;
    E_ary(i) = terr;
end
errmean = mean(E_ary);
errstd = std(E_ary);
fprintf(fid, 'Testing error mean: %.2f, std: %.2f.\n', errmean, errstd);
fprintf(fid, 'Error rate mean: %f, std: %f.\n', errmean/(samples-partrains), errstd/(samples-partrains));
if ~strcmp(fsrt.method, 'NO')
    comean = mean(com_ary);
    comstd = std(com_ary);
    fprintf(fid, 'Components mean: %.2f, std: %.2f.\n', comean, comstd);
end

lend = cputime;
ltime = lend - lstart;
fprintf(fid, 'time: %f.\n', ltime);

if fid ~= 1
    fclose(fid);
end



















