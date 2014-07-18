function [T, U, gK, mc, com] = mykpls(data, label, types, com, ratio, ker, para)
% Kernel Partial Least Squares.
if nargin < 6
    ker = 'linear';
    para = [];
end
% Initial kernel matrix.
K = kf(ker, para, data, data);
mc = mean(K);
K = kcenter(K, mc); gK = K;
% Centering b;
codemat = codegen(types, 'ova');
if types < 3
    b = label;
    b = b - mean(b);
    b0 = b;
else
    samples = size(data, 1);
    for i = 1:length(label)
        b(i, :) = codemat(label(i), :);
    end
    b = b - repmat(mean(b), [samples 1]);
    bmat = b*b';
    bmat0 = bmat;
end
% PLS components
if ~isempty(com) && com == -1 || ~isempty(ratio)
    com = Inf;
end
for i = 1:com
    if types < 3
        t = K*b; t = t/norm(t);
        U(:, i) = b*(b'*t);
    else
        [V D] = eig(K*bmat);
        D = diag(D);
        [dummy idx] = max(D);
        t = V(:, idx);
        if ~isreal(t)
            com = i - 1;
            break;
        end
        U(:, i) = bmat*t;
    end
    T(:, i) = t;
    if types < 3
        b = b - t*(t'*b);
        lab_rat = 1 - sum(b.^2)/sum(b0.^2); % Ratio of label variance captured. 
    else
        bmat = bmat - t*(t'*bmat) - bmat*t*t' + t*((t'*bmat)*t)*t';
        lab_rat = 1 - sum(sum(bmat.^2))/sum(sum(bmat0.^2)); % Ratio of label variance captured. 
    end
    K = K - t*(t'*K) - K*t*t' + t*((t'*K)*t)*t';
    var_rat = 1 - sum(sum(K.^2))/sum(sum(gK.^2));   % Ratio of variable variance captured
    if com == Inf && lab_rat > ratio && var_rat > 0.9
        com = i;
        break;
    end
end










