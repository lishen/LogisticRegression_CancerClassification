function grank = trank(train_data, train_label, conlvl)
% Student T-test to select genes based on the P-value.

[samples vars] = size(train_data);
if exist('conlvl') && ~isempty(conlvl)
    crival = tinv(conlvl, samples - 2);
end
gscore = zeros(vars, 1);

idx1 = find(train_label == 1);
idx2 = find(train_label == -1);
n1 = length(idx1);
n2 = length(idx2);
for i = 1:vars
    m1 = mean(train_data(idx1, i));
    m2 = mean(train_data(idx2, i));
    var1 = var(train_data(idx1, i));
    var2 = var(train_data(idx2, i));
    gscore(i) = abs(m1-m2)/sqrt(var1/n1+var2/n2);
end
if exist('conlvl') && ~isempty(conlvl)
    grank = find(gscore >= crival);
else
    [gscore grank] = sort(gscore);
    gscore = flipud(gscore);
    grank = flipud(grank);
end










