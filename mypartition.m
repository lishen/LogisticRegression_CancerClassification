function [train_data, train_label, test_data, test_label] = mypartition(data, label, trains)
% Partition into training and testing data. Equal distribution. 
samples = size(data, 1);
ratio = trains/samples;
if length(find(label == -1))    % Binary.
    pos_idx = find(label == 1); neg_idx = find(label == -1);
    poss = length(pos_idx); negs = length(neg_idx);
    pos_idx = pos_idx(randperm(poss)); neg_idx = neg_idx(randperm(negs));
    train_poss = ceil(poss*ratio);
    train_negs = floor(negs*ratio);
    train_idx = [pos_idx(1:train_poss); neg_idx(1:train_negs)];
    test_idx = [pos_idx(train_poss+1:poss); neg_idx(train_negs+1:negs)];
else    % Multiclass.
    i = 1;
    train_idx = []; test_idx = [];
    while 1
        cls_idx = find(label == i);
        ncls = length(cls_idx);
        if ncls == 0
            break;
        end
        cls_idx = cls_idx(randperm(ncls));
        trains = ceil(ncls*ratio);
        train_idx = [train_idx; cls_idx(1:trains)];
        test_idx = [test_idx; cls_idx(trains+1:ncls)];
        i = i + 1;
    end
end
train_data = data(train_idx, :);
test_data = data(test_idx, :);
train_label = label(train_idx);
test_label = label(test_idx);



