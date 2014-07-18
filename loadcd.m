function datasrt = loadcd(datasrt)
%--------------------Training Dataset------------------------------%
datasrt.train_label = zeros(datasrt.trains, 1);
datasrt.train_data = zeros(datasrt.trains, datasrt.genes);
fid = fopen([datasrt.path '/' datasrt.file '_train.data']);

for i = 1:datasrt.trains
    datasrt.train_data(i, :) = fscanf(fid, '%f,');
    label = fgets(fid);
    label = strtok(label);
    if strcmpi(label, datasrt.label{1}) == 1
        datasrt.train_label(i) = 1;
    elseif strcmpi(label, datasrt.label{2}) == 1
        datasrt.train_label(i) = -1;
    end
end
fclose(fid);

%-----------------------Test Dataset-------------------------------------%
datasrt.test_label = zeros(datasrt.tests, 1);
datasrt.test_data = zeros(datasrt.tests, datasrt.genes);
fid = fopen([datasrt.path '/' datasrt.file '_test.data']);

for i = 1:datasrt.tests
    datasrt.test_data(i, :) = fscanf(fid, '%f,');
    label = fgets(fid);
    label = strtok(label);
    if strcmpi(label, datasrt.label{1}) == 1
        datasrt.test_label(i) = 1;
    elseif strcmpi(label, datasrt.label{2}) == 1
        datasrt.test_label(i) = -1;
    end
end
fclose(fid);






