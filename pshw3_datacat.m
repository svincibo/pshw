clear all; close all; clc;

% Next steps:
% Separate based on condition, letters or numbers.
% Get missing data filled in, Emily Wagner.
% Integrate behavioral measures at t1 to see if any predict who learns and
% who does not learn to draw letters or numbers.

rootdir = '/Volumes/Seagate/project-preschool-handwriting';

% remove = [321 416 421];

%% Dataset 1: Behavioral data.
d1 = readtable(fullfile(rootdir, 'supportfiles', 'rawdata', 'results_standard_assessments.csv'));
d1 = sortrows(d1, "ID");
if exist('remove')
idxtemp = d1.ID ~= remove; idxtemp2 = idxtemp(:, 1).*idxtemp(:, 2).*idxtemp(:, 3); idx = find(idxtemp2);
d1 = d1(idx, :);
end

%% Dataset 2: handwriting images.
dtemp1 = readtable(fullfile(rootdir, 'supportfiles', 'mturkdata', 'letters-all-randomized.csv'));
dtemp1.Var3 = repmat({'letter'}, [size(dtemp1, 1) 1]);
dtemp2 = readtable(fullfile(rootdir, 'supportfiles', 'mturkdata', 'numbers-all-randomized.csv'));
dtemp2.Var3 = repmat({'digit'}, [size(dtemp2, 1) 1]);
dtemp2.Symbol = cellstr(string(dtemp2.Symbol));
d2 = cat(1, dtemp1, dtemp2);
d2 = d2(:, [1:3 end]);
clear dtemp1 dtemp2;

%% Dataset 3: mturk classifications.
mturk_contents = dir(fullfile(rootdir, 'supportfiles', 'mturkdata', 'mturk_batches'));
% Remove the '.' and '..' files.
mturk_contents = mturk_contents(arrayfun(@(x) x.name(1), mturk_contents) ~= '.');
% Keep only names that are digit folders.
% mturk_contents = mturk_contents(arrayfun(@(x) x.name(1), mturk_contents) == 'd');
for i = 1:size(mturk_contents, 1)

    % Read in data for this mturk batch.
    d3_temp = readtable(fullfile(mturk_contents(i).folder, mturk_contents(i).name));

    % Add a column for symbol type.
    if contains(mturk_contents(i).name, 'digit')

        d3_temp(:, end+1) = repmat({'digit'}, [size(d3_temp, 1) 1]);
        d3_temp.Answer_category_label = cellstr(string(d3_temp.Answer_category_label));

    elseif contains(mturk_contents(i).name, 'letter')

        d3_temp(:, end+1) = repmat({'letter'}, [size(d3_temp, 1) 1]);

    end

    if i == 1

        d3 = d3_temp;

    else

        d3 = cat(1, d3, d3_temp);

    end

    clear d3_temp;

end
d3.Properties.VariableNames{end} = 'SymbolType';

% Add the symbol name provided in d2 to each row of d3 using image_url as the link.
for i = 1:size(d3, 1)

    currentimage = d3.Input_image_url{i};
    idx = find(contains(d2.image_url, currentimage));

    currentsymbol = d2.Symbol{idx};
    if i == 1

        d3(i, end+1) = {currentsymbol};

    else

        d3(i, end) = {currentsymbol};
    end

    clear currentsymbol currentimage idx;

end
d3.Properties.VariableNames{end} = 'Symbol';

% Create accuracy variable.
d3.trial_accuracy = strcmp(d3.Answer_category_label, d3.Symbol);

% Collapse across trials for each symbol instance by averaging accuracy across the 4 trials.
hitid = unique(d3.HITId);
for i = 1:size(hitid, 1)

    % Get index for all rows associated with this hit, i.e., this image.
    currenthitid = hitid(i);
    idx = find(contains(d3.HITId, currenthitid));

    % Get the image url for this hit, i.e., this image.
    currentimageurl = unique(d3.Input_image_url(idx));
    d4_image(i) = currentimageurl;

    % Get subID associated with this hit, i.e., this image.
    temp = currentimageurl{:};
    d4_subID(i) = {temp(9:11)};

    % Get age of this subject.
    idx3 = find(d1.ID == str2num(temp(9:11))); clear temp;
    d4_age(i) = {d1.age(idx3)};

    % Get age group of this subject.
    d4_agegroup(i) = {d1.grp_age(idx3)};

    % Get sex of this subject. F = 1, M = 2
    tempsex = d1.sex{idx3};
    if strcmp(tempsex, 'F' )
        d4_sex(i) = {'1'};
    elseif strcmp(tempsex, 'M' )
        d4_sex(i) = {'2'};
    end
    clear tempsex;

    % Get behavioral data for this subject.
    d4_vmi(i) = {d1.VMI_1(idx3)};
    d4_vp(i) = {d1.VP_1(idx3)};
    d4_mc(i) = {d1.MC_1(idx3)};
    d4_larue(i) = {d1.LaRue_1(idx3)};
    d4_lk(i) = {d1.L_K__1(idx3)};
    d4_rhyme(i) = {d1.Rhyme_1(idx3)};
    d4_sound(i) = {d1.Sound_1(idx3)};
    d4_vdis(i) = {d1.V_Dis_1(idx3)};
    d4_name(i) = {d1.Name_1(idx3)};
    d4_cat(i) = {d1.Cat__1(idx3)};
    d4_rec(i) = {d1.Rec__1(idx3)};
    d4_lkk(i) = {d1.L_K__1_1(idx3)};
    clear idx3;

    % Get symbol type associated with this hit, i.e., this image.
    d4_symboltype(i) = unique(d3.SymbolType(idx));

    % Get the symbol associated with this hit, i.e., this image.
    d4_symbol(i) = unique(d3.Symbol(idx));

    % Get the week associated with this hit, i.e., this image.
    temp = currentimageurl{:};
    idx2 = strfind(temp, 'week');
    d4_week(i) = num2cell(str2num(temp(idx2+4)));
    clear idx2;

    % Get the trial associated with this hit, i.e., this image.
    idx2 = strfind(temp, '.png');
    temp2 = temp(idx2-1);
    % fixing annoyances of the png file names, sometimes has + at the end before .png
    if strcmp(temp2, '+') || strcmp(temp2, '$')
        d4_trial(i) = num2cell(str2num(temp(idx2-2)));
    elseif strcmp(temp2, 'k')
        d4_trial(i) = num2cell(str2num(temp(idx2-7)));
    else
        d4_trial(i) = num2cell(str2num(temp(idx2-1)));
    end
    clear temp temp2 idx2;

    % Get the average amount of time spent on this hit, i.e., this image, across workerids.
    d4_time(i) = num2cell(nanmean(d3.WorkTimeInSeconds(idx)));

    % Get the average accuracy of labeling this hit, i.e., this image, across workerids.
    d4_accuracy(i) = num2cell(nanmean(d3.trial_accuracy(idx)));

    % Calculate confusability score. 
    c = d3.Answer_category_label(idx); target = d4_symbol(i);
    d4_answer1(i) = c(1); d4_answer2(i) = c(2); d4_answer3(i) = c(3); d4_answer4(i) = c(4); d4_answer5(i) = c(5);
%     if all(contains(c, target)) %all labels are correct
%         d4_confusability(i) = num2cell(((length(unique(c))-1)/length(c)));
%     elseif any(contains(c, target)) % at least one label is correct
%         d4_confusability(i) = num2cell(((length(unique(c))-1)/length(c)));
%     else % no labels are correct
        d4_confusability(i) = num2cell(((length(unique(c)))/length(c)));
%     end

    clear currenthitid currentimageurl c;

end

% Normalize confusability metric so that it ranges from 0 to 1 like
% accuracy.
d4_confusability = num2cell((cell2mat(d4_confusability)-min(cell2mat(d4_confusability)))/(max(cell2mat(d4_confusability))-min(cell2mat(d4_confusability))));

d4 = array2table([d4_subID', d4_age', d4_agegroup', d4_sex', d4_image', d4_symboltype', d4_symbol', d4_week', d4_trial', d4_time', d4_accuracy', d4_confusability', ...
    d4_answer1', d4_answer2', d4_answer3', d4_answer4', d4_answer5', ...
    d4_vmi', d4_vp', d4_mc', d4_larue', d4_lk', d4_rhyme', d4_sound', d4_vdis', d4_name', d4_cat', d4_rec', d4_lkk']);
d4.Properties.VariableNames{1} = 'subid';
d4.Properties.VariableNames{2} = 'age';
d4.Properties.VariableNames{3} = 'agegroup';
d4.Properties.VariableNames{4} = 'sex';

d4.Properties.VariableNames{5} = 'image';
d4.Properties.VariableNames{6} = 'symboltype';
d4.Properties.VariableNames{7} = 'symbol';
d4.Properties.VariableNames{8} = 'week';
d4.Properties.VariableNames{9} = 'trial';

d4.Properties.VariableNames{10} = 'time';
d4.Properties.VariableNames{11} = 'legibility';
d4.Properties.VariableNames{12} = 'confusability';

d4.Properties.VariableNames{13} = 'label1';
d4.Properties.VariableNames{14} = 'label2';
d4.Properties.VariableNames{15} = 'label3';
d4.Properties.VariableNames{16} = 'label4';
d4.Properties.VariableNames{17} = 'label5';

d4.Properties.VariableNames{18} = 'vmi';
d4.Properties.VariableNames{19} = 'vp';
d4.Properties.VariableNames{20} = 'mc';
d4.Properties.VariableNames{21} = 'larue';
d4.Properties.VariableNames{22} = 'lk';
d4.Properties.VariableNames{23} = 'rhyme';
d4.Properties.VariableNames{24} = 'sound';
d4.Properties.VariableNames{25} = 'vdis';
d4.Properties.VariableNames{26} = 'name';
d4.Properties.VariableNames{27} = 'cat';
d4.Properties.VariableNames{28} = 'rec';
d4.Properties.VariableNames{29} = 'lkk';

% clear d4_subID d4_age d4_agegroup d4_image d4_symboltype d4_symbol d4_week d4_time d4_accuracy;
% clear ;
d4 = sortrows(d4, ["subid", "symbol", "week", "trial"]);

if exist('remove')
idxtemp = d4.subid ~= remove; idxtemp2 = idxtemp(:, 1).*idxtemp(:, 2).*idxtemp(:, 3); idx = find(idxtemp2);
d4 = d4(idx, :);
end

% Write out the mturk data, d4, as csv file.
n = length(unique(d4.subid));
filename = sprintf('pshw_mturkdata_behdata_n%s_%s', num2str(n), datestr(now,'yyyymmdd'));
writetable(d4, fullfile(rootdir, 'supportfiles', [filename '.csv']));

