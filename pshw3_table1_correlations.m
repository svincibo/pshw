clear all; close all; clc;
figcount = 0;
rootdir = '/Volumes/Seagate/project-preschool-handwriting';
d = readtable(fullfile(rootdir, 'supportfiles', 'pshw_mturkdata_behdata_n37_20230929.csv'));
figcount = 0;

data = [];
alphabet = unique(d.symbol);
subids = unique(d.subid);
weeks = unique(d.week);
for s = 1:length(subids)

    idx_subid = find(d.subid == subids(s));
    %     disp(subids(s))
    count = 0;

    for a = 1:length(alphabet)

        idx_symbol = find(strcmp(d.symbol, alphabet{a}));
        %         wcount = 0;

        for w = 1:length(weeks)

            count = count + 1;
            %             wcount = wcount + 1;

            idx_week = find(d.week == w);

            idx = intersect(intersect(idx_subid, idx_symbol), idx_week);

            %             % Ensure that participants and symbols with perfect scores on day 1 are effectively removed.
            %             if w == 1 && nanmean(d.legibility(idx)) == 1
            %                 data(s, count) = NaN;
            %             elseif w ~= 1 && isnan(data(s, count-1))
            %                 data(s, count) = NaN;
            %             else
            data_legibility(s, count) = nanmean(d.legibility(idx));
            data_confusability(s, count) = nanmean(d.confusability(idx));

            %             end
            datalabel(count) = {strcat(alphabet{a}, num2str(w))};

            clear idx_week;

        end

        clear idx_symbol;

    end

    behlets(s, :) = d(idx_subid(1), [1:2 18:end]);

    clear idx_subid;

end
data = cat(2, subids, data_legibility, data_confusability);
data(all(isnan(data(:, 2:end)), 2), :) = [];
for i = 1:length(data(:, 1))
    idx(i) = find(table2array(behlets(:, 1)) == data(i, 1));
end
letters = [behlets(idx, :) array2table(data(:, 2:end))];
letters.Properties.VariableNames = [behlets.Properties.VariableNames strcat(datalabel, 'leg') strcat(datalabel, 'con')];
clear data datalabel;

% Estimate beta  and y-intercept for each letter and each child.
for a = 1:length(alphabet)

    idx = find(contains(letters.Properties.VariableNames, 'leg') & contains(letters.Properties.VariableNames, alphabet{a}));

    for s = 1:length(letters.subid)

        if length(find(isnan(table2array(letters(s, idx))))) >= 2

            tempbeta(s) = NaN;
            tempintercept(s) = NaN;

        else

            x = 1:length(idx); %repmat(1, [1 length(idx)])]';
            y = table2array(letters(s, idx));
            mdl = fitlm(x, y, 'linear');
            tempbeta(s) = round(table2array(mdl.Coefficients(2, 1)), 2);

            % Select baseline as the y-intercept from the fitted line.
            tempintercept(s) = round(table2array(mdl.Coefficients(1, 1)), 2);

            % Select baseline as the day one accuracy.
            %         tempintercept(s) = y(1); % This option looks visually the same but is not tractable for the correlation calculation due to the 1s.

        end

        clear x y mdl;

    end

    beta(:, a) = tempbeta;
    intercept(:, a) = tempintercept;

    clear idx tempbeta tempintercept;

end
beta = array2table([letters.subid beta], 'VariableNames', [{'subid2'}; strcat(alphabet, 'beta')]);
meanbeta = array2table(nanmean(table2array(beta(:, 2:end)), 2), 'VariableNames', {'meanbeta'});
intercept = array2table([letters.subid intercept], 'VariableNames', [{'subid3'}; strcat(alphabet, 'int')]);
meanint = array2table(nanmean(table2array(intercept(:, 2:end)), 2), 'VariableNames', {'meanint'});

% Construct visual-motor skill composite variable.
vmskill = array2table(mean([letters.vmi letters.vp letters.mc], 2), 'VariableNames', {'vmskill'});

% Construct literacy composite variable.
literacy = array2table(mean([letters.larue letters.vdis], 2), 'VariableNames', {'literacy'});

% Construct language/phonological composite variable.
language = array2table(mean([letters.rhyme letters.sound], 2), 'VariableNames', {'language'});

% Add beta and intercept to the letters table. The letters table already includes beh data.
letters = [letters meanbeta beta(:, 2:end) meanint intercept(:, 2:end) vmskill literacy language];

% Zscore all variables, except age.
temp = table2array(letters(:, 3:end));
z = (temp-nanmean(temp, 1))./nanstd(temp, [], 1); 
zletters = [letters(:, 1:2) array2table(z)]; 
zletters.Properties.VariableNames = letters.Properties.VariableNames;
clear z temp;

%% Calculate correlation matrix.
X = table2array(zletters(:, [3:6 10 8:9 521:523 2]));
[rval, pval] = corrcoef(X);

% Plot correlation matrix.
figure(1); imagesc(rval);
xticks(1:length(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2])));
xticklabels(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2]));
yticks(1:length(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2])));
yticklabels(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2]));
colorbar;

% Export tables.
out = array2table(rval);
out.Properties.VariableNames = zletters.Properties.VariableNames([3:6 10 8:9 521:523 2]);
writetable(out, fullfile(rootdir, 'output', 'table1_rval.csv'));
clear out;

out = array2table(pval);
out.Properties.VariableNames = zletters.Properties.VariableNames([3:6 10 8:9 521:523 2]);
writetable(out, fullfile(rootdir, 'output', 'table1_pval.csv'));
clear out;

%% Calculate correlation matrix controlling for age.
X = table2array(zletters(:, [3:6 10 8:9 521:523 2]));
[rval, pval] = partialcorr(X(:, 1:end-1), X(:, end));

% Plot correlation matrix controlling for age.
figure(1); imagesc(rval);
xticks(1:length(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2])));
xticklabels(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2]));
yticks(1:length(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2])));
yticklabels(zletters.Properties.VariableNames([3:6 10 8:9 521:523 2]));
colorbar;

% Export tables.
out = array2table(rval);
out.Properties.VariableNames = zletters.Properties.VariableNames([3:6 10 8:9 521:523]);
writetable(out, fullfile(rootdir, 'output', 'table1_rval_agecontrolled.csv'));
clear out;

out = array2table(pval);
out.Properties.VariableNames = zletters.Properties.VariableNames([3:6 10 8:9 521:523]);
writetable(out, fullfile(rootdir, 'output', 'table1_pval_agecontrolled.csv'));
clear out;