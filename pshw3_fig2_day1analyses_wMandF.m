clear all; close all; clc;

rootdir = '/Volumes/Seagate/project-preschool-handwriting';
d = readtable(fullfile(rootdir, 'supportfiles', 'pshw_mturkdata_behdata_n37_20230929.csv'));
figcount = 0;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0;
alphablend = .8;
yticklength = 0;
xticklength = 0.02;

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

            % All.
            idx = intersect(intersect(idx_subid, idx_symbol), idx_week);
            data_legibility_a(s, count) = nanmean(d.legibility(idx));
            data_confusability_a(s, count) = nanmean(d.confusability(idx));

            % Female.
            idx_female = find(d.sex==1);
            idx = intersect(intersect(intersect(idx_subid, idx_symbol), idx_week), idx_female);
            data_legibility_f(s, count) = nanmean(d.legibility(idx));
            data_confusability_f(s, count) = nanmean(d.confusability(idx));

            % Male. 
            idx_male = find(d.sex==2);
            idx = intersect(intersect(intersect(idx_subid, idx_symbol), idx_week), idx_male);
            data_legibility_m(s, count) = nanmean(d.legibility(idx));
            data_confusability_m(s, count) = nanmean(d.confusability(idx));

            %             end
            datalabel(count) = {strcat(alphabet{a}, num2str(w))};

            clear idx idx_female idx_male idx_week;

        end

        clear idx_symbol;

    end

    behlets(s, :) = d(idx_subid(1), [1:2 18:end]);

    clear idx_subid;

end
data = cat(2, subids, data_legibility_a, data_confusability_a, data_legibility_f, data_confusability_f, data_legibility_m, data_confusability_m);
data(all(isnan(data(:, 2:end)), 2), :) = [];
for i = 1:length(data(:, 1))
    idx(i) = find(table2array(behlets(:, 1)) == data(i, 1));
end
letters = [behlets(idx, :) array2table(data(:, 2:end))];
letters.Properties.VariableNames = [behlets.Properties.VariableNames strcat(datalabel, 'lega') strcat(datalabel, 'cona') ...
    strcat(datalabel, 'legf') strcat(datalabel, 'conf') strcat(datalabel, 'legm') strcat(datalabel, 'conm')];
clear data datalabel;
n = size(letters, 1);
% % Estimate beta  and y-intercept for each letter and each child.
% for a = 1:length(alphabet)
% 
%     idx = find(contains(letters.Properties.VariableNames, 'acc') & contains(letters.Properties.VariableNames, alphabet{a}));
% 
%     for s = 1:length(letters.subid)
% 
%         if length(find(isnan(table2array(letters(s, idx))))) >= 2
% 
%             tempbeta(s) = NaN;
%             tempintercept(s) = NaN;
% 
%         else
% 
%             x = [1:length(idx)]; %repmat(1, [1 length(idx)])]';
%             y = table2array(letters(s, idx));
%             mdl = fitlm(x, y, 'linear');
%             tempbeta(s) = round(table2array(mdl.Coefficients(2, 1)), 2);
% 
%             % Select baseline as the y-intercept from the fitted line.
%             tempintercept(s) = round(table2array(mdl.Coefficients(1, 1)), 2);
% 
%             % Select baseline as the day one accuracy.
%             %         tempintercept(s) = y(1); % This option looks visually the same but is not tractable for the correlation calculation due to the 1s.
% 
%         end
% 
%         clear x y mdl;
% 
%     end
% 
%     beta(:, a) = tempbeta;
%     intercept(:, a) = tempintercept;
% 
%     clear idx tempbeta tempintercept;
% 
% end
% beta = array2table([letters.subid beta], 'VariableNames', [{'subid2'}; strcat(alphabet, 'beta')]);
% meanbeta = array2table(nanmean(table2array(beta(:, 2:end)), 2), 'VariableNames', {'meanbeta'});
% intercept = array2table([letters.subid intercept], 'VariableNames', [{'subid3'}; strcat(alphabet, 'int')]);
% meanint = array2table(nanmean(table2array(intercept(:, 2:end)), 2), 'VariableNames', {'meanint'});
% 
% % Add beta and intercept to the letters table. The letters table already includes beh data.
% letters = [letters meanbeta beta(:, 2:end) meanint intercept(:, 2:end)];
% 
% % % Zscore all variables, except age.
% % temp = table2array(letters(:, 3:end));
% % z = (temp-nanmean(temp, 1))./nanstd(temp, [], 1); 
% % zletters = [letters(:, 1:2) array2table(z)]; 
% % zletters.Properties.VariableNames = letters.Properties.VariableNames;
% % clear z temp;

%% Get day 1 accuracy for all.
day1acc_idx = find(contains(letters.Properties.VariableNames, '1lega'));
day1acc = letters(:, day1acc_idx);
day1acc.Properties.VariableNames = extractBefore(day1acc.Properties.VariableNames, '1lega');
acc_a = nanmean(table2array(day1acc), 1); n = (20+17)/2; na = n;
%length(find(~isnan(table2array(day1acc(:, 1))))); na=37; %note: 20 letters, 17 digits
accstd_a = 1.96*nanstd(table2array(day1acc), [], 1)/sqrt(n);

day1acc_idx = find(contains(letters.Properties.VariableNames, '1legf'));
day1acc = letters(:, day1acc_idx);
day1acc.Properties.VariableNames = extractBefore(day1acc.Properties.VariableNames, '1legf');
acc_f = nanmean(table2array(day1acc), 1); n = (7+13)/2; nf = n;
%length(find(~isnan(table2array(day1acc(:, 1))))); nf=20; %note: 7 letters, 13 digits
accstd_f = 1.96*nanstd(table2array(day1acc), [], 1)/sqrt(n);

day1acc_idx = find(contains(letters.Properties.VariableNames, '1legm'));
day1acc = letters(:, day1acc_idx);
day1acc.Properties.VariableNames = extractBefore(day1acc.Properties.VariableNames, '1legm');
acc_m = nanmean(table2array(day1acc), 1); n = (13 + 4)/2; nm = n;
%n = length(find(~isnan(table2array(day1acc(:, 1))))); nm=17; %note: 13 letters, 4 digits
accstd_m = 1.96*nanstd(table2array(day1acc), [], 1)/sqrt(n);

[h,p,ci,stats] = ttest2(acc_f, acc_m)
disp(['ALL: Mean legibility is ' num2str(mean(acc_a)) ', with SD = ' num2str(std(acc_a)) ', and n = ' num2str(na) '.'])
disp(['FEMALE: Mean legibility is ' num2str(mean(acc_f)) ', with SD = ' num2str(std(acc_f)) ', and n = ' num2str(nf) '.'])
disp(['MALE: Mean legibility is ' num2str(mean(acc_m)) ', with SD = ' num2str(std(acc_m)) ', and n = ' num2str(nm) '.'])

%% Get day 1 confusability.
day1con_idx = find(contains(letters.Properties.VariableNames, '1cona'));
day1con = letters(:, day1con_idx);
day1con.Properties.VariableNames = extractBefore(day1con.Properties.VariableNames, '1cona');
con_a = nanmean(table2array(day1con), 1); n = (20+17)/2; na = n;
%n = length(find(~isnan(table2array(day1con(:, 1))))); na=37; %note: 20 letters, 17 digits
constd_a = 1.96*nanstd(table2array(day1con), [], 1)/sqrt(n);

day1con_idx = find(contains(letters.Properties.VariableNames, '1conf'));
day1con = letters(:, day1con_idx);
day1con.Properties.VariableNames = extractBefore(day1con.Properties.VariableNames, '1conf');
con_f = nanmean(table2array(day1con), 1); n = (7+13)/2; nf = n;
%n = length(find(~isnan(table2array(day1con(:, 1))))); nf=20; %note: 7 letters, 13 digits
constd_f = 1.96*nanstd(table2array(day1con), [], 1)/sqrt(n);

day1con_idx = find(contains(letters.Properties.VariableNames, '1conm'));
day1con = letters(:, day1con_idx);
day1con.Properties.VariableNames = extractBefore(day1con.Properties.VariableNames, '1conm');
con_m = nanmean(table2array(day1con), 1); n = (13+4)/2; nm = n;
%n = length(find(~isnan(table2array(day1con(:, 1))))); nm=17; %note: 13 letters, 4 digits
constd_m = 1.96*nanstd(table2array(day1con), [], 1)/sqrt(n);

[h,p,ci,stats] = ttest2(con_f, con_m)
disp(['ALL: Mean confusability is ' num2str(mean(con_a)) ', with SD = ' num2str(std(con_a)) ', and n = ' num2str(na) '.'])
disp(['FEMALE: Mean confusability is ' num2str(mean(con_f)) ', with SD = ' num2str(std(con_f)) ', and n = ' num2str(nf) '.'])
disp(['MALE: Mean confusability is ' num2str(mean(con_m)) ', with SD = ' num2str(std(con_m)) ', and n = ' num2str(nm) '.'])

%% Plot Figure 1c: Relationship between accuracy and confusability.
figcount = figcount + 1; figure(figcount); hold on;
% scatter(acc, con, 'w'); xlabel('accuracy'); ylabel('confusability'); 
% ylim([0 1]); xlim([0 1]);
T=text(con_a, acc_a, day1con.Properties.VariableNames);
for t = 1:length(T); T(t).FontName = "Arial"; T(t).FontSize=18;  end
[r, p] = corrcoef(acc_a, con_a)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 0.25];
xax.TickValues = [0 0.08 0.16 0.24];
xax.TickLabels = {'0', '0.08', '0.16', '0.24'};
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0.25 1];
yax.TickValues = [0.25 0.50 0.75 1.00];
yax.TickLabels = {'0.25', '0.50', '0.75', '1.00'};
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.FontName = fontname;
yax.FontSize = fontsize;
% yax.FontAngle = fontangle;

% general
a = gca;
box off

a.XLabel.String = 'Confusability';
a.XLabel.FontSize = fontsize;

a.YLabel.String = 'Legibility';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['plot_figure2c_n=' num2str(n)]), '-dpng')

%% Plot Figure 1a. Rank ordered accuracy.
figcount = figcount + 1; figure(figcount); hold on;
 
[~, idx] = sort(acc_a, 'ascend');
acc_a = acc_a(idx); acc_f = acc_f(idx); acc_m = acc_m(idx);
accstd_a = accstd_a(idx); accstd_f = accstd_f(idx); accstd_m = accstd_m(idx);
accsym = day1acc.Properties.VariableNames(idx); clear idx;

% Plot female.
for p = 1:length(acc_f)
      plot([acc_f(p) - abs(accstd_f(p)) acc_f(p) + abs(accstd_f(p))], [p p], 'Color', 'k')
end
scatter(acc_f, 1:length(acc_f), 100, 'Filled', 'MarkerFaceColor', [178 34 34]/255, ...
    'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 1); 

% Plot male.
for p = 1:length(acc_m)
      plot([acc_m(p) - abs(accstd_m(p)) acc_m(p) + abs(accstd_m(p))], [p p], 'Color', 'k')
end
scatter(acc_m, 1:length(acc_m),100, 'Filled', 'MarkerFaceColor', [70 130 180]/255, ...
    'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 1); 

% Plot all.
for p = 1:length(acc_a)
      plot([acc_a(p) - abs(accstd_a(p)) acc_a(p) + abs(accstd_a(p))], [p p], 'Color', 'k')
end
scatter(acc_a, 1:length(acc_a), 100, 'Filled', 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'w'); 

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.25 1];
xax.TickValues = [0.25 0.50 0.75 1.00];
xax.TickLabels = {'0.25', '0.50', '0.75', '1.00'};
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(acc_a)+1];
yax.TickValues = 1:length(acc_a);
yax.TickDirection = 'out';
yax.TickLabels = accsym;
yax.FontName = fontname;
yax.FontSize = fontsize/2;

% general
a = gca;
box off

a.XLabel.String = 'Legibility';
a.XLabel.FontSize = fontsize;

a.YLabel.String = 'Symbol';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['plot_figure2a_n=' num2str(n)]), '-dpng')

hold off;

%% Plot Figure 1b. Rank ordered confusability.
figcount = figcount + 1; figure(figcount); hold on;

[~, idx] = sort(con, 'descend');
con = con(idx);
constd = constd(idx);
consym = day1con.Properties.VariableNames(idx); clear idx;

% scatter(con, 1:length(con), 100, 'Filled'); hold on;
for p = 1:length(con)
      plot([con(p) - abs(constd(p)) con(p) + abs(constd(p))], [p p], 'Color', 'k')
end
scatter(con, 1:length(con), 100, 'Filled', 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'w'); 

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 0.25];
xax.TickValues = [0 0.08 0.16 0.24];
xax.TickLabels = {'0', '0.08', '0.16', '0.24'};
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(con)+1];
yax.TickValues = 1:length(con);
yax.TickDirection = 'out';
yax.TickLabels = consym;
yax.FontName = fontname;
yax.FontSize = fontsize/2;

% general
a = gca;
box off

a.XLabel.String = 'Confusability';
a.XLabel.FontSize = fontsize;

a.YLabel.String = 'Symbol';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['plot_figure2b_n=' num2str(n)]), '-dpng')


