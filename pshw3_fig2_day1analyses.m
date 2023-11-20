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

            idx = intersect(intersect(idx_subid, idx_symbol), idx_week);

            %             % Ensure that participants and symbols with perfect scores on day 1 are effectively removed.
            %             if w == 1 && nanmean(d.accuracy(idx)) == 1
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

    behlets(s, :) = d(idx_subid(1), [1:2 16:end]);

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

%% Get day 1 accuracy.
day1acc_idx = find(contains(letters.Properties.VariableNames, '1leg'));
day1acc = letters(:, day1acc_idx);
day1acc.Properties.VariableNames = extractBefore(day1acc.Properties.VariableNames, '1leg');

acc = nanmean(table2array(day1acc), 1);
accstd = 1.96*nanstd(table2array(day1acc), [], 1)/sqrt(n);

%% Get day 1 confusability.
day1con_idx = find(contains(letters.Properties.VariableNames, '1con'));
day1con = letters(:, day1con_idx);
day1con.Properties.VariableNames = extractBefore(day1con.Properties.VariableNames, '1con');

con = nanmean(table2array(day1con), 1);
constd = 1.96*nanstd(table2array(day1con), [], 1)/sqrt(n);

%% Plot Figure 1c: Relationship between accuracy and confusability.
figcount = figcount + 1; figure(figcount); hold on;
% scatter(acc, con, 'w'); xlabel('accuracy'); ylabel('confusability'); 
% ylim([0 1]); xlim([0 1]);
T=text(con, acc, day1con.Properties.VariableNames);
for t = 1:length(T); T(t).FontName = "Arial"; T(t).FontSize=18;  end
[r, p] = corrcoef(acc, con)

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
 
[~, idx] = sort(acc, 'ascend');
acc = acc(idx);
accstd = accstd(idx);
accsym = day1acc.Properties.VariableNames(idx); clear idx;

% scatter(acc, 1:length(acc), 100, 'Filled'); hold on;
for p = 1:length(acc)
      plot([acc(p) - abs(accstd(p)) acc(p) + abs(accstd(p))], [p p], 'Color', 'k')
end
scatter(acc, 1:length(acc), 100, 'Filled', 'MarkerFaceColor', 'k', ...
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
yax.Limits = [0 length(acc)+1];
yax.TickValues = 1:length(acc);
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


