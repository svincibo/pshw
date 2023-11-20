clear all; close all; clc;

% https://www.qualtrics.com/support/stats-iq/analyses/regression-guides/interpreting-residual-plots-improve-regression/#ExaminingPredictedResidual

figcount = 0;
rootdir = '/Volumes/Seagate/project-preschool-handwriting';
d = readtable(fullfile(rootdir, 'supportfiles', 'pshw_mturkdata_behdata_n37_20230929.csv'));
figcount = 0;
symbolset = unique(d.symbol);

% d = d(yidx, :);

% Construct literacy composite variable.
literacy = array2table(mean([d.larue d.vdis], 2), 'VariableNames', {'literacy'});

% Construct language/phonological composite variable.
language = array2table(mean([d.rhyme d.sound], 2), 'VariableNames', {'language'});

% Construct visual-motor skill composite variable.
vmskill = array2table(mean([d.vmi d.vp d.mc], 2), 'VariableNames', {'vmskill'});

% Add new composite variables to d.
d = [d literacy language vmskill];

%% Re-oganize data to average over trial.
dnew = table();
subids = unique(d.subid); n = length(subids);

weeks = unique(d.week);
count = 0;
for s = 1:length(subids)

    idx_subid = find(d.subid == subids(s));
    %     disp(subids(s))

    for a = 1:length(symbolset)

        idx_symbol = find(strcmp(d.symbol, symbolset{a}));
        %         wcount = 0;

        for w = 1:length(weeks)

            idx_week = find(d.week == w);

            idx = intersect(intersect(idx_subid, idx_symbol), idx_week);

            if ~isempty(idx) % then this subject is not a digit/letter subject so skip

                count = count + 1;

                % Select all data for the first trial for this symbol/participant/week.
                dnew(count, :) = d(idx(1), :);

                % Replace the trial 1 legibility with the average legibility over the 4 trials.
                legibilitytemp(count) = nanmean(d.legibility(idx));

                % Replace the trial 1 confusability with the average confusability over the 4 trials.
                confusabilitytemp(count) = nanmean(d.confusability(idx));

            end

            clear idx_week;

        end

        clear idx_symbol;

    end

    clear idx_subid;

end

% Replace the trial legibility data with the averaged legibility data.
dnew.legibility = legibilitytemp'; clear legibilitytemp;

% Replace the trial confusability data with the averaged confusability data.
dnew.confusability = confusabilitytemp'; clear confusabilitytemp;

% Remove the time, trials, and the labels columns, as they are no longer accurate/needed.
dnew.trial = []; dnew.time = [];
dnew.label1 = []; dnew.label2 = []; dnew.label3 = []; dnew.label4 = []; dnew.label5 = [];

%% Estimate beta and y-intercept for each symbol and each child.
count = 0;
for s = 1:length(subids)

    idx_subid = find(dnew.subid == subids(s));

    for a = 1:length(symbolset)

        idx_symbol = find(strcmp(dnew.symbol, symbolset{a}));

        idx = intersect(idx_subid, idx_symbol);

        if ~isempty(idx)

                count = count + 1;

                % Select all data for the first trial for this symbol/participant/week.
                dnewest(count, :) = dnew(idx(1), :);

                % Select tempdata for linear fitting, calculation of learning slope for this subject and this symbol.
                tempdata = dnew(idx, :);

                % Learning slopes, but only if they have at least two practice days and did not receive perfect scores on all days.
                if length(~isnan(tempdata.legibility)) >= 2 && sum(tempdata.legibility) < length(tempdata.legibility) && sum(tempdata.confusability) > 0

                    % Learning slope for legibility.
                    x = tempdata.week;
                    y = tempdata.legibility; %(tempdata.legibility-nanmean(tempdata.legibility))./nanstd(tempdata.legibility, [], 1);
                    mdl = fitlm(x, y, 'linear');
                    legbeta(count) = round(table2array(mdl.Coefficients(2, 1)), 2);
                    legint(count) = round(table2array(mdl.Coefficients(1, 1)), 2);

                    % Learning slope for confusability.
                    x = tempdata.week;
                    y = tempdata.confusability; %(tempdata.confusability - nanmean(tempdata.confusability))./nanstd(tempdata.confusability, [], 1);
                    mdl = fitlm(x, y, 'linear');
                    conbeta(count) = round(table2array(mdl.Coefficients(2, 1)), 2);
                    conint(count) = round(table2array(mdl.Coefficients(1, 1)), 2);

                else

                    legbeta(count) = NaN;
                    legint(count) = NaN;
                    conbeta(count) = NaN;
                    conint(count) = NaN;

                end

                % Go ahead and grab the day 1 scores for use later. Just easy to grab here.
                idx_week1 = find(tempdata.week == 1);
                if ~isempty(idx_week1)
                    
                    legday1(count) = tempdata.legibility(idx_week1);
                    conday1(count) = tempdata.confusability(idx_week1);

                else

                    legday1(count) = NaN;
                    conday1(count) = NaN;

                end

                clear tempdata idx_week1;

        end

        clear idx idx_symbol;

    end

    clear idx_subid;

end

% Append beta_legibility, int_legibility, beta_confusability, and int_confusability.
legbeta = array2table(legbeta', 'VariableNames', {'legbeta'});
legint = array2table(legint', 'VariableNames', {'legint'});
conbeta = array2table(conbeta', 'VariableNames', {'conbeta'});
conint = array2table(conint', 'VariableNames', {'conint'});
legday1 = array2table(legday1', 'VariableNames', {'legday1'});
conday1 = array2table(conday1', 'VariableNames', {'conday1'});

dnewest = [dnewest legbeta legint legday1 conbeta conint conday1];

% Remove week, legibility, and confusability columns, as they are no longer accurate/relevant.
dnewest.week = []; dnewest.legibility = []; dnewest.confusability= []; 

% Remove image and symboltype for ease later on; image and symboltpe are cells.
dnewest.image = []; dnewest.symboltype = [];

% We predicted that the more accurately a participant wrote a symbol on day one, the less gains we would see in their symbol copying over the 6-week period. 
% In other words, we expect that those with the lowest scores on day one had the greatest opportunity for growth. Results supported our hypothesis.
% X = dnewest(find(~isnan(dnewest.legbeta)), :);
% modelspec = 'legbeta ~ legint + (1|subid) + (1|symbol)';
% mdl = fitlme(X, modelspec); mdl
% [r, p] = corrcoef(dnewest.legbeta(~isnan(dnewest.legbeta)), dnewest.legint(~isnan(dnewest.legint)))
% [r, p] = corrcoef(dnewest.conbeta(~isnan(dnewest.conbeta)), dnewest.conint(~isnan(dnewest.conint)))
% Note: this makes more sense to do after we have collapsed across participants (below).

%% Get at table that lists the mean beta and std for each symbol across participants. Do this for legibility and confusability.

for a = 1:length(symbolset)

    idx = find(strcmp(dnewest.symbol, symbolset{a}));

    legbetamean(a) = nanmean(dnewest.legbeta(idx)); legbeta95ci(a) = 1.96*nanstd(dnewest.legbeta(idx), [], 1)/sqrt(n);
    legintmean(a) = nanmean(dnewest.legint(idx)); legint95ci(a) = 1.96*nanstd(dnewest.legint(idx), [], 1)/sqrt(n);

    conbetamean(a) = nanmean(dnewest.conbeta(idx)); conbeta95ci(a) = 1.96*nanstd(dnewest.conbeta(idx), [], 1)/sqrt(n);
    conintmean(a) = nanmean(dnewest.conint(idx)); conint95ci(a) = 1.96*nanstd(dnewest.conint(idx), [], 1)/sqrt(n);

    legday1here(a) = nanmean(dnewest.legday1(idx)); legday195ci(a) = 1.96*nanstd(dnewest.legday1(idx), [], 1)/sqrt(n);
    conday1here(a) = nanmean(dnewest.conday1(idx)); conday195ci(a) = 1.96*nanstd(dnewest.conday1(idx), [], 1)/sqrt(n);

end
t = array2table([legbetamean', legbeta95ci', legbetamean'-legbeta95ci', legbetamean'+legbeta95ci', ...
    legintmean', legint95ci', legintmean'-legint95ci', legintmean'+legint95ci', ...
    legday1here', legday195ci', ...
    conbetamean', conbeta95ci', conbetamean'-conbeta95ci', conbetamean'+conbeta95ci',...
    conintmean', conint95ci', conintmean'-conint95ci', conintmean'+conint95ci',...
    conday1here',conday195ci'], ...
    'VariableNames', {'legbetamean', 'legbeta95ci', 'legbetalow', 'legbetahigh', ...
    'legintmean', 'legint95ci', 'legintlow', 'leginthigh',...
    'legday1', 'legday195ci', ...
    'conbetamean', 'conbeta95ci', 'conbetalow', 'conbetahigh', ...
    'conintmean', 'conint95ci', 'conintlow', 'coninthigh',...
    'conday1', 'conday195ci'});
t = [array2table(symbolset, 'VariableNames', {'symbol'}) t];

writetable(t, fullfile(rootdir, 'supportfiles', 'table4_betas_ints.csv'))

%% Grab correlation between day1 scores and the int estimates to validate using the int estimates.
[r, p] = corrcoef(legintmean, legday1here);
[r, p] = corrcoef(conintmean, conday1here);

%% Are the betas for legibility highly correlated with the betas for confusability?
[r, p] = corrcoef(t.legbetamean, t.conbetamean)
% yes, r = -0.49, p = 0.003, so just use legibility going forward
% the betas for confusability are also more tightly clustered around each
% other, so the betas for legibility are likely more sensitive

%% Plot relationship between meanbeta and meanint for legibility and confusability.
figcount = figcount + 1; figure(figcount); hold on;
T=text(t.legintmean, t.legbetamean, t.symbol);
for t2 = 1:length(T); T(t2).FontName = "Arial"; T(t2).FontSize=18;  end
[r, p] = corrcoef(t.legintmean, t.legbetamean)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.25 1];
xax.TickValues = [0.25 0.50 0.75 1.00];
xax.TickLabels = {'0.25', '0.50', '0.75', '1.00'};
xax.TickDirection = 'out';
xax.TickLength = [0 0];
xax.FontName = 'Arial';
xax.FontSize = 16;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [-0.025 0.05];
yax.TickValues = [-0.025 0 0.025 0.05];
yax.TickLabels = {'-0.025', '0', '0.025', '0.05'};
yax.TickDirection = 'out';
yax.TickLength = [0 0];
yax.FontName = 'Arial';
yax.FontSize = 16;
% yax.FontAngle = fontangle;

% general
a = gca;
box off

a.XLabel.String = 'Baseline (intercept)';
a.XLabel.FontSize = 16;

a.YLabel.String = 'Learning slope (beta)';
a.YLabel.FontSize = 16;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1]);

print(fullfile(rootdir, 'plots', ['plot_figure5c_legibility_n=' num2str(n)]), '-dpng')

figcount = figcount + 1; figure(figcount); hold on;
T=text(t.conintmean, t.conbetamean, t.symbol);
for t2 = 1:length(T); T(t2).FontName = "Arial"; T(t2).FontSize=18;  end
[r, p] = corrcoef(t.conintmean, t.conbetamean)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 0.30];
xax.TickValues = [0.25 0.50 0.75 1.00];
xax.TickLabels = {'0.25', '0.50', '0.75', '1.00'};
xax.TickDirection = 'out';
xax.TickLength = [0 0];
xax.FontName = 'Arial';
xax.FontSize = 16;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [-0.025 0.025];
yax.TickValues = [-0.025 0 0.025 0.05];
yax.TickLabels = {'-0.025', '0', '0.025', '0.05'};
yax.TickDirection = 'out';
yax.TickLength = [0 0];
yax.FontName = 'Arial';
yax.FontSize = 16/2;
% yax.FontAngle = fontangle;

% general
a = gca;
box off

a.XLabel.String = 'Baseline (intercept)';
a.XLabel.FontSize = 16;

a.YLabel.String = 'Learning slope (beta)';
a.YLabel.FontSize = 16;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1]);

print(fullfile(rootdir, 'plots', ['plot_figure5c_confusability_n=' num2str(n)]), '-dpng')

%% Plot Figure 1a. Rank ordered legibility.
figcount = figcount + 1; figure(figcount); hold on;
 
[~, idx] = sort(t.legbetamean, 'ascend');
acc = t.legbetamean(idx);
accstd = t.legbetastd(idx);
accsym = t.symbol(idx); clear idx;

% scatter(acc, 1:length(acc), 100, 'Filled'); hold on;
for p = 1:length(acc)
      plot([acc(p) - abs(accstd(p)) acc(p) + abs(accstd(p))], [p p], 'Color', 'k')
end
scatter(acc, 1:length(acc), 100, 'Filled', 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'w'); 

% xaxis
xax = get(gca,'xaxis');
xax.Limits = [-0.025 0.05];
xax.TickValues = [-0.025 0 0.025 0.05];
xax.TickLabels = {'-0.025', '0', '0.025', '0.05'};
xax.TickDirection = 'out';
xax.TickLength = [0 0];
xax.FontName = 'Arial';
xax.FontSize = 16;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(acc)+1];
yax.TickValues = 1:length(acc);
yax.TickDirection = 'out';
yax.TickLabels = accsym;
yax.FontName = 'Arial';
yax.FontSize = 16/2;

% general
a = gca;
box off

a.XLabel.String = 'Legibility learning slope (beta)';
a.XLabel.FontSize = 16;

a.YLabel.String = 'Symbol';
a.YLabel.FontSize = 16;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1])

plot([0 0], [-5 65], 'k');

print(fullfile(rootdir, 'plots', ['plot_figure5a_legibility_n=' num2str(n)]), '-dpng')

hold off;

%% Plot Figure 1b. Rank ordered confusability.
figcount = figcount + 1; figure(figcount); hold on;

[~, idx] = sort(conbetamean, 'descend');
con = conbetamean(idx);
constd = conbetastd(idx);
consym = t.symbol(idx); clear idx;

% scatter(con, 1:length(con), 100, 'Filled'); hold on;
for p = 1:length(con)
      plot([con(p) - abs(constd(p)) con(p) + abs(constd(p))], [p p], 'Color', 'k')
end
scatter(con, 1:length(con), 100, 'Filled', 'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'w'); 

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [-0.04 0.04];
xax.TickValues = [-0.04 -0.02 0 0.02 0.04];
xax.TickLabels = {'-0.04', '-0.02', '0', '0.02', '0.04'};
xax.TickDirection = 'out';
xax.TickLength = [0 0];
xax.FontName = 'Arial';
xax.FontSize = 16;
% xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(con)+1];
yax.TickValues = 1:length(con);
yax.TickDirection = 'out';
yax.TickLabels = consym;
yax.FontName = 'Arial';
yax.FontSize = 16/2;

% general
a = gca;
box off

a.XLabel.String = 'Confusability learning slope (beta)';
a.XLabel.FontSize = 16;

a.YLabel.String = 'Symbol';
a.YLabel.FontSize = 16;
a.YLabel.FontAngle = 'italic';
pbaspect([1 1 1]);

plot([0 0], [-5 65], 'k');

print(fullfile(rootdir, 'plots', ['plot_figure5b_confusability_n=' num2str(n)]), '-dpng')

