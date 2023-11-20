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

% Grab the mean for vmskill, literacy, and phonological awareness for Table 1.
disp(['The mean of the Visual-motor Skill composite scores is ' num2str(mean(dnewest.vmskill)) ', with a SD = ' num2str(std(dnewest.vmskill)) '.']);
disp(['The mean of the Literacy Skill composite scores is ' num2str(mean(dnewest.literacy)) ', with a SD = ' num2str(std(dnewest.literacy)) '.']);
disp(['The mean of the Phonological Skill composite scores is ' num2str(mean(dnewest.language)) ', with a SD = ' num2str(std(dnewest.language)) '.']);

%% Model selection and fitting: literacy + language + vmskill.
% linear mixed-effects modeling starting notion:
% fixed effects as predictor variables: standard assessments
% random effects portion corresponds to 1 + subid, because the intercept is included by default
% subid is the grouping variable
% symbol-1|subid term to account for symbol effects within each subid
% include a –1 to eliminate the intercept from the second random-effect

% Ensure that sex and agegroup are set as categorical variables.
dnewest.sex = categorical(dnewest.sex);
dnewest.agegroup = categorical(dnewest.agegroup);

% Be sure to z-score the variables of interest. 
dnewest.legbeta = (dnewest.legbeta - nanmean(dnewest.legbeta))/nanstd(dnewest.legbeta); 
dnewest.legint = (dnewest.legint - nanmean(dnewest.legint))/nanstd(dnewest.legint);
dnewest.vmi = (dnewest.vmi - nanmean(dnewest.vmi))/nanstd(dnewest.vmi); 
dnewest.vp = (dnewest.vp - nanmean(dnewest.vp))/nanstd(dnewest.vp); 
dnewest.mc = (dnewest.mc - nanmean(dnewest.mc))/nanstd(dnewest.mc); 
dnewest.larue = (dnewest.larue - nanmean(dnewest.larue))/nanstd(dnewest.larue); 
dnewest.lk = (dnewest.lk - nanmean(dnewest.lk))/nanstd(dnewest.lk); 
dnewest.rec = (dnewest.rec - nanmean(dnewest.rec))/nanstd(dnewest.rec); 
dnewest.vdis = (dnewest.vdis - nanmean(dnewest.vdis))/nanstd(dnewest.vdis); 
dnewest.sound = (dnewest.sound - nanmean(dnewest.sound))/nanstd(dnewest.sound); 
dnewest.vmskill = (dnewest.vmskill - nanmean(dnewest.vmskill))/nanstd(dnewest.vmskill);
dnewest.literacy = (dnewest.literacy - nanmean(dnewest.literacy))/nanstd(dnewest.literacy); 
dnewest.language = (dnewest.language - nanmean(dnewest.language))/nanstd(dnewest.language); 

% Check to see if age should be included as a covariate of no interest.
X = dnewest(~isnan(dnewest.legbeta), :);

modelspec = 'legbeta ~ agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp
% sex, and age*sex are significant predictors and will be controlled for 

%% Check for multicollinearity among dependent variables before moving forward with the multiple linear mixed-effects modeling.

% Calculate correlation and significance. 
[rval, pval] = corrcoef([dnewest.vmi dnewest.vp dnewest.mc ...
    dnewest.larue dnewest.vdis dnewest.rhyme dnewest.sound ...
    dnewest.vmskill dnewest.literacy dnewest.language ...
    dnewest.age]);

% Plot.
figure(1); imagesc(rval); colorbar; caxis([-1 1]);
labels = {'vmi', 'vp', 'mc', ...
    'larue', 'vdis', 'rhyme', 'sound', ...
    'vmskill', 'literacy', 'language', ...
    'age'};
xticks(1:length(labels)); xticklabels(labels);
yticks(1:length(labels)); yticklabels(labels);

% Export tables.
out = array2table(rval);
out.Properties.VariableNames = labels;
writetable(out, fullfile(rootdir, 'output', 'table1_rval.csv'));
clear out;

out = array2table(pval);
out.Properties.VariableNames = labels;
writetable(out, fullfile(rootdir, 'output', 'table1_pval.csv'));
clear out;

% Calculate correlation and significance, controlling for age. 
[rval, pval] = partialcorr([dnewest.vmi dnewest.vp dnewest.mc ...
    dnewest.larue dnewest.vdis dnewest.rhyme dnewest.sound ...
    dnewest.vmskill dnewest.literacy dnewest.language], ...
    dnewest.age);

% Plot.
figure(1); imagesc(rval); colorbar; caxis([-1 1]);
labels = {'vmi', 'vp', 'mc', ...
    'larue', 'vdis', 'rhyme', 'sound', ...
    'vmskill', 'literacy', 'language'};
xticks(1:length(labels)); xticklabels(labels);
yticks(1:length(labels)); yticklabels(labels);

% Export tables.
out = array2table(rval);
out.Properties.VariableNames = labels;
writetable(out, fullfile(rootdir, 'output', 'table1_rval_agecontrolled.csv'));
clear out;

out = array2table(pval);
out.Properties.VariableNames = labels;
writetable(out, fullfile(rootdir, 'output', 'table1_pval_agecontrolled.csv'));
clear out;





% high correlation among predictors, so concerned about multicolinearity.




%% Analysis 1, legibility: 
%% Perform the planned mixed-model linear regression to demonstrate that inital legibility/confusability predicts change in legibility/confusability.
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme1 = fitlme(X, 'legbeta ~ legint + agegroup*sex + (1|symbol) + (1|subid)'); lme1

figure(1); plotResiduals(lme1, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[Pa, Fa, DF1a, DF2a] = coefTest(lme1)
% test significance of the overall model
% significant, F(4, 534) = 57.69, p = 1.85e-40, R2 = 

anova(lme1,'dfmethod','satterthwaite')
% gets p-values for inference of fixed effects across levels of the random effect

% Linear model with only agegroup = 1 (younger children).
X = dnewest(~isnan(dnewest.legbeta), :);
X = X(find(double(X.agegroup) == 1), :);

removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme1 = fitlme(X, 'legbeta ~ legint + sex + (1|symbol) + (1|subid)'); lme1

[Pa, Fa, DF1a, DF2a] = coefTest(lme1)
anova(lme1,'dfmethod','satterthwaite')

% independent samples t-test
[H,P,CI,STATS] = ttest2(X.legbeta(find(double(X.sex)==1), :), X.legbeta(find(double(X.sex)==2), :))

X = dnewest(~isnan(dnewest.legbeta), :);
X = X(find(double(X.agegroup) == 1), :);
mfy = mean(X.legbeta(find(double(X.sex)==1))); sdfy = std(X.legbeta(find(double(X.sex)==1))); cify = 1.96*(sdfy/sqrt(length(X.legbeta(find(double(X.sex)==1)))));
mmy = mean(X.legbeta(find(double(X.sex)==2))); sdmy = std(X.legbeta(find(double(X.sex)==2))); cimy = 1.96*(sdmy/sqrt(length(X.legbeta(find(double(X.sex)==2)))));

clear X;

% Linear model with only agegroup = 2 (older children).
X = dnewest(~isnan(dnewest.legbeta), :);
X = X(find(double(X.agegroup) == 2), :);

removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme1 = fitlme(X, 'legbeta ~ legint + sex + (1|symbol) + (1|subid)'); lme1

[Pa, Fa, DF1a, DF2a] = coefTest(lme1)
anova(lme1,'dfmethod','satterthwaite')

% independent samples t-test
[H,P,CI,STATS] = ttest2(X.legbeta(find(double(X.sex)==1), :), X.legbeta(find(double(X.sex)==2), :))

X = dnewest(~isnan(dnewest.legbeta), :);
X = X(find(double(X.agegroup) == 2), :);
mfo = mean(X.legbeta(find(double(X.sex)==1))); sdfo = std(X.legbeta(find(double(X.sex)==1))); cifo = 1.96*(sdfo/sqrt(length(X.legbeta(find(double(X.sex)==1)))));
mmo = mean(X.legbeta(find(double(X.sex)==2))); sdmo = std(X.legbeta(find(double(X.sex)==2))); cimo = 1.96*(sdmo/sqrt(length(X.legbeta(find(double(X.sex)==2)))));

clear X;

% Plot result: For younger children, female legibility slopes were greater
% than male legibility slopes. For older children, there was no significant
% difference in legibilty slopes between females and males.
a = [1 2 3 4]; 
b = [mfy mmy mfo mmo];
c = [cify cimy cifo cimo];

figure(1); hold on;
ba = bar(a, b);
for i = 1:length(a); plot([a(i) a(i)], [b(i)-c(i) b(i)+c(i)], 'k-', 'LineWidth', 1.5);
end
ba.FaceColor = 'w';
ba.EdgeColor = 'k';
ba.BarWidth = 0.5;
ba.LineWidth = 1.5;
ba.FaceAlpha = 0.4;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = 0.5; xlim_hi = 4.5;
ylim_lo = -0.5; ylim_hi = 0.5;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = a;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = {'Female', 'Male', 'Female', 'Male'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = 'Younger                               Older';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [-0.50 -0.25 0 0.25 0.50];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'-0.50', '-0.25', '0', '0.25', '0.50'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = 'Legibility Score';

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
yax.Label.String = {'Slope for Legibility Score'; '(over 6 weeks)'};

print(fullfile(rootdir, 'plots', 'plot_figure4a_agesexinteraction'), '-dpng')

%% Analysis 1, confusability: 
%% Perform the planned mixed-model linear regression to demonstrate that inital legibility/confusability predicts change in legibility/confusability.
X = dnewest(~isnan(dnewest.conbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'conbeta ~ conint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme1 = fitlme(X, 'conbeta ~ conint + agegroup*sex + (1|symbol) + (1|subid)'); lme1

figure(1); plotResiduals(lme1, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[Pa, Fa, DF1a, DF2a] = coefTest(lme1)
% test significance of the overall model
% significant, F(4, 534) = 57.69, p = 1.85e-40, R2 = 

anova(lme1,'dfmethod','satterthwaite')

% Linear model with only agegroup = 1 (younger children).
X = dnewest(~isnan(dnewest.legbeta), :);
X = X(find(double(X.agegroup) == 1), :);

removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme1 = fitlme(X, 'conbeta ~ conint + sex + (1|symbol) + (1|subid)'); lme1

[Pa, Fa, DF1a, DF2a] = coefTest(lme1)
anova(lme1,'dfmethod','satterthwaite')

% independent samples t-test
[H,P,CI,STATS] = ttest2(X.conbeta(find(double(X.sex)==1), :), X.conbeta(find(double(X.sex)==2), :))

X = dnewest(~isnan(dnewest.conbeta), :);
X = X(find(double(X.agegroup) == 1), :);
mfy = mean(X.conbeta(find(double(X.sex)==1))); sdfy = std(X.conbeta(find(double(X.sex)==1))); cify = 1.96*(sdfy/sqrt(length(X.conbeta(find(double(X.sex)==1)))));
mmy = mean(X.conbeta(find(double(X.sex)==2))); sdmy = std(X.conbeta(find(double(X.sex)==2))); cimy = 1.96*(sdmy/sqrt(length(X.conbeta(find(double(X.sex)==2)))));

clear X;

% Linear model with only agegroup = 2 (older children).
X = dnewest(~isnan(dnewest.legbeta), :);
X = X(find(double(X.agegroup) == 2), :);

removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme1 = fitlme(X, 'conbeta ~ conint + sex + (1|symbol) + (1|subid)'); lme1

[Pa, Fa, DF1a, DF2a] = coefTest(lme1)
anova(lme1,'dfmethod','satterthwaite')

% independent samples t-test
[H,P,CI,STATS] = ttest2(X.conbeta(find(double(X.sex)==1), :), X.conbeta(find(double(X.sex)==2), :))

X = dnewest(~isnan(dnewest.conbeta), :);
X = X(find(double(X.agegroup) == 2), :);
mfo = mean(X.conbeta(find(double(X.sex)==1))); sdfo = std(X.conbeta(find(double(X.sex)==1))); cifo = 1.96*(sdfo/sqrt(length(X.conbeta(find(double(X.sex)==1)))));
mmo = mean(X.conbeta(find(double(X.sex)==2))); sdmo = std(X.conbeta(find(double(X.sex)==2))); cimo = 1.96*(sdmo/sqrt(length(X.conbeta(find(double(X.sex)==2)))));

clear X;

% Plot result: For younger children, female legibility slopes were greater
% than male legibility slopes. For older children, there was no significant
% difference in legibilty slopes between females and males.
a = [1 2 3 4]; 
b = [mfy mmy mfo mmo];
c = [cify cimy cifo cimo];

figure(1); hold on;
ba = bar(a, b);
for i = 1:length(a); plot([a(i) a(i)], [b(i)-c(i) b(i)+c(i)], 'k-', 'LineWidth', 1.5);
end
ba.FaceColor = 'w';
ba.EdgeColor = 'k';
ba.BarWidth = 0.5;
ba.LineWidth = 1.5;
ba.FaceAlpha = 0.4;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = 0.5; xlim_hi = 4.5;
ylim_lo = -0.01; ylim_hi = 0.02;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = a;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = {'Female', 'Male', 'Female', 'Male'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = 'Younger                               Older';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [-0.01 0 0.01 0.02];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'-0.01', '0', '0.01', '0.02'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = {'Slope for Confusability Score'; '(over 6 weeks)'};

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
% title('Interaction between Age and Sex for Confusability');

print(fullfile(rootdir, 'plots', 'plot_figure4b_agesexinteraction'), '-dpng')

