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

% high correlation among predictors, so concerned about multicolinearity, checks below.

%% Assess need for multi-level model.
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp

% Add residuals to the data table.
X = [array2table(residuals(lmetemp), 'VariableNames', {'residuals'}) X];
clear lmetemp;

% Compare literacy and vmskill in predicting beta, accounting for intercept base and random factors.
lme_intonly = fitlme(X, 'residuals ~ 1'); lme_intonly

lme_intrandsub = fitlme(X, 'residuals ~ 1 + (1|subid)'); lme_intrandsub
compare(lme_intonly, lme_intrandsub)

lme_intrandsym = fitlme(X, 'residuals ~ 1 + (1|symbol)'); lme_intrandsym
compare(lme_intonly, lme_intrandsym)

lme_intrandsubsym = fitlme(X, 'residuals ~ 1 + (1|subid) + (1|symbol)'); lme_intrandsubsym
compare(lme_intonly, lme_intrandsubsym)
compare(lme_intrandsub, lme_intrandsubsym)
compare(lme_intrandsym, lme_intrandsubsym)

%% Step 1: VM Skill predicting learning
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp

% Add residuals to the data table.
X = [array2table(residuals(lmetemp), 'VariableNames', {'residuals'}) X];
clear lmetemp;

%[an aside to get male and female means
idx_yf = find(double(X.agegroup) == 1 & double(X.sex) == 1);
legslope_yf = X.residuals(idx_yf);
subid_yf = unique(X.subid(idx_yf));

for i = 1:length(subid_yf)

    idx_yf2 = find(subid_yf == subid_yf(i));
    legout_yf(i) = mean(legslope_yf(idx_yf2));

end

idx_ym = find(double(X.agegroup) == 1 & double(X.sex) == 2);
legslope_ym = X.residuals(idx_ym);
subid_ym = unique(X.subid(idx_ym));

for i = 1:length(subid_ym)

    idx_ym2 = find(subid_ym == subid_ym(i));
    legout_ym(i) = mean(legslope_yf(idx_ym2));

end
subout = [subid_yf; subid_ym]; legout = [legout_yf'; legout_ym'];
dout = array2table([subout legout], "VariableNames", {'subid', 'legibility'});
writetable(dout, fullfile(rootdir, 'output', 'youngmaleandfemale_legibility_sex.csv'));

% aside ending here]

% Compare literacy and vmskill in predicting beta, accounting for intercept base and random factors.
modelspec = 'residuals ~ vmskill';
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme_vmskill = fitlme(X, 'residuals ~ vmskill + (1|subid) + (1|symbol)'); lme_vmskill
% fixed: intercept is significant, p = 0.0001, vmskill is significant, p = 0.0002
% random: subid and symbol are significant

% figure(1); plotResiduals(lme_vmskill, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[Pa, Fa, DF1a, DF2a] = coefTest(lme_vmskill)
% test significance of the overall model
% significant, F = 13.91, p = 2.14e-4

anova(lme_vmskill,'dfmethod','satterthwaite')
% gets p-values for inference of fixed effects across levels of the random effect
% vmskill is significant

%% Step 2: Does adding literacy skill and phonological skill to the VMSKILL model improve the fit? 
%% and which of the three are most predictive of learning?
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp

% Add residuals to the data table.
X = [array2table(residuals(lmetemp), 'VariableNames', {'residuals'}) X];
clear lmetemp;

% Compare literacy and vmskill in predicting beta, accounting for intercept base and random factors.
modelspec = 'residuals ~ vmskill';
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme_vmskill_litskill_phonskill = fitlme(X, 'residuals ~ vmskill + literacy + language + (1|subid) + (1|symbol)'); lme_vmskill_litskill_phonskill
% fixed: intercept is significant, p = 0.0001, vmskill is significant, p = 0.0002
% random: subid and symbol are significant

figure(1); plotResiduals(lme_vmskill_litskill_phonskill, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[Pa, Fa, DF1a, DF2a] = coefTest(lme_vmskill_litskill_phonskill)
% test significance of the overall model
% significant, F = 13.91, p = 2.14e-4

anova(lme_vmskill_litskill_phonskill,'dfmethod','satterthwaite')
% gets p-values for inference of fixed effects across levels of the random effect
% vmskill is significant

%% Step 3: Perform the planned mixed-model linear regression to test relationship between vmskill, literacy, and symbol production learning.
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp

% Add residuals to the data table.
X = [array2table(residuals(lmetemp), 'VariableNames', {'residuals'}) X];
clear lmetemp;

% Compare literacy and vmskill in predicting beta, accounting for intercept base and random factors.
modelspec = 'residuals ~ vmskill*literacy*language'; %vmskill*literacy';
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme2 = fitlme(X, 'residuals ~ vmskill*literacy*language + (1|symbol) + (1|subid)'); lme2

% fixed: intercept is significant, p = 0.005, vmskill is significant, p =
% 0.04, literacy is not significant, p = 0.39, literacy*vmskill is
% significant, p = 5.6e-8
% random: subid and symbol are significant

figure(1); plotResiduals(lme2, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[Pa, Fa, DF1a, DF2a] = coefTest(lme2)
% test significance of the overall model
% significant, F = 8.47, p = 2.42e-4

anova(lme2,'dfmethod','satterthwaite')
% gets p-values for inference of fixed effects across levels of the random effect

% [Pa, Fa, DF1a, DF2a] = coefTest(lme2, [0,1,1, 0], 0 ,'dfmethod', 'satterthwaite') % should test if literacy + vmskill = 0, i.e., [0 1 -1] contrast
% Compare literacy and vmskill to see which is the better predictor
% vmskill is a stronger predictor than literacy, F(1, 36) = 8.40, p = 0.0063. 

%% ANALYSIS 4: Interrogate the VM*LIT*PHON interaction.
%% PLOT: VM vs RES

% Define variable of most interest: vmskill.
% x = X.vmskill;

% Define first interacting variable: phonological skill.
z = X.language;
zh = mean(z) + std(z);
zl = mean(z) - std(z);

% Define second interacting varaible: literacy skill.
w = X.literacy;
wh = mean(w) + std(w);
wl = mean(w) - std(w);

% Get coefficients from regression.
cint = double(lme2.Coefficients(1,2)); %intercept
cx = double(lme2.Coefficients(4,2)); %vmskill
cz = double(lme2.Coefficients(3,2)); %phonology
cw = double(lme2.Coefficients(2,2)); %literacy
cxz = double(lme2.Coefficients(7,2));%vmskill, phon
cxw = double(lme2.Coefficients(6,2));%vmskill, lit
czw = double(lme2.Coefficients(5,2));%phon, lit
cxzw = double(lme2.Coefficients(8,2));%vmskill, phon, lit

% Create slopes.
zhwh = cx  + zh*cxz + wh*cxw + zh*zh*cxzw;
zhwl = cx + zh*cxz + wl*cxw + zh*wl*cxzw;
zlwh = cx + zl*cxz + wh*cxw + zl*wh*cxzw;
zlwl = cx + zl*cxz + wl*cxw + zl*wl*cxzw;

% Create intercepts.
izhwh = cint + zh*cz + wh*cw + zh*wh*czw;
izhwl = cint + zh*cz + wl*cw + zh*wl*czw;
izlwh = cint + zl*cz + wh*cw + zl*wh*czw;
izlwl = cint + zl*cz + wl*cw + zl*wl*czw;

% Create predicted y-values.
x0 = min(X.vmskill):0.005:max(X.vmskill);
HHy = izhwh + x0.*zhwh; %high phon, high lit
HLy = izhwl + x0.*zhwl; %high phon, low lit
LHy = izlwh + x0.*zlwh; %low phon, high lit
LLy = izlwl + x0.*zlwl; %low phon, low lit

% Plot.
close all;
figure(2); hold on;
s = scatter(X.vmskill, X.residuals);
% s.CData = [0 0 0]/255;
s.Marker = 'o';
s.MarkerEdgeColor = [255 255 255]/255;
s.MarkerFaceColor = [150 150 150]/255;
s.MarkerFaceAlpha = .5;
s.SizeData = 36;

p4 = plot(x0, LLy); %low phon, low lit, black
p4.LineWidth = 3;
p4.LineStyle = '-';
p4.Color = [0 0 0]/255;

p2 = plot(x0, HLy, 'g'); %high phon, low lit, black
p2.LineWidth = 3;
p2.LineStyle = '--';
p2.Color = [0 0 0]/255;

p3 = plot(x0, LHy); %low phon, high lit, gray
p3.LineWidth = 3;
p3.LineStyle = '-.';
p3.Color = [0 0 0]/255;

p1 = plot(x0, HHy); %high phon, high lit, gray
p1.LineWidth = 3;
p1.LineStyle = ':';
p1.Color = [0 0 0]/255;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = -2; xlim_hi = 2.5;
ylim_lo = -1.5; ylim_hi = 1.5;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1.5 -1 -.5 0 .5 1 1.5 2 2.5];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = {'-2.0', '-1.5', '-1.0', '-0.5', '0', '0.5', '1.0', '1.5', '2.0', '2.5'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = 'Visual-motor Skill Composite Score (z-scored)';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [-1.5 -1.0 -0.5 0 .5 1 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'-1.5', '-1.0', '-0.5', '0', '0.5', '1.0', '1.5'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = {'Slope for Legibility Score'; '(over 6 weeks, z-scored, adjusted)'};

l = legend({'', 'Low literacy, Low phonological', 'Low literacy, High phonological',  ...
    'High literacy, Low phonological', 'High literacy, High phonological'});
l.FontName = 'Arial';
l.FontSize = 12;
l.Location = 'northeast';
l.Orientation = 'vertical';

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
% title('Interaction between Age and Sex for Confusability');

print(fullfile(rootdir, 'plots', 'plot_figure5_vmlitphoninteraction_vmres'), '-dpng')

%% ANALYSIS 4: Interrogate the VM*LIT*PHON interaction.
%% PLOT: LIT vs RES

% Define variable of most interest: vmskill.
% x = X.literacy;

% Define first interacting variable: phonological skill.
z = X.language;
zh = mean(z) + std(z);
zl = mean(z) - std(z);

% Define second interacting varaible: literacy skill.
w = X.vmskill;
wh = mean(w) + std(w);
wl = mean(w) - std(w);

% Get coefficients from regression.
cint = double(lme2.Coefficients(1,2)); %intercept
cx = double(lme2.Coefficients(2,2)); %literacy
cz = double(lme2.Coefficients(3,2)); %phonology
cw = double(lme2.Coefficients(4,2)); %vmskill
cxz = double(lme2.Coefficients(5,2));%literacy, phon 
cxw = double(lme2.Coefficients(6,2));%lit, vmskill
czw = double(lme2.Coefficients(7,2));%phon, vmskill
cxzw = double(lme2.Coefficients(8,2));%vmskill, phon, lit

% Create slopes.
zhwh = cx  + zh*cxz + wh*cxw + zh*zh*cxzw;
zhwl = cx + zh*cxz + wl*cxw + zh*wl*cxzw;
zlwh = cx + zl*cxz + wh*cxw + zl*wh*cxzw;
zlwl = cx + zl*cxz + wl*cxw + zl*wl*cxzw;

% Create intercepts.
izhwh = cint + zh*cz + wh*cw + zh*wh*czw;
izhwl = cint + zh*cz + wl*cw + zh*wl*czw;
izlwh = cint + zl*cz + wh*cw + zl*wh*czw;
izlwl = cint + zl*cz + wl*cw + zl*wl*czw;

% Create predicted y-values.
x0 = min(X.literacy):0.005:max(X.literacy);
HHy = izhwh + x0.*zhwh; %high phon, high vmskill
HLy = izhwl + x0.*zhwl; %high phon, low vmskill
LHy = izlwh + x0.*zlwh; %low phon, high vmskill
LLy = izlwl + x0.*zlwl; %low phon, low vmskill

% Plot.
figure(3); hold on;
s = scatter(X.literacy, X.residuals);
% s.CData = [0 0 0]/255;
s.Marker = 'o';
s.MarkerEdgeColor = [255 255 255]/255;
s.MarkerFaceColor = [150 150 150]/255;
s.MarkerFaceAlpha = .5;
s.SizeData = 36;

p4 = plot(x0, LLy); %low phon, low vmskill, black
p4.LineWidth = 3;
p4.LineStyle = '-';
p4.Color = [0 0 0]/255;

p2 = plot(x0, HLy, 'g'); %high phon, low vmskill, black
p2.LineWidth = 3;
p2.LineStyle = '--';
p2.Color = [0 0 0]/255;

p3 = plot(x0, LHy); %low phon, high vmskill, gray
p3.LineWidth = 3;
p3.LineStyle = '-.';
p3.Color = [0 0 0]/255;

p1 = plot(x0, HHy); %high phon, high vmskill, gray
p1.LineWidth = 3;
p1.LineStyle = ':';
p1.Color = [0 0 0]/255;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = -2; xlim_hi = 2.5;
ylim_lo = -1.5; ylim_hi = 1.5;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1.5 -1 -.5 0 .5 1 1.5 2 2.5];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = {'-2.0', '-1.5', '-1.0', '-0.5', '0', '0.5', '1.0', '1.5', '2.0', '2.5'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = 'Literacy Composite Score (z-scored)';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [-1.5 -1.0 -0.5 0 .5 1 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'-1.5', '-1.0', '-0.5', '0', '0.5', '1.0', '1.5'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = {'Slope for Legibility Score'; '(over 6 weeks, z-scored, adjusted)'};

l = legend({'', 'Low vmskill, Low phonological', 'Low vmskill, High phonological',  ...
    'High vmskill, Low phonological', 'High vmskill, High phonological'});
l.FontName = 'Arial';
l.FontSize = 12;
l.Location = 'northeast';
l.Orientation = 'vertical';

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
% title('Interaction between Age and Sex for Confusability');

print(fullfile(rootdir, 'plots', 'plot_figure5_vmlitphoninteraction_litres'), '-dpng')


%% ANALYSIS 4: Interrogate the VM*LIT*PHON interaction.
%% PLOT: PHON vs RES

% Define variable of most interest: vmskill.
% x = X.language;

% Define first interacting variable: phonological skill.
z = X.literacy;
zh = mean(z) + std(z);
zl = mean(z) - std(z);

% Define second interacting varaible: literacy skill.
w = X.vmskill;
wh = mean(w) + std(w);
wl = mean(w) - std(w);

% Get coefficients from regression.
cint = double(lme2.Coefficients(1,2)); %intercept
cx = double(lme2.Coefficients(3,2)); %phonology
cz = double(lme2.Coefficients(2,2)); %literacy
cw = double(lme2.Coefficients(4,2)); %vmskill
cxz = double(lme2.Coefficients(5,2));%literacy, phon 
cxw = double(lme2.Coefficients(7,2));%phon, vmskill
czw = double(lme2.Coefficients(6,2));%lit, vmskill
cxzw = double(lme2.Coefficients(8,2));%vmskill, phon, lit

% Create slopes.
zhwh = cx  + zh*cxz + wh*cxw + zh*zh*cxzw;
zhwl = cx + zh*cxz + wl*cxw + zh*wl*cxzw;
zlwh = cx + zl*cxz + wh*cxw + zl*wh*cxzw;
zlwl = cx + zl*cxz + wl*cxw + zl*wl*cxzw;

% Create intercepts.
izhwh = cint + zh*cz + wh*cw + zh*wh*czw;
izhwl = cint + zh*cz + wl*cw + zh*wl*czw;
izlwh = cint + zl*cz + wh*cw + zl*wh*czw;
izlwl = cint + zl*cz + wl*cw + zl*wl*czw;

% Create predicted y-values.
x0 = min(X.language):0.005:max(X.language);
HHy = izhwh + x0.*zhwh; %high phon, high vmskill
HLy = izhwl + x0.*zhwl; %high phon, low vmskill
LHy = izlwh + x0.*zlwh; %low phon, high vmskill
LLy = izlwl + x0.*zlwl; %low phon, low vmskill

% Plot.
figure(4); hold on;
s = scatter(X.language, X.residuals);
% s.CData = [0 0 0]/255;
s.Marker = 'o';
s.MarkerEdgeColor = [255 255 255]/255;
s.MarkerFaceColor = [150 150 150]/255;
s.MarkerFaceAlpha = .5;
s.SizeData = 36;

p4 = plot(x0, LLy); %low lit, low vmskill, black
p4.LineWidth = 3;
p4.LineStyle = '-';
p4.Color = [0 0 0]/255;

p2 = plot(x0, HLy, 'g'); %high lit, low vmskill, black
p2.LineWidth = 3;
p2.LineStyle = '--';
p2.Color = [0 0 0]/255;

p3 = plot(x0, LHy); %low lit, high vmskill, gray
p3.LineWidth = 3;
p3.LineStyle = '-.';
p3.Color = [0 0 0]/255;

p1 = plot(x0, HHy); %high lit, high vmskill, gray
p1.LineWidth = 3;
p1.LineStyle = ':';
p1.Color = [0 0 0]/255;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = -2; xlim_hi = 2.5;
ylim_lo = -1.5; ylim_hi = 1.5;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [-2 -1.5 -1 -.5 0 .5 1 1.5 2 2.5];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = {'-2.0', '-1.5', '-1.0', '-0.5', '0', '0.5', '1.0', '1.5', '2.0', '2.5'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = 'Phonology Composite Score (z-scored)';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [-1.5 -1.0 -0.5 0 .5 1 1.5];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'-1.5', '-1.0', '-0.5', '0', '0.5', '1.0', '1.5'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = {'Slope for Legibility Score'; '(over 6 weeks, z-scored, adjusted)'};

l = legend({'', 'Low vmskill, Low literacy', 'Low vmskill, High literacy',  ...
    'High vmskill, Low literacy', 'High vmskill, High literacy'});
l.FontName = 'Arial';
l.FontSize = 12;
l.Location = 'northeast';
l.Orientation = 'vertical';

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
% title('Interaction between Age and Sex for Confusability');

print(fullfile(rootdir, 'plots', 'plot_figure5_vmlitphoninteraction_phonres'), '-dpng')

%% Be careful about multi-collinearity: 
%% Perform individual regressions for each predictor to confirm that it matches with the full model (re: multicolinearity).

% LITERACY
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp

% Add residuals to the data table.
X = [array2table(residuals(lmetemp), 'VariableNames', {'residuals'}) X];
clear lmetemp;

% Compare literacy and vmskill in predicting beta, accounting for intercept base and random factors.
modelspec = 'residuals ~ literacy';
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme_literacy = fitlme(X, 'residuals ~ literacy + (1|subid) + (1|symbol)'); lme_literacy
% fixed: intercept is significant, p = 0.0008, literacy is significant, p = 0.00003
% random: subid is significant, ci for symbol cannot be calculated so
% likely overparameterized

figure(1); plotResiduals(lme_literacy, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[p,F,r] = coefTest(lme_literacy)
% test significance of the overall model
% significant, F = 12.96, p = 3.51e-4

anova(lme_literacy,'dfmethod','satterthwaite')
% gets p-values for inference of fixed effects across levels of the random effect
% literacy is significant

% PHONOLOGY
X = dnewest(~isnan(dnewest.legbeta), :);

% Control for the relationship between starting point (interceptbase) and learning rate (beta).
modelspec = 'legbeta ~ legint + agegroup*sex'; %doesn't make sense to remove at this point
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lmetemp = fitlme(X, modelspec); lmetemp

% Add residuals to the data table.
X = [array2table(residuals(lmetemp), 'VariableNames', {'residuals'}) X];
clear lmetemp;

% Compare literacy and vmskill in predicting beta, accounting for intercept base and random factors.
modelspec = 'residuals ~ language';
removeidx = findoutliers_lm(X, modelspec);
X(removeidx, :) = [];
lme_language = fitlme(X, 'residuals ~ language + (1|subid) + (1|symbol)'); lme_language
% fixed: intercept is significant, p = 0.0001, vmskill is significant, p = 0.0002
% random: subid and symbol are significant

figure(1); plotResiduals(lme_language, 'fitted') 
% has a small artifact in residuals, but generally good and ok to use

[p,F,r] = coefTest(lme_language)
% test significance of the overall model
% significant, F = 13.91, p = 2.14e-4

anova(lme_language,'dfmethod','satterthwaite')
% gets p-values for inference of fixed effects across levels of the random effect
% vmskill is significant

%% Other random checks.
% Validation: VMI should correlate with age.
scatter(letters.vmi, letters.age)
modelspec = 'vmi ~ age';
removeidx = findoutliers_lm(letters, modelspec);
letters(removeidx, :) = [];
lmtemp = fitlm(letters, modelspec); lmtemp

% Does VMI correlate with meanbeta?
scatter(letters.vmi, letters.meanbeta)
modelspec = 'vmi ~ age + meanbeta';
removeidx = findoutliers_lm(letters, modelspec);
letters(removeidx, :) = [];
lmtemp = fitlm(letters, modelspec); lmtemp

% Does VMI correlate with meanbeta?
scatter(letters.age, letters.meanbeta)
modelspec = 'meanbeta ~ age';
removeidx = findoutliers_lm(letters, modelspec);
letters(removeidx, :) = [];
lmtemp = fitlm(letters, modelspec); lmtemp

% Does VMI correlate with meanbeta?
scatter(letters.meanint, letters.meanbeta)
modelspec = 'meanbeta ~ meanint';
removeidx = findoutliers_lm(letters, modelspec);
letters(removeidx, :) = [];
lmtemp = fitlm(letters, modelspec); lmtemp

% Does VMI correlate with meanbeta?
modelspec = 'meanbeta ~ vmi + meanint + age';
removeidx = findoutliers_lm(letters, modelspec);
letters(removeidx, :) = [];
lmtemp = fitlm(letters, modelspec); lmtemp