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

%% Get day 1 confusability matrix.
figcount = figcount + 1; figure(figcount); 
dlet = d;
truelabels = cat(1, dlet.symbol, dlet.symbol, dlet.symbol, dlet.symbol, dlet.symbol);
predictedlabels = cat(1, dlet.label1, dlet.label2, dlet.label3, dlet.label4, dlet.label5);
cm = confusionchart(truelabels, predictedlabels);
cm.RowSummary = 'row-normalized';
cm.XLabel = 'Symbol Label Assigned by MTurk Worker';
cm.YLabel = 'Symbol Copied by the Participant';
sortClasses(cm,'descending-diagonal');
print(fullfile(rootdir, 'plots', ['plot_figure3_lettersandsymbols_n=' num2str(n)]), '-dpng')

figcount = figcount + 1; figure(figcount); %letters
dlet = d(find(strcmp(d.symboltype, 'letter')), :);
truelabels = cat(1, dlet.symbol, dlet.symbol, dlet.symbol, dlet.symbol, dlet.symbol);
predictedlabels = cat(1, dlet.label1, dlet.label2, dlet.label3, dlet.label4, dlet.label5);
cm = confusionchart(truelabels, predictedlabels);
cm.RowSummary = 'row-normalized';
sortClasses(cm,'descending-diagonal');
print(fullfile(rootdir, 'plots', ['plot_figure3_letters_n=' num2str(n)]), '-dpng')

figcount = figcount + 1; figure(figcount); %digits
ddig = d(find(strcmp(d.symboltype, 'digit')), :);
truelabels = cat(1, ddig.symbol, ddig.symbol, ddig.symbol, ddig.symbol, ddig.symbol);
predictedlabels = cat(1, ddig.label1, ddig.label2, ddig.label3, ddig.label4, ddig.label5);
cm = confusionchart(truelabels, predictedlabels);
cm.RowSummary = 'row-normalized';
sortClasses(cm,'descending-diagonal');
print(fullfile(rootdir, 'plots', ['plot_figure3_digits_n=' num2str(n)]), '-dpng')

% 
% %% Get day 1 accuracy for only letters.
% day1accletters = day1con(:, 11:36); % letters only
% keep = ~all(isnan(table2array(day1accletters)), 2); % letters only
% day1accletters = day1accletters(keep, :);
% 
% X = table2array(day1accletters);
% 
% [rho, pmat] = corrcoef(X, 'rows', 'pairwise');
% 
% % D = squareform(pdist(rho));
% 
% % Calculate eigenvalues.
% [Y,eigvals] = cmdscale(1-rho);
% 
% % Calculate normalized eigenvalues.
% neigvals = eigvals/max(abs(eigvals));
% 
% % Display eigenvalues and normalized eigenvalues.
% disp('eigenvals       eigvals/max(abs(eigvals))')
% disp([eigvals neigvals]);
% 
% % % Calculate error.
% % maxerr = max(abs(D - pdist(Y(:,1))))/max(D);
% % disp(['Max relative error for 1D: ' num2str(maxerr) '.']);
% % 
% % maxerr = max(abs(D - pdist(Y(:,1:2))))/max(D);
% % disp(['Max relative error for 2D: ' num2str(maxerr) '.']);
% % 
% % maxerr = max(abs(D - pdist(Y(:,1:3))))/max(D) ;
% % disp(['Max relative error for 3D: ' num2str(maxerr) '.']);
% % 
% % maxerr = max(abs(D - pdist(Y)))/max(D) ;
% % disp(['Max relative error for full reconstruction: ' num2str(maxerr) '.']);
% 
% figcount = figcount + 1; figure(figcount); 
% bar(eigvals)
% 
% 
% n_dim = 2; % can only handle n_dim = 2 at this point (and likely that is all that is needed for this project)
% n_clust = 4; 
% n_rep = 10000;
% % Perform kmeans clustering on the predicted data.
% [clust_idx, C, sumd, D] = kmeans(Y(:, 1:n_dim), n_clust, 'Replicates', n_rep, 'Options', statset('Display', 'final')); % kmeans uses squared Euclidean distances by default
% 
% % Make a mesh; compute the distance from each centroid, C, to points on the grid (procedure found here: https://www.mathworks.com/help/stats/kmeans.html).
% xlimlo = -1; ylimlo = -1;
% xlimhi = 1; ylimhi = 1;
% 
% x = xlimlo:0.01:xlimhi; 
% y = ylimlo:0.01:ylimhi;
% [xG,yG] = meshgrid(x, y);
% XGrid = [xG(:), yG(:)]; % Defines a fine grid on the plot
% idx2Region = kmeans(XGrid, n_clust, 'MaxIter', 1, 'Start', C);
% 
% % Calculate silhouette values.
% s = silhouette(Y(:, 1:n_dim), clust_idx, 'sqEuclidean');
% 
% 
% figcount = figcount + 1; figure(figcount); 
% darkgray = [155 155 155]/255; gray = [200 200 200]/255; lightgray = [100 100 100]/255; grayish = [10 10 10]/255;
% markersize = 20; 
% if n_clust == 4
% 
%     % Add grid.
%     g = gscatter(XGrid(:,1), XGrid(:,2), idx2Region, [darkgray; gray; lightgray; grayish], '.');
%     g(1).MarkerSize = markersize; g(2).MarkerSize = markersize; g(3).MarkerSize = markersize; g(4).MarkerSize = markersize;
% 
% elseif n_clust == 3
%     
%     % Add grid.
%     g = gscatter(XGrid(:,1), XGrid(:,2), idx2Region, [darkgray; gray; lightgray], '.');
%     g(1).MarkerSize = markersize; g(2).MarkerSize = markersize; g(3).MarkerSize = markersize;
%     
% elseif n_clust == 2
%     
%     % Add grid.
%     g = gscatter(XGrid(:,1), XGrid(:,2), idx2Region, [gray; lightgray], '.');
%     g(1).MarkerSize = markersize; g(2).MarkerSize = markersize;
%     
% end
% hold on;
% ylim([-.5 .5]); xlim([-.5 .5]);
% plot(Y(:, 1), Y(:, 2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', markersize, ...
%     'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'w');
% text(Y(:, 1)+0.002, Y(:, 2)+0.002, day1accletters.Properties.VariableNames);

