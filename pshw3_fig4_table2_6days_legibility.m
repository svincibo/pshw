clear all; close all; clc;

rootdir = '/Volumes/Seagate/project-preschool-handwriting';
d = readtable(fullfile(rootdir, 'supportfiles', 'pshw_mturkdata_behdata_n37_20230929.csv'));
figcount = 0;
alphabet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', ...
    'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', ...
    'W', 'X', 'Y', 'Z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

% Calculate means and standard deviations for each symbol.
day1idx = find(d.week == 1); day2idx = find(d.week == 2); day3idx = find(d.week == 3); 
day4idx = find(d.week == 4); day5idx = find(d.week == 5); day6idx = find(d.week == 6); 
yidx = find(d.agegroup == 1); oidx = find(d.agegroup == 2);

nall = length(unique(d.subid)); 
ny = length(unique(d.subid(find(d.agegroup == 1))));
no = length(unique(d.subid(find(d.agegroup == 2))));

for i = 1:length(alphabet)

    idx_symbol = find(strcmp(d.symbol, alphabet{i}));

    day1(i) = nanmean(d.legibility(intersect(idx_symbol, day1idx)));
    day1std(i) = nanstd(d.legibility(intersect(idx_symbol, day1idx)));
    day2(i) = nanmean(d.legibility(intersect(idx_symbol, day2idx)));
    day2std(i) = nanstd(d.legibility(intersect(idx_symbol, day2idx)));
    day3(i) = nanmean(d.legibility(intersect(idx_symbol, day3idx)));
    day3std(i) = nanstd(d.legibility(intersect(idx_symbol, day3idx)));
    day4(i) = nanmean(d.legibility(intersect(idx_symbol, day4idx)));
    day4std(i) = nanstd(d.legibility(intersect(idx_symbol, day4idx)));
    day5(i) = nanmean(d.legibility(intersect(idx_symbol, day5idx)));
    day5std(i) = nanstd(d.legibility(intersect(idx_symbol, day5idx)));
    day6(i) = nanmean(d.legibility(intersect(idx_symbol, day6idx)));
    day6std(i) = nanstd(d.legibility(intersect(idx_symbol, day6idx)));

    day1y(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, yidx), day1idx))); 
    day1ystd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, yidx), day1idx))); 
    day2y(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, yidx), day2idx)));
    day2ystd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, yidx), day2idx)));
    day3y(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, yidx), day3idx)));
    day3ystd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, yidx), day3idx)));
    day4y(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, yidx), day4idx)));
    day4ystd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, yidx), day4idx)));
    day5y(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, yidx), day5idx)));
    day5ystd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, yidx), day5idx)));
    day6y(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, yidx), day6idx)));
    day6ystd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, yidx), day6idx)));

    day1o(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, oidx), day1idx))); 
    day1ostd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, oidx), day1idx)));
    day2o(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, oidx), day2idx)));
    day2ostd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, oidx), day2idx)));
    day3o(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, oidx), day3idx)));
    day3ostd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, oidx), day3idx)));
    day4o(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, oidx), day4idx)));
    day4ostd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, oidx), day4idx)));
    day5o(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, oidx), day5idx)));
    day5ostd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, oidx), day5idx)));
    day6o(i) = nanmean(d.legibility(intersect(intersect(idx_symbol, oidx), day6idx)));
    day6ostd(i) = nanstd(d.legibility(intersect(intersect(idx_symbol, oidx), day6idx)));

    clear idx_symbol;

end

% Bar plot, w/confidence intervals: Letters and numbers.
figcount = figcount + 1; figure(figcount); hold on;

subplot(1, 2, 1); hold on; % letters
keep = 1:26;

xtemp = repmat(1, [1 length(day1(keep))]); ytemp = day1(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(2, [1 length(day2(keep))]); ytemp = day2(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(3, [1 length(day3(keep))]); ytemp = day3(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(4, [1 length(day4(keep))]); ytemp = day4(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(5, [1 length(day5(keep))]); ytemp = day5(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(6, [1 length(day6(keep))]); ytemp = day6(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;

xtemp = repmat(8, [1 length(day1y(keep))]); ytemp = day1y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(9, [1 length(day2y(keep))]); ytemp = day2y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(10, [1 length(day3y(keep))]); ytemp = day3y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(11, [1 length(day4y(keep))]); ytemp = day4y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(12, [1 length(day5y(keep))]); ytemp = day5y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(13, [1 length(day6y(keep))]); ytemp = day6y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;

xtemp = repmat(15, [1 length(day1o(keep))]); ytemp = day1o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(16, [1 length(day2o(keep))]); ytemp = day2o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(17, [1 length(day3o(keep))]); ytemp = day3o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(18, [1 length(day4o(keep))]); ytemp = day4o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(19, [1 length(day5o(keep))]); ytemp = day5o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(20, [1 length(day6o(keep))]); ytemp = day6o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;

a = [1 2 3 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20]; 
b = [mean(day1(keep)) mean(day2(keep)) mean(day3(keep)) mean(day4(keep)) mean(day5(keep)) mean(day6(keep)) ...
    mean(day1y(keep)) mean(day2y(keep)) mean(day3y(keep)) mean(day4y(keep)) mean(day5y(keep)) mean(day6y(keep)) ...
    mean(day1o(keep)) mean(day2o(keep)) mean(day3o(keep)) mean(day4o(keep)) mean(day5o(keep)) mean(day6o(keep))];
c = 1.96*[std(day1(keep))/sqrt(length(day1(keep))) std(day2(keep))/sqrt(length(day2(keep))) std(day3(keep))/sqrt(length(day3(keep))) std(day4(keep))/sqrt(length(day4(keep))) std(day5(keep))/sqrt(length(day5(keep))) std(day6(keep))/sqrt(length(day6(keep))) ...
std(day1y(keep))/sqrt(length(day1y(keep))) std(day2y(keep))/sqrt(length(day2y(keep))) std(day3y(keep))/sqrt(length(day3y(keep))) std(day4y(keep))/sqrt(length(day4y(keep))) std(day5y(keep))/sqrt(length(day5y(keep))) std(day6y(keep))/sqrt(length(day6y(keep))) ...
std(day1(keep))/sqrt(length(day1o(keep))) std(day2o(keep))/sqrt(length(day2o(keep))) std(day3o(keep))/sqrt(length(day3o(keep))) std(day4o(keep))/sqrt(length(day4o(keep))) std(day5o(keep))/sqrt(length(day5o(keep))) std(day6o(keep))/sqrt(length(day6o(keep)))];
ba = bar(a, b); 
for i = 1:length(a); plot([a(i) a(i)], [b(i)-c(i) b(i)+c(i)], 'k-', 'LineWidth', 1.5); 
end
ba.FaceColor = 'w';
ba.EdgeColor = 'k';
ba.BarWidth = 0.75;
ba.LineWidth = 1.5;
ba.FaceAlpha = 0.4;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = 0; xlim_hi = 21;
ylim_lo = 0.25; ylim_hi = 1.00;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = a;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = {'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = '    All                    Younger                 Older';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0.25 0.50 0.75 1.00];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'0.25', '0.50', '0.75', '1.00'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = 'Legibility Score';

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
title('Letters');

subplot(1, 2, 2); hold on; % numbers
keep = 27:36;

xtemp = repmat(1, [1 length(day1(keep))]); ytemp = day1(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(2, [1 length(day2(keep))]); ytemp = day2(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(3, [1 length(day3(keep))]); ytemp = day3(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(4, [1 length(day4(keep))]); ytemp = day4(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(5, [1 length(day5(keep))]); ytemp = day5(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(6, [1 length(day6(keep))]); ytemp = day6(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;

xtemp = repmat(8, [1 length(day1y(keep))]); ytemp = day1y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(9, [1 length(day2y(keep))]); ytemp = day2y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(10, [1 length(day3y(keep))]); ytemp = day3y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(11, [1 length(day4y(keep))]); ytemp = day4y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(12, [1 length(day5y(keep))]); ytemp = day5y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(13, [1 length(day6y(keep))]); ytemp = day6y(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;

xtemp = repmat(15, [1 length(day1o(keep))]); ytemp = day1o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(16, [1 length(day2o(keep))]); ytemp = day2o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(17, [1 length(day3o(keep))]); ytemp = day3o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(18, [1 length(day4o(keep))]); ytemp = day4o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(19, [1 length(day5o(keep))]); ytemp = day5o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;
xtemp = repmat(20, [1 length(day6o(keep))]); ytemp = day6o(keep); s=scatter(xtemp, ytemp, 'o', 'SizeData', 20, 'CData', [0 0 0]/255, 'MarkerEdgeAlpha', 0.50);
clear xtemp ytemp s;

a = [1 2 3 4 5 6 8 9 10 11 12 13 15 16 17 18 19 20]; 
b = [mean(day1(keep)) mean(day2(keep)) mean(day3(keep)) mean(day4(keep)) mean(day5(keep)) mean(day6(keep)) ...
    mean(day1y(keep)) mean(day2y(keep)) mean(day3y(keep)) mean(day4y(keep)) mean(day5y(keep)) mean(day6y(keep)) ...
    mean(day1o(keep)) mean(day2o(keep)) mean(day3o(keep)) mean(day4o(keep)) mean(day5o(keep)) mean(day6o(keep))];
c = 1.96*[std(day1(keep))/sqrt(length(day1(keep))) std(day2(keep))/sqrt(length(day2(keep))) std(day3(keep))/sqrt(length(day3(keep))) std(day4(keep))/sqrt(length(day4(keep))) std(day5(keep))/sqrt(length(day5(keep))) std(day6(keep))/sqrt(length(day6(keep))) ...
std(day1y(keep))/sqrt(length(day1y(keep))) std(day2y(keep))/sqrt(length(day2y(keep))) std(day3y(keep))/sqrt(length(day3y(keep))) std(day4y(keep))/sqrt(length(day4y(keep))) std(day5y(keep))/sqrt(length(day5y(keep))) std(day6y(keep))/sqrt(length(day6y(keep))) ...
std(day1(keep))/sqrt(length(day1o(keep))) std(day2o(keep))/sqrt(length(day2o(keep))) std(day3o(keep))/sqrt(length(day3o(keep))) std(day4o(keep))/sqrt(length(day4o(keep))) std(day5o(keep))/sqrt(length(day5o(keep))) std(day6o(keep))/sqrt(length(day6o(keep)))];
ba = bar(a, b); 
for i = 1:length(a); plot([a(i) a(i)], [b(i)-c(i) b(i)+c(i)], 'k-', 'LineWidth', 1.5); 
end
ba.FaceColor = 'w';
ba.EdgeColor = 'k';
ba.BarWidth = 0.75;
ba.LineWidth = 1.5;
ba.FaceAlpha = 0.4;

fontsize = 16;
fontname = 'Arial';
xticklength = 0;

xlim_lo = 0; xlim_hi = 21;
ylim_lo = 0.40; ylim_hi = 1.00;

xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = a;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = {'Pre', 'Post', 'Pre', 'Post', 'Pre', 'Post'};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.Label.String = '    All                    Younger                 Older';

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0.40 0.60 0.80 1.00];
yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'0.40', '0.60', '0.80', '1.00'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.Label.FontAngle= 'italic';
% yax.TickLabelRotation = 0;
yax.Label.String = 'Legibility Score';

g = gca;
g.TickDir = 'out';
box off; pbaspect([1 1 1]);
title('Digits');

print(fullfile(rootdir, 'plots', ['plot_bar_6days_legibility_n=' num2str(nall)]), '-dpng')

% Export table.
out = array2table([day1' day1std' day2' day2std' day3' day3std' day4' day4std' day5' day5std' day6' day6std' ...
    day1y' day1ystd' day2y' day2ystd' day3y' day3ystd' day4y' day4ystd' day5y' day5ystd' day6y' day6ystd' ...
    day1o' day1ostd' day2o' day2ostd' day3o' day3ostd' day4o' day4ostd' day5o' day5ostd' day6o' day6ostd']);
out = [array2table(alphabet', 'VariableNames', {'alphabet'}) out];
out.Properties.VariableNames = {'letter', 'day1', 'day1std', 'day2', 'day2std', 'day3', 'day3std', 'day4', 'day4std', 'day5', 'day5std', 'day6', 'day6std', ...
    'dayy1', 'day1ystd', 'day2y', 'day2ystd', 'day3y', 'day3ystd', 'day4y', 'day4ystd', 'day5y', 'day5ystd', 'day6y', 'day6ystd', ...
    'day1o', 'day1ostd', 'day2o', 'day2ostd', 'day3o', 'day3ostd', 'day4o', 'day4ostd', 'day5o', 'day5ostd', 'day6o', 'day6ostd'};
writetable(out, fullfile(rootdir, 'supportfiles', 'table2_6days_legibility.csv'))

