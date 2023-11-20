function [outliers] = findoutliers_lm(d, modelspec)

% Returns the INDICDES of outliers in a table (d) based on a specific model (modelspec). modelspec is a string. 

% modelspec = 'assoc ~ hip*age + sex';

% Fit robust regression model.
mdlr = fitlm(d, modelspec, 'RobustOpts', 'on');

outliers = [];
% Examine model residuals: boxplot of raw residuals.
figure(1); f = figure('visible', 'off');
m = mdlr.Residuals.Raw;
e = eps(max(m(:)));
boxplot(m)
% ylabel('Raw Residuals')
% Suppress figure display.
set(gcf,'Visible','off');              
set(0,'DefaultFigureVisible','off');

% Get indices of the outliers.
h = flipud(findobj(gcf,'tag','Outliers')); % flip order of handles
for jj = 1 : length( h )
    x =  get( h(jj), 'XData' );
    y =  get( h(jj), 'YData' );
    for ii = 1 : length( x )
        if not( isnan( x(ii) ) )
            ix = find( abs( m(:,jj)-y(ii) ) < e );
            outliers = cat(1, outliers, ix);
            %                 text( x(ii), y(ii), sprintf( '\\leftarrowY%02d', ix ) )
        end
    end
end

% Examine robust weights: boxplot of robust weights.
figure(2); f = figure('visible', 'off');
m = mdlr.Robust.Weights;
e = eps(max(m(:)));
boxplot(m);
% ylabel('Robust Beta-Weights')
% Suppress figure display.
set(gcf,'Visible','off');              
set(0,'DefaultFigureVisible','off');

% Get indices of the outliers.
h = flipud(findobj(gcf,'tag','Outliers')); % flip order of handles
for jj = 1 : length( h )
    x =  get( h(jj), 'XData' );
    y =  get( h(jj), 'YData' );
    for ii = 1 : length( x )
        if not( isnan( x(ii) ) )
            ix = find( abs( m(:,jj)-y(ii) ) < e );
            outliers = cat(1, outliers, ix);
            %                 text( x(ii), y(ii), sprintf( '\\leftarrowY%02d', ix ) )
        end
    end
end

outliers = sort(unique(outliers));

f = figure('visible', 'on');

close all;
