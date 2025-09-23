close all;
clear;

% Parameters
n = 100;  % Sample size
sens_range_se = linspace(0.6, 1.0, 101);  % Sensitivity range for SE plots
spec_range_se = linspace(0.6, 1.0, 101);  % Specificity range for SE plots

% Create mesh grid
[S, F] = meshgrid(sens_range_se, 1 - spec_range_se);  % F = 1 - specificity

% Function to compute Standard Error surface
% compute_se = @(y_true, s, f) sqrt((y_true * (1 - y_true) / n) ./ (s - f).^2);

% Function to compute Standard Error surface in terms of p
compute_se = @(p, s, f) sqrt((f + (s - f) .* p) .* (1 - f - (s - f) .* p) ./ (n * (s - f).^2));


% Compute SE for y_true = 0.5 and y_true = 0.9
SE_50 = 100*compute_se(0.5, S, F); 
SE_90 = 100*compute_se(0.9, S, F); 

vmin = min([SE_50(:); SE_90(:)]);  % Global minimum across both datasets
vmax = max([SE_90(:); SE_50(:)]);  % Global maximum across both datasets

% Plot settings
colormap_choice = colormap(brewermap([],"YlGnBu"));
fz = 24; % Font size
annot_line_color = 0.2*[1, 1, 1]; % Grey color for annotations
circle_size = 100; % Marker size for circles

% Plot combined figures using tiledlayout
%figure('Position', [100, 100, 1300, 500]);
%t = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');

figure('Position', [100, 100, 1200, 500]);
t = tiledlayout(1, 2, 'TileSpacing', 'none', 'Padding', 'none');

annot_positions1 = [compute_sens_spec(0.07,0.5,100), compute_sens_spec(0.1,0.5,100)];
annot_positions2 = [compute_sens_spec(0.07,0.9,100), compute_sens_spec(0.1,0.9,100)];

% --- Figure 2A --
nexttile;
imagesc([0.6, 1], [0.6, 1], SE_50);
clim([vmin, vmax]);  % Set fixed color limits
set(gca, 'YDir', 'normal', 'FontSize', fz);
colormap(colormap_choice);
% colorbar;
hold on;
% grid on;

% Add contour lines
contour(sens_range_se, spec_range_se, SE_50, 100*[0.06, 0.07, 0.08, 0.1, 0.15, 0.2], ...
    '--', 'Color', 'black', 'LineWidth', 0.5, 'ShowText', 'on');

% Add diagonal line
plot([0.6, 1], [0.6, 1], ':', 'Color', [179, 27, 27]/256, 'LineWidth', 1);

% Add grey annotation lines and circles
for pos = annot_positions1
    plot([min(sens_range_se),pos],[pos,pos], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
    plot([pos,pos],[pos,min(spec_range_se)], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
    plot(pos, pos, 'o', 'MarkerEdgeColor', annot_line_color, 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');
end

se_req = compute_sensitivity(0.070, 0.50, 100, 0.9995);
plot([min(sens_range_se),se_req],[0.9995,0.9995], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
plot([se_req,se_req],[0.9995,min(spec_range_se)], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
plot(se_req, 0.9995, 'o', 'MarkerEdgeColor', annot_line_color, 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');


hold off;
axis square;
xlabel('Sensitivity (%)', 'FontSize', fz+10);
ylabel('Specificity (%)', 'FontSize', fz+10);
title('Standard Error ($\hat{p} = 50\%$)', 'Interpreter', 'latex', 'FontSize', fz, 'FontWeight', 'normal');
text(-0.22, 1.02, 'A', 'Units', 'normalized', 'FontSize', 36, 'FontWeight', 'normal');
xticks(0.6:0.1:1);
xticklabels(string(xticks() * 100));
yticks(0.6:0.1:1);
yticklabels(string(yticks() * 100));
xlim([0.595, 1.005]);  % Add margin on the left
ylim([0.595, 1.005]);


% --- Figure 2B ---
nexttile;
imagesc([0.6, 1], [0.6, 1], SE_90);
clim([vmin, vmax]);  % Set fixed color limits
set(gca, 'YDir', 'normal', 'FontSize', fz);
colormap(colormap_choice);
% colorbar;
hold on;
% grid on;

% Add contour lines
contour(sens_range_se, spec_range_se, SE_90, 100*[0.04, 0.06, 0.07, 0.08, 0.1, 0.15,0.20], ...
    '--', 'Color', 'black', 'LineWidth', 0.5, 'ShowText', 'on');

% Add diagonal line
plot([0.6, 1], [0.6, 1], ':', 'Color', [179, 27, 27]/256, 'LineWidth', 1);

% Add grey annotation lines and circles
for pos = annot_positions2
    plot([min(sens_range_se),pos],[pos,pos], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
    plot([pos,pos],[pos,min(spec_range_se)], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
    plot(pos, pos, 'o', 'MarkerEdgeColor', annot_line_color, 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');
    % text(pos+0.01,pos+0.01,['(',num2str(round(10000*pos)/100),',',num2str(round(10000*pos)/100),')'])
end

se_req = compute_sensitivity(0.070, 0.90, 100, 0.9995);
plot([min(sens_range_se),se_req],[0.9995,0.9995], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
plot([se_req,se_req],[0.9995,min(spec_range_se)], ':', 'Color', annot_line_color, 'LineWidth', 1.5);
plot(se_req, 0.9995, 'o', 'MarkerEdgeColor', annot_line_color, 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'none');

hold off;
axis square;
xlabel('Sensitivity (%)', 'FontSize', fz+10);
ylabel('Specificity (%)', 'FontSize', fz+10);
title('Standard Error ($\hat{p} = 90\%$)', 'Interpreter', 'latex', 'FontSize', fz, 'FontWeight', 'normal');
text(-0.22, 1.02, 'B', 'Units', 'normalized', 'FontSize', 36, 'FontWeight', 'normal');

xticks(0.6:0.1:1);
xticklabels(string(xticks() * 100));
yticks(0.6:0.1:1);
yticklabels(string(yticks() * 100));
xlim([0.595, 1.005]);  % Add margin on the left
ylim([0.595, 1.005]);


cb = colorbar('Location', 'eastoutside');  % Attach to the current figure

% Get current colorbar position
cbPos = cb.Position;

% Dynamically adjust width (e.g., increase by 50%)
scaleFactor = 1;  % Increase width by 50%
cbPos(3) = cbPos(3) * scaleFactor;  % Adjust width

% Optionally, slightly move it to the right to maintain spacing from plots
cbPos(1) = cbPos(1) + 0.03;  % Shift colorbar horizontally

% Apply the new position
cb.Position = cbPos;

% Customize the colorbar appearance
cb.Label.String = 'Standard Error (%)';
cb.Label.Rotation = 90;
cb.Label.FontSize = fz;
cb.Label.FontName = 'Times';


% Set font to Times for all elements
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times');

% Save combined figure (uncomment if needed)
exportgraphics(gcf, 'Figure2_Combined.png', 'Resolution', 300);
exportgraphics(gcf, 'Figure2_Combined.pdf','BackgroundColor','white','ContentType','vector');
