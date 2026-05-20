function createPolarPlot_diff(data, se, groups, types, config, titleStr, tractType, ylims, includeControls)
    % CREATEPOLARPLOT_DIFF Creates a polar plot using Cartesian coordinates with different y-limits for each axis
    % and displays mean values for each group directly on the plot.
    %
    % Usage:
    %   createPolarPlot_diff(data, se, groups, types, config, titleStr, tractType, ylims, includeControls)
    %
    % Inputs:
    %   data            - Matrix of values to plot (rows: groups, columns: types)
    %   se              - Matrix of standard errors (same size as data)
    %   groups          - Cell array of group names (clinical measures)
    %   types           - Cell array of type names (e.g., {'RRMS', 'SPMS', 'PPMS'})
    %   config          - Configuration structure containing msTypeColors
    %   titleStr        - Title string for the plot
    %   tractType       - String indicating tract type ('dwi' or 'template')
    %   ylims           - Matrix of y-limits [min max] for each group (rows: groups, columns: [min max])
    %   includeControls - Boolean indicating whether to include controls in legend
    %
    % Example:
    %   createPolarPlot_diff(data, se, groups, {'RRMS', 'SPMS', 'PPMS'}, config, ...
    %                   'Left Hemisphere', 'dwi', ylims, false)
    
    % Input validation
    if nargin < 9
        includeControls = false;
    end
    
    % Define fixed angles for each measure (in degrees)
    % 0°: x9HPTND, 60°: DMT, 120°: MSPro, 180°: x9HPTD, 240°: EDSS, 300°: T25FW
    theta_deg = [0, 60, 120, 180, 240, 300];
    theta = deg2rad(theta_deg);
    
    % Create a regular axes
    ax = gca;
    hold on;
    
    % Set aspect ratio to equal and remove axes
    axis equal;
    axis off;
    
    % Draw background grid circles
    numCircles = 5; % Number of concentric circles to draw
    radii = linspace(0, 1, numCircles+1);
    radii = radii(2:end); % Remove the center point
    
    % Create full circle for theta
    thetaCircle = linspace(0, 2*pi, 100);
    
    % Draw concentric circles (grid lines)
    for r = radii
        x = r * cos(thetaCircle);
        y = r * sin(thetaCircle);
        plot(x, y, 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.0);
    end
    
    % Draw radial lines from center to each group
    for i = 1:length(theta)
        plot([0, cos(theta(i))], [0, sin(theta(i))], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.0);
    end
    
    % Plot each type with its designated color
    for i = 1:length(types)
        type = types{i};
        color = config.msTypeColors.(type) / 255; % Convert RGB from 0-255 to 0-1
        
        % Extract data for this type
        dataForType = data(:,i);
        
        % Create coordinates for mean line
        x = zeros(length(theta), 1);
        y = zeros(length(theta), 1);
        
        % Create coordinates for error bounds if SE provided
        if ~isempty(se)
            seForType = se(:,i);
            x_upper = zeros(length(theta), 1);
            y_upper = zeros(length(theta), 1);
            x_lower = zeros(length(theta), 1);
            y_lower = zeros(length(theta), 1);
        end
        
        for j = 1:length(theta)
            % Normalize the value to [0,1] based on its ylims
            r = (dataForType(j) - ylims(j,1)) / (ylims(j,2) - ylims(j,1));
            x(j) = r * cos(theta(j));
            y(j) = r * sin(theta(j));
            
            % Error bounds if SE provided
            if ~isempty(se)
                r_upper = min(1, (dataForType(j) + seForType(j) - ylims(j,1)) / (ylims(j,2) - ylims(j,1)));
                r_lower = max(0, (dataForType(j) - seForType(j) - ylims(j,1)) / (ylims(j,2) - ylims(j,1)));
                
                x_upper(j) = r_upper * cos(theta(j));
                y_upper(j) = r_upper * sin(theta(j));
                x_lower(j) = r_lower * cos(theta(j));
                y_lower(j) = r_lower * sin(theta(j));
            end
        end
        
        % Add closing points
        x = [x; x(1)];
        y = [y; y(1)];
        
        if ~isempty(se)
            x_upper = [x_upper; x_upper(1)];
            y_upper = [y_upper; y_upper(1)];
            x_lower = [x_lower; x_lower(1)];
            y_lower = [y_lower; y_lower(1)];
            
            % Create shaded error area
            x_patch = [x_upper; flipud(x_lower)];
            y_patch = [y_upper; flipud(y_lower)];
            patch(x_patch, y_patch, color', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        
        % Plot mean line and markers
        plot(x, y, 'Color', color, 'LineWidth', 2, 'DisplayName', type, ...
             'Marker', 'o', 'MarkerFaceColor', color, 'MarkerSize', 8, ...
             'MarkerIndices', 1:length(theta));
    end
    
    % Add group labels and mean values around the circle
    for i = 1:length(groups)
        % Position labels slightly outside the circle
        labelR = 1.15;
        x = labelR * cos(theta(i));
        y = labelR * sin(theta(i));
        
        % Adjust text alignment based on position in the circle
        if abs(x) < 0.1 % near vertical axis
            halign = 'center';
        elseif x > 0
            halign = 'left';
        else
            halign = 'right';
        end
        
        if abs(y) < 0.1 % near horizontal axis
            valign = 'middle';
        elseif y > 0
            valign = 'bottom';
        else
            valign = 'top';
        end
        
        % Create label with measure name and mean values for each group
        label = groups{i};
        for j = 1:length(types)
            mean_val = data(i,j);
            label = sprintf('%s\n%s: %.1f', label, types{j}, mean_val);
        end
        text(x, y, label, 'HorizontalAlignment', halign, 'VerticalAlignment', valign, 'FontWeight', 'bold');
    end
    
    % Add legend with only one entry per type
    % Get the current legend objects
    h = findobj(gca, 'Type', 'line');
    % Only use the last entries for the legend
    if includeControls
        legend(h(1:4), 'Location', 'northeastoutside');
    else
        legend(h(1:3), 'Location', 'northeastoutside');
    end
    
    % Add title
    title(titleStr);
    
    % Set axis limits
    axisLimit = 1.2;
    xlim([-axisLimit, axisLimit]);
    ylim([-axisLimit, axisLimit]);
end 