function createPolarPlot(data, se, groups, types, config, titleStr, tractType, maxR, includeControls, varargin)
    % CREATEPOLARPLOT Creates a polar plot using Cartesian coordinates
    %
    % Usage:
    %   createPolarPlot(data, se, groups, types, config, titleStr, tractType, maxR, includeControls)
    %
    % Inputs:
    %   data            - Matrix of values to plot (rows: groups, columns: types)
    %   se              - Matrix of standard errors (same size as data)
    %                     Use [] if no standard errors should be plotted
    %   groups          - Cell array of group names
    %   types           - Cell array of type names (e.g., {'RRMS', 'SPMS', 'PPMS'})
    %   config          - Configuration structure containing msTypeColors
    %   titleStr        - Title string for the plot
    %   tractType       - String indicating tract type ('dwi' or 'template')
    %   maxR            - Maximum radius for plot (single value) or [min, max] radius range
    %   includeControls - Boolean indicating whether to include controls in legend
    %
    % Example:
    %   createPolarPlot(data, se, groups, {'RRMS', 'SPMS', 'PPMS'}, config, ...
    %                   'Left Hemisphere', 'dwi', [], false)
    %   createPolarPlot(data, se, groups, {'RRMS', 'SPMS', 'PPMS'}, config, ...
    %                   'Left Hemisphere', 'dwi', [0.5 2.0], false)
    
    % Input validation

    if nargin < 9
        includeControls = false;
    end
    if nargin < 8
        maxR = [];
    end

    % Parse optional log scale flag
    if ~isempty(varargin)
        logScale = logical(varargin{1});
    else
        logScale = false;
    end

    % Transformation helper
    logTransform = @(x) log10(x + 1); % avoids log(0)

    if logScale
        dataPlot = logTransform(data);
        if ~isempty(se)
            sePlot = logTransform(data + se) - logTransform(data);
        else
            sePlot = [];
        end
    else
        dataPlot = data;
        sePlot = se;
    end
    
    % Convert data to polar coordinates with 45-degree rotation
    theta = linspace(0, 2*pi, length(groups)+1);
    theta = theta(1:end-1) + pi/4; % Add 45-degree rotation
    
    % Create a regular axes
    ax = gca;
    hold on;
    
    % Set aspect ratio to equal and remove axes
    axis equal;
    axis off;
    
    % Handle maxR input - can be single value or [min, max] range
    if isempty(maxR)
        if ~isempty(se)
            maxROrig = max(max(data(:) + se(:))) * 1.1; % 10% margin
        else
            maxROrig = max(data(:)) * 1.1;
        end
        if maxROrig == 0
            maxROrig = 1; % Prevent zero radius
        end
        minROrig = 0; % Default minimum radius
    elseif length(maxR) == 1
        % Single value - use as maximum radius, minimum is 0
        maxROrig = maxR;
        minROrig = 0;
    elseif length(maxR) == 2
        % Two values - use as [min, max] range
        minROrig = min(maxR);
        maxROrig = max(maxR);
    else
        error('maxR must be either a single value or a two-element array [min, max]');
    end

    % For plotting, we always want to use 0 as center and scale the data
    % Data normalization: map [minROrig, maxROrig] to [0, maxRPlot]
    dataRange = maxROrig - minROrig;
    if dataRange == 0
        dataRange = 1; % Prevent division by zero
    end
    
    if logScale
        % Apply log transformation first, then normalize
        dataLogTransformed = logTransform(dataPlot);
        minROrigLog = logTransform(minROrig);
        maxROrigLog = logTransform(maxROrig);
        dataRangeLog = maxROrigLog - minROrigLog;
        
        dataNormalized = (dataLogTransformed - minROrigLog) / dataRangeLog;
        if ~isempty(sePlot)
            seLogTransformed = logTransform(dataPlot + sePlot) - dataLogTransformed;
            seNormalized = seLogTransformed / dataRangeLog;
        else
            seNormalized = [];
        end
        maxRPlot = 1; % Normalized scale (0 to 1)
    else
        % Normalize data to [0, maxRPlot] range
        dataNormalized = (dataPlot - minROrig) / dataRange;
        if ~isempty(sePlot)
            % For SE, we need to scale it proportionally
            seNormalized = sePlot / dataRange;
        else
            seNormalized = [];
        end
        maxRPlot = 1; % Normalized scale (0 to 1)
    end
    
    minRPlot = 0; % Always use 0 as center for plotting
    
    % Draw background grid circles
    numCircles = 5; % Number of concentric circles to draw
    radiiOrig = linspace(minROrig, maxROrig, numCircles+1);
    radiiOrig = radiiOrig(2:end); % Remove the innermost point to avoid duplicate with minROrig

    % Convert original radii to normalized plotting positions
    if logScale
        radiiOrigLog = logTransform(radiiOrig);
        radiiPlot = (radiiOrigLog - minROrigLog) / dataRangeLog * maxRPlot;
    else
        radiiPlot = (radiiOrig - minROrig) / dataRange * maxRPlot;
    end
    
    % Create full circle for theta
    thetaCircle = linspace(0, 2*pi, 100);
    
    % Draw concentric circles (grid lines) - make more obvious
    for idxR = 1:numel(radiiPlot)
        rPlot = radiiPlot(idxR);
        rLabel = radiiOrig(idxR);
        x = rPlot * cos(thetaCircle);
        y = rPlot * sin(thetaCircle);
        plot(x, y, 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.0);
        
        % Optionally add radial labels using original scale values
        text(rPlot * 1.05, 0, num2str(rLabel, '%.1f'), 'FontSize', 8, 'HorizontalAlignment', 'left');
    end
    
    % Draw radial lines from center to maximum radius - make more obvious
    for i = 1:length(theta)
        plot([0, maxRPlot * cos(theta(i))], [0, maxRPlot * sin(theta(i))], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.0);
    end
    
    % Plot each type with its designated color
    for i = 1:length(types)
        type = types{i};
        color = config.msTypeColors.(type) / 255; % Convert RGB from 0-255 to 0-1
        
        % Extract data for this type
        dataForType = dataNormalized(:,i);
        
        % Create coordinates for mean line
        x = zeros(length(theta), 1);
        y = zeros(length(theta), 1);
        
        % Create coordinates for error bounds if SE provided
        if ~isempty(seNormalized)
            seForType = seNormalized(:,i);
            x_upper = zeros(length(theta), 1);
            y_upper = zeros(length(theta), 1);
            x_lower = zeros(length(theta), 1);
            y_lower = zeros(length(theta), 1);
        end
        
        for j = 1:length(theta)
            % Mean coordinates
            r = dataForType(j);
            x(j) = r * cos(theta(j));
            y(j) = r * sin(theta(j));
            
            % Error bounds if SE provided
            if ~isempty(seNormalized)
                r_upper = r + seForType(j);
                r_lower = max(0, r - seForType(j)); % Ensure non-negative (center is 0)
                
                x_upper(j) = r_upper * cos(theta(j));
                y_upper(j) = r_upper * sin(theta(j));
                x_lower(j) = r_lower * cos(theta(j));
                y_lower(j) = r_lower * sin(theta(j));
            end
        end
        
        % Add closing points
        x = [x; x(1)];
        y = [y; y(1)];
        
        if ~isempty(seNormalized)
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
        plot(x, y, 'Color', color, 'LineWidth', 2, 'DisplayName', type);
        %         % Plot mean line and markers
        % plot(x, y, 'Color', color, 'LineWidth', 2, 'DisplayName', type, ...
        %      'Marker', 'o', 'MarkerFaceColor', color, 'MarkerSize', 8, ...
        %      'MarkerIndices', 1:length(theta));
    end
    
    % Add group labels around the circle
    for i = 1:length(groups)
        % Position labels slightly outside the maximum radius
        labelR = maxRPlot * 1.15;
        x = labelR * cos(theta(i));
        y = labelR * sin(theta(i));
        
        % Adjust text alignment based on position in the circle
        if abs(x) < 0.1 * maxRPlot % near vertical axis
            halign = 'center';
        elseif x > 0
            halign = 'left';
        else
            halign = 'right';
        end
        
        if abs(y) < 0.1 * maxRPlot % near horizontal axis
            valign = 'middle';
        elseif y > 0
            valign = 'bottom';
        else
            valign = 'top';
        end
        
        text(x, y, groups{i}, 'HorizontalAlignment', halign, 'VerticalAlignment', valign, 'FontWeight', 'bold');
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
    axisLimit = maxRPlot * 1.2;
    xlim([-axisLimit, axisLimit]);
    ylim([-axisLimit, axisLimit]);
end 