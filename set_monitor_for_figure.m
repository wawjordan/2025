function set_monitor_for_figure(fh,mon)

% Get pixel position of monitors.
% fh = figure('units', 'pixels');
fh.set('units','pixels');
MP = get(0, 'MonitorPositions');
N = size(MP, 1);

% Might want to set an initial position this to some reasonable location
% in the event of the window being "Restored Down".
newPosition = MP(1,:);

% Set position of fh to be within a secondary monitor, if possible.
if size(MP, 1) == 1
    % Single monitor -- do nothing.
    
else
    for i = 1:mon
        % Multiple monitors - shift to the Nth monitor.
        newPosition(1) = newPosition(1) + MP(i,1);
    end
end
fh.set('Position', newPosition, 'units', 'normalized');
% fh.WindowState = 'maximized'; % Maximize with respect to current monitor.

end