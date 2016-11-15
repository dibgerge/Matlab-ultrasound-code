%% Input arguments: 
%   c: A cell array of the velocities. One cell per mode. 
%   fstartInd: The frequency index for each mode. 
%   f: The frequency vector. 
%   props: a structure for drawing the lines
function drawVelocityDispersion(c, fInd, f, props)

% check the number of arguments
if nargin < 3 || nargin > 4 
    error('wrong number of arguments.');
end

% set default values of props
if nargin == 3
    props = struct;    
end

% check validity of c 
if ~iscell(c) ||  ~isnumeric([c{:}])
    error('c should be cells of numeric arrays.');
end

% check validity of startInd
if ~isnumeric([fInd{:}]) || length(fInd) ~= length(c)
    error('fInd should be a numeric array with same length as c.');
end

% Check if  f is valid and have correct values
if ~isnumeric(f) || ~isempty(find(f<0,1))
    error('Invalid values in the frequencies f.');
end

if ~isempty(find([fInd{:}] > length(f),1))
    error('fInd values could not be greater than the length of f.');
end

% check validity of props 
if ~isstruct(props) 
    error('props should be a structure.');
end

% set default value sof props 
if ~isfield(props, 'spec')
    props.spec = 'r';
end

if ~isfield(props, 'linewidth') 
    props.linewidth = 3;
elseif ~isnumeric(props.linewidth) || ~isscalar(props.linewidth) || props.linewidth < 1
    error('props.linewidth has an invalid value.');
end

if ~isfield(props, 'xlabel')
    props.xlabel = 'Frequency (kHz)';
end

if ~isfield(props, 'prefix')
    props.prefix = 'S';
end
    

oldhold = ishold;

%set(gcf, 'Position', [2200, 100, 1300, 800]);

set(gca, 'MinorGridLineStyle', 'none', 'XGrid', 'off', ...
    'YGrid', 'off', 'ZGrid', 'off', 'Box', 'off', 'FontSize', 40, ...
    'FontWeight', 'normal')


for i=1:length(c)
    hold on;
    % do not plot curves with less than 3 points
    if length(c{i}) < 3
        continue;
    end
    pl = plot(f(fInd{i})/1e6, 1e-3*c{i}, props.spec, 'linewidth', props.linewidth);
    if ~isempty(props.prefix)
        text(f(fInd{i}(3))/1e6, 1e-3*c{i}(3)+0.30, [props.prefix '_{' num2str(i-1) '}'], ...
            'FontWeight', 'normal', 'FontSize', 40, 'Color', get(pl, 'Color'));
    end
end


% xlabel(props.xlabel)
xlabel('Frequency (MHz)')
ylabel('Velocity (km/s)')
% xlim([f(1)/1e3 f(end)/1e3])

if ~oldhold
    hold off;
end

end

