%%% Probe the Cauchy cwt for its properties 
%%% First reveal the mother function (cwt of a point)

%%%% How many Scales is enough, for the Cauchy? 

% Let's analyze pure periodic striped patterns
% Bigger array needed here
arraySize2 = 512;
% Create the x and y coordinate arrays
x = 1:arraySize2;
y = 1:arraySize2;
[X, Y] = meshgrid(x, y);

% Logarithmic 
Scales = 10.^(1:.025:3) /4. ;
% Equal spaced
    % Scales = 2:2:50;

wavenums = 2.^(1:1:5);

% Initialize the 2D array with zeros
y = zeros(size(Scales,2), size(wavenums,2));

% Fill y array with results 
for iwave = 1:size( wavenums,2 )  % 2.^(1:1:5) %4:8:80

    halfwave = arraySize2./wavenums(iwave)./2;

% Create the sinusoidal pattern
    stripes = sin(2 * pi * wavenums(iwave) * X / arraySize2);

    Angles = 0;

    cwtCauchy = cwtft2(stripes,wavelet="cauchy",scales=Scales, angles=Angles);
    spec = squeeze( cwtCauchy.cfs );

    % Display the image
    figure(wavenums(iwave));
    subplot(211)
    contour(stripes); colorbar;
    subplot(212)
    bar( Scales, squeeze(abs(spec(256,256,:)))); 
    title(string(wavenums(iwave))+': half wavelength is '+string(halfwave)); 
    xlabel('Scales');

    y(:,iwave) = squeeze(abs(spec(256,256,:)));
    y(:,iwave) = y(:,iwave) *2.0 ./ transpose(Scales);

end % wavenumber loop

% Quick plot 
% bar(Scales,y)





%%%% FANCY PLOT 

data = y;  % 81 rows, 5 columns

% Define colors for each case (customize as needed)
colors = [
    0.0000    0.4470    0.7410  % Blue
    0.8500    0.3250    0.0980  % Red
    0.9290    0.6940    0.1250  % Yellow
    0.4940    0.1840    0.5560  % Purple
    0.4660    0.6740    0.1880   % Green
    ];

% Calculate bar width and offsets
numCases = size(data, 2);
barWidth = 2; % Adjust as needed
groupWidth = barWidth * numCases; % Total width of a group of bars
xCenters = Scales; % - groupWidth/2; % + barWidth/2; % Centers of the bar groups

% Create the bar plot
figure(1);
hold on; % Keep all bars on the same plot

for i = 1:numCases
    % Plot the bars for the current case
    bar(xCenters, data(:, i), 20, 'FaceColor', colors(i,:), 'EdgeColor','k'); % Use specified color
end

hold off;

% Customize the plot
xlabel('Scale S');  % Or your x-axis label
ylabel('abs(spec) *2/S');             % Or your y-axis label
title('Cauchy Scale spectra of unit monochromatic waves');

% Set x-axis ticks and labels (optional)
% xticks(1:size(data,1)); % Show ticks for each data point
xlim([0 250]); 

% Add a legend
legend({'128', '64', '32', '16', '8'}, 'Location', 'best'); % Customize labels

% Improve appearance (optional)
% grid on;
%xlim([0 size(data,1)+1]); % Set axis limits to avoid cutting off bars
%ylim([min(min(data)) max(max(data))]); % Set axis limits to avoid cutting off bars

title('Cauchy Scales spectrum for monochromatic waves, by half-wavelength')
