%%% Probe the Cauchy cwt for its properties 
%%% First reveal the mother function (cwt of a point)

% Define the size of the array
arraySize = 300;
image = zeros(arraySize, arraySize);

% Set a single pixel to 1 (center of the image for best visualization)
image(arraySize/2, arraySize/2) = 1;

% Pick a NW-SE angle and Scale for visual pleasure
Angles = pi/4 ;
Scales = 20 ;

cwtCauchy = cwtft2(image,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwtCauchy.cfs );

% Create the isosurface plot of real, imaginary
figure(1);
subplot(131)
mesh( real(spec) ); title('a) Cauchy real part'); axis off;
lighting phong
camlight

subplot(132)
mesh( imag(spec) ); title('b) Cauchy imaginary part'); axis off;
lighting phong
camlight


%%%% How many angles is enough, for the Cauchy case? 
% Let's try it for a single line pattern 
subplot(133); 
image(arraySize/2, :) = 1;
Scales = 10;
Angles = 0:pi/12:pi ;

cwtCauchy = cwtft2(image,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwtCauchy.cfs );

plot( squeeze(real(spec(150,150,:))) );
bar( Angles*180/pi, squeeze(real(spec(150,150,:))) ); 
title('c) 12 angle bins'); xlabel('degrees')



%%%% How many Scales is enough, for the Cauchy? 

%%%%%% PDF file with many figures 
% Create a PDF file name
%fails pdfFileName = 'Cauchy_scales_wavenumbers.pdf';
% Create the PDF file with existing plot
%fails print(gcf, '-dpdf', pdfFileName);


% Let's analyze pure periodic striped patterns
% Bigger array needed here
arraySize2 = 512;
% Create the x and y coordinate arrays
x = 1:arraySize2;
y = 1:arraySize2;
[X, Y] = meshgrid(x, y);

for wavenumber = 2.^(1:1:5) %4:8:80

    %wavenumber = 4; % Adjust for more or fewer stripes
    halfwave = arraySize2/wavenumber/2;

% Create the sinusoidal pattern
    stripes = sin(2 * pi * wavenumber * X / arraySize2);

% Display the image
    figure(wavenumber);
    subplot(211)
    contour(stripes); colorbar;

% Logarithmic 
    Scales = 10.^(1:.025:3) /4. ;
% Equal spaced
    % Scales = 2:2:50;

    Angles = 0;

    cwtCauchy = cwtft2(stripes,wavelet="cauchy",scales=Scales, angles=Angles);
    spec = squeeze( cwtCauchy.cfs );

    subplot(212)
    bar( Scales, squeeze(abs(spec(150,150,:))) ); 
    title(string(wavenumber)+': half wavelength is '+string(halfwave)); xlabel('Scales')

    % bar plot with overlays, doesn't work
    % figure(2); hold off; bar( Scales, squeeze(abs(spec(150,150,:))) ); 

    % Append to the existing PDF file on subsequent iterations
    % fails print(gcf, '-dpdf', '-append', pdfFileName);

end % wavenumber loop

%fails close(gcf); % Close the current figure to prevent it from being displayed
%fails disp(['PDF file saved as: ', pdfFileName]);
