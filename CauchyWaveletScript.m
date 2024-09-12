% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");

video = read(v);

% red channel brightness only 
red = squeeze( video(:,:,1,:) );
% size(red) 1080, 1920, 115

%%% SELECT ONE FRAME FOR SPATIAL ANALYSIS
% First argument is y increasing downward, second is x from left
% frame50=red(:,:,50); % for vonkarmans case
frame50=red(100:700,150:880,3); % for synthetic case

% QUICK SHOW 
figure(1)
imshow(frame50); colorbar;


%%% Wavelet transform 

% a set of 7 Angles from 0 to pi (N, NNW, WNW, W, WSW, SSW, )
Angles = 0:pi/12:pi ;

% a LOGARITHMIC set of 10 Scales
% Scales = [2,5,10,20,40,80,160,320,640,1000] % 2,5,10 too small, 1000 too big
% Better range 
Scales = 10.^(1:.1:1.9) ;

cwtCauchy = cwtft2(frame50,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwtCauchy.cfs );

% SIZE: size Y=1080, X=1920, 10 SCALES, 8 ANGLES 



% PLOT ANNOTATED IMAGE

% function image_with_wavelet_overlay(spec, Scales, scale, angle)
figure(1)
image_with_wavelet_overlay(frame50, spec, Scales, 4, 4)
title('scale4, angle4')
hold off 

% Let's hunt for easter eggs of wavelet power, without knowing location. 
% Compute area-averaged power by angle and scale: away from the edges 
% mean() averages over the first dimension, so two of those will make

power = abs(spec) .^2;
innerpower = squeeze(  mean(mean( power(100:500,100:600, :,:) ))  );


% There's a mean increase of power with scale, normalize it away
figure(2)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(Scales, meanbyscale); title('mean power by scale')

% Normalize by that mean increase with scale, call it anglespec:
anglespec = innerpower .* 0;     % right sized container for angle spectrum
for isc = 1:10   %size( transpose(Scales) )
    anglespec(:,isc) = squeeze(innerpower(:,isc)) ./ transpose(meanbyscale);
end 

% The angle spectrum 
figure(3)
pcolor(anglespec); colorbar(); 
title('areameanpower/meanbyscale')
hold off 


% Prominent is 2,9
figure(4)
image_with_wavelet_overlay(frame50, spec, Scales, 2, 9)
title('scale 2, angle 9')
hold off 

% another is 8,5
figure(5)
image_with_wavelet_overlay(frame50, spec, Scales, 8, 5)
title('scale 8, angle 5')
hold off 