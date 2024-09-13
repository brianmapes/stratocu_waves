% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
%red = squeeze( video(:,:,1,:) );
%size(red)  % 1080, 1920, 115
% frame50=red(:,:,50); % for vonkarmans case

% v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");
%v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
%video = read(v);
%red = squeeze( video(:,:,1,:) ); % red channel brightness only 
%frame50=red(100:700,150:880,3); % for synthetic case

v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze( video(:,:,1,:) ); % red channel brightness only 
frame50 = red(1:1000,600:1450,3); % for real data image

% QUICK SHOW 
figure(1)
imshow(frame50); colorbar;

%%% Wavelet transform 

% a set of 7 Angles from 0 to pi (N, NNW, WNW, W, WSW, SSW, )
Angles = 0:pi/24:pi ;

% a LOGARITHMIC set of 10 Scales
% Scales = [2,5,10,20,40,80,160,320,640,1000] % 2,5,10 too small, 1000 too big
% Better range 
Scales = 10.^(1:.05:1.9) ;

cwtCauchy = cwtft2(frame50,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwtCauchy.cfs );


% Let's hunt for easter eggs of wavelet power, without knowing location. 
% Compute area-averaged power by angle and scale: away from the edges 
% mean() averages over the first dimension, so two of those will make

power = abs(spec) .^2;

% for synthetic data (has axes in image area, must clip them) 
% innerpower = squeeze(  mean(mean( power(100:500,100:600, :,:) ))  );

% for actual data (avoid edge effects --> some "inner" box)
innerpower = squeeze(  mean(mean( power(50:800,200:800, :,:) ))  );


% There's a mean increase of power with scale, normalize it away
figure(2)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(Scales, meanbyscale); title('mean power by scale'); xlabel('scale (pixels)')

% Normalize by that mean increase with scale, call it anglespec:
anglespec = innerpower .* 0;     % right sized container for angle spectrum
for isc = 1:24   %size( transpose(Scales) )
    anglespec(:,isc) = squeeze(innerpower(:,isc)) ./ transpose(meanbyscale);
end 

% The angle spectrum 
figure(3)
pcolor(anglespec); colorbar(); 
xlabel('Angle index'); ylabel('Scale index')
title('areameanpower/meanbyscale')
hold off 

figure(4)
pcolor(Angles*180/pi, Scales, anglespec); colorbar(); 
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('areameanpower/meanbyscale')
hold off 


figure(7)
image_with_wavelet_overlay(frame50, spec, Scales, 9,8); title('scale 9 angle 8')
figure(8)
image_with_wavelet_overlay(frame50, spec, Scales, 15,14); title('scale 15 angle 14')
