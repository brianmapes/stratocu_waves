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
Angles = 0:pi/12:pi

% a LOGARITHMIC set of 10 Scales
% Scales = [2,5,10,20,40,80,160,320,640,1000] % 2,5,10 too small, 1000 too big
% Better range 
Scales = 10.^(1.2:.2:3);

cwtCauchy = cwtft2(frame50,wavelet="cauchy",scales=Scales, angles=Angles);

spec = squeeze( cwtCauchy.cfs );

% RESULT: size 1080,1920, 10 SCALES, 8 ANGLES


% Amplitude seems to grow roughly linearly with scale. 
% figure(2); plot( Scales,squeeze(abs(spec(500,500,:,7))) )

% function image_with_wavelet_overlay(spec, Scales, scale, angle)

image_with_wavelet_overlay(frame50, spec, Scales, 3, 6);

power = abs(spec) .^2;