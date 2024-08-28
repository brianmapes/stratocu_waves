% Import video 
v = VideoReader("Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
video = read(v);

% red channel brightness only 
red = squeeze( video(:,:,1,:) );
% size(red) 1080, 1920, 115

frame50=red(:,:,50);
figure(1)
imshow(frame50); colorbar;



%%% Wavelet transform 

% a set of 8 Angles from 0 to pi 
Angles = 0:pi/8:pi-pi/8

% a set of 10 Scales
Scales = 5:5:50

cwtCauchy = cwtft2(frame50,wavelet="cauchy",scales=Scales, angles=Angles);

spec = squeeze( cwtCauchy.cfs );
size(spec) 

% RESULT: size 1080,1920, 10 SCALES, 8 ANGLES


% Overlay wavelet power on image 
figure(1)
imshow(frame50); colorbar; axis on

hold on

posLevels = 10:20:90; negLevels = -90:20:-10;

% Real part is crests and trofs, imag is gradients, abs is a magnitude map 
contour( real(spec(:,:,10,7)), LevelList=posLevels,EdgeColor='red' );
contour( real(spec(:,:,10,7)), LevelList=negLevels,EdgeColor='blue' );

% "power" is amplitude abs() squared 
contour( abs(spec(:,:,10,7)).^2 );