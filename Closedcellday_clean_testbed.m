% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
%red = squeeze( video(:,:,1,:) );
%size(red)  % 1080, 1920, 115
% frame50=red(:,:,50); % for vonkarmans case

% v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");
%v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
%video = read(v);
%red = squeeze( video(:,:,1,:) ); % red channel brightness only 
%frame =red(100:700,150:880,3); % for synthetic case
%frame2=red(100:700,150:880,4); DT = 30*60; % 4th image in video, 1/2 h later

v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze( video(:,:,1,:) ); % red channel brightness only 

% Full-res Frames. fx,fy,fframe
fframe = red(1:1000,600:1450,3) ; % for real data image, 3rd image in video
fframe2= red(1:1000,600:1450,4) ; DT = 30*60; % 4th image in video, 1/2 h later
% coordinate arrays at full resolution
fy = 1:size(fframe,1);
fx = 1:size(fframe,2);


% resize images (arrays) to make computations faster: x,y,frame  
% shrinkfactor=1;
shrinkfactor=5; invshrinkfactor = 0.2;
frame = imresize(fframe , invshrinkfactor);
frame2= imresize(fframe2, invshrinkfactor);
% coordinate arrays at low resolution
y = (1:size(frame,1) ) *shrinkfactor;
x = (1:size(frame,2) ) *shrinkfactor;


% QUICK SHOW of fframe contoured over by frame-frame2
figure(1)
image(fx,fy,fframe); colorbar; axis on; hold on; 
contour(x,y,frame2-frame);
title('Frame 1, and contours of Frame2-Frame1')


%%% Wavelet transform inputs (scales and angles) 
% 24 Angles from 0 to pi 
NSCALES = 24;
Angles = 0:pi/NSCALES:pi ;

% a POWERS OF 10 set of Scales, shrinkfactor resizes them
% (~feature size in pixels)
Scales = 10.^(1:.05:1.9);
Scales = Scales/shrinkfactor; 


%%% Call wavelet spectrum 

%%% Ready to test with Morlet 
% wavelet_name = "morlet"; params = [2.0, 3.0, 1.0]
% wavelet_struct = {'name': wavelet_name, 'param': params}
%spec = cwtft2(frame,wavelet="wavelet_name",wavelet_struct,
% scales=Scales,angles=Angles);

% Cauchy one-liner for frame 
cwt  = cwtft2(frame ,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwt.cfs );

% Second image frame2 
cwt2 = cwtft2(frame2,wavelet="cauchy",scales=Scales, angles=Angles);
spec2= squeeze( cwt2.cfs );

% Compute the power spectra and the cross spectrum
power  = abs(spec) .^2;
power2 = abs(spec2).^2;
xspec  = spec2 .* conj(spec);


% Let's hunt for peaks of wavelet power, without knowing location. 
% Power is the squared amplitude, a nonlinear function before averaging
% Compute area-averaged power and xspec, result is fn of angle and scale
% keep away from edges where the unpadded, unwindowed FFT has artifacts

% mean() averages over the first dimension, so two of those is spatial avg
% specify "inner" to be the area a max(Scales) buffer away from edges 

buffer = round( max(Scales) );  % don't use Scales 1/2 the domain or bigger!

innerpower = squeeze( mean(mean( power(buffer:size(power,1)-buffer, ...
                                       buffer:size(power,2)-buffer, :,:) )));
innerpower2= squeeze( mean(mean(power2(buffer:size(power,1)-buffer, ...
                                       buffer:size(power,2)-buffer, :,:) )));

% There's a mean increase of power with scale, divide by it. 
% Why is it there? Study later: is abs(spec) an amplitude in the 
% units of the data? Or is it an INTEGRAL which inherits a factor of scale?
% If that were known, we could clear this up mathematically. 
% For now, we emphasize ANGULAR peaks (what we mean by waves: anisotropy) 
% So just divide by the mean, as usual for power (like the F test in stats)

figure(2)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(Scales, meanbyscale); title('mean power by scale'); xlabel('scale (pixels)')

% Normalize by that mean increase with scale, call it anglespec:
anglespec = innerpower .* 0;     % right sized container for angle spectrum
for isc = 1:NSCALES
    anglespec(:,isc) = squeeze(innerpower(:,isc)) ./ transpose(meanbyscale);
end 

% The angle spectrum: let's have a quick look in index space 
figure(3)
imagesc(anglespec); colorbar(); 
xlabel('Angle index'); ylabel('Scale index')
title('areameanpower/meanbyscale')

% Labeled: this should be a polar plot, using angle as the azimuth 
figure(4)
pcolor(Angles*180/pi, Scales, anglespec); colorbar(); 
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('areameanpower/meanbyscale (learn to do polar plot)')



% Grab the peaks where anglespec exceeds a threshold strength
[row,col] = find( imregionalmax(anglespec) & anglespec>1 )

% illustrate them with calls to image_with_wavelet_overlay()
for ipeak = 1:size(row)
    figure(10+ipeak);
    isc = row(ipeak); ian = col(ipeak);
    contours = 2; 
    image_with_wavelet_overlay(frame, spec, Scales, isc,ian, contours); 
    title('scale, angle: '+string( Scales(isc))+' pixels, ' ...
          +string(Angles(ian)*180./pi)+'deg');
end


% Make coherence-squared and angle spectra from spec .* conj(spec2)
% If spec = R1 + iQ1 with Q meaning quadrature,
% xspec = spec .* conj(spec2) = (R1*R2)+(Q1*Q2) + i (R1Q2+R2Q1)

xpowr = abs(xspec);
%power  = abs(spec) .^2;
%power2 = abs(spec2).^2;

% the area-averaged complex spectrum is innerxspec
innerxspec = squeeze( mean(mean( xspec(buffer:size(xspec,1)-buffer, ...
                                       buffer:size(xspec,2)-buffer, :,:) )));

% coh2 is its magnitude (normalized by abs(spec) and abs(spec2))
inner_coh2 = abs(innerxspec) ./ sqrt(innerpower) ./ sqrt(innerpower2);

% inner_angle is atan(imag/real)  
inner_angl = angle(innerxspec);


figure(5)
imagesc(Angles*180/pi, Scales, inner_coh2); colorbar(); 
xlabel('Angle (deg)'); ylabel('Scale (pixels, approx.)')
title('Squared coherence between images')

figure(6)
pcolor(Angles*180/pi, Scales, inner_angl); colorbar(); hold on
contour(Angles*180/pi, Scales, inner_coh2); colorbar(); 
xlabel('Angle (deg)'); ylabel('Scale (pixels, approx.)')
title('Phase angle and coh2, averaged over Inner area')


