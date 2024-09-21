% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
%red = squeeze( video(:,:,1,:) );
%size(red)  % 1080, 1920, 115
% fframe =red(:,:,50); % for vonkarmans case
% fframe2=red(:,:,60); % for vonkarmans case

v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze( video(:,:,1,:) ); % red channel brightness only 
% Full-res Frames. fx,fy,fframe
fframe = red(1:1000,600:1450,3) ; % 3rd image in video
fframe2= red(1:1000,600:1450,4) ; % 4th image in video, 1/2 h later

% Synthetic data, must strip off axes it was plotted with 
%v = VideoReader("DATA/syn_scwave_warp_mod.mp4");    % weak brightness mod
%v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");  % x5 brightness modul.
%video = read(v);
%red = squeeze( video(:,:,1,:) ); % red channel brightness only 
%fframe =red(125:700,200:880,3); % for synthetic case, must strip off axes
%fframe2=red(125:700,200:880,4); % 4th image in video, 1/2 h later



% coordinate arrays at full resolution
fy = 1:size(fframe,1);
fx = 1:size(fframe,2);
DT = "half-hour";

% resize images (arrays) to make computations faster: x,y,frame  
% shrinkfactor=1;
shrinkfactor=5; invshrinkfactor = 0.2;
frame = imresize(fframe , invshrinkfactor);
frame2= imresize(fframe2, invshrinkfactor);
% coordinate arrays at low resolution
y = (1:size(frame,1) ) *shrinkfactor;
x = (1:size(frame,2) ) *shrinkfactor;


close all % close any figures from last run 


% QUICK SHOW of fframe contoured over by frame2
figure(1)
image(fx,fy,fframe); colormap('gray'); colorbar; axis on; hold on; 
levels = cast(max(fframe(:)),'single') *(20:30)/30.;
contour(x,y,frame2, LevelList=levels, EdgeColor='blue');
title('Frame 1, and contours of Frame2 (low res)')


%%% Wavelet transform inputs (scales and angles) 
% 24 Angles from 0 to pi 
NANGLES = 24;
Angles = 0:pi/NANGLES:pi ;

% a POWERS OF 10 set of Scales, shrinkfactor resizes them
% (~feature size in pixels)
Scales = 10.^(1:.05:1.9); % largest is just a skinch too large, small too
Scales = 10.^(1.05:.05:1.8); % try this trim  

Scales = Scales/shrinkfactor; 
NSCALES = size(Scales,2);
fScales = Scales*shrinkfactor; 

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
% specify "inner" to be the area a 1-Scale length buffer away from edges 

innerpower = squeeze( mean(mean(power))) .*0; % RIGHT SHAPED CONTAINERS
innerpower2=innerpower; innerxspec=innerpower; 

% specify "inner" as 1-Scale length buffer from edges, SCALE BY SCALE
for isc = 1:NSCALES
    buffer = round(Scales(isc)*3);  % for Scales too big, NaN problems appear
    innerpower(isc,:) = squeeze( mean(mean( power(buffer:size(power,1)-buffer, ...
                                           buffer:size(power,2)-buffer, isc,:) )));
    innerpower2(isc,:)= squeeze( mean(mean(power2(buffer:size(power,1)-buffer, ...
                                           buffer:size(power,2)-buffer, isc,:) )));
    innerxspec(isc,:) = squeeze( mean(mean( xspec(buffer:size(xspec,1)-buffer, ...
                                           buffer:size(xspec,2)-buffer, isc,:) )));
end 

% There's a mean increase of power with scale, divide by it. 
% Why is it there? Study later: is abs(spec) an amplitude in the 
% units of the data? Or is it an INTEGRAL which inherits a factor of scale?
% If that were known, we could clear this up mathematically. 
% For now, we emphasize ANGULAR peaks (what we mean by waves: anisotropy) 
% So just divide by the mean, as usual for power (like the F test in stats)

figure(2)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(fScales, meanbyscale); title('mean power by scale'); 
xlabel('fScale (pixels in full-res image)')

% At each angle, normalize by that angle-mean increase with scale, 
% call the result scaleanglespec:

scaleanglespec = innerpower .* 0;     % right sized container for scl-ang 
scaleanglespec2 = innerpower .* 0;    % right sized container again
for iangle = 1:NANGLES
    scaleanglespec(:,iangle) =  squeeze(innerpower (:,iangle)) ./ ...
        transpose(meanbyscale);
    scaleanglespec2(:,iangle) = squeeze(innerpower2(:,iangle)) ./ ...
        transpose(meanbyscale);
end 

% Make coherence-squared and angle spectra from spec .* conj(spec2)
% Everything is an inner area mean (far from boundaries) 
% If spec = R1 + iQ1 with Q meaning quadrature,
% xspec = spec .* conj(spec2) = (R1*R2)+(Q1*Q2) + i (R1Q2+R2Q1)

% coh2 is its magnitude (normalized by abs(spec) and abs(spec2))
inner_coh2 = abs(innerxspec) ./ sqrt(innerpower) ./ sqrt(innerpower2);

% inner_phasediff is atan(imag/real), phase differ from frame to frame2
inner_phasediff = angle(innerxspec);



% The angle spectrum: let's have a quick look, first in index space 
figure(3)
subplot(131)
pcolor(Angles*180/pi, fScales, scaleanglespec); colorbar(); hold on; 
contour(Angles*180/pi, fScales, imregionalmax(scaleanglespec) & scaleanglespec>1)
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('areameanpower/meanbyscale and its peaks frame1')

% Labeled space: this should be a polar plot, using angle as the azimuth 
subplot(132)
pcolor(Angles*180/pi, fScales, scaleanglespec2); colorbar(); hold on; 
contour(Angles*180/pi, fScales, imregionalmax(scaleanglespec2) & scaleanglespec2>1)
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('areameanpower/meanbyscale and peaks frame2')

% coh2 overlaid by power2
subplot(133)
pcolor(Angles*180/pi, fScales, inner_coh2); colorbar()
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('coherence squared')


%%% For every peak in the spectrum averaged over whole image area, 
%%% Let's locate an area where it is prominent, and set the contour 
%%% interval to make that area clear on an annotated image. 

% SCALE AND ANGLE of all spectral peaks where scaleanglespec also
% exceeds a threshold (like 1, the ratio of power to anglemean power)
[row,col] = find( imregionalmax(scaleanglespec) & scaleanglespec>1 );


% Loop over these spectral peaks. For each, find its peak activity on map. 
 for ipeak = 1:size(row,1)
    isc = row(ipeak)
    ian = col(ipeak);

%%% Assign a phase speed. Note this applies to whole INNERMEAN area!
%%% Could redo the inner_phasediff calculation AREA BY AREA if needed  
% Phase speed (pixels per image timegap) is phasediff/pi * Scale
% Scale appears to be the size of a ridge or trough (pi of angle)
% Then, express it in original image pixel units 
    speed = abs( inner_phasediff(isc,ian) )/pi * fScales(isc);

% print header about this spectral peak
    "Scale-Angle spectrum peak " + string(ipeak)
    "at Scale (pixels), Angle (deg): "+...
        string(fScales(isc))+', '+string(Angles(ian)*180/pi)
    "Its mean phase speed is "+string(speed) + " pixels per DT "+DT

% Where is this **spectral** peak spatially, in image space? Where are
% locations on the map where there is a **spatial** peak of this power
    powermap = abs(spec(:,:,isc,ian)) .^2;
    powerpeaks = powermap( find( imregionalmax(powermap) )); % values 

% find x and y locations of all peaks of power at this scale,angle
    [ypeaks,xpeaks] = find( imregionalmax(powermap) );
    fxpeaks = xpeaks * shrinkfactor; fypeaks = ypeaks * shrinkfactor;

% Find the BIGGEST peak that is NOT within one fScale of the boundary

% First, find ones that are clean of the boundary:
    TOOCLOSE = fScales(isc);
    FARENOUGH = ...
        fxpeaks>          TOOCLOSE & fypeaks>          TOOCLOSE & ...
        (max(fx)-fxpeaks)>TOOCLOSE & (max(fy)-fypeaks)>TOOCLOSE ;
    clean = find( FARENOUGH );

% Select the BIGGEST POWER one, to use ypeak,xpeak for a plot & stats
% If there is not a "clean" biggest, use biggest 
    biggest = find( powerpeaks == max(powerpeaks) );
    cbiggest = find( powerpeaks(clean) == max(powerpeaks(clean)) );
    ypeak = ypeaks(clean(cbiggest));
    xpeak = xpeaks(clean(cbiggest));

% if none is clean, show biggest even if too near boundary (or 
% (or should I continue loop discarding this one?)
    if(size(clean,1) == 0)
        "NO PEAKS ARE AWAY FROM EDGES"
        ypeak = ypeaks(biggest);
        xpeak = xpeaks(biggest);
    end

    fxpeak = xpeak * shrinkfactor; 
    fypeak = ypeak * shrinkfactor;

% OLD DUMB Find the peak closest to the center of the map 
%    centerdistance = (fypeaks - mean(y)).^2 + (fxpeaks - mean(x)).^2;
%    centermost = find( centerdistance==min(centerdistance) );


% for selected peak ypeak,xpeak 
% assign amplitudes for setting contour levels on the annotated plot 
    ampli = abs( spec (ypeak,xpeak,isc,ian) );
    ampli2= abs( spec2(ypeak,xpeak,isc,ian) );

% BAIL OUT STEP not to make a plot if amplitude is too small 
    if(ampli < 1)
        'Amplitude '+string(ampli)+' too small, skipping it'
        continue
    end


% print information 
    "Spatial peaks  are at of its power: x=" + ...
        string(fxpeaks) + ', y=' + string(fypeaks)
    "Biggest one is "+ ...
        string(fxpeak) + ' , ' + string(fypeak)
    "Amplitude 1,2 = " + string( ampli ) + ', ' + string(ampli2)
    "Speed " +string(speed) + ' pixels per ' + DT

 % Annotate the full frame image with this peak,scaling contours for ampli 
 % contours are -10 to 10 / clevfactor so make clev=ampli/3 or so 
    figure(10+ipeak);
    clev = ampli/3; 

    image_with_wavelet_overlay(fframe, spec,x,y, Scales, isc,ian, clev); 

    title('scale, angle: '+string(fScales(isc))+' pixels, ' ...
          +string(Angles(ian)*180./pi)+'deg, BestInner peak @x,y = ' + ...
           string(fxpeak) + ' , ' + string(fypeak)...
          +' clev is '+string(clev) )

end



% Coherence squared is in [0,1]
%figure(5)
%pcolor(Angles*180/pi, fScales, inner_coh2); colorbar(); 
%xlabel('Angle (deg)'); ylabel('Scale (pixels, approx.)')
%title('Squared coherence between images')

% Phase angle is in [-pi,pi], multiply by scale to get a speed 
%figure(6)
%pcolor(Angles*180/pi, fScales, inner_phasediff); colorbar(); hold on
%contour(Angles*180/pi, fScales, inner_coh2,'black', ...
%    LevelList=(1:80)/100.); 
%colorbar(); 
%xlabel('Angle (deg)'); ylabel('Scale (pixels, approx.)')
%title('Phase angle masked by low coh2, averaged over Inner area')


