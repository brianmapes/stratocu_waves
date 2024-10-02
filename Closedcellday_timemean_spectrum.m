% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
%red = squeeze( video(:,:,1,:) );
%size(red)  % 1080, 1920, 115
% fframe =red(:,:,50); % for vonkarmans case
% fframe2=red(:,:,60); % for vonkarmans case

v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze( video(:,:,1,:) ); % red channel brightness only 
% Full-res Frames. fx,fy,fframe 600:1450 in x to avoid broken area?
raw = red(1:1000,1:1450,:) ; % all B&W images in video

% Synthetic data, must strip off axes it was plotted with 
%v = VideoReader("DATA/syn_scwave_warp_mod.mp4");    % weak brightness mod
%v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");  % x5 brightness modul.
%video = read(v);
%red = squeeze( video(:,:,1,:) ); % red channel brightness only 
%fframe =red(125:700,200:880,3); % for synthetic case, must strip off axes
%fframe2=red(125:700,200:880,4); % 4th image in video, 1/2 h later


%%% Pre-process the data. Compress dynamic range of each frame, 
% Subtract the mean, detrend in x and y, window the values (mute the edges).
fframes = raw; % make a copy of the right size 
for it = 1:size(raw,3)
    fframes(:,:,it) = preprocess_img(raw(:,:,it));
end

% coordinate arrays at full resolution
fy = 1:size(fframes,1);
fx = 1:size(fframes,2);
DT = "half-hour";

% resize images (arrays) to make computations faster: x,y,frame  
% shrinkfactor=1;
shrinkfactor=5; invshrinkfactor = 0.2;
frames = imresize(fframes , invshrinkfactor);

% coordinate arrays at low resolution
y = (1:size(frames,1) ) *shrinkfactor;
x = (1:size(frames,2) ) *shrinkfactor;

close all % close any figures from last run 


%%% Wavelet transform inputs (scales and angles) 
% 24 Angles from 0 to pi 
NANGLES = 24;
Angles = 0:pi/NANGLES:pi ;

% a POWERS OF 10 set of Scales, shrinkfactor resizes them
% (Scale ~ feature size in pixels, empirically)
%Scales = 10.^(1:.05:1.9); % largest is just a skinch too large, small too
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
cwt  = cwtft2(frames ,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwt.cfs );

% Compute the power spectrum and its TIME AVERAGE 
powpow  = abs(spec) .^2;
power = squeeze( mean(powpow,3) );  % time average! 

% Let's hunt for peaks of wavelet power, without knowing location. 
% Compute area-averaged power and xspec, result is fn of angle and scale
% keep away from edges where the unpadded, unwindowed FFT has artifacts

% mean() averages over the first dimension, so two of those is spatial avg
% specify "inner" to be the area a 1-Scale length buffer away from edges 

%innerpower  = squeeze( mean(mean(power ))) .*0; % RIGHT-SHAPED CONTAINER
innerpowpow = squeeze( mean(mean(powpow))) .*0; % time dependent 

% specify "inner" as 1-Scale length buffer from edges, SCALE BY SCALE
for isc = 1:NSCALES
    buffer = round(Scales(isc)*3);  % for Scales too big, NaN problems appear
    for it = 1:size(fframes,3)
        innerpowpow(it,isc,:) = squeeze( mean(mean( ...
                  powpow(buffer:size(powpow,1)-buffer, ...
                         buffer:size(powpow,2)-buffer, it,isc,:) )));
    end
end

innerpower = squeeze(mean(innerpowpow,1)); % time mean 

% There's a mean increase of power with scale, divide by it. 
% Why is it there? Study later: is abs(spec) an amplitude in the 
% units of the data? Or is it an INTEGRAL which inherits a factor of scale?
% If that were known, we could clear this up mathematically. 
% For now, we emphasize ANGULAR peaks (what we mean by waves: anisotropy) 
% So just divide by the mean, as usual for power (like the F test in stats)

figure(3)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(fScales, meanbyscale); title('mean power by scale'); 
xlabel('fScale (pixels in full-res image)')

% At each angle, normalize by that angle-mean increase with scale, 
% call the result scaleanglespec:

scaleanglespect = innerpowpow .* 0;   % right sized container time-scl-ang 
for iangle = 1:NANGLES
    for it = 1:size(fframes,3)
        scaleanglespect(it,:,iangle) = ...
            innerpowpow (it,:,iangle) ./ meanbyscale;
    end
end 
scaleanglespec = squeeze(mean(scaleanglespect,1));   % time mean
% The angle spectrum: let's have a quick look, first in index space 
figure(4)

isosurface(scaleanglespect,1)

%pcolor(Angles*180/pi, fScales, scaleanglespec); colorbar(); hold on; 
%contour(Angles*180/pi, fScales, imregionalmax(scaleanglespec) & scaleanglespec>1)
%xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
%title('time-areameanpower/meanbyscale and its peaks')

%%% For every peak in the spectrum averaged over whole image area, 
%%% Let's locate an area where it is prominent, and set the contour 
%%% interval to make that area clear on an annotated image. 

% SCALE AND ANGLE of all spectral peaks where scaleanglespec also
% exceeds a threshold (like 1, the ratio of power to anglemean power)
[row,col] = find( imregionalmax(scaleanglespec) & scaleanglespec>1 );


% Loop over time-mean spectral peaks. For each, find peak activity on map. 
 for ipeak = 1:size(row,1)
    isc = row(ipeak);
    ian = col(ipeak);

% print header about this spectral peak
    "Scale-Angle spectrum peak " + string(ipeak)
    "at Scale (pixels), Angle (deg): "+...
        string(fScales(isc))+', '+string(Angles(ian)*180/pi)

    % Where is this **spectral** peak spatially, in image space? Where are
% locations on the map where there is a **spatial** peak of this power
    powerpeaks = power( find( imregionalmax(power) )); % values 

% find x,y locations of all peaks of time-mean power at this scale,angle
    [ypeaks,xpeaks] = find( imregionalmax(power) );
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

% for selected peak ypeak,xpeak 
% assign amplitudes for setting contour levels on the annotated plot 
%    ampli = abs( spec (ypeak,xpeak,isc,ian) );
%    ampli2= abs( spec2(ypeak,xpeak,isc,ian) );
% BAIL OUT STEP not to make a plot if amplitude is too small 
%    if(ampli < 0)
%        'Amplitude '+string(ampli)+' too small, skipping it'
%        continue
%    end

% Loop over all frames, making a quick look of each 
    for it = 1:size(fframes,3)

 % Annotate the full frame image with this peak,scaling contours for ampli 
 % contours are -10 to 10 / clevfactor so make clev=ampli/3 or so 
        figure(100*ipeak + it);
        specnow = squeeze( spec(:,:,it,:,:) );
        ampli = abs( specnow(ypeak,xpeak,isc,ian) );
        clev = ampli/5; 

        image_with_wavelet_overlay(squeeze(raw(:,:,it)), ...
                           specnow,x,y, Scales, isc,ian, clev); 

        title('scale, angle: '+string(fScales(isc))+' pixels, ' ...
          +string(Angles(ian)*180./pi)+'deg, BestInner peak @x,y = ' + ...
           string(fxpeak) + ' , ' + string(fypeak)...
          +' clev is '+string(clev) )
    end % end it loop over times

end

