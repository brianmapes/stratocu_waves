function output = preprocess_img(raw)

% Pre-process sat cloud image data before wavelet projection/diagnosis 
dat = cast(raw,'double');

%%% First, compress dynamic range on the dark end (we seek waves in cloud) 
% Avoid zeros for log and sqrt functions
whiterer = cast(max(dat(:)),'double') * 1.05;

logdarkness   = log ( whiterer - dat );
sqrt_darkness = sqrt( whiterer - dat );

%%%% Subtract the mean and detrend and window the data (mute the edges).
% Detrend removes the mean and linear gradient (like sunshine brightness)
% If A is a matrix, then detrend operates on each column separately, 
% subtracting each trend from the corresponding column of A.
% Since "columns" are privileged, apply both to img and transpose(img)

% Do it for -log darkness and -sqrt darkness
det = detrend(-logdarkness);
detlog = transpose( detrend( transpose(det) ));

det = detrend(-sqrt_darkness);
detsqrt = transpose( detrend( transpose(det) ));

%%% Mute the values (which are now statistically centered on 0) at edges 
% This may be called a "windowing" like Hanning or Hamming 
% or "vignette" in 2D image space 

window = cast( raw*0.0 + 1.0, 'double');  % right sized array of ones 
% make the ones go to zero near the edges... somehow...


% The output of all the above 
output = detlog .* window;

end