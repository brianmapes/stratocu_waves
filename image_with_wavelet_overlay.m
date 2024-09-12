
function output = image_with_wavelet_overlay(img,spec, Scales, scale, angle)
    % Overlay wavelet power on image 
    imshow(img); colorbar; axis on
    
    hold on

    posLevels = 1:2:9;
    negLevels = -9:2:-1;

    % Adjust contour levels by Scale as factor (for real/imag), factor^2 (for power)
    factor = Scales(scale);
    
    % Real part is crests and trofs, imag is gradients, abs is a magnitude map 
    contour( real(spec(:,:,scale,angle)), LevelList=posLevels*factor, EdgeColor='red' );
    contour( real(spec(:,:,scale,angle)), LevelList=negLevels*factor, EdgeColor='blue' );
    
    % "power" is amplitude abs() squared 
    contour( (abs(spec(:,:,scale,angle))*factor) .^2 );

    % legend
end

