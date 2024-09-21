
function output = image_with_wavelet_overlay(img,spec,x,y, Scales, scale, angle,clevfactor)
    % Overlay wavelet power on image 
    image(img); colormap(gray); colorbar; axis on

    hold on

    posLevels = (1:2:9)   /clevfactor;
    negLevels = (-9:2:-1) /clevfactor;

    % Adjust contour levels by Scale as factor (for real/imag), factor^2 (for power)
    sfactor = Scales(scale);
    
    % Real part is crests and trofs, imag is gradients, abs is a magnitude map 
    contour( x,y,real(spec(:,:,scale,angle)), LevelList=posLevels*sfactor, EdgeColor='red' );
    contour( x,y,real(spec(:,:,scale,angle)), LevelList=negLevels*sfactor, EdgeColor='blue' );
    
    % "power" is amplitude abs() squared 
    contour( x,y,(abs(spec(:,:,scale,angle))*sfactor) .^2, EdgeColor='white' );

    % legend
end

