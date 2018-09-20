function p = GammaPower (parameters, x)
%This is a gamma power function. Helpful in fitting luminance measurements from displays for gamma correction.
%Of the two parameters, gamma is the exponent and A is a constant, usually around 1.
%In the case of display calibration, the dependent variable x would be the raw luminance (RGB) values sent to the screen.
%The gamma exponents for each RGB channel are what you need to compute the inverse gamma and thereby gamma 'correct' your screen.
%R Maloney, April 2015

Gam = parameters(1); 
A = parameters(2); 

p = A .* (x .^ Gam );

end

