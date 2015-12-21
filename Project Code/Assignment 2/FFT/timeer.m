x = rand(2^12,1);
clc
fprintf('Using the slow direct method (Matrix product): ');
tic
slowft(x);
toc


fprintf('\n\nUsing the recursive method (radix-2): ');
tic
symft(x);
toc


fprintf('\n\nUsing the inbuilt method (MATLAB-FFTW): ');
tic
fft(x);
toc


