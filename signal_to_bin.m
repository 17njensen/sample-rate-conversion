function audio2bin(infile,outfile)
if(nargin == 1)
    outfile = ['test_',infile(1:max(strfind(infile, '.'))), 'bin'];
end
fprintf('Input file = %s\n', infile);
fprintf('Output file = %s\n', outfile);

%fft
[x,fs] = audioread(infile);
% %for our signals in lab 4
% %sound(x,fs);
% time = [0:(length(x)-1)]/fs;
% t1 = 1;
% t2 = 1.01;
% i1 = round(t1*fs);
% i2 = round(t2*fs);
% nfft = 2^20;
% freq = ([0:nfft-1]/nfft - 0.5)*fs;
% X = fft(x(i1:i2),nfft);
% fmax = max(X(:))
% %X = lowpass(X,500,fs);
% %spectrogram
% [a,fs1]= audioread(infile);
% nfft2 = 2^12;
% overlap = round(0.8*nfft2);
% window = hamming(nfft2);
% fprintf('sample rate = %d\n', fs);
fid = fopen(infile,'rb');
ndim = fread(fid,1,'int');
nchan = fread(fid,1,'int');
dim0 = fread(fid,1,'int');
dim1 = fread(fid,1,'int');
dim2 = fread(fid,1,'int');

%file writing
x = x.';
sx = size(x);
y = [1 size(x) fs 0];
fid = fopen(outfile, 'wb');
fwrite(fid,[1 size(x) fs 0], 'int');
fwrite(fid,x(:),'float');
fclose(fid);

% FFT plot
%subplot(2,2,[3 4]);
plot(freq,20*log10(abs(fftshift(X))));
xlabel('frequency [Hz]','FontSize',10);
ylabel('magnitude [dB]','FontSize',10);
title('Frequency plot of note.wav - Problem 2','FontSize', 10);
xlim([0 500]);
grid on;

whos
return;