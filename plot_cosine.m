%fid = fopen('test_out_16.bin');
%fid = fopen('test_out_8.bin');
%fid = fopen('test_out_4.bin');
%fid = fopen('ghost_out.bin');
%fid = fopen('poly_ghost.bin');
%fid = fopen('poly_16.bin');
%fid = fopen('poly_8.bin');
fid = fopen('poly_4.bin');
%fid = fopen('ghost.bin');
ndim  = fread(fid,1,'int');
nchan = fread(fid,1,'int');
dim0  = fread(fid,1,'int');
dim1  = fread(fid,1,'int');
dim2  = fread(fid,1,'int');
[x_out,fs_out] = fread(fid,inf,'float');
nfft = 2^9;
nvals = (0:nfft-1)/nfft;
X_out = fft(x_out,nfft);
%plot(nvals,(abs(fftshift(X_out))));
fprintf('fs = %d\n',fs_out);
fprintf('ndim = %d  nchan = %d  dim0 = %d, dim1 = %d, dim2 = %d\n',ndim,nchan,dim0,dim1,dim2);
%fid1 = fopen('ghost.bin');
fid1 = fopen('test_out_8.bin');
ndim1  = fread(fid1,1,'int');
nchan1 = fread(fid1,1,'int');
dim01  = fread(fid1,1,'int');
dim11  = fread(fid1,1,'int');
dim21  = fread(fid1,1,'int');
[x_in,fs_in] = fread(fid1,inf,'float');

fid2 = fopen('test_out_4.bin');
ndim2  = fread(fid2,1,'int');
nchan2 = fread(fid2,1,'int');
dim02  = fread(fid2,1,'int');
dim12  = fread(fid2,1,'int');
dim22  = fread(fid2,1,'int');
[x_a,fs_a] = fread(fid2,inf,'float');
fprintf('fs_input = %d\n',fs_in);
Lh = 154;
%nfft = 10240000;
N = 100*Lh;
n = 0:1:N-1;
nvals = (0:nfft-1)/nfft;
%fo = 1/16; 
%fo = 1/8; 
fo = 1/4;
xc = cos(2*pi*fo*n);
%plot(nvals,abs(fft(xc,nfft)))
%hold on;
fout = (2/3)*(1/16);
%fout = (2/3)*(1/8);
%fout = (2/3)*(1/4);
xout = cos(2*pi*fout*n);
%plot(nvals,abs(fft(xout,nfft)))
%legend('original', 'result')
%sound(x_out,fs_out)

fclose(fid);
plot(nvals,20*log10(abs(fftshift(X_out))));
%hold on;
%plot(nvals,20*log10(abs(fftshift(fft(x_in,nfft)))));
%hold on;
%plot(nvals,20*log10(abs(fftshift(fft(x_a, nfft)))));
hold on;
plot(nvals,20*log10(abs(fftshift(fft(xc,nfft)))));
%stem(x)
%plot(x_out)
%plot(fftshift(abs(x_out)));
%plot(abs(fftshift(x_out)));
%hold on;
%plot(x_in)
%plot(abs(fftshift(x_in)));
%hold on;
%plot(x_a)
%plot(abs(fftshift(x_a)));
%plot(x_in)
%hold on;
%stem(xc)
%plot(xc)
%hold on;
%stem(xout)
%plot(xout)
xlabel('Frequency(Hz)', 'FontSize', 10);
ylabel('Amplitude(dB)','FontSize', 10);
title('Cosine - Polyphase','FontSize', 10);
legend('Output', 'Input of 1/4')
whos
