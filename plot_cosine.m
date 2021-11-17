fid = fopen('test_out_16.bin');
%fid = fopen('test_out_8.bin');
%fid = fopen('test_out_4.bin');
ndim  = fread(fid,1,'int');
nchan = fread(fid,1,'int');
dim0  = fread(fid,1,'int');
dim1  = fread(fid,1,'int');
dim2  = fread(fid,1,'int');
[x,fs] = fread(fid,inf,'float');

Lh = 479;
N = 100*Lh;
n = 0:1:N-1;
fo = 1/16; 
%fo = 1/8; 
%fo = 1/4;
xc = cos(2*pi*fo*n);
fout = (2/3)*(1/16);
%fout = (2/3)*(1/8);
%fout = (2/3)*(1/4);
xout = cos(2*pi*fout*n);


fclose(fid);
%stem(x)
plot(x)
hold on;
%stem(xc)
plot(xc)
hold on;
%stem(xout)
plot(xout)

legend('code', '1/16','actual')
