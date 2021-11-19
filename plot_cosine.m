%fid = fopen('test_out_16.bin');
%fid = fopen('test_out_8.bin');
%fid = fopen('test_out_4.bin');
%fid = fopen('ghost_out.bin');
fid = fopen('poly_ghost.bin');
ndim  = fread(fid,1,'int');
nchan = fread(fid,1,'int');
dim0  = fread(fid,1,'int');
dim1  = fread(fid,1,'int');
dim2  = fread(fid,1,'int');
[x_out,fs_out] = fread(fid,inf,'float');
fprintf('fs = %d\n',fs_out);
fprintf('ndim = %d  nchan = %d  dim0 = %d, dim1 = %d, dim2 = %d\n',ndim,nchan,dim0,dim1,dim2);
fid1 = fopen('ghost.bin');
ndim1  = fread(fid1,1,'int');
nchan1 = fread(fid1,1,'int');
dim01  = fread(fid1,1,'int');
dim11  = fread(fid1,1,'int');
dim21  = fread(fid1,1,'int');
[x_in,fs_in] = fread(fid1,inf,'float');
fprintf('fs_input = %d\n',fs_in);
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
sound(x_out,fs_out)

fclose(fid);
%stem(x)
plot(x_out)
hold on;
plot(x_in)
%hold on;
%stem(xc)
%plot(xc)
%hold on;
%stem(xout)
%plot(xout)

%legend('code', '1/16','actual')
whos
