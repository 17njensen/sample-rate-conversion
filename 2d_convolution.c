#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "image2.txt"
//*******************************************************************************************//
//  Nik Jensen
//  Digital Signal Image Processing
//  A02250195
//*******************************************************************************************//
//pgm file format (MxN matrix = row X col)
//P5 - file type
//Width (n)
//Height (m) 
//Max num val in image
//cycles from top row down to bottom row (of size height) while each row is of size width
//******************************************************************************************//
/*
Pesudocode
    Read in input file
        Fetch image type
            Print warning if not correct image type
        Fetch image width and height
        Fetch max num val
        Loop i=0, i<height-1
            Loop j=0, j<width-1
                increment i, storing each bit into the [j]th row's [i] column
    Convolve with h
        determine size of h
        formulate array of size (width + hm -1)x(height + hn -1) for zero pad 
        convolve:
            Loop x,y = 0, x,y<(z_height,z_width)
                Loop m=0, m<(z_height)
                    Loop n=0, n<(z_width)
                        z_array[x][y] = f(m,n)*g(x-m,y-n)
    Output results
        Print image type
        Print image width and height
            z_width, z_height
        print max num val
        Print z_array results
    close allocated arrays
    close in and out files


*/

/* https://www.geeksforgeeks.org/how-to-read-a-pgmb-format-image-in-c/
*/
typedef struct{
    char *in_type;
    int in_width;
    int in_height;
    int in_max;
} file_header;

void main(int argc, char** argv) {
    FILE* infile;
    FILE* outfile;
    FILE* impfile;
    unsigned char buffer[1096];
    if (NULL == (infile = fopen(argv[1], "rb"))) { //error check and open file
        printf("error: Cannot open image.pgm file for input.\n");
        return;
    }
    if (NULL == (outfile = fopen(argv[2], "wb"))) { //error check and open file
        printf("error: Cannot open output file for writing.\n");
        return;
    }
    if (NULL == (impfile = fopen(argv[3], "rb"))) { //error check and open file
        printf("error: Cannot open impulse file for reading.\n");
        return;
    }
    char line[1096];
    char line1[1096];

    // rotate matrix
    // for(int i = 0; i < 3; i++){
    //     for(int j = 0; j < 3; j++){
    //         // int temp = lpf[i][j];
    //         // lpf[i][j] = lpf[3-i-1][3-j-1];
    //         // lpf[3-i-1][3-j-1] = temp;

    //         // int temp = H[i][j];
    //         // H[i][j] = H[5-i-1][5-j-1];
    //         // H[5-i-1][5-j-1] = temp;

    //         int temp = S1[i][j];
    //         S1[i][j] = S1[3-i-1][3-j-1];
    //         S1[3-i-1][3-j-1] = temp;
    //         int temp1 = S2[i][j];
    //         S2[i][j] = S2[3-i-1][3-j-1];
    //         S2[3-i-1][3-j-1] = temp1;
    //     }
    // }
    //fetch image header
    file_header h0, h_1;
    for (int i = 0; i <= 21; i++) {
        if (fgets(line, sizeof(line), infile)) { //fetch one line at a time for header
            if (line[0] == '#') { //if a comment exists in the header
                continue;
            }
            if (i == 0) { //file type
                printf("in here\n");
                h0.in_type = "P5";
                continue;
            }
            if (i == 20) {
                for (int j = 0; j < 2; j++) {
                    char* delim = strtok(line, " ");
                    if (j == 0) {
                        h0.in_width = atoi(delim);
                        delim = strtok(NULL, " ");
                    }
                    else if (j == 1) {
                        h0.in_height = atoi(delim);
                        delim = strtok(NULL, " ");
                    }
                }
                continue;
            }
            if (i == 21) {
                h0.in_max = atoi(line);
                continue;
            }
        }
    }
    for (int i = 0; i <= 21; i++) {
        if (fgets(line1, sizeof(line1), impfile)) { //fetch one line at a time for header
            if (line1[0] == '#') { //if a comment exists in the header
                continue;
            }
            if (i == 0) { //file type
                printf("in here\n");
                h_1.in_type = "P5";
                continue;
            }
            if (i == 20) {
                for (int j = 0; j < 2; j++) {
                    char* delim = strtok(line1, " ");
                    if (j == 0) {
                        h_1.in_width = atoi(delim);
                        delim = strtok(NULL, " ");
                    }
                    else if (j == 1) {
                        h_1.in_height = atoi(delim);
                        delim = strtok(NULL, " ");
                    }
                }
                continue;
            }
            if (i == 21) {
                h_1.in_max = atoi(line1);
                continue;
            }
        }
    }
    printf("\tfile type = %s\n", h_1.in_type);
    printf("\theight = %d\n", h_1.in_height);
    printf("\twidth = %d\n", h_1.in_width);
    printf("\tMax num val in image = %d\n", h_1.in_max);
    int Rh = 160; //length of impulse response
    int Ch = 165;
    int Rx = h0.in_width; // Width of input image (n)
    int Cx = h0.in_height; //Height of input image (m)
    int Ry = Rx + Rh - 1; //length of convolution result
    int Cy = Cx + Ch - 1;
    int Rz = Rx + 2*(Rh-1); //length with zero padding
    int Cz = Cx + 2*(Ch-1); //height with zero padding
    printf("Rh = %d, Rx = %d, Ry = %d, Rz = %d\n", Rh, Rx, Ry, Rz);

    unsigned char* x = (unsigned char*)calloc(sizeof(unsigned char), Rz*Cz);
    unsigned char* h = (unsigned char*)calloc(sizeof(unsigned char), Rh*Ch);
    float h1[165][160];
    // float h2[160][160];
    float* h2 = (float*)calloc(sizeof(float), Rh*Ch);
    float* x1 = (float*)calloc(sizeof(float), Rz*Cz);
    float* x_sq = (float*)calloc(sizeof(float), Rz*Cz);
    float* y = (float*)calloc(sizeof(float), Ry*Cy);
    float* y1 = (float*)calloc(sizeof(float), Ry*Cy);
    float* y2 = (float*)calloc(sizeof(float), Ry*Cy);
    unsigned char* y0 = (unsigned char*)calloc(sizeof(unsigned char), Ry*Cy);
    // int a = Rz*Ch + Rh+1; //start Ch rows down to include zero padding
    int a = Ch-1;
    while (!feof(infile)){
        //read in each line with each pixel of size unsigned char (8 bytes)
        fread(x + (a*Cz + Ch-1), sizeof(unsigned char), Rx, infile); //input file
        a++;
    }
    while (!feof(impfile)){
        //read in each line with each pixel of size unsigned char (8 bytes)
        fread(h, sizeof(unsigned char), Rh, impfile); //impulse
    }
    for(int k = 0; k < Rz*Cz; k++){ //convert x to float for convolution and square for part 4
        x1[k] = (float)(x[k]*x[k]);//have to square each pixel that is read in. Possibly read in, square, then convert to float
    }
    printf("h[10] = %u\n",h[10]);
    for(int m = 0; m < Rh; m++){
        for(int n = 0; n < Ch; n++){
            h1[m][n] = 1.0;
        }
    }
    for(int k = 0; k < Rh*Ch; k++){
        h2[k] = (float)h[k];
    }

    printf("h2[10] = %f\n",h2[10]);
    printf("done reading\n");
    float tmp0,tmp1;
    int i;
    printf("Start Convolution\n");
    //for part 1, 2, 3
    for (int m = 0; m < Ry; m++){ //go to next row
        printf("m = %d\n", m);
        for (int n = 0; n < Cy; n++){ //go through columns of one row 
        printf("n = %d\n", n);
            //convolution with impulse
            for (tmp0 = 0.0,tmp1 = 0.0,i=0; i<Rh; i++){ //convolves for one point of n inside row m
                for (int j = 0; j<Ch; j++){
                    // tmp0 += h[i][j] * x1[(m+i)*Cz + (n+j)]; //for convolution
                    tmp0 += h1[i][j] * x1[(m+i)*Cz + (n+j)];
                    tmp1 += h2[m*Ch + n] * x1[(m+i)*Cz + (n+j)];
                }
            }
            y1[m*Cy+n] = tmp0;
            y2[m*Cy+n] = tmp1;
            y[m*Cy+n] = ((y2[m*Cy+n])/(y2[m*Cy+n]));
        }
    }
    printf("finished convolution\n");
    for(int k = 0; k < Ry*Cy; k++){//convert output to unsigned char
        y0[k] = (unsigned char)y[k];
    }

    for(int i = 0; i < Ry; i++){ //check if values are higher or lower than bounds provided by file
        for(int j = 0; j < Cy; j++){
            if(y0[i*Cy+j] > h0.in_max){
                y0[i*Cy + j] = h0.in_max;
            }
            else if(y0[i*Cy + j] < 0){
                y0[i*Cy + j] = 0;
            }
        }
    }
    printf("finish and write\n");
    // fprintf(outfile, "%s\n%d\n%d\n%d\n", h0.in_type, Rz, Cz, h0.in_max);
    // fwrite(x_sq, sizeof(unsigned char), Rz*Cz, outfile);
    fprintf(outfile, "%s\n%d\n%d\n%d\n", h0.in_type, Ry, Cy, h_1.in_max);
    fwrite(y0, sizeof(unsigned char), Ry*Cy, outfile);
    fclose(infile);
    fclose(outfile);
    fclose(impfile);
    printf("done\n");
}
