#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "filter_3.h"
#include "cos_16.h"
#include "cos_8.h"
#include "cos_4.h"
#include <time.h>

/*
This is for the second computer assignment in DSIP. 
Sample rate of input is 11,025 Hz
Period = T = 1/11,025 Hz = 9.07x10^-5

3 Sections: 
1. Upsample by L
2. Convolve with a lowpass
3. Downsample by M 
Lh = 154
L/M = 3/2

Upsampling: 
insert L-1 zeros into the input signal
Period becomes T/L
Signal size becomes 3*Fx.d0

Low Pass Filter info:
->  L > M so our cutoff frequency is pi/L = pi/3
->  cutoff is Fc = 1837.5
->  Has a gain of L = 3
->  Sampling freq is Fs*(3/2) = 11025 *(3/2) = 16,537.5 Hz

execution time comparison:
Standard: 0.07800 sec OR 7.8ms
Polyphase: 0.015 sec OR 1.5ms

https://en.wikipedia.org/wiki/Downsampling_(signal_processing)
https://en.wikipedia.org/wiki/Upsampling
*/

#define L 3
#define M 2

typedef struct 
{ 
 int ndim; //number of dimensions 
 int nchan; //number of channels 
 int d0;  //length of the first dimension 
 int d1;  //length of second dimension or sample rate if audio 
 int d2;  //length of third dimension 
} dsp_file_header; 

typedef struct
{
    float* r0;
    float* r1;
    float* r2;
} r_filters;

/**
 * @brief Upsampler
 * @param l L - our sampling constant
 * @param x is our input
 * @param xl is our output
 * @param Lzl is the length of our output (with padding)
 * @param Lx is our input length
 * @param Lh is our impulse length
 */
float* upsample(int l, float* x, float* xl, int Lzl, int Lx, int Lh){
    printf("Start Upsampling\n");
    int j = 0;
    int i;
    if(&Lh == NULL){ //in case of no zero padded input
        for(i = 0; i < (Lzl); i++){
            if(i%(L) == 0){
                // printf("i = %d, j = %d\n",i,j);
                xl[i] = x[j];
                j++;
            }
        }
    }
    else{ //zero padded input
        for(i = Lh-1; i < (Lzl-Lh-1); i++){
        // for(int i = Lh-1, j = 0; i < (Lzl-Lh-1), j < Lx; i++, j++){ 
            // printf("i = %d\tj = %d\n",i,j);
            if(i%(L) == 0){
                // printf("i = %d, j = %d\n",i,j);
                xl[i] = x[j];
                j++;
            }
        }
    }
    // printf("%d\n",sizeof(xl)); //i = 144175      j = 47899
    return xl;
}

/**
 * @brief Downsampler
 * @param m M - the downsampler constant
 * @param x is the input
 * @param y is the output
 * @param Lx is the length of the input
 * @param Ly is the length of our output
 */
float* downsample(int m, float* x, float* y, int Lx, int Ly, int Lh){
    printf("Start Downsampling\n");
    int j = 0;
    if(&Lh == NULL){ //in case of no zero padded input
        for(int i = 0; i < Lx; i++){
            if(i%(m) == 0){
                // printf("i = %d, j = %d\n",i,j);
                y[j] = x[i];
                j++;
            }
        }
    }
    else{ //zero padded input
        for(int i = Lh-1; i < (Lx-Lh-1); i++){
            if(i%(m) == 0){
                // printf("i = %d, j = %d\n",i,j);
                y[j] = x[i];
                j++;
            }
        }
    }
    return y;
}

/**
 * @brief Convolution function
 * @param x is input
 * @param h is what x gets convovled with
 * @param y is the output
 * @param Lx is len of input
 * @param Lh is len of what gets convolved
 * @param ly is the zero padded output length
 */
float* conv(float* x, float* h, float* y, int Lx, int Lh, int Ly){
    printf("Start Convolution\n");
    int j;
    for (int i = 0; i < Ly; i++) { 
        for (j = 0; j < Lh; j++) { 
            y[i] += h[j] * x[i + j]; //multiply and accumulate (MAC) 
        } 
    }
    return y;
}

/**
 * @brief Delay function. Delays the whole input by delay_val
 * @param input is input
 * @param output is output 
 * @param in_length is input length
 * @param out_length is output length
 * @param delay_val is how much the input gets delayed by
 */
float* shift(float* input, float* output, int in_length, int out_length, int delay_val){
    for(int i = 0; i < in_length; i++){
        if(i+delay_val < out_length && i+delay_val > 0){
            output[i+delay_val] = input[i]; //get a move on!
        }
    }
    return output;
}

/**
 * @brief fetches r functions for polyphase
 * @param m M
 * @param r0 r0 input
 * @param r1 r1 input
 * @param r2 r2 input
 * @param h impulse
 * @param Lh impulse length
 */
float* fetch_r_into_3(float* r0, float* r1, float* r2, float* h, int Lh){
    int j = 0, k = 0, n = 0;
    for(int i = 0; i < Lh; i++){
        if(i%3 == 0){ //n = 0, 3, 6, ...
            r2[j] = h[i];
            j++;
        }
        else if(i%3 == 1){ //n = 1, 4, 7, ...
            r1[k] = h[i];
            k++;
        }
        else if(i%3 == 2){ //n = 2, 5, 8, ...
            r0[n] = h[i];
            n++;
        }
    }
    return r0,r1,r2;
}
/**
 * @brief Splits R into two Rs
 * @param input input to be split
 * @param r0 first output
 * @param r1 second output
 * @param length length of the input
 */
float* split_r_in_2(float* input, float* r0, float* r1, int length){
    int j = 0;
    int k = 0;
    for(int i = 0; i < length; i++){
        if(i%2 == 0){ //n = 0, 2, 4, ...
            r1[j] = input[i];
            j++;
        }
        else if(i%2 == 1){ //n = 1, 3, 5, ...
            r0[k] = input[i];
            k++;
        }
    }
    return r0, r1;
}

void main(int argc, char** argv){
    //read in file
    FILE* fx;
    FILE* fy;
    if (NULL == (fx = fopen(argv[1], "rb"))) { //error check and open file
        printf("error: Cannot open input file.\n");
        return;
    }
    if (NULL == (fy = fopen(argv[2], "wb"))) { //error check and open file
        printf("error: Cannot open output file for writing.\n");
        return;
    }
    //grab headers of each file 
    dsp_file_header h0, h1, ho; 
    fread(&h0, sizeof(dsp_file_header), 1, fx); 
    memcpy(&ho, &h0, sizeof(dsp_file_header));
    int a = L;
    int b = M;
    float fs = (float)h0.d1;
    float fs_out = ((float)a/(float)b)*fs;
    printf("%d, %d\n",a,b);
    printf("%f\n",fs_out);
    ho.d1 = (int)16537.5;
    fwrite(&ho, sizeof(dsp_file_header), 1, fy); 
    printf("ndim = %d, nchan = %d, d0 = %d, d1 = %d, d2 = %d\n", h0.ndim, h0.nchan, h0.d0, h0.d1, h0.d2);

    int Lh = sizeof(h)/sizeof(h[0]); //coefficients from filter.h
    
    // int Lx = h0.d0;
    // int Lx = sizeof(x_16)/sizeof(x_16[0]);
    // int Lx = sizeof(x_8)/sizeof(x_8[0]);
    int Lx = sizeof(x_4)/sizeof(x_4[0]);
    int Lxl = h0.d0 * L;
    int Lz = Lxl + 2*(Lh-1); //length output for upsampled
    int Lv = Lxl + (Lh-1); //length of filtered upsampled signal 
    int Ly = Lv / M; //final output length
    // printf("Lh = %d, Lx = %d, Ly = %d, Lv = %d, Lzl = %d\n", Lh, Lx, Ly, Lv, Lz); 
    // float* x = calloc(sizeof(float), Lx); //allocate data for file store 
    // float* v = calloc(sizeof(float), Lz);
    // float* y = calloc(sizeof(float), Ly); 
    // float* xl = calloc(sizeof(float), Lz); 
    
    // while (!feof(fx)) { 
    //     fread(x, sizeof(float), Lx, fx);
    // } 
    printf("Create r statements\n");
    int Lrk = Lh/L;
    int Lrkk = Lrk/M;
    float* r0 = calloc(sizeof(float), Lrk);
    float* r1 = calloc(sizeof(float), Lrk);
    float* r2 = calloc(sizeof(float), Lrk);
    r0,r1,r2 = fetch_r_into_3(r0, r1, r2, h, Lh);
    float* r00 = calloc(sizeof(float), Lrkk);
    float* r01 = calloc(sizeof(float), Lrkk);
    printf("SPLIT r statements\n");
    r00, r01 = split_r_in_2(r0, r00, r01, Lrk);
    float* r10 = calloc(sizeof(float), Lrkk);
    float* r11 = calloc(sizeof(float), Lrkk);
    printf("SPLIT r statements\n");
    r10, r11 = split_r_in_2(r1, r10, r11, Lrk);
    float* r20 = calloc(sizeof(float), Lrkk);
    float* r21 = calloc(sizeof(float), Lrkk);
    printf("SPLIT r statements\n");
    r20, r21 = split_r_in_2(r2, r20, r21, Lrk);
    free(r0);
    free(r1);
    free(r2);

    int Ld_pad = (Lx/2)+2*(Lh-1);
    int Lvk = (Lx/2)+(Lh-1); 
    clock_t begin = clock(); 
    //1.)
    float* x0_1 = calloc(sizeof(float), Lx);
    float* x0_2 = calloc(sizeof(float), Lx);
    //x0_1 = shift(x, x0_1, Lx, Lx, -2);
    // x0_1 = shift(x_16, x0_1, Lx, Lx, -2);
    // x0_1 = shift(x_8, x0_1, Lx, Lx, -2);
    x0_1 = shift(x_4, x0_1, Lx, Lx, -2);
    x0_2 = shift(x0_1, x0_2, Lx, Lx, 1);
    float* xd0_1 = calloc(sizeof(float), Ld_pad);
    float* xd0_2 = calloc(sizeof(float), Ld_pad);
    xd0_1 = downsample(M, x0_1, xd0_1, Lx, Ld_pad, Lh);
    xd0_2 = downsample(M, x0_2, xd0_2, Lx, Ld_pad, Lh);
    float* v0_1 = calloc(sizeof(float), Lvk);
    float* v0_2 = calloc(sizeof(float), Lvk);
    v0_1 = conv(xd0_1, r00, v0_1, Ld_pad, Lrkk, Lvk);
    v0_2 = conv(xd0_2, r01, v0_2, Ld_pad, Lrkk, Lvk);
    float* v_0 = calloc(sizeof(float), Lvk);
    for(int i = 0; i < Lvk; i++){
        v_0[i] = v0_1[i] + v0_2[i];
    }
    float* yu_0 = calloc(sizeof(float), Lvk*3);
    yu_0 = upsample(L, v_0, yu_0, 3*Lvk, Lvk, NULL);
    float* y_0 = calloc(sizeof(float), 3*Lvk);
    y_0 = shift(yu_0, y_0, 3*Lvk, 3*Lvk, 1);
    free(x0_1);
    free(x0_2);
    free(xd0_1);
    free(xd0_2);
    free(v0_1);
    free(v0_2);
    free(v_0);
    free(yu_0);
    free(r00);
    free(r01);

    //2.)
    float* x1_1 = calloc(sizeof(float), Lx);
    float* x1_2 = calloc(sizeof(float), Lx);
    // x1_1 = shift(x, x1_1, Lx, Lx, -2);
    // x1_1 = shift(x_16, x1_1, Lx, Lx, -2);
    // x1_1 = shift(x_8, x1_1, Lx, Lx, -2);
    x1_1 = shift(x_4, x1_1, Lx, Lx, -2);
    x1_2 = shift(x1_1, x1_2, Lx, Lx, 1);
    float* xd1_1 = calloc(sizeof(float), Ld_pad);
    float* xd1_2 = calloc(sizeof(float), Ld_pad);
    xd1_1 = downsample(M, x1_1, xd1_1, Lx, Ld_pad, Lh);
    xd1_2 = downsample(M, x1_2, xd1_2, Lx, Ld_pad, Lh);
    float* v1_1 = calloc(sizeof(float), Lvk);
    float* v1_2 = calloc(sizeof(float), Lvk);
    v1_1 = conv(xd1_1, r10, v1_1, Ld_pad, Lrkk, Lvk);
    v1_2 = conv(xd1_2, r11, v1_2, Ld_pad, Lrkk, Lvk);
    float* v_1 = calloc(sizeof(float), Lvk);
    for(int i = 0; i < Lvk; i++){
        v_1[i] = v1_1[i] + v1_2[i];
    }
    float* y_1 = calloc(sizeof(float), Lvk*3);
    y_1 = upsample(L, v_1, y_1, 3*Lvk, Lvk, NULL);
    free(x1_1);
    free(x1_2);
    free(xd1_1);
    free(xd1_2);
    free(v1_1);
    free(v1_2);
    free(v_1);
    free(r10);
    free(r11);
    


    // //3.)
    float* x2_1 = calloc(sizeof(float), Lx);
    float* x2_2 = calloc(sizeof(float), Lx);
    // x2_1 = shift(x, x2_1, Lx, Lx, -2);
    // x2_1 = shift(x_16, x2_1, Lx, Lx, -2);
    // x2_1 = shift(x_8, x2_1, Lx, Lx, -2);
    x2_1 = shift(x_4, x2_1, Lx, Lx, -2);
    x2_2 = shift(x2_1, x2_2, Lx, Lx, 1);
    float* xd2_1 = calloc(sizeof(float), Ld_pad);
    float* xd2_2 = calloc(sizeof(float), Ld_pad);
    xd2_1 = downsample(M, x2_1, xd2_1, Lx, Ld_pad, Lh);
    xd2_2 = downsample(M, x2_2, xd2_2, Lx, Ld_pad, Lh);
    float* v2_1 = calloc(sizeof(float), Lvk);
    float* v2_2 = calloc(sizeof(float), Lvk);
    v2_1 = conv(xd2_1, r20, v2_1, Ld_pad, Lrkk, Lvk);
    v2_2 = conv(xd2_2, r21, v2_2, Ld_pad, Lrkk, Lvk);
    float* v_2 = calloc(sizeof(float), Lvk);
    for(int i = 0; i < Lvk; i++){
        v_2[i] = v2_1[i] + v2_2[i];
    }
    float* y_2 = calloc(sizeof(float), Lvk*3);
    y_2 = upsample(L, v_2, y_2, 3*Lvk, Lvk, NULL);
    free(x2_1);
    free(x2_2);
    free(xd2_1);
    free(xd2_2);
    free(v2_1);
    free(v2_2);
    free(v_2);
    free(r20);
    free(r21);


    // //finish section
    float* y_01 = calloc(sizeof(float), 3*Lvk);
    for(int i = 0; i < 3*Lvk; i++){
        y_01[i] = y_0[i] + y_1[i];
    }
    float* y01 = calloc(sizeof(float), 3*Lvk);
    y01 = shift(y_01, y01, 3*Lvk, 3*Lvk, 1);
    float* y = calloc(sizeof(float), 3*Lvk);
    for(int i = 0; i < 3*Lvk; i++){
        y[i] = y_2[i] + y01[i];
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
    printf("Time of standard form = %f\n", time_spent);
    fwrite(y, sizeof(float), 3*Lvk, fy);
    printf("Lh = %d\n", Lh);
    // //output to file
    fclose(fx); 
    fclose(fy); 
}