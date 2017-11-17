# video-quality-metrics
This program is dedicated to measuring such objective video quality metrics as PSNR and mean SSIM (MSSIM). It's main feature is a time-efficient implementation optimized with CPU vector instructions and multicore parallelism. The program has an option to print metric values for all frames in a sequential order as well as average results for entire video. The input videos must be in the \*.y4m raw format.

Both quality metrics are calculated only for luminance component of the video frames assuming the original videos are in the YC<sub>b</sub>C<sub>r</sub> colour space.

### PSNR
The PSNR metric is calculated as:

![PSNR formula](https://latex.codecogs.com/gif.latex?PSNR%20%3D%2010%20%5Ccdot%20%5Ctextrm%7Blog%7D_%7B10%7D%20%5Cleft%28%20%5Cfrac%7BL%5E2%7D%7B%20%5Cdisplaystyle%7B%5Cfrac%7B1%7D%7BN%7D%20%5Ccdot%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20%5Cleft%28x_i%20-%20y_i%20%5Cright%29%5E2%7D%20%7D%20%5Cright%29%2C)

where _N_ is a number of pixels in a frame; _L_ = 255 is a dynamic range of pixel brightness; x<sub>i</sub> and y<sub>i</sub> are values of the corresponding frame pixels.

### Mean SSIM
MSSIM is implemented as an arithmetical average of the SSIM metric in 11x11 pixel windows according to the original paper:
> Zhou Wang, A. C. Bovik, H. R. Sheikh and E. P. Simoncelli, 
> "Image quality assessment: from error visibility to structural similarity," 
> in IEEE Transactions on Image Processing, vol. 13, no. 4, pp. 600-612, April 2004.

The SSIM in the 11x11 window is calculated as follows:

![SSIM formula](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20%26%20SSIM%20%3D%20%5Cfrac%7B%282%5Cmu_x%5Cmu_y%20&plus;%20C_1%29%282%5Csigma_%7Bxy%7D%20&plus;%20C_2%29%7D%7B%28%5Cmu_x%5E2%20&plus;%20%5Cmu_y%5E2%20&plus;%20C_1%29%28%5Csigma_x%5E2%20&plus;%20%5Csigma_y%5E2%20&plus;%20C_2%29%7D%20%5C%5C%20%26%20%5Cmu_x%20%3D%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20w_i%20x_i%20%5C%5C%20%26%20%5Csigma_x%20%3D%20%5Csqrt%7B%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20w_i%20%28x_i%20-%20y_i%29%5E2%7D%20%5C%5C%20%26%20%5Csigma_%7Bxy%7D%20%3D%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20w_i%20%28x_i%20-%20%5Cmu_x%29%20%28y_i%20-%20%5Cmu_y%29%2C%20%5Cend%7Balign*%7D)

where _N_ = 121 is the number of pixels in the window; C<sub>1</sub> = (0.01 * L)<sup>2</sup>; C<sub>2</sub>=(0.03 * L)<sup>2</sup>; _L_ = 255 is a dynamic range; x<sub>i</sub> and y<sub>i</sub> are corresponding pixel values; w<sub>i</sub> represents a value from the 11x11 matrix of samples from the 2D Gaussian with σ = 1.5:

| -5 | -4 | -3 | -2 | -1 | center | +1 | +2 | +3 | +4 | +5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|0.000001 | 0.000008 | 0.000037 | 0.000112 | 0.000219 | 0.000274 | 0.000219 | 0.000112 | 0.000037 | 0.000008 | 0.000001|
|0.000008 | 0.000058 | 0.000274 | 0.000831 | 0.001619 | 0.002021 | 0.001619 | 0.000831 | 0.000274 | 0.000058 | 0.000008|
|0.000037 | 0.000274 | 0.001296 | 0.003937 | 0.007668 | 0.009577 | 0.007668 | 0.003937 | 0.001296 | 0.000274 | 0.000037|
|0.000112 | 0.000831 | 0.003937 | 0.011960 | 0.023294 | 0.029091 | 0.023294 | 0.011960 | 0.003937 | 0.000831 | 0.000112|
|0.000219 | 0.001619 | 0.007668 | 0.023294 | 0.045371 | 0.056662 | 0.045371 | 0.023294 | 0.007668 | 0.001619 | 0.000219|
|0.000274 | 0.002021 | 0.009577 | 0.029091 | 0.056662 | 0.070762 | 0.056662 | 0.029091 | 0.009577 | 0.002021 | 0.000274|
|0.000219 | 0.001619 | 0.007668 | 0.023294 | 0.045371 | 0.056662 | 0.045371 | 0.023294 | 0.007668 | 0.001619 | 0.000219|
|0.000112 | 0.000831 | 0.003937 | 0.011960 | 0.023294 | 0.029091 | 0.023294 | 0.011960 | 0.003937 | 0.000831 | 0.000112|
|0.000037 | 0.000274 | 0.001296 | 0.003937 | 0.007668 | 0.009577 | 0.007668 | 0.003937 | 0.001296 | 0.000274 | 0.000037|
|0.000008 | 0.000058 | 0.000274 | 0.000831 | 0.001619 | 0.002021 | 0.001619 | 0.000831 | 0.000274 | 0.000058 | 0.000008|
|0.000001 | 0.000008 | 0.000037 | 0.000112 | 0.000219 | 0.000274 | 0.000219 | 0.000112 | 0.000037 | 0.000008 | 0.000001|

### Compiling and running
The program uses Intel AVX/AVX2 instructions and therefore requires an appropriate CPU to run.

It also needs Qt library for parsing command line arguments and time measurement. In order to build the application use Qt Creator or type in Linux command line: `qmake video-quality-metrics.pro && make`.

Run application using the following syntax:
```
video-quality-metrics <video_file_1> <video_file_2> [-nopsnr|-nomssim] [-noprogress] [-notime] [-print-per-frame] [-threads <int>]
```

### Example output
```
y4m-psnr-mssim-fp32 original.y4m decoded.y4m -nopsnr -notime -print-per-frame -threads 6
Frame  Y-MSSIM
  1   0.9615719
  2   0.9586021
  3   0.9563300
  4   0.9555718
  5   0.9553111
[...]
average 0.9575674
```
