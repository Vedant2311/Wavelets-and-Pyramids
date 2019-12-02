# Image-Processing-using-Wavelets

## Haar Wavelets Denoising

Implemented the Haar wavelet transform for images, which recursively decomposes an image into approximation and detail coefficients. The output is in the form of a single 2D array containing the coefficients at all scales. Also it is necessary to ensure that the transform is orthogonal, i.e. the sum of squared values is the same before and after the transform; this requires scaling the coefficients by  √2 at each step. Also implemented the inverse Haar transform, which reconstructs the original image from the wavelet coefficients

Multiresolution representations such as pyramids and wavelets can be used for noise removal. For natural images, most of the detail coefficients are close to zero (due to smooth regions) while a few are quite large (due to edges). Therefore, suppressing low-amplitude detail coefficients while retaining high-amplitude ones is likely to reduce noise in smooth regions without significantly affecting edges and other important image features. Two common strategies to do so are hard thresholding and soft thresholding

Thus, here the haar transform of the image is obtained, then thresholded, and then it's inverse haar transform is computed to obtain the denoised image. Normally, it is observed that soft thresholding works better than hard thresholding, but still it is dependant on the user to set the type of thresholding and the thresholding value

## Haar Wavelets Compression

The fact that most detail coefficients are close to zero also suggests a simple method for lossy image compression: simply retain only the top k% wavelet coefficients, and replace the rest with zeroes. Then the wavelet representation will contain very long runs of zeroes, which can be compressed easily using run length encoding. Thus, a function is implemented to take an Image and a compression parameter **k** which writes the compressed version of the image into a file of a binary format. A decompression function is also implemented which will read the compressed file and produce a lossy reconstruction of the original image

## Lifting Scheme Denoising and Compression

Also, implemented the Lifting scheme wavelet transform (http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=E1704C91F41A9502C98F620C0C6A36FC?doi=10.1.1.47.4142&rep=rep1&type=pdf), which is a really fast method (As compared to the Haar transform) and produces the output which is nearly the same as obtained for the Haar transform. 
