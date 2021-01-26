# Image-Blending-Using-Pyramids

Implemented the Gaussian and Laplacian Pyramid as described by (http://persci.mit.edu/pub_pdfs/spline83.pdf). The output pyramid is the list of smaller and smaller images. Also, implemented a method to reconstruct the original image from the Laplacian Pyramid as well as did soft and hard thresholding to obtain a denoised image from the Laplacian Pyramid. 

Now, Image pyramids can be used to seamlessly blend two images. To do so, you will need three images of the same size: the two images a and b to be blended, and a region image r defining which regions the two images should contribute to. Construct the Laplacian pyramids of a and b, blend their corresponding levels using the Gaussian pyramid of r, and reconstruct the blended image. This method is implemented for greyscale as well as colored images
