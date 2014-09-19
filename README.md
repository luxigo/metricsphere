Copied from http://www.cs.jhu.edu/~misha/Code/MetricSphere/

CODE DESCRIPTION
This distribution provides an implementation of our streaming multigrid Poisson solver for stitching, smoothing, and sharpening of spherical images represented by the commonly used equirectangular map:
Stitching: Given an input image consisting of composited images, and given a mask which assigns the same color ID to pixels coming from the same source image, our code outputs the stitched image whose gradients are a best-fit to the gradients of the composite image, with forced zero gradients across seam boundaries.
Smoothing/Sharpening: Given an input image, a gradient modulation term, and a pixel fidelity term, our code outputs the best-fit image whose gradients match the modulated gradients of the input image and whose pixel values match the pixel values of the original image. Formally, if F is the original image, α is the fidelity term, and β is the gradient modulation term, the output image G is the image minimizing the sum of square norms:

EXECUTABLE ARGUMENTS
StreamingSphericalSolver:
--in <input spherical image>
This string is the the name of the image file that is to be processed. (Currently supported file-types include PNG, JPEG, BMP, WDP, TIFF, and our tiled image format iGrid).
--out <output spherical image>
This string is the the name of the image file to which the output will be written.
[--labels <input mask image>]
This optional string is the name of the image file serving as the mask for stitching. (Since the values of the mask are used to determine if adjacent pixels in the composite image come from the same source, the mask should not be compressed using lossy compression. Similarly, in representing the composited pixels, be wary of using JPEG compression. Even at 100% quality, it can blur out the seams between images, so that setting the seam-crossing gradient to zero is no longer sufficient.)
[--stencilIO <input/output integral table>]
This optional string is the the name of the file from/to which the table of integrals used for computing the coefficients of the system matrices are read/written.
[--vCycles <number of v-cycles>]
This optional integer specifies the number of V-cycles that should be performed. (Default value is 1.)
[--iters <number of Gauss-Seidel iterations>]
This optional integer specifies the number of Gauss-Seidel iterations that should be performed per pass. (Default value is 5.)
[--highPrecision]
By default, the system is solved using low-precision floating point values (single precision for the arithmetics, half precision for storage). If this argument is specified both the arithmetics and the I/O are done with double precision floating point values.
[--minMGRes <coarsest resolution of the multigrid solver>]
This optional integer specifies the coarsest resolution at which Gauss-Seidel relaxation should be performed by the multigrid solver. Once the solver reaches this resolution, it solves using a conjugate-gradients solver. (Default value is 64.)
[--inCoreRes <in-core/out-of-core transition resolution>]
This optional integer specifies the transition at which the solver should switch from an out-of-core to an in-core solver. (Default value is 1024.)
[--quality <JPEG compression quality>]
This optional integer, in the range [0,100], specifies the compression quality that should be used if the output image is in JPEG format. (Default value is 100.)
[--iWeight <pixel fidelity term>]
If the system is solving the Poisson equation to perform image smoothing or sharpening, this value specifies the fidelity term α.
[--gScale <gradient modulation term>]
If the system is solving the Poisson equation to perform image smoothing or sharpening, this value specifies the gradient modulation β.
[--hdrLabels]
By default, image pixels are read in at 8 bits per pixel. When using either PNG or TIFF output, enabling this flag will force reading in of the label file at 16 bits per pixel.
[--unknownIndex <RGB mask>]
When stitching images, this triplet of integer values can be used to identify the regions of the mask corresponding to missing data. When specified, unknown pixel values are set to black and gradients across the boundaries of unknown data are maintained, to avoid a harmonic fill-in.
[--verbose]
If this optional argument is specified, the solver will output the magnitude of the residuals at the different levels.
[--progress]
If this optional argument is specified, the solver will show the progress of the solver through the restriction and prolongation phases.
DETAILS
There are several important points to keep in mind regarding the implementation of the code:
The multigrid solver assumes that the width and height of the input are powers-of-two and that the width is twice the height. If the input image(s) do not have this property, the executable performs on-the-fly upsampling to the next power-of-two image, and then down-samples the solution before outputting the image to a file.
The definition of the solver requires the computation of several integrals at high precision. This calculation is costly but is independent of the image. To this end, the --stencilIO parameter supports the input/output of the integral table. If this parameter is specified, the executable will check if an integral table at sufficiently high resolution exists at the location described by the argument of --stencilIO, and only if it does not will the integrals be computed. Then, if the the table was computed (i.e. either the table did not exist or it was not of sufficiently high resolution) the integral table will be written out to the argument of the --stencilIO. Note that due to the multigrid structure, the low-resolution table is always contained within the higher-resolution table so over-writing does with the higher-res table does not lose information.
To facilitate the execution, we provide an integral table that is suitable for processing images up to a resolution of 65,536 x 32,768.

