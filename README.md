# SemiGlobalMathingImplementation
This is matlab implementation of disparity map generation from stereo images with semi global matching algorithm.

'Heiko Hirschmuller. Stereo processing by semiglobal matching and mutual information. 
Pattern Analysis and Machine Intelligence, IEEE Transactions on, 30(2):328ñ341, 2008'

The main function is How2RunSGMWrapper.m in which all the parameters are setts.
SGMWrapper.m is called from How2RunSGMWrapper.m and does the actual processing.

## Inputs

![LeftImage](https://github.com/kobybibas/SemiGlobalMathingImplementation/blob/master/image_left.png)
![RightImage](https://github.com/kobybibas/SemiGlobalMathingImplementation/blob/master/image_right.png)

## Output

![DisparityMap](https://github.com/kobybibas/SemiGlobalMathingImplementation/blob/master/disparity_map.tiff)
