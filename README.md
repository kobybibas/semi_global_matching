# SemiGlobalMathingImplementation
This is matlab implementation of disparity map generation from stereo images with semi global matching algorithm.

```
@article{hirschmuller2007stereo,
  title={Stereo processing by semiglobal matching and mutual information},
  author={Hirschmuller, Heiko},
  journal={IEEE Transactions on pattern analysis and machine intelligence},
  volume={30},
  number={2},
  pages={328--341},
  year={2007},
  publisher={IEEE}
}
```

1. The main function is How2RunSGMWrapper.m in which all the parameters are setts.
2. SGMWrapper.m is called from How2RunSGMWrapper.m and does the actual processing.

## Inputs

![LeftImage](https://github.com/kobybibas/SemiGlobalMathingImplementation/blob/master/image_left.png)
![RightImage](https://github.com/kobybibas/SemiGlobalMathingImplementation/blob/master/image_right.png)

## Output

![DisparityMap](https://github.com/kobybibas/SemiGlobalMathingImplementation/blob/master/disparity_map.tiff)
