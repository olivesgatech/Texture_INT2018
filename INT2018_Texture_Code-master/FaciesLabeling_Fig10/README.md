

## A comparative study of texture attributes for characterizing subsurface structures in seismic volumes 

###  About: 

This code was used to generate the results in Figure 10 of the paper: 

```
Z. Long, Y. Alaudah, M. Qureshi, Y. Hu, Z. Wang, M. Alfarraj, G. AlRegib, A. Amin, M. Deriche, S. Al-Dharrab, and H. Di, â€œ A comparative study of texture attributes for characterizing subsurface structures in seismic volumes,â€ Interpretation, vol. 6, no. 4, pp. T1055-T1066, Nov. 2018.
```

If you use this code, please consider citing the paper. 

### How to run: 

The code performs facies labeling to generate results in Figure 10 of the paper. The facies data are included in the subdirectory of `facies`. 

Before running the code, you will need to download and install the VLFeat library from http://www.vlfeat.org/.

The code also uses some programs under its parent directory (i.e., `getmapping.m` and programs associated with the descriptors used such as GLCM, M-CLBP, and LRI), so you need to add those related directories to your MATLAB path. Then, all you need to do is to run `main.m`.

