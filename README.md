# Fast 3D Image Moments using the DPM algorithm
Accurate calculation of raw 3D image moments with O(L+N+M) multiplications for grayscale volumes. See article [here](https://arxiv.org/abs/2012.08099).

### Project
The project contains an implementation of the algorithm. For comparison purposes it also contains implementations for other approaches.

### The Discrete Projection Moment algorithm
The DPM algorithm reduces this problem from 3D moments of an LxMXM array, to 1D moments of several projection line integrals. The original volume can be projected vertically, horizontally, at other orientations, and then summed along those axes to provide a set of 2D images. The 2D images can be further projected to 1D line integrals.. The 3D moments are linear conbinations of 2D moments - and in turn the 2D moments become linear combinations of 1D moments of these line integrals and there is no loss of information. This reduces the number of multiplications from O(L.M.N) to O(L+M+N). See the timings on an Intel Core i5 7th Gen below. .

### Results


![Timings](https://github.com/wild-ig/dpm_moments/raw/main/comparison.png)