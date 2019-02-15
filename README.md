u-track-PMMS implements a stochastic smoothing framework for piecewise-stationary dynamic model in multiple particle tracking.

![alt text](https://raw.githubusercontent.com/proudot/u-track-PMMS/master/img/illu.png)


# Usage 

The approach is implemented on to of the flexible u-track framework in Matlab. Please see the scriptUsage folder for examples. 

Tested with Matlab 2013-2016b. 

# Manuscript 

One of the major challenges in multiple particle tracking is the capture of extremely heterogeneous movements of objects in crowded scenes. The presence of numerous assignment candidates in the expected range of particle motion makes the tracking ambiguous and induces false positives. Lowering the ambiguity by reducing the search range, on the other hand, is not an option, as this would increase the rate of false negatives. We propose here a piecewise-stationary motion model (PMM) for the particle transport along an iterative smoother that exploits recursive tracking in multiple rounds in forward and backward temporal directions. By fusing past and future information, our method, termed PMMS, can recover fast transitions from freely or confined diffusive to directed motions with linear time complexity. To avoid false positives, we complemented recursive tracking with a robust inline estimator of the search radius for assignment (a.k.a. gating), where past and future information are exploited using only two frames at each optimization step. We demonstrate the improvement of our technique on simulated data especially the impact of density, variation in frame to frame displacements, and motion switching probability. We evaluated our technique on the 2D particle tracking challenge dataset published by Chenouard et al. in 2014. Using high SNR to focus on motion modeling challenges, we show superior performance at high particle density. On biological applications, our algorithm allows us to quantify the extremely small percentage of motor-driven movements of fluorescent particles along microtubules in a dense field of unbound, diffusing particles. We also show with virus imaging that our algorithm can cope with a strong reduction in recording frame rate while keeping the same performance relative to methods relying on fast sampling.

Roudot, P., L. Ding, K. Jaqaman, C. Kervrann, and G. Danuser. “Piecewise-Stationary Motion Modeling and Iterative Smoothing to Track Heterogeneous Particle Motions in Dense Environments.” IEEE Transactions on Image Processing 26, no. 11 (November 2017): 5395–5410. https://doi.org/10.1109/TIP.2017.2707803.
