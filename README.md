# Tensor robust PCA (TRPCA) and tensor completion based on tensor nuclear norm under linear transform

### Introduction

In [1], we propose a new tensor nuclear norm induced by the tensor-tensor product (t-product) [2] and apply it to tensor robust PCA (TRPCA) with exact recovery guarantee in theory. We also apply this tensor nuclear norm for tensor completion and provide the theoretical recovery guarantee in [3]. The original t-product [2] uses the discrete Fourier transform and uses the fast Fourier transform (FFT) for efficient computing. It is further generlaized to the t-product under arbitrary invertible linear transform in [4]. 

In [5], we show that if the linear transform satisfies <a href="https://www.codecogs.com/eqnedit.php?latex=L^\top&space;L&space;=&space;LL^\top=\ell&space;I" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L^\top&space;L&space;=&space;LL^\top=\ell&space;I" title="L^\top L = LL^\top=\ell I" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=\ell>0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\ell>0" title="\ell>0" /></a>, then we can define a more general tensor nuclear norm induced by the t-product under this linear transform. We then apply this generalized tensor nuclear norm for tensor completion [5] and tensor robust PCA (TRPCA) [6] under linear transform and provide the exact recovery guarantee in theory. Intuitively, our works, [5] and [6], can be regarded as extensions of our former work [3] and [1], respectively. Howver, we restrict the linear transform <a href="https://www.codecogs.com/eqnedit.php?latex=L" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L" title="L" /></a> to be a real matrix in [5][6] whie [3][1] uses the complex discrete Fourier transform. This makes some key differences in the proofs of the recovery guarantee.

We provide the codes for the following models:
<ol>    
<li><b> Tensor completion based on tensor nuclear norm under linear transform [5]</b><br/>
    
<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{\mathcal{X}}&space;\|\mathcal{X}\|_*,&space;\&space;\text{s.t.}&space;\&space;\mathcal{P}_{\Omega}(\mathcal{X})=\mathcal{P}_{\Omega}(\mathcal{M})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{\mathcal{X}}&space;\|\mathcal{X}\|_*,&space;\&space;\text{s.t.}&space;\&space;\mathcal{P}_{\Omega}(\mathcal{X})=\mathcal{P}_{\Omega}(\mathcal{M})" title="\min_{\mathcal{X}} \|\mathcal{X}\|_*, \ \text{s.t.} \ \mathcal{P}_{\Omega}(\mathcal{X})=\mathcal{P}_{\Omega}(\mathcal{M})" /></a>

  
<li><b> Tensor robust PCA (TRPCA) based on tensor nuclear norm under linear transform [6]</b><br/>
  
<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{\mathcal{L},&space;\mathcal{E}}&space;\|\mathcal{L}\|_*&plus;\lambda\|\mathcal{E}\|_1,&space;\text{&space;s.t.&space;}&space;\mathcal{X}&space;=&space;\mathcal{L}&plus;\mathcal{E}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{\mathcal{L},&space;\mathcal{E}}&space;\|\mathcal{L}\|_*&plus;\lambda\|\mathcal{E}\|_1,&space;\text{&space;s.t.&space;}&space;\mathcal{X}&space;=&space;\mathcal{L}&plus;\mathcal{E}" title="\min_{\mathcal{L}, \mathcal{E}} \|\mathcal{L}\|_*+\lambda\|\mathcal{E}\|_1, \text{ s.t. } \mathcal{X} = \mathcal{L}+\mathcal{E}" /></a>
</ol>


### Related Toolboxes
<ul>
  <li> <a href="https://github.com/canyilu/tproduct" class="textlink">Tensor-Tensor Product Toolbox</a></li>       
  <li> <a href="https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA" class="textlink">Tensor robust principal component analysis </a></li>  
  <li> <a href="https://github.com/canyilu/tensor-completion-tensor-recovery" class="textlink">Tensor completion by tensor nuclear norm minimization</a></li>
  <li> <a href="https://github.com/canyilu/LibADMM" class="textlink">A Library of ADMM for Sparse and Low-rank Optimization </a></li>
</ul>


### References
<ol>
<li> Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin, Shuicheng Yan. Tensor Robust Principal Component Analysis with A New Tensor Nuclear Norm. TPAMI. 2019

<li> Misha E.Kilmer, Carla D. Martin. Factorization strategies for third-order tensors. Linear Algebra and its Applications, 435(3):641–658, 2011

<li> Canyi Lu, Jiashi Feng, Zhouchen Lin, Shuicheng Yan. Exact Low Tubal Rank Tensor Recovery from Gaussian Measurements. International Joint Conference on Artificial Intelligence (IJCAI). 2018
  
<li> Eric Kernfeld, Misha Kilmer, Shuchin Aeron. Tensor-Tensor Products with Invertible Linear Transforms. Linear Algebra and its Applications. 2015.

<li> Canyi Lu, Xi Peng, Yunchao Wei. Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible Linear Transforms. IEEE International Conference on Computer Vision and Pattern Recognition (CVPR), 2019

<li> Canyi Lu, Pan Zhou. Exact Recovery of Tensor Robust Principal Component Analysis under Linear Transforms. arXiv preprint arXiv:1907.08288. 2019

<li> Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin, Shuicheng Yan. Tensor Robust Principal Component Analysis: Exact Recovery of Corrupted Low-Rank Tensors via Convex Optimization. CVPR, 2016
</ol>
