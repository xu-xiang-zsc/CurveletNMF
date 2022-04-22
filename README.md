# CurveletNMF

Hyperspectral unmixing(HU) is an efficient wayto extract component information from mixed pixels in remotely sensed imagery. Nonnegative matrix factorization (NMF) based unmixing methods have been widely used due to their ability to extract endmembers (pure spectral signatures) and their corresponding fractional abundances in  simultaneous fashion. In this article, we propose a new sparse-constrained NMF method for HU purposes. Unlike most sparse regularizers, imposed on abundances (vectors
or matrix) directly, our method imposes sparse constraints on a transformed abundance domain. It is based on the assumption that sparsity, when applied on a well-designed transform domain, leads to sparser representations than those in the corresponding source domain (e.g., natural images are approximately sparse in
a wavelet domain). In this regard, we specifically explore sparsity on a curvelet transformed domain of abundances. Moreover, we consider the Chambolle–Pock algorithm to solve the involved optimization model, so as to obtain a fast and stable solution. 

Please cite:
[1] X. Xu, J. Li, S. Li, and A. Plaza. Curvelet Transform Domain-Based Sparse Nonnegative Matrix Factorization for Hyperspectral Unmixing. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2020, 13: 4908-4924.
