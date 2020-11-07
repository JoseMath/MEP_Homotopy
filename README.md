This repository contains code for the paper "Fiber product homotopy method for multiparameter eigenvalue problems".


Our homotopy method tracks several copies of the eigenvalue λ; they should all converge to the same value if our method performs correctly.

The λi1 coordinates of all eight paths for t ∈ [0.95, 1] and i = 1, 2, 3. | The λi1 coordinates of the path highlighted in the left panel for t ∈ [0, 1] and i = 1, 2, 3.
:-------------------------:|:-------------------------:
![1](https://raw.githubusercontent.com/JoseMath/MEP_Homotopy/master/img/eightEndPaths.png)  |  ![](https://raw.githubusercontent.com/JoseMath/MEP_Homotopy/master/img/path1.png)



- `MATLAB`, this folder contains our MATLAB implementation of the fiber product homotopy method, including the parallelization version. To compare with the Delta method, you will need to install the MATLAB package `MultiEigPar` first.
- `Bertini`, this folder contains Macaulay2 implementation for running both the fiber product homotopy method and the diagonal coefficient homotopy method on Bertini.


All our experiments are done on [RCC Midway](https://rcc.uchicago.edu/) with 28 cpu cores and 56G memory.
