# Soft Quantization using Entropic Regularization
The quantization problem aims to find the best possible approximation of probability measures on R𝑑
using finite, discrete measures. The Wasserstein distance is a typical choice to measure the quality of the
approximation.

This contribution investigates the properties and robustness of the entropy-regularized quantization
problem, which relaxes the standard quantization problem. The proposed approximation technique
naturally adopts the softmin function, which is well known for its robustness in terms of theoretical and
practicability standpoints. Moreover, we use the entropy-regularized Wasserstein distance to evaluate
the quality of the soft quantization problem’s approximation, and we implement a stochastic gradient
approach to achieve the optimal solutions. The control parameter in our proposed method allows for the
adjustment of the optimization problem’s difficulty level, providing significant advantages when dealing
with exceptionally challenging problems of interest. As well, this contribution empirically illustrates the
performance of the method in various expositions.


When you are using this code, please cite the paper.

<a id="1">[1]</a> Rajmadan Lakshmanan and Alois Pichler. (2023). [Soft Quantization using Entropic Regularization](https://www.mdpi.com/1099-4300/25/10/1435). 

This paper also comprehensively explains the Soft Quantization using Entropic Regularization (Not yet updated).


## Directory structure

| File/Folder   | Purpose                                                                                   |
| ------------- |-------------------------------------------------------------------------------------------|   
| src           | Soft Quantization algorithm from Section 4.1.1 of [[1]](#1) |
| Graphs and images        |  Graphs and images of the numerical experiments.               |


