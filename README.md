This repository contains the R code used to implement the simulation studies in the paper:

[**Non-parametric Bayesian inference via loss functions under model misspecification**](https://arxiv.org/abs/2103.04086)  by Yu Luo, David A. Stephens, Daniel J. Graham, Emma J. McCoy

The code is organized by simulation studies corresponding to the examples presented in the paper, together with additional experiments assessing the validity and calibration of the Gibbs posterior.

## Validity check
The ```Validity check/``` folder contains code for calibration experiments for checking the validity of the proposed generalized Bayesian loss-based method, corresponding to Section 3.3 in the paper.

## Calibration
The ```Calibration/``` folder contains code for calibration experiments for quantile regression in Section 3.4

## Example 1
Implements Example 1 in Section 6.1.

 + Simulates data under the Example 1 data-generating process.
 + Computes posterior distributions using the proposed Bayesian loss function.

## Example 2
Implements Example 1 in Section 6.2.

 + Simulates data under the Example 2 data-generating process.
 + Computes posterior distributions under high-dimensional confounding..


## Example 3
Implements Example 1 in Section 6.3.

 + Simulates data under the Example 3 data-generating process.
 + Compares the proposed method with flexible modeling approaches.
