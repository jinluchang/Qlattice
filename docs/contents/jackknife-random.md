# A Jackknife-bootstrap hybrid resampling method

> Luchang Jin
> 2025/12/29

$$
\def\ba#1\ea{\begin{align}#1\end{align}}
\def\nn{\nonumber}
\def\ra{\rangle}
\def\la{\langle}
\def\bra{\big\rangle}
\def\bla{\big\langle}
\def\Bra{\Big\rangle}
\def\Bla{\Big\langle}
\def\ud{\mathrm{d}}
\nn
$$

The method is proposed in the following paper by Chien-Fu Jeff Wu (吳建福).

https://projecteuclid.org/journals/annals-of-statistics/volume-14/issue-4/Jackknife-Bootstrap-and-Other-Resampling-Methods-in-Regression-Analysis/10.1214/aos/1176350142.full

```
@article{10.1214/aos/1176350142,
author = {C. F. J. Wu},
title = {{Jackknife, Bootstrap and Other Resampling Methods in Regression Analysis}},
volume = {14},
journal = {The Annals of Statistics},
number = {4},
publisher = {Institute of Mathematical Statistics},
pages = {1261 -- 1295},
keywords = {$M$-regression, balanced residuals, bias reduction, bias-robustness, bootstrap, Fieller's linterval, generalized linear models, jackknife percentile, Linear regression, Nonlinear regression, representation of the least squares estimator, variable jackknife, Weighted jackknife},
year = {1986},
doi = {10.1214/aos/1176350142},
URL = {https://doi.org/10.1214/aos/1176350142}
}
```

Below, we concisely describe our implementation of the method in the context of lattice QCD calculations. Let $C_j$ be the initial data, and $j$ is the index of the configuration. For example, $C_j$ can be a correlation function measured on the configuration $j$. For a particular $j$, $C_j$ can be one number or a set of numbers. The average of the data is:

$$
\ba
C_\text{avg} = \frac{1}{N} \sum_{j} C_j
,
\\
\ea
$$

where the summation ranges over all available configurations, $N$ is the total number of available configurations for this observable $C$.

We intend to define the Jackknife-bootstrap hybrid (J-B hybrid) samples to fluctuate around $C_\text{avg}$ similar to how $C_\text{avg}$ fluctuate around the true expectation value of $C$. The total number of J-B hybrid samples is $N_\text{rs}$. Similar to the standard bootstrap procedure, this number is adjustable. The definition of the J-B hybrid samples is

$$
\ba
\overline{C}_{i}
=
C_\text{avg}
+
\sum_{j}
\frac{r_{i,j}}{\sqrt{N(N-1)}}
(C_j - C_\text{avg})
,
\\
\ea
$$

where $i$ is the resampling sample index that ranges from $1$ to $N_\text{rs}$. The random weights $r_{i,j}$ follow the standard normal distribution with

$$
\ba
\mathrm{E}(r_{i,j}) =& 0
,
\\
\mathrm{E}(r_{i,j}^2) =& 1
.
\\
\ea
$$

The random numbers $r_{i,j}$ with different $i$ or $j$ indices are statistically independent. Note that the $j$ index should uniquely label the configuration, including both the ID for the ensemble and the trajectory number of the configuration within the ensemble.

After the J-B hybrid samples are obtained, we can calculate the estimation of the central value and the statistical error of observable $O$.

$$
\ba
O_\text{avg} =& O(C_\text{avg})
,
\\
O_\text{err} =& \sqrt{\frac{1}{N_\text{rs}}\sum_{i=1}^{N_\text{rs}} (O(\overline{C}_i) - O_\text{avg})^2}
.
\\
\ea
$$

## Blocking

To deal with possible correlation between the data from different configurations, we need to introduce blocking. In the J-B hybrid resampling method, we implement the blocking procedure as follow.

We introduce the blocking function acting on the J-B hybrid sample index $i$ and the configuration index $j$

$$
\ba
b(i,j).
\ea
$$

The function should return unique label for the block that the configuration $j$ belongs to. Note that the blocking schemes can be different for different J-B hybrid sample index ($i$). Typically, we should keep the blocking size the same. However, we may choose different the blocking boundaries for different J-B hybrid samples. The number of configurations within a block is denoted as

$$
\ba
N_{b(i,j)}.
\ea
$$

With blocking, the definition of the average is the same as before,

$$
\ba
C_\text{avg} = \frac{1}{N} \sum_{j} C_j
.
\\
\ea
$$

The definition of the J-B hybrid samples is slightly altered as

$$
\ba
\overline{C}_{i}
=
C_\text{avg}
+
\sum_{j}
\frac{r_{i,b(i,j)}}{\sqrt{N(N-N_{b(i,j)})}}
(C_j - C_\text{avg})
.
\\
\ea
$$

Note that the random weights $r_{i,b(i,j)}$ depends on the label of the block ($b(i,j)$), instead of the index of the configuration ($j$).

The remaining procedures are the same as before.
