# lingamr: LiNGAM Algorithms for Causal Discovery

R implementation of LiNGAM (Linear Non-Gaussian Acyclic Model)
algorithms for causal discovery, following Shimizu et al. (2011)
<https://www.jmlr.org/papers/v12/shimizu11a.html>. Based on the 'Python'
implementation by Ikeuchi et al. (2023)
<https://github.com/cdt15/lingam>. The VAR-LiNGAM residual diagnostics
are inspired by the 'VARLiNGAM' R code of Moneta et al.
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>.

## Details

tutorial:
[`vignette("lingamr", package = "lingamr")`](https://morimotoosamu.github.io/lingamr/articles/lingamr.md)

## References

Shimizu, S., et al. (2011). DirectLiNGAM: A direct method for learning a
linear non-Gaussian structural equation model. *Journal of Machine
Learning Research*, 12, 1225-1248.

Hyvärinen, A., Zhang, K., Shimizu, S., & Hoyer, P. O. (2010). Estimation
of a structural vector autoregression model using non-Gaussianity.
*Journal of Machine Learning Research*, 11, 1709-1731.

Moneta, A., Entner, D., Hoyer, P. O., & Coad, A. (2013). Causal
inference by independent component analysis: Theory and applications.
*Oxford Bulletin of Economics and Statistics*, 75(5), 705-730. VARLiNGAM
R code:
<https://sites.google.com/site/dorisentner/publications/VARLiNGAM>

Shimizu, S. (2012). Joint estimation of linear non-Gaussian acyclic
models. *Neurocomputing*, 81, 104-107. (MultiGroup Direct LiNGAM)

Tashiro, T., Shimizu, S., Hyvärinen, A., & Washio, T. (2014).
ParceLiNGAM: A causal ordering method robust against latent confounders.
*Neural Computation*, 26(1), 57-83.

Maeda, T. N., & Shimizu, S. (2020). RCD: Repetitive causal discovery of
linear non-Gaussian acyclic models with latent confounders. *AISTATS
2020*, PMLR 108, 735-745.

Wang, Y. S., & Drton, M. (2020). High-dimensional causal discovery under
non-Gaussianity. *Biometrika*, 107(1), 41-59.

Zeng, Y., Shimizu, S., Matsui, H., & Sun, F. (2022). Causal discovery
for linear mixed data. *Proceedings of the First Conference on Causal
Learning and Reasoning (CLeaR 2022)*, PMLR 177, 994-1009.

Python implementation (DirectLiNGAM, VAR-LiNGAM, MultiGroup,
ParceLiNGAM, RCD, LiM, HighDim): cdt15/lingam,
<https://github.com/cdt15/lingam>

## See also

Useful links:

- <https://github.com/morimotoosamu/lingamr>

- <https://morimotoosamu.github.io/lingamr/>

- Report bugs at <https://github.com/morimotoosamu/lingamr/issues>

## Author

**Maintainer**: Osamu Morimoto <galactic.supermarket@gmail.com>
\[copyright holder\]

Authors:

- Osamu Morimoto <galactic.supermarket@gmail.com> \[copyright holder\]

Other contributors:

- T. Ikeuchi \[copyright holder\]

- G. Haraoka \[copyright holder\]

- M. Ide \[copyright holder\]

- W. Kurebayashi \[copyright holder\]

- S. Shimizu \[copyright holder\]
