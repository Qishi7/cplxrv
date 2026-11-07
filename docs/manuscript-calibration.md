# Manuscript calibration

The primary simulation comparison uses real `X`, proper errors (`relation = 0`),
and `a_s = b_s = 0.5` for each error variance. Specify the following `hyper`
lists explicitly; package defaults are unchanged.

| Prior | Paired complex fit (`model = "complex"`) | Each scalar fit (`model = "real"`) |
| --- | --- | --- |
| Student-t | `list(a_s=0.5, b_s=0.5, nu=1)` | `list(a_s=0.5, b_s=0.5, a_tau=1, b_tau=1)` |
| Lasso | `list(a_s=0.5, b_s=0.5, r=0.5, delta=0.5)` | `list(a_s=0.5, b_s=0.5, r=0.5, delta=0.8966996574)` |
| GDP | `list(a_s=0.5, b_s=0.5, alpha=1, eta=1)` | `list(a_s=0.5, b_s=0.5, alpha=1, eta=1.3403844863)` |

For example, using `y` and real `X` from the generic demo, the calibrated
independent scalar GDP fits use:

```r
paper_gdp <- list(a_s=0.5, b_s=0.5, alpha=1, eta=1.3403844863)
fit_R <- cplx_shrink(Re(y), X, prior="gdp", model="real", hyper=paper_gdp)
fit_I <- cplx_shrink(Im(y), X, prior="gdp", model="real", hyper=paper_gdp)
```

M-ReIm requires two independent fits. M-Scalar-Pooled instead requires an
explicit block-diagonal real design and stacked response, sharing one error
variance and, for lasso, one global rate. A single ordinary fitting call does
not construct these simulation strategies. These parameter examples do not
supply the paper's datasets, seed mapping, four-chain schedules, diagnostics,
scoring, or full ablation and fMRI workflows. For the H1 sensitivity only,
the scalar lasso uses `delta=0.5` and scalar GDP uses `eta=1`.

