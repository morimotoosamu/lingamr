# Acyclicity constraint h(W) and the matrix power term used in its gradient (Yu et al. 2019 formulation, as used by the Python source)

`E = M^(d-1)` is computed by the sequential product for small `d`
(bit-identical to the Python source's loop) and by binary exponentiation
for larger `d`, where the O(d^4) sequential product starts to dominate;
the switch only reorders floating-point multiplications.

## Usage

``` r
lim_h(W, d)
```
