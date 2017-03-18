# CPD model

Python implemetation of the Chemical Percolation Devolatilization (**CPD**) model.
The original model was developed in Fortran at BYU: http://www.et.byu.edu/~tom/devolatilization/CPD%20model.html

The model requires **PKP** as dependencies (ask for having access).

![Comparison between original CPD (dashed lines) and new implementation (solid lines)](./cpd.png)

## Validation

Check [validation](./notebook/validation.ipynb)

* Illinois coal shows differences. Recheck code!

## TODO

* Add gas concentration (done small differences)
* Add nitrogen formation
* Add drop tube mode
* Test more (done check notebooks)
* Speed-up run (working)
* Use Numba or Cython to speed up
  - Check why Illinois shows different results with Numba!
