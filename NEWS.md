# MMEFM 0.0.0.9000

* Implement main-effect and common-component estimation with supplied or automatically selected ranks, Van Loan or Procrustes alignment, and Van Loan or ALS refitting.
* Implement manuscript eigenvalue-ratio rank selection, with indexwise cross-group eigenvalue maxima for global ranks and full consecutive-ratio search ranges.
* Add structured `mmefm_fit` results with print, summary, fitted-component, and residual methods.
* Add global common-factor existence detection and group screening using a shared circular-shift bootstrap cutoff.
* Add cross-platform PSOCK bootstrap execution with matching serial/parallel random inputs, preserved caller RNG state, and worker cleanup.
* Add simulation with heterogeneous groups, Gaussian and raw Student-t innovations, piecewise memory, optional global ownership sets, and optional component storage.
* Correct combined spatial-rank feasibility to the IC1-centered dimensions `p_m - 1` and `q_m - 1`, allowing equality without rank repair.
