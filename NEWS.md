# MMEFM 0.1.0

* Simulate grouped matrix-valued time series with heterogeneous dimensions, Gaussian or raw Student-t innovations, and optional model-component storage.
* Generate piecewise ARFIMA series with stable AR(p) dynamics and shared or series-specific fractional memory.
* Estimate global and local main effects and common components at supplied or automatically selected ranks, with Van Loan alignment and refitting by default and optional Procrustes alignment and ALS refitting.
* Select the eight loading rank classes using unperturbed eigenvalue ratios, indexwise cross-group maxima for global spectra, and full consecutive-ratio search ranges. Check combined ranks against the IC1-centered spatial dimensions without rank repair.
* Return structured model fits with print and summary methods, fitted-component reconstruction, and final residuals.
* Detect global common factors and screen participating groups using a shared circular-shift bootstrap cutoff, with optional PSOCK execution and matching serial/parallel random inputs.
