# focus-series
Tools to analyse the focus-series of the Fraunhofer telescope at the Wendelstein observatory.

## CLI Scripts
```bash
focus.py [-h] -t temperature -zd zenith-distance
```
Standalone script implementing the surface interpolant of the focus function.
Evaluates the spline at the given temperature and zenith distance.
Simple linear extrapolation outside the splines' convex hull.
```
tsi_mutual_info_map.py [-h] [-f file] [-n number-of-parameters] [--plot-path plot-path] [-s size]
psf_mutual_info_map.py
psf_tsi_mutual_info_map.py
```
These take files containing mutual information scores as generated by `tsi_mutual_info.py`, `psf_mutual_info.py`
and `psf_tsi_mutual_info.py` and plot a heatmap of the correlation matrix. By providing a maximum number of parameters
per axis, the plot will be divided into several files accordingly. While dividing the plot reduces the overview, this 
can be advantageous for printed display of a large parameter space.
```
analyse_wind.py [-h] [-p parameter] [--psf-path psf-path] [--meteo-path meteo-path] 
                [--plot-path plot-path] [--size-multiplier size-multiplier]
```
Uses the preprocessed meteological data as generated by `preprocess_meteo.py` and correlates the wind velocity and direction
to a given psf parameter. The data is then reduced to a polar heatmap and ploted.