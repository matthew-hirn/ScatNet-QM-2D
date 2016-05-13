# ScatNet-QM-2D
Wavelet scattering regression of quantum chemical energies for planar molecules.

This software computes two-dimensional wavelet scattering transforms of planar molecules and corresponding molecular energy regressions, as described in the paper "**Wavelet scattering regression of quantum chemical energies**," by Matthew Hirn, Stéphane Mallat, and Nicolas Poilvert.

## Installation
The folder ScatNet-QM-2D can be saved to any location. In Matlab, add to the path the top level folder ScatNet-QM-2D and include all subfolders.

## Overview

Within the top level folder are two main folders:
- **scatnet_light**: A modified version of the ScatNet Light package (https://github.com/edouardoyallon/ScatNetLight), originally developed by Edouard Oyallon with contributions by Matthew Hirn.
- **scat_qm**: The Scat QM plug-in to the ScatNet Light package used to compute quantum molecular energy regressions.

Most users will only need to work with a few files in the scat_qm folder.

## The Data Sets

Two data sets are included in the package. They are located in **scat_qm → data → data_sets**. The two data sets are:
- **qm7_2d.mat**: These are the 454 planar molecules from the QM7 data set (http://quantum-machine.org/datasets/).
- **qm_2d.mat**: This is a new data base, consisting of 4357 planar organic molecules. More information on its contents and how it was generated can be found in the qm_2d_README.rst file located in the same folder.

## Main Scripts

The main scripts are in **scat_qm → scripts → main**. There are six scripts, broken into three types:
- **qm_compute_xxx_2d.m**: Computes the dictionary coefficients of the specified type xxx (either **fourier** or **scat**). Note that wavelet invariant coefficients are contained within the scattering invariant coefficient computation.
- **qm_regression_xxx_ols_only.m**: Loads an already computed dictionary from qm_compute_xxx_2d.m, and computes a greedy regression using orthogonal least square using 5 specified folds. No bagging is performed, and the optimal number of regression terms is not estimated from the training data. These scripts generate the data used in Figure 4 in "Wavelet scattering regression of quantum chemical energies."
- **qm_regression_xxx.m**: Loads an already computed dictionary from qm_compute_xxx_2d.m, and computes a greedy regression using orthogonal least square using 5 specified folds. Bagging is utilized and the optimal number of regression terms is estimated from the training data. These scripts generate the data in Table 1 in "Wavelet scattering regression of quantum chemical energies."

## Displaying Figures

Scripts for generating figures from "Wavelet scattering regression of quantum chemical energies" can be found in the folder **scat_qm → scripts → display**.
- Figure 3 can be generated using **qm_scat_display_first_layer.m**.
- Figure 4 from can be generated using **qm_plot_rmse_dictionaries.m** and the outputs from qm_regression_xxx_ols_only.m.
- Figures 5, 6, and 7 can be generated using first **qm_compute_scat_avg_weights.m** and then **qm_display_scat_avg_weights.m**.

## Additional Documentation

Additional documentation for the ScatNet-QM-2D package can be found in the **ScatNet-QM-2D_manual.pdf** file.

## Authors

This package was developed and written by Matthew Hirn, Stéphane Mallat, and Nicolas Poilvert. Correspondence should be sent to Matthew Hirn at mhirn@msu.edu.

## License

Copyright 2016 Matthew Hirn, Stéphane Mallat, Nicolas Poilvert

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
