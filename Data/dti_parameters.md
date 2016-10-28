## ndmg diffusion parameters

#### Registration
| Step                          | Package | Function        | Parameters                | Description |
|-------------------------------|---------|-----------------|---------------------------|-------------|
| DTI De-noising                | FSL     | eddy_correct    | None specified (defaults) |  |
| MPRAGE Skull-stripping        | FSL     | bet             | `-B `                     | reduction of image bias and neck voxels |
| DTI align to MPRAGE           | FSL     | epi_reg         | None specified (defaults) | |
| MPRAGE align to MNI 1mm       | FSL     | flirt           | `-cost mutualinfo`;  `-bins 256`; `-dof 12`; `-searchrx -180 180`; `-searchry -180 180`; `-searchrz -180 180` | evaluated with mutual information, using 256 bin histograms of intensity, a 12 degree of freedom model, and searching the entire physical space for possible alignments |
| Apply MPRAGE transform to DTI | FSL     | flirt           | `-interp trilinear`;  `-applyxfm` | applying the transforms computed with trilinear interpolation |
| Resample aligned DTI          | nilearn  | resample_img   | `interpolation="nearest"` | resampling with nearest-neighbour interpolation |


#### Diffusion processing
| Step                          | Package | Function        | Parameters                | Description |
|-------------------------------|---------|-----------------|---------------------------|-------------|
| Tensor Fitting                | dipy    | TensorModel     | None specified (defaults) | |
| Tractography                  | dipy    | EuDX            | `a_low=0.1`               | FA stopping threshold for fiber tracking |
