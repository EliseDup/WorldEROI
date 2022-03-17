## Acknowledging authorship
In the academic spirit of collaboration, the source code should be appropriately acknowledged in the resulting scientific disseminations.

You are welcome to report any bugs related to the code to the following:
 elise.dupont1@gmail.com

## Content

## How to run the model
Download the source files from git from a terminal with the following command (or directly via the website) :

`git clone https://github.com/EliseDup/WorldEROI`

The main function of the model is world_grid.world_grid_eroi() that return a pandas dataframe with a row per grid cell.
It contains all the information to estimate the wind and solar potential in every geographical location: estimated energy outputs [EJ/year], energy inputs [EJ/year] and associated EROI [-]

## References
* [Dupont et al., 2017] Dupont E., Koppelaar R. and Jeanmart H., Global available wind energy with physical and energy return on investment constraints, Applied Energy 209 (2018) 322–338

* [Dupont et al., 2019] Dupont E., Koppelaar R. and Jeanmart H., Global available solar energy under physical and energy
return on investment constraints, Applied Energy 257 (2020) 113968