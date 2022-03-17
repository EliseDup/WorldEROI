## Acknowledging authorship
In the academic spirit of collaboration, the source code should be appropriately acknowledged in the resulting scientific disseminations.

You are welcome to report any bugs related to the code to the following:
 elise.dupont1@gmail.com

## Content
The main function of the model is world_grid.world_grid_eroi() that returns a pandas.dataframe with a row per grid cell.

The dataframe contains all the information to estimate wind and solar potential in every geographical location.

Input data from external databases are described in XLS file: *data/Legend_Wind_Solar_Data.xls*
* Latitude
* Longitude
* GHI : global horizontal irradiation [kWh/m^2/day]
* DNI : direct normal irradiation [kWh/m^2/day]
* ...

Then the estimated energy produced and energy inputs are computed within the *world_grid.world_grid_eroi()* method for every renewable technology considered.

* wind_onshore_e, wind_offshore_e, pv_e, csp_e = estimated energy outputs [EJ/year]
* wind_onshore_e_in, ... : energy inputs [EJ/year]
* wind_onshore_eroi, ... : associated EROI [-]
* ...

## How to run the model
Download the source files from git from a terminal with the following command (or directly via the website) :

`git clone https://github.com/EliseDup/WorldEROI`

* Download the dataframe from a python command line. 

First install the dependencies :
  
`pip install -r requirements.txt`

Then launch python and load the database with the world_grid_eroi() method :

`>>> import world_grid`

`>>> df = world_grid.world_grid_eroi()`


Ignore the 'ParserWarning' and 'DtypeWarning'.

That is it, you should have all the date in your dataframe *df*.


* Or ... Do more with the model with a Python IDE like PyCharm :-)

## References
* [Dupont et al., 2017] Dupont E., Koppelaar R. and Jeanmart H., Global available wind energy with physical and energy return on investment constraints, Applied Energy 209 (2018) 322–338

* [Dupont et al., 2019] Dupont E., Koppelaar R. and Jeanmart H., Global available solar energy under physical and energy
return on investment constraints, Applied Energy 257 (2020) 113968