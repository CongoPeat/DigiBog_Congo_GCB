# DigiBog_Congo
Files for the version of DigiBog_Congo used in Young et al. (2023) published in Global Change Biology (Simulating carbon accumulation and loss in the central Congo peatlands). The model is set up to run in 1-D. These files were used in the simulations of the CEN-17.4 core near Ekologouma in the Central Congo Basin. The inputs here produce the main simulation results in the paper. Sensitivity tests can be performed using the parameter modifications described in the paper. See [Garcin et al. (2022)](https://www.nature.com/articles/s41586-022-05389-3) for the background to these simulations. DigiBog models are written in Fortran and must be compiled using a suitable compiler (e.g., gfortran, ifort) - see here https://fortran-lang.org if unsure about how to proceed, or contact the lead author of the paper to disuss the model. Simulations were run on Linux and Mac PCs.

The files are:

## Program files
```
DigiBog_Main_Congo.f90
DigiBog_Hydro_lumped.f90
DigiBog_new_layer.f90
global_def.f90
```
## Input files
```
010_DigiBog_CP_IN_information.txt
020_DigiBog_CP_IN_net_rain.txt
030_DigiBog_CP_IN_temp.txt
035_DigiBog_CP_IN_co2_ppm.txt
040_DigiBog_CP_IN_column_status.txt
041_DigiBog_CP_IN_column_status.txt
050_DigiBog_CP_IN_baltitude.txt
055_DigiBog_CP_IN_pft_params.comment (describes values in file 055*.txt)
055_DigiBog_CP_IN_pft_params.txt
```
## Bash wrapper scripts for makefiles (optional)
```
DB_run_make.sh
DB_make_debug.sh
```

## Makefiles (optional)
```
Make_DB
Make_DB_debug
```

If compiling the model using these scripts type ` :$ . DB_run_make.sh` in the terminal. 

## Main model run output files (** = files used in results)
```
060_DigiBog_CP_OUT_layer_mass.txt
**070_DigiBog_CP_OUT_column_height.txt**
**080_DigiBog_CP_OUT_wt_height.txt**
090_DigiBog_CP_OUT_transmiss.txt
100_DigiBog_CP_OUT_mass_area.txt
**110_DigiBog_CP_OUT_wt_depth.txt**
120_DigiBog_CP_OUT_layer_wet_prop.txt
130_DigiBog_CP_OUT_col_t_step.txt
140_DigiBog_CP_OUT_t_step_sum.txt
160_DigiBog_CP_OUT_wt_depth_summer.txt
170_DigiBog_CP_OUT_run_time.txt
180_DigiBog_CP_OUT_layers_in_columns.txt
190_DigiBog_CP_OUT_rem_mass.txt
200_DigiBog_CP_OUT_init_mass.txt
210_DigiBog_CP_OUT_layer_heights.txt
**220_DigiBog_CP_OUT_layer_age.txt**
```


