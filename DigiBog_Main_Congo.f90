program DIGIBOG_CONGO

  !   DigiBog_Congo. a model to simulate peat accumulation in interfluvial
  !   domed peatlands in the Central Congo Basin.
  !   Copyright (c) 2023  the authors (see below)

  !   This program is free software: you can redistribute it and/or modify
  !   it under the terms of the gnu General Public License as published by
  !   the Free Software Foundation, either version 3 of the License, or
  !   (at your option) any later version.

  !   This program is distributed in the hope that it will be useful,
  !   but without any warranty; without even the implied warranty of
  !   merchantability or fitness for a particular purpose.  See the
  !   gnu General Public License for more details.

  !   You should have received a copy of the gnu General Public License
  !   along with this program.  If not, see https://www.gnu.org/licenses/.


!New version as of 05/04/2022. Modifications tracked via the DigiBog
!GitHub repository - https://github.com/DigiBog/DigiBog_Congo.git

! -----------------------------------------------------------------------------
! Section 1.0 Program header
! -----------------------------------------------------------------------------

  !Description (different model versions included here).
  !a model to simulate water tables  and peat accumulation / wastage in
  !peatlands and other shallow aquifers using the Boussinesq equation.
  !The model uses centimetres and seconds.
  !Details of how the model works may be obtained from the code owners
  !(email addresses below).

  !Summary of code configuration
  !Net rainfall = annual, monthly or weekly
  !Temperature = annually
  !Layers = lumped

  !Current code owners
  !Paul j. Morris, Andy j. Baird*, Lisa r. Belyea, Dylan m. Young, Pete Gill
  !*School of Geography
  !University of Leeds
  !Leeds
  !LS2 9JT
  !p.j.morris@reading.ac.uk
  !a.j.baird@leeds.ac.uk
  !l.belyea@qmul.ac.uk
  !d.m.young@leeds.ac.uk


  !Modification history of code
  !Programmer           Date           Modifications
  !============================================================================
  !Andrew j. Baird      08/04/2005     Original 1.5-d code
  !----------------------------------------------------------------------------
  !Paul j. Morris       24/02/2006     Conversion to 2.5-d and Fortran
  !                                    Standards implemented
  !----------------------------------------------------------------------------
  !Paul j. Morris       27/02/2006     Testing, minor corrections
  !----------------------------------------------------------------------------
  !Paul j. Morris       05/03/2006     Subroutines 'column_activation' and
  !                                    'steady_state_check' written. Former
  !                                    no longer exists (removed 18/07/2013).
  !----------------------------------------------------------------------------
  !Paul j. Morris       19/06/2006     Testing completed
  !----------------------------------------------------------------------------
  !Paul j. Morris       20/03/2007     Code cleaning
  !----------------------------------------------------------------------------
  !Paul j. Morris       24/04/2007     Further cleaning, including removal of
  !                                    replicated spatial step reference in
  !                                    'move_water' subroutine
  !----------------------------------------------------------------------------
  !Paul j. Morris       09/05/2007     Above-ground storage facilitated in
  !                                    2.5-d version
  !----------------------------------------------------------------------------
  !Paul j. Morris       02/07/2008     Final code cleaning, full annotation
  !----------------------------------------------------------------------------
  !Andy j. Baird        18/07/2013     Addition of variable net rainfall read
  !                                    in from file. Change to how boundary
  !                                    conditions specified; zero-flow
  !                                    Neumann condition also now allowed for,
  !                                    as are internal boundaries. Change to
  !                                    how output written to file (can now
  !                                    write to file multiple times before end
  !                                    of run).
  !----------------------------------------------------------------------------
  !Andy j. Baird        27/08/2013     Code cleaning and de-bugging of version
  !                                    from 18/07/2013
  !----------------------------------------------------------------------------
  !Dylan m. Young       20/05/14       Combine 1D and Hydro models for use in
  !                                    blanket bog simulation of management
  !                                    impacts. For testing, set fixed Dirichlet
  !                                    condition of column elevation * 0.5.
  !                                    Added memory allocation for x_flux and
  !                                    y_flux in wat_k_mean (Hydro_procs) for
  !                                    de-bugging purposes.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       11/06/14       Reorganisation of code as suggested by
  !                                    ajb and agreed with ajb and pjm.
  !                                    Removal of de-bugging code and resetting
  !                                    x_flux and y_flux to local arrays.
  !                                    Removal of steady-state check.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       18/06/14       Addition of variable Dirichlet condition
  !                                    water-table read in from the parameter
  !                                    file.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       13/07/14       Added timestep multiplier for calculated
  !                                    timestep for stability improvements
  !-----------------------------------------------------------------------------
  !Dylan m. Young       01/01/15       Added mineral layer for sloping model
  !                                    stability. Can be used for both raised
  !                                    and blanket bog models.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       17/03/15       Weather now read-in on a weekly basis
  !                                    (can be used for stable or variable
  !                                    weather). Mean annual temp. is used for
  !                                    new litter production.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       20/03/15       Recalcitrance added to anoxic
  !                                    decomposition.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       28/04/15       Added write-out for k profile and maximum
  !                                    timestep calc.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       06/12/16       Updated recalcitrance for anoxic
  !                                    decomp based on Morris et. al.,2015
  !                                    Rainfall varies weekly, temperature is
  !                                    set on an annual basis
  !-----------------------------------------------------------------------------
  !Pete Gill            04/2017        Added layer lumping routine to reduce
  !                                    model run time.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       06/07/17       Updated recalcitrance to include oxic and
  !                                    anoxic
  !                                    decomp based on Morris et. al.,2015
  !-----------------------------------------------------------------------------
  !Dylan m. Young       13/05/20       New version including bug fixes in the
  !                                    lumping routine for calculating the
  !                                    layer remaining mass correctly and to
  !                                    add the recalculation of transmissivity
  !                                    after layers have been lumped. This
  !                                    version supersedes all previous versions
  !                                    of the ecohydrological model code.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       27/08/20       New version for CongoPeat. Includes one
  !                                    Tree pft with three litter fractions
  !                                    leaves, wood and roots. The functions are
  !                                    from HPMTrop.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       11/11/21       All code modifications now documented on
  !                                    GitHub. Only major revisions will be
  !                                    noted here.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       05/04/22       Modified to include labile and
  !                                    recalcitrant pools in each litter
  !                                    fraction. The new parameters are read in
  !                                    from a new input file in 2D format.
  !-----------------------------------------------------------------------------
  !Code description.
  !Language: Fortran 95.
  !Software standards: Andrews, p. (1998). Unified Model Documentation Paper No.
  !3, Software Standards for the Unified Model: Fortran and Unix, Version 7.2,
  !Met Office, uk, 17 pp.
  !See also https://fortran-lang.org/learn/best_practices/style_guide/
  !-----------------------------------------------------------------------------

  !Modules used:
  use hydro_procedures_transmax
  use new_litter_procedures

  implicit none

  !-----------------------------------------------------------------------------
  !Format output table for remaining mass (y,z, rem mass, leaves mass, wood
  !mass, roots mass, layer thickness, column height. And layer thickness (x,y,
  !z, leaves thickness, wood thickness, roots thickness, grass thickness, total
  !layer thickness, and column height)
  1100 format (1X,I3,I3,I5,I5,4(20F20.10))
  !-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Section 2.0 Definitions; Declarations
!-------------------------------------------------------------------------------

! Single-valued variable (and model parameter) declarations.

  integer :: alloc_error,        &  !Memory allocation error flag
             t_extent,           &  !Number of simulated model years
             t_add,              &  !Additional years after spin-up
             output_interval,    &  !Interval (years) for writing to output
                                    !files
             sub_year_counter,   &  !Counter for use in main model loop
             sub_year_counter_su,&  !Counter for use in main model loop
                                    !for summer water tables
             output_counter,     &  !Output counter for use in main time loop
             year_counter,       &  !Counter for total number of model years
             sub_annual_counter, &  !Counter for total weeks or months
             new_net_rain,       &  !Flag for reading in new net rain
             x_extent,           &  !Model grid extent in x direction (number)
             y_extent,           &  !Model grid extent in y direction (number)
             z_extent,           &  !Spatial index counter
             x,                  &  !Spatial index counter
             y,                  &  !Spatial index counter
             z,                  &  !Spatial index counter
             i,                  &  !Index counter
             j,                  &  !Index counter
             k,                  &  !for exitstat
             pft,                &  !pft index counter
             pft_frac,           &  !pft fraction
             peat_top,           &  !marker for layer at the peatland surface
             litter_frac,        &  !litter fractions index counter
             recal_start,        &  !Start of range for recalcitrance prop.n
             recal_end,          &  !End of range for recalcitrance prop.n
             dens_start,         &  !Start of range for density parameters
             dens_end,           &  !End of range for density parameters
             new_m_start,        &  !Start of range for new mass values
             new_m_end,          &  !End of range for new mass values
             pft_prop,           &  !line in file 055_* for pft fraction
             net_rain_res,       &  !Net rainfall temporal resolution
             leaves,             &  !litter fraction number
             wood,               &  !litter fraction number
             roots,              &  !litter fraction number
             n_pft,              &  !number of PFTs
             n_fractions,        &  !number of litter fractions
             n_pools,            &  !number of pools per fraction
             n_pft_params,       &  !no. of parameters for each pft
             n_decay_params,     &  !no. of decay parameters for each pft
             spin_up,            &  !Running a spin up 1 = yes
             layer_exclusion,    &  !Don't lump these layers (numbered from top
                                    !of a column)
             lump_index

  real(kind=q) :: base_temp,           &  !deg. Celsius
                  Q10_oxic,            &  !Factor (dimensionless)
                  Q10_anoxic,          &  !Factor (dimensionless)
                  porosity,            &  !Drainable porosity (cm^3 cm^-3 --
                                          !a proportion)
                  k_param_a,           &  !The a parameter in the k equation
                                          !(cm yr^-1)
                  k_param_b,           &  !The b parameter in the k equation
                                          !(dimensionless)
                  timestep,            &  !yr
                  t_min,               &  !min timestep
                  t_max,               &  !max timestep
                  t_min_years,         &  !No of years before timestep increases
                  t_max_years,         &  !Period over which timestep increases
                  t_rate,              &  !Rate of increase in timestep
                  transmax,            &  !Maximum transmissivity for timestep
                  t_step_multi,        &  !timestep multiplication factor
                  t_step_sum_netr,     &  !sum of sub annual timesteps
                  t_step_sum_su,       &  !sum of summer timesteps
                  mean_t_step,         &  !Average annual timesteps
                  annual_tsteps,       &  !Number of annual timesteps
                  net_rain_tstep,      &  !Number of weekly timesteps
                  rainfall,            &  !Net rainfall (cm yr^-1)
                  temperature,         &  !deg. c
                  co2_ppm_ref,         &   !annual CO2 concentration (ppm)
                  co2_ppm_curr,        &  !annual CO2 concentration (ppm)
                  co2_ppm_chge,        &  !change in annual CO2 ppm
                  beta_fact,           &  !chge in npp per 100 ppm CO2
                  spatial_step,        &  !Model horizontal increment
                                          !(x & y) (cm)
                  pond_depth,          &  !Depth of surface ponding (cm)
                  time_start,          &  !Seconds
                  time_finish,         &  !Seconds
                  time_elapsed,        &  !Seconds
                  weight_z,            &  !For weighted avg. lumping
                  weight_z_minus,      &  !For weighted avg. lumping
                  threshold,           &  !Lump layers < this height (cm)
                  mineral,             &  !mineral soil thickness (cm)
                  pond_olf,            &  !pond parameter for olf
                  height_olf,          &  !peat height parameter for olf
                  max_root,            &  !maximum rootig depth (cm)
                  root_depth              !rooting depth (cm)


  !Climate intervals
  character(len=1) :: weather_error       !error flag
  character(len=1) :: type_error          !error flag

  !Input and output files
  character(len=40) :: data_file_name_010, & !in Model run information file
                      data_file_name_020, &  !in net rainfall values
                      data_file_name_030, &  !in temperature values
                      data_file_name_035, &  !in  CO2 concentration
                      data_file_name_040, &  !in column activation status file
                      data_file_name_050, &  !in base altitude input file
                      data_file_name_055, &  !in pft parameters in 2D array
                      data_file_name_060, &  !out Layer mass remaining output
                      data_file_name_070, &  !out Column height output file
                      data_file_name_080, &  !out Water-table height output file
                      data_file_name_090, &  !out Transmissivity output file
                      data_file_name_100, &  !out Mass per area output file
                      data_file_name_110, &  !out Mean annual column water table
                      data_file_name_120, &  !out Layer wet proportion file
                      data_file_name_130, &  !out timestep
                      data_file_name_140, &  !out timestep sum
                      data_file_name_160, &  !out Summer water table
                      data_file_name_170, &  !out Run Time
                      data_file_name_180, &  !out Layers per column
                      data_file_name_190, &  !out Layer mass remaining,
                                             !elevation, and thickness
                      data_file_name_200, &  !out initial layer mass
                      data_file_name_210, &  !out final layer height
                      data_file_name_220     !out layer age

!Prevent output files being overwritten if the simulation is a spin-up read run
character (len=10) :: file_status, file_pos

!For the 3-d files
character (len=40) :: data_file_name_510, & !Wet proportion out
                      data_file_name_520, & !Layer storage out
                      data_file_name_530, & !Layer_age out
                      data_file_name_540, & !Layers in layers out
                      data_file_name_550, & !Layer attributes out
                      data_file_name_560, & !Layer mass all out
                      data_file_name_570, & !pft layer mass
                      data_file_name_580, & !Water table height out
                      data_file_name_590, & !Number of layers out
                      data_file_name_600    !Transmissivity out

  character(len = 255) :: cwd

  real(kind=q), allocatable, dimension(:) :: fract_recal, & !Recalcitrance prop
                                                            !of litter fractions
                                             new_layer_mass !Mass of litter
                                                            !g cm^-2

  real(kind=q), allocatable, dimension(:,:) :: base_altitude,     & !Above datum
                                               water_change,      &
                                               water_table,       & !Above base
                                               wk_mean,           & !Depth-av. k
                                               col_wt_depth,      & !cm
                                               col_wt_sum,        & !cm
                                               col_wt_sum_su,     & !cm
                                               col_wt_depth_av,   & !cm
                                               col_wt_depth_av_su,& !cm
                                               col_mass_per_area, & !g cm^-2
                                               ann_prod,          & !g cm^-2
                                               density,           & !g cm^-3
                                               pft_params,        & !various
                                               pft_params_ref,    & !various
                                               decay_params         !prop. yr-1


  integer, allocatable, dimension(:,:) :: no_layers, &  !Number of layers (x,y)
                                          root_count    !Number of layers (x,y)
                                                         !within rooting depth

  logical, allocatable, dimension(:,:,:) :: root_mask !Layers within rooting
                                                      !depth (true or false)

  integer, allocatable, dimension(:,:,:) :: layers_in_layer


  real(kind=q), allocatable, dimension(:,:,:) :: wet_proportion, & !proportion
                                                 layer_age, & !years
                                                 bg_root_mass !g cm^2

  !layer_attributes stores layer thickness, k (cm yr^-1) and s (cm^3 cm^-3)
  !transmissivity stores layer elevation above base and transmissivity
  !layer_storage stores layer elevation (cm) above base and each layer's water
  !capacity as volume per unit area; i.e. expressed as a depth (cm)
  !layer_mass stores layer current (g cm^-2), initial (g cm^-2)
  !and remaining mass (proportion).
  real(kind=q), allocatable, dimension(:,:,:,:) :: layer_attributes, &
                                                   transmissivity, &
                                                   layer_storage, &
                                                   layer_mass
  !Layer mass for all pfts
  !Stores the values for all pfts, their litter fractions (e.g. leaves) and the
  !labile and recalcitrant pools of material. Mass (g cm^-2), thickness (cm) and
  !mass remaining (proportion) are stored for each litter fraction as
  !1 - current mass labile, 2 - current mass recalcitrant,
  !3 - original mass labile, 4 - original mass recalcitrant,
  !5 - thickness (whole fraction), 6 - remaining mass (whole fraction).
  !i.e. pft_layer_mass(x, y, z, pft, fraction, 1:6)
  real(kind=q), allocatable, dimension(:,:,:,:,:,:) :: pft_layer_mass

  !activation_status indicates whether/how column participates in simulation
  character(8), allocatable, dimension(:,:) :: activation_status

!-------------------------------------------------------------------------------
! Section 3.0 Data Input and output; Memory Management
!-------------------------------------------------------------------------------
  !Name input and output data files.

  !Inputs
  data_file_name_010 = "010_DigiBog_CP_IN_information.txt"
  data_file_name_020 = "020_DigiBog_CP_IN_net_rain.txt"
  data_file_name_030 = "030_DigiBog_CP_IN_temp.txt"
  data_file_name_035 = "035_DigiBog_CP_IN_co2_ppm.txt"
  data_file_name_040 = "040_DigiBog_CP_IN_column_status.txt"
  data_file_name_050 = "050_DigiBog_CP_IN_baltitude.txt"
  data_file_name_055 = "055_DigiBog_CP_IN_pft_params.txt"

! Open files for input data.
  open(unit=010, file=data_file_name_010, status="old")
  open(unit=020, file=data_file_name_020, status="old")
  open(unit=030, file=data_file_name_030, status="old")
  open(unit=035, file=data_file_name_035, status="old")
  open(unit=040, file=data_file_name_040, status="old")
  open(unit=050, file=data_file_name_050, status="old")
  open(unit=055, file=data_file_name_055, status="old")

  !Read data from information input data file.
  read (010, *) n_pft
  read (010, *) n_fractions
  read (010, *) n_pools
  read (010, *) n_pft_params
  read (010, *) recal_start
  read (010, *) recal_end
  read (010, *) dens_start
  read (010, *) dens_end
  read (010, *) new_m_start
  read (010, *) new_m_end
  read (010, *) pft_prop
  read (010, *) base_temp
  read (010, *) Q10_oxic
  read (010, *) Q10_anoxic
  read (010, *) porosity
  read (010, *) k_param_a
  read (010, *) k_param_b
  read (010, *) t_extent
  read (010, *) annual_tsteps
  read (010, *) output_interval
  read (010, *) x_extent
  read (010, *) y_extent
  read (010, *) spatial_step
  read (010, *) pond_depth
  read (010, *) max_root
  read (010, *) t_step_multi
  read (010, *) layer_exclusion
  read (010, *) threshold
  read (010, *) mineral
  read (010, *) net_rain_res
  read (010, *) t_min
  read (010, *) t_max
  read (010, *) t_min_years
  read (010, *) t_max_years
  read (010, *) pond_olf
  read (010, *) height_olf
  read (010, *) co2_ppm_ref
  read (010, *) beta_fact
  read (010, *) spin_up
  read (010, *) t_add

  !Output files
  if (spin_up == 2) then
    file_status = "old"
    file_pos = "append"
  else
    file_status = "replace"
    file_pos = "asis"
  end if

  data_file_name_060 = "060_DigiBog_CP_OUT_layer_mass.txt"
  data_file_name_070 = "070_DigiBog_CP_OUT_column_height.txt"
  data_file_name_080 = "080_DigiBog_CP_OUT_wt_height.txt"
  data_file_name_090 = "090_DigiBog_CP_OUT_transmiss.txt"
  data_file_name_100 = "100_DigiBog_CP_OUT_mass_area.txt"
  data_file_name_110 = "110_DigiBog_CP_OUT_wt_depth.txt"
  data_file_name_120 = "120_DigiBog_CP_OUT_layer_wet_prop.txt"
  data_file_name_130 = "130_DigiBog_CP_OUT_col_t_step.txt"
  data_file_name_140 = "140_DigiBog_CP_OUT_t_step_sum.txt"
  data_file_name_160 = "160_DigiBog_CP_OUT_wt_depth_summer.txt"
  data_file_name_170 = "170_DigiBog_CP_OUT_run_time.txt"
  data_file_name_180 = "180_DigiBog_CP_OUT_layers_in_columns.txt"
  data_file_name_190 = "190_DigiBog_CP_OUT_rem_mass.txt"
  data_file_name_200 = "200_DigiBog_CP_OUT_init_mass.txt"
  data_file_name_210 = "210_DigiBog_CP_OUT_layer_heights.txt"
  data_file_name_220 = "220_DigiBog_CP_OUT_layer_age.txt"

  open(unit=060, file=data_file_name_060, status=file_status, &
    position=file_pos)
  open(unit=070, file=data_file_name_070, status=file_status, &
    position=file_pos)
  open(unit=080, file=data_file_name_080, status=file_status, &
    position=file_pos)
  open(unit=090, file=data_file_name_090, status=file_status, &
    position=file_pos)
  open(unit=100, file=data_file_name_100, status=file_status, &
    position=file_pos)
  open(unit=110, file=data_file_name_110, status=file_status, &
    position=file_pos)
  open(unit=120, file=data_file_name_120, status=file_status, &
    position=file_pos)
  open(unit=130, file=data_file_name_130, status=file_status, &
    position=file_pos)
  open(unit=140, file=data_file_name_140, status=file_status, &
    position=file_pos)
  open(unit=160, file=data_file_name_160, status=file_status, &
    position=file_pos)
  open(unit=170, file=data_file_name_170, status=file_status, &
    position=file_pos)
  open(unit=180, file=data_file_name_180, status=file_status, &
    position=file_pos)
  open(unit=190, file=data_file_name_190, status=file_status, &
    position=file_pos)
  open(unit=200, file=data_file_name_200, status=file_status, &
    position=file_pos)
  open(unit=210, file=data_file_name_210, status=file_status, &
    position=file_pos)
  open(unit=220, file=data_file_name_220, status="replace", &
    position=file_pos)

  !issue #9
  !Set up output files if the simulation is a spin-up run.
  if (spin_up == 1) then
    data_file_name_510 = "510_DigiBog_CP_OUT_wet_prop.dat"
    data_file_name_520 = "520_DigiBog_CP_OUT_layer_store.dat"
    data_file_name_530 = "530_DigiBog_CP_OUT_layer_age.dat"
    data_file_name_540 = "540_DigiBog_CP_OUT_lyrs_in_lyrs.dat"
    data_file_name_550 = "550_DigiBog_CP_OUT_lyr_attr.dat"
    data_file_name_560 = "560_DigiBog_CP_OUT_lyr_mass.dat"
    data_file_name_570 = "570_DigiBog_CP_OUT_pft_lyr_mass.dat"
    data_file_name_580 = "580_DigiBog_CP_OUT_water_table.dat"
    data_file_name_590 = "590_DigiBog_CP_OUT_no_layers.dat"
    data_file_name_600 = "600_DigiBog_CP_OUT_transmiss.dat"
  end if

  !Binary files for 3-d version
  if (spin_up == 2) then
  !Copy the binary files for a 3-d simulation to the inputs dir.
  !Delete any binary files in the input dir and then copy the news ones there.
  !This forces the user to run a spin up to generate the 1-d inputs.
    !Check the directory exists; if not, create it.
    call EXECUTE_COMMAND_LINE( &
    "if [ ! -d inputs_3d ]; then mkdir inputs_3d; fi", &
      wait = .true.)
    !Delete any old files in the inputs dir and copy the new ones
    call EXECUTE_COMMAND_LINE("rm -f inputs_3d/[5,6]*0*.dat 2>/dev/null && mv [5,6]*0*.dat inputs_3d/ 2>/dev/null", &
        wait = .true., exitstat = k)
    if (k /= 0) then
      write (*, *) ""
      write (*, *) "** No output files to move to input_3d/ **"
      stop
    end if
  endif

  !Calculate total number of parameters for all pft fractions
  n_pft_params = (n_pft_params * n_fractions) + 1 !pft proportion of new mass

  !Is this run a spin-up?
  if (spin_up == 1) then
    print *
    write (*, *) "This is a spin-up write run (spin_up = 1)"
  else if (spin_up == 2) then
    print *
    write (*, *) "This is a spin-up read run (spin_up = 2)"
  else
    print *
    write (*, *) "This is a standard model run (no spin up)"
  end if

  write (*, '(A31)') "Is this  correct (y/n)?"
  read *, type_error
    select case ( type_error )
  case ("n")
    write (*, *) "incorrect simulation type"
    stop
  case default
    write (*, *) "simulation type okay - simulation continues"
  end select
  write (*, *) " "

  !check if the climate resolution is correct
  write (*,*) "Net rain resolution =",net_rain_res, "(1=year, 2=month, 3=week)"

  write (*, '(A31)') "Is this  correct (y/n)?"
  read *, weather_error
  select case ( weather_error )
    case ("n")
      write (*, *) "Incorrect interval for climate input"
        stop
    case default
      write (*, *) "Climate resolution okay - simulation continues"
  end select

  !Set the extent of layers
  !note: for spin-up read simulations.
  !This would be a problem without lumping because the array wouldn't
  !accomodate the new layers added when running a spin-up read.
  z_extent = t_extent + 2 !Ponding layer + added mineral + peat layers

  !Allocate memory to arrays.
  !Columns --------------------------------------------------------------------

  allocate(col_wt_depth(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for col_wt_depth"
    stop
  end if

  allocate(col_wt_depth_av(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for col_wt_depth_av"
    stop
  end if

  allocate(col_wt_depth_av_su(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for col_wt_depth_av"
    stop
  end if

  allocate(col_wt_sum(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for col_wt_sum"
    stop
  end if

  allocate(col_wt_sum_su(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for col_wt_sum_su"
    stop
  end if

  allocate(col_mass_per_area(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for col_mass_per_area"
    stop
  end if

  allocate(base_altitude(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for base_altitude"
    stop
  end if

  allocate(ann_prod(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for ann_prod"
    stop
  end if


  !Hydrology ------------------------------------------------------------------

  allocate(water_change(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for water_change"
    stop
  end if

  allocate(water_table(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for water_table"
    stop
  end if

  allocate(wk_mean(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for wk_mean"
    stop
  end if

  allocate(activation_status(x_extent, y_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for activation_status"
    stop
  end if

  allocate(wet_proportion(x_extent, y_extent, z_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for wet_proportion"
    stop
  end if

  allocate(transmissivity(x_extent, y_extent, z_extent, 2), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for transmissivity"
    stop
  end if


  !Layers and layer lumping ----------------------------------------------------

  allocate(layer_storage(x_extent, y_extent, z_extent, 2), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for layer_storage"
    stop
  end if

  allocate(no_layers(x_extent,y_extent), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not alocate space for no_layers"
    stop
  end if

  allocate(layer_age(x_extent, y_extent, z_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for wet_proportion"
    stop
  end if

  allocate(layers_in_layer(x_extent, y_extent, z_extent), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for wet_proportion"
    stop
  end if

  allocate(layer_attributes(x_extent, y_extent, z_extent, 3), &
                                                             stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for layer_attributes"
    stop
  end if


  !Below ground roots ----------------------------------------------------------

  allocate(root_count(x_extent,y_extent), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not alocate space for root_count"
    stop
  end if

  allocate(bg_root_mass(x_extent,y_extent,n_pft), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not alocate space for bg_root_mass"
    stop
  end if

  !The count of layers for rooting depth ignores the pond
  allocate(root_mask(x_extent, y_extent, z_extent - 1), stat=alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not allocate space for root_mask"
    stop
  end if


  !Plant litter fractions-------------------------------------------------------

  !Recalcitrance proportions for each litter fraction
  allocate(fract_recal(n_fractions), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for fract_recal"
    stop
  end if

  !New mass values for each litter fraction in each pft
  allocate(new_layer_mass(n_fractions), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for new_layer_mass"
    stop
  end if

  !Dry bulk density for each litter fraction in each pft
  allocate(density(n_pft, n_fractions), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for density"
    stop
  end if

  !Total layer mass array
  allocate(layer_mass(x_extent, y_extent, z_extent, 3), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for layer_mass"
    stop
  end if

  !Mass of all pfts
  allocate(pft_layer_mass(x_extent, y_extent, z_extent, n_pft, n_fractions, 6),&
          stat = alloc_error)
  if (alloc_error /= 0) then
    write (*, *) "Model could not allocate space for pft_layer_mass"
    stop
  end if

  !All pft parameters for labile and recalcitrant pools and bulk density for
  !each fraction in each pft
  allocate(pft_params(n_pft, n_pft_params), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not alocate space for pft_params"
    stop
  end if
  !Save original (reference parameters for CO2 concentration calcs)
  allocate(pft_params_ref(n_pft, n_pft_params), stat = alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not alocate space for pft_params_ref"
    stop
  end if

  !There are four decay parameters per litter fraction (2 pools x 2 parameters)
  !the 2 parameters are ox and anox.
  allocate(decay_params(n_pft, n_fractions * n_pools * 2), &
          stat = alloc_error)
  if (alloc_error /= 0) then
    write (*,*) "Model could not alocate space for decay_params"
    stop
  end if

!-------------------------------------------------------------------------------
! Section 4.0 Initialisation
!-------------------------------------------------------------------------------

  !Initialise all single-valued variables and arrays.

  !Single-valued variables.
  year_counter = 1
  sub_annual_counter = 0
  annual_tsteps = 1.0
  net_rain_tstep = 0
  timestep = t_min / real(525600)
  t_step_sum_netr = 0.0
  t_step_sum_su = 0.0
  mean_t_step = 0.0
  sub_year_counter = 0
  sub_year_counter_su = 0
  output_counter = 0
  x = 1
  y = 1
  z = 1
  i = 1
  j = 1
  new_net_rain = 1
  weight_z = 0
  weight_z_minus = 0
  pft = 1
  pft_frac = 1
  litter_frac = 1
  leaves = 1
  wood = 2
  roots = 3
  peat_top = 0
  n_decay_params = n_fractions * n_pools

  !Arrays.
  no_layers = 0
  layer_mass = 0.0
  pft_layer_mass = 0.0
  wet_proportion = 0.0
  col_wt_depth = 0.0
  col_wt_sum = 0.0
  col_wt_sum_su = 0.0
  col_wt_depth_av = 0.0
  col_wt_depth_av_su = 0.0
  col_mass_per_area = 0.0
  layer_attributes = 0.0
  water_change = 1.0
  wk_mean = 1.0
  transmissivity = 1.0
  layer_storage = 1.0
  water_table = 0.0
  base_altitude = 0.0
  layer_age = 0.0
  layers_in_layer = 0
  root_count = 0
  bg_root_mass = 0.0
  root_mask = .false.
  ann_prod = 0.0
  decay_params = 0.0
  fract_recal = 0.0
  density = 0.0
  new_layer_mass = 0.0

  !Read data from files to rank-two arrays for active columns
  do x = 1, x_extent
    do y = 1, y_extent
      read (040, *) activation_status(x, y)
      read (050, *) base_altitude(x, y)
    end do
  end do

  !Read data for pft decay parameters into array. Because the input is in
  !column format, it needs to be read in as an implicit do loop without comments
  !The pft_params file is organised as follows:
  !oxic decay labile for each fraction, anoxic decay labile for each fraction
  !oxic decay recalcitrant for each fraction, anoxic decay recalcitrant for each
  !fraction. Proportion of recalcitrant material for each fraction. Density for
  !each fraction. New mass for each fraction. Proportion of pft in forest.
  !This equals (n_pft_params (from the 010 file) * n_fractions) + 1 (pft split)
  !note: The number of columns in the file must equal n_pft.

  read (055, *) ((pft_params(i, j), i = 1, n_pft), &
          j = 1, n_pft_params)

  !Copy the pft params arrary for use in CO2 concentrations sub-routine
  pft_params_ref = pft_params

  !Bulk density values for each fraction in each pft
  density(:, :) = pft_params(:, dens_start:dens_end)

  ! Initialise the properties for a layer of mineral soil.
  do x = 1, x_extent
    do y = 1, y_extent
      if (activation_status(x, y) == "on") then
        !Specify layer thickness
        layer_attributes(x, y, 1, 1) = mineral
        !Calculate layer 1 k
        layer_attributes(x, y, 1, 2) = 3153.6 !Layer k = 10-4cm s-1
        !Layer porosity
        layer_attributes(x, y, 1, 3) =  0.2
        !Column elevation for layer 1 = layer thickness of layer 1
        transmissivity(x, y, 1, 1) = layer_attributes(x, y, 1, 1)
        !Calculate layer 1 transmissivity = layer 1 thickness * layer 1 k
        transmissivity(x, y, 1, 2) =  transmissivity(x, y, 1, 1) &
          * layer_attributes(x, y, 1, 2)
        !Set wet proportion (saturated)
        wet_proportion(x, y, 1) = 1.0
        !Initial column properties
        !Assign mass per area = layer 1 mass
        !col_mass_per_area(x, y) = layer_mass(x, y, 1, 1)
        !Set the water table to the top of layer 1 for active columns
        water_table(x, y) = transmissivity(x, y, 1, 1) * wet_proportion(x, y, 1)
      end if
    end do
  end do

  !Skip these routines if the simulation is not a spin-up write.
  if (spin_up /= 2) then
    !Assign three layers to active columns - mineral, z1 and pond
    do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
          no_layers(x, y) = 3
        end if
      end do
    end do

    !Initialise pond (z = 3) properties.
    layer_attributes(:x_extent, :y_extent, 3, 1) = pond_depth !thickness
    layer_attributes(:x_extent, :y_extent, 3, 2) = 0.0 !k
    layer_attributes(:x_extent, :y_extent, 3, 3) = 1.0 !s

    !Calculate the peat properties for z2 litter fractions (leaves, wood, roots)
    ! Initialise properties of first layer of peat profile for "on" columns.
    ! pft_params(pft,recal_start:recal_end) holds the proportion of
    ! recalcitrant material in the litter fractions where pft is the number
    ! of the pft.

    !Array of new litter mass

    !Modify the mass to be added to z2 according to the CO2 concentration
    call co2_conc(pft_params, beta_fact, co2_ppm_ref, new_m_start, &
        new_m_end, pft_params_ref, co2_ppm_curr)

    do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
          !Calculations and assignments for columns with active status.
          !Calculate layer 2 mass and thickness
          !For each pft
          do pft = 1, n_pft
            !Set the recalcitrant proportions for the pft
            fract_recal = pft_params(pft,recal_start:recal_end)
            !Set the new litter mass for each fraction
            new_layer_mass = pft_params(pft,new_m_start:new_m_end) * &
            pft_params(pft, pft_prop)
            !For each litter fraction
            do litter_frac = 1, n_fractions
              !Current mass and original mass of the two pools
              !Labile pool
              pft_layer_mass(x, y, 2, pft, litter_frac, (/1,3/)) = &
              new_layer_mass(litter_frac) * (1 - fract_recal(litter_frac))
              !Recalcitrant pool
              pft_layer_mass(x, y, 2, pft, litter_frac, (/2,4/)) = &
              new_layer_mass(litter_frac) * fract_recal(litter_frac)
              !Set thickness of each fraction
              pft_layer_mass(x, y, 2, pft, litter_frac, 5) = &
              sum(pft_layer_mass(x, y, 2, pft, litter_frac, (/1,2/))) / &
              density(pft,litter_frac)
              !Set remaining mass of each fraction (= 100%)
              pft_layer_mass(x, y, 2, pft, litter_frac, 6) = 1.0
            end do
          end do

          !Calculate total layer 2 mass and set total original mass
          !Add the current mass for labile and recalcitrant pools for each
          !fraction, and set both current and original mass for total mass
          !array.
          layer_mass(x, y, 2, 1:2) = &
          sum(pft_layer_mass(x, y, 2, :, :, (/1,2/)))

          !Set z2 remaining mass to 100%
          layer_mass(x, y, 2, 3) = 1.0

          !Calculate total layer 2 thickness
          layer_attributes(x, y, 2, 1) = &
          sum(pft_layer_mass(x, y, 2, :, :, 5))

          !Calculate layer 2 k
          !Equation below from Morris et al. (2011) - see main model loop.
          layer_attributes(x,y, 2, 2) =  k_param_a * (exp (k_param_b))

          !Layer porosity
          layer_attributes(x, y, 2, 3) =  porosity

          !Column elevation for layer 2 = layer thickness of layer 1 + layer 2
          transmissivity(x, y, 2, 1) = transmissivity(x, y, 1, 1) + &
          layer_attributes(x, y, 2, 1)

          !Calculate layer 2 transmissivity = layer 1 thickness * layer 1 k
          transmissivity(x, y, 2, 2) =  transmissivity(x, y, 1, 1) &
          + (layer_attributes(x, y , 2 ,1) &
          * layer_attributes(x, y, 2, 2))
          !Set wet proportion (saturated)
          wet_proportion(x, y, 2) = 1.0

          !Initial column properties
          !Assign mass per area = layer 2 mass
          col_mass_per_area(x, y) = layer_mass(x, y, 2, 1)

          !Set the water table to the top of layer 2 for active columns
          water_table(x, y) = transmissivity(x, y, 2, 1)

          !Set the age of the first peat layer (i.e. z = 2) as the oldest
          layer_age(x, y, 2) = t_extent
          layers_in_layer(x, y, 2) = 1

        end if
      end do
    end do
  end if

  !If this simulation is a spin-up read, then read in the *.dat files
  if (spin_up == 2) then
     !Read in the binary files
    call getcwd(cwd)
    open(unit = 510, &
           file =trim (cwd)//'/inputs_3d/510_DigiBog_CP_OUT_wet_prop.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 520, &
           file =trim (cwd)//'/inputs_3d/520_DigiBog_CP_OUT_layer_store.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 530, &
             file =trim (cwd)//'/inputs_3d/530_DigiBog_CP_OUT_layer_age.dat', &
             status = 'old', form = "unformatted", action = "read" &
           )
    open(unit = 540, &
           file =trim (cwd)//'/inputs_3d/540_DigiBog_CP_OUT_lyrs_in_lyrs.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 550, &
           file =trim (cwd)//'/inputs_3d/550_DigiBog_CP_OUT_lyr_attr.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 560, &
           file =trim (cwd)//'/inputs_3d/560_DigiBog_CP_OUT_lyr_mass.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 570, &
           file =trim (cwd)//'/inputs_3d/570_DigiBog_CP_OUT_pft_lyr_mass.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 580, &
           file =trim (cwd)//'/inputs_3d/580_DigiBog_CP_OUT_water_table.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 590, &
           file =trim (cwd)//'/inputs_3d/590_DigiBog_CP_OUT_no_layers.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )
    open(unit = 600, &
           file =trim (cwd)//'/inputs_3d/600_DigiBog_CP_OUT_transmiss.dat', &
           status = 'old', form = "unformatted", action = "read" &
         )

    !Note:this currently only works for 1-d to 1-d or 1-d to 2-d, not 2-d to 2-d
    read (510) wet_proportion(2, 2, 2:)
    read (520) layer_storage(2, 2, 2:, :)
    read (530) layer_age(2, 2, 2:)
    read (540) layers_in_layer(2, 2, 2:)
    read (550) layer_attributes(2, 2, 2:, :)
    read (560) layer_mass(2, 2, 2:, :)
    read (570) pft_layer_mass(2, 2, 2:, :, :, :)
    read (580) water_table(2, 2)
    read (590) no_layers(2, 2)
    read (600) transmissivity(2, 2, 2:, :)

    !Copy the 1D column and layers to all other x, y locations
    do y = 2 , y_extent - 1
      do x = 2, x_extent - 1
        wet_proportion(x, y, 2:) = wet_proportion(2, 2, 2:)
        layer_storage(x, y, 2:, :) = layer_storage(2, 2, 2:, :)
        !Set the layer age to take into account the additional time
        layer_age(x, y, 2:) = layer_age(2, 2, 2:) + t_add
        layers_in_layer(x, y, 2:) = layers_in_layer(2, 2, 2:)
        layer_attributes(x, y, 2:, :) =layer_attributes(2, 2, 2:, :)
        layer_mass(x, y, 2:, :) = layer_mass(2, 2, 2:, :)
        pft_layer_mass(x, y, 2:, :, :, :) = &
          pft_layer_mass(2, 2, 2:, :, :, :)
        transmissivity(x, y, 2:, :) = transmissivity(2, 2, 2:, :)
        water_table(x, y) = water_table(2, 2)
        no_layers(x, y) = no_layers(2, 2)
      end do
    end do

    !Set the year counter
    year_counter = t_extent + 1

    !Add a new layer (would usually be added when model initialises for year 1)
    !Calculate the root depth and belowground root mass
    call rooting(x_extent, y_extent, no_layers, transmissivity, max_root, &
      root_depth, mineral, root_count, root_mask, bg_root_mass, n_pft, &
      new_layer_mass, pft_params, new_m_start, new_m_end, pft_prop, &
      activation_status, roots)

    !Add a new layer to the peat surface
    call add_layer(x_extent, y_extent, no_layers, layer_attributes, &
        pond_depth, activation_status)

    !Modify the mass of the new top layer according to the CO2 concentration
    call co2_conc(pft_params, beta_fact, co2_ppm_ref, new_m_start, &
        new_m_end, pft_params_ref, co2_ppm_curr)

    !Add above- and belowground litter
    call new_layer(x_extent, y_extent, no_layers, layer_attributes, &
        activation_status, root_depth, root_count, layer_age, t_extent, t_add, &
        layer_mass, pft_layer_mass, n_pft, fract_recal, recal_start, &
        recal_end, new_m_start, new_m_end, pft_prop, leaves, wood, roots, &
        new_layer_mass, year_counter, layers_in_layer, col_wt_depth_av, &
        n_fractions, n_pools, density, bg_root_mass, transmissivity, porosity, &
        col_mass_per_area, k_param_a, k_param_b, pft_params)

    !Don't reset the time step to the min value
    timestep = t_max / real(525600)
  end if

  !Set water-table height for Dirichlet condition. x1 and x_extent are
  !"neu" or "off" cells as is y_extent
  !For the first CongoPeat simulations, the following isn't relevant (we use
  !a 1D model with no subsurface flow between columns).
  Diri_wt: do x = 2, (x_extent -1)
    do y = 1, y_extent
      if (activation_status(x, y) == "diri") then
         !Set to height of the peat surface
         water_table(x, y) =  &
           transmissivity(x, y + 1, (no_layers(x, y + 1) - 1), 1)
         !Or set a value independent of the next active column e.g.
         !water_table(x, y) = 5
      end if
    end do
  end do Diri_wt

  !Determine the read frequency for net rainfall (annual (1), monthly (2),
  !weekly (3))
  !Annual
  if (net_rain_res == 1) then
          net_rain_tstep = annual_tsteps
  !Monthly
  else if (net_rain_res == 2) then
          net_rain_tstep = annual_tsteps / real(12)
  !Weekly
  else if (net_rain_res == 3) then
          net_rain_tstep = annual_tsteps / real(52)
  end if

  !Timestep management. If a variable timestep is not used, one that increases
  !linearly can be used to speed up the model. The rate of increase and min and
  !max values are currently ascertained through trial and error.
  !Set the rate of timestep increase between t_min and t_max
  t_rate = (t_max - t_min) / (t_max_years - t_min_years)

!-------------------------------------------------------------------------------
! Section 5.0 Main Calculations; Time Management
!-------------------------------------------------------------------------------
  call cpu_time(time_start)

  !1.Commence main annual loop -------------------------------------------------
  Annual_loop: do

    !Weather
    !If annual weather, read in net rainfall
    if (net_rain_res == 1) then
        read (020, *) rainfall
    end if

    !Read annual temperature
    read (030, *) temperature

    !Calculate the difference in CO2 ppm between the current and reference
    !values
    co2_ppm_chge = co2_ppm_curr - co2_ppm_ref

    !Display model clock on screen in units of years.
    write (*, *) "Model year ", year_counter

    !read in the oxic and anoxic parameters for the litter fractions
    !Calculate decay rate variables for current year.

    !New routines for calculating decay parameters for the litter pools
    !Oxic and anoxic decay rates for labile pools in each parameter.
    !These are the first 6 values from each pft
    do i = 1, n_pft
      do j = 1, n_fractions * n_pools
        decay_params(i, j) = pft_params(i ,j) &
          * Q10_oxic**((temperature - base_temp) / 10.0)
      end do
    end do

    !Oxic and anoxic decay rates for recalcitrant pools in each parameter.
    !These are the second set of values (n_fractions x n_pools) from each pft
    do i = 1, n_pft
      do j = (n_fractions * n_pools) + 1, n_fractions * n_pools * 2
       decay_params(i, j) = pft_params(i ,j) &
           * Q10_anoxic**((temperature - base_temp) / 10.0)
      end do
    end do

    !2.Start sub annual time steps -------------------------------------------
    Sub_annual_loop: do

      !2.1 Net rainfall and counters
      !Check if a new sub-annual time step has been signalled with new_year.
      !If so read in the next values for net rainfall
      if (net_rain_res /= 1) then
        if (new_net_rain == 1) then
          !Read next net rainfall
          read (020, *) rainfall
          !Reset time counters for next trip through sub annual loop
          !Summer water tables are June, July, August which doesn't correspond
          !to the wet seasons in central Congo.
          t_step_sum_netr = 0.0
          t_step_sum_su = 0.0
          new_net_rain = 0
        end if
      end if

      !Count the number of sub-annual loops (≡ annual_tsteps)
      !the count is used in the average water-table calc in the annual loop
      sub_year_counter = sub_year_counter + 1
      !Summer water tables (nothern hemisphere). Change these monthly values
      !to output the required time series (e.g. if there are two wet seasons
      !for example)
      if (net_rain_res /= 1) then
        if (net_rain_res == 3) then
          if (sub_annual_counter >= 22 .and. sub_annual_counter <= 35) then
              sub_year_counter_su = sub_year_counter_su + 1
          end if
        else if (net_rain_res == 2) then
          if (sub_annual_counter == 6 .or. sub_annual_counter == 7 .or. &
             sub_annual_counter == 8) then
              sub_year_counter_su = sub_year_counter_su + 1
          end if
        end if
      end if

      !Re-initialise column properties
      col_mass_per_area = 0.0

      !This code only  applies if a variable time step is used.
      !Set the timestep for this iteration to be equal to the relevant timescale
      !and reset transmax to 0.0
      transmax = 0.0

      !2.2 Calculate water storage and depth averaged k
      !Calculate amount of water (expressed as a depth) which may be stored in
      !each layer within each column.
      call layer_water_depth(x_extent, y_extent, no_layers, layer_attributes, &
                        layer_storage, activation_status)

      !Calculate the depth-averaged k below the water table for each column
      call wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)

      !=========================================================================
      !The model now uses a linearly increasing time step which is then
      !fixed after a specified time - see the parameters file.

      !2.3 Set maximum timesteps depending on elapsed model years

      !Calculate the time step
      if (year_counter > t_min_years + t_max_years) then
          timestep = t_max / real(525600)
        else if (year_counter <=  t_min_years) then
          timestep = t_min / real(525600)
        else
          timestep = (((year_counter - t_min_years) * t_rate) + t_min) / &
                        real(525600)
      end if

      !Keep the sum of sub-annual timesteps equal to the net rain interval time
      !step to prevent the over-accumulation of time in the model time extent
      if (timestep + t_step_sum_netr > net_rain_tstep) then
        timestep = net_rain_tstep - t_step_sum_netr
        !For sub-annual weather, set the new net rain read flag
        if (net_rain_res /= 1) then
          new_net_rain = 1
        end if
      !Increment the sub annual counter
      sub_annual_counter = sub_annual_counter + 1
      end if

      !Sum up the value of each timestep to keep track of sub-annual looping
      t_step_sum_netr = t_step_sum_netr + timestep
      !Summer water tables.
      if (net_rain_res /= 1) then
        if (net_rain_res == 3) then
          if (sub_annual_counter >= 22 .and. sub_annual_counter <= 35) then
            t_step_sum_su = t_step_sum_su + timestep
          end if
        else if (net_rain_res == 2) then
          if (sub_annual_counter == 6 .or. sub_annual_counter == 7 .or. &
                  sub_annual_counter == 8) then
              t_step_sum_su = t_step_sum_su + timestep
          end if
        end if
      end if

      !2.4 Move water and update the water-table position
      !Calculate the amount of water (expressed as a depth) that moves between
      !each column. The flow law is based on the Boussinesq equation.
      call move_water(x_extent, y_extent, timestep, spatial_step, rainfall, &
                      base_altitude, water_change, water_table, wk_mean, &
                      activation_status, &
                      pond_olf, height_olf, mineral, no_layers, layer_storage)

      !Update position (height) of the water table in each column
      call water_table_update(x_extent, y_extent, no_layers, water_change, &
                              water_table, layer_attributes, layer_storage, &
                              activation_status)

      !Calculate new water-table depth for each active column and use to
      !calculate mean water-table depth for year. The mean value is used to
      !calculate the new annual layer litter production.
      Water_table_depth: do x = 1, x_extent
        do y = 1, y_extent
          if (activation_status(x, y) == "on") then
            col_wt_depth(x, y) = transmissivity(x, y, (no_layers(x,y) - 1), 1) &
                               - water_table(x, y)
            col_wt_sum(x, y)  = col_wt_sum(x, y) + col_wt_depth(x, y)
            !Summer water tables
            if (net_rain_res /= 1) then
              if (net_rain_res == 3) then
                if (sub_annual_counter >= 22 &
                        .and. sub_annual_counter <= 35) then
                  col_wt_sum_su(x, y) = col_wt_sum_su(x, y) + col_wt_depth(x, y)
                end if
              else if (net_rain_res == 2) then
                if (sub_annual_counter == 6 .or. sub_annual_counter == 7 .or. &
                        sub_annual_counter == 8) then
                  col_wt_sum_su(x, y) = col_wt_sum_su(x, y) + col_wt_depth(x, y)
                end if
              end if
            end if
          end if
        end do
      end do Water_table_depth

      !2.5 Decompose peat in each layer using info. from previous timestep.
      !Use proportionate decay model. Amount of peat lost depends on whether
      !peat above or below the water table. Scan upwards through all model
      !layers. Use new layer properties to update column properties.
      Decompose: do x = 1, x_extent
        do y = 1, y_extent
          if (activation_status(x, y) == "on") then
            !Decompose peat layers only
            do z = 2, (no_layers(x,y) - 1)
              !Update layer mass for each litter pool
              !applied only to oxic and anoxic decomposition
              !For each pft
              do pft = 1, n_pft
                !For each litter fraction
                do litter_frac = 1, n_fractions
                  !Labile pool ====
                  !Oxic decay
                  pft_layer_mass(x, y, z, pft, litter_frac, 1) = &
                    (pft_layer_mass(x, y, z, pft, litter_frac, 1) * &
                      (1 - wet_proportion(x, y, z)) * &
                        (exp (- decay_params(pft, litter_frac) * timestep))) + &
                  !Add anoxic decay
                  (pft_layer_mass(x, y, z, pft, litter_frac, 1) * &
                    wet_proportion(x, y, z) * &
                    (exp (- decay_params(pft, litter_frac + n_fractions)  * &
                          timestep)))

                  !Recalcitrant pool ====
                  !Oxic decay
                  pft_layer_mass(x, y, z, pft, litter_frac, 2) = &
                    (pft_layer_mass(x, y, z, pft, litter_frac, 2) * &
                      (1 - wet_proportion(x, y, z)) * &
                        (exp (- decay_params(pft, litter_frac + &
                        (n_fractions * n_pools)) * &
                        timestep))) + &
                  !Add anoxic decay
                  (pft_layer_mass(x, y, z, pft, litter_frac, 2) * &
                    wet_proportion(x, y, z) * &
                    (exp (- decay_params(pft,litter_frac + n_fractions + &
                            (n_fractions * n_pools))  * timestep)))

                  !Set thickness of each fraction
                  pft_layer_mass(x, y, z, pft, litter_frac, 5) = &
                    sum(pft_layer_mass(x, y, z, pft, litter_frac, (/1,2/))) / &
                      density(pft,litter_frac)

                  !Set remaining mass of each fraction
                  pft_layer_mass(x, y, z, pft, litter_frac, 6) = &
                    sum(pft_layer_mass(x, y, z, pft, litter_frac, (/1,2/))) / &
                      sum(pft_layer_mass(x, y, z, pft, litter_frac, (/3,4/)))
                end do
              end do

              !Sum of layer thicknesses (i.e. all fractions in all pfts)
              layer_attributes(x, y, z, 1) = &
                sum(pft_layer_mass(x, y, z, :, :, 5))

              !Calculate layer total mass
              layer_mass(x, y, z, 1) = &
                sum(pft_layer_mass(x, y, z, :, :, (/1,2/)))

              !Calculate layer total remaining mass
              layer_mass(x, y, z, 3) = layer_mass(x, y, z, 1) / &
                layer_mass(x, y, z, 2)

              !Current layer mass for addition to mass per area array
              col_mass_per_area(x, y) = col_mass_per_area(x, y) + &
                                      layer_mass(x, y, z, 1)
            end do
          end if
        end do
      end do Decompose ! End of decomposition loop

      !Recalculate layer k. This loop has no effect in these 1-D simulations
      !For peat layers (>= z2) only
      Layer_k: do x = 1, x_extent
        do y = 1, y_extent
          if (activation_status(x, y) == "on") then
            do z = 2, (no_layers(x,y) - 1)
            layer_attributes(x, y, z, 2) = k_param_a * (exp (k_param_b &
                                           * layer_mass(x, y, z, 3)))
            end do
          end if
        end do
      end do Layer_k

      !Calculate transmissivity profile for each column
      call trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                        transmissivity, activation_status)

      !Recalculate layer wet proportion
      Layer_wet_prop: do x = 1, x_extent
        do y = 1, y_extent
          if (activation_status(x, y) == "on") then
            do z = 1, (no_layers(x,y) - 1)
              !Wet proportion
              if (z == 1) then
                if (water_table(x, y) >= transmissivity(x, y, z, 1)) then
                    wet_proportion(x, y, 1) = 1.0
                  else
                    wet_proportion(x, y, 1) =  water_table(x, y) / &
                                               layer_attributes(x, y, 1, 1)
                  end if
                else if (water_table(x, y) >= transmissivity(x, y, z, 1)) then
                  wet_proportion(x, y, z) = 1.0
                else if (water_table(x, y) <= transmissivity(x, y, z - 1, 1)) &
                  then
                    wet_proportion(x, y, z) = 0.0
                else
                  wet_proportion(x, y, z) = (water_table(x, y) -&
                  transmissivity(x, y, z - 1, 1)) / layer_attributes(x, y, z, 1)
              end if
            end do
          end if
        end do
      end do Layer_wet_prop

      !2.6 Check if sub annual total has been reached
      if (net_rain_res == 3 .and. sub_annual_counter == 52) then
        exit
      else if (net_rain_res == 2 .and. sub_annual_counter == 12) then
        exit
      else if (net_rain_res == 1 .and. sub_annual_counter == 1) then
        !Reset time step sum for annual weather
        t_step_sum_netr = 0
        exit
      end if

    !End sub-annual loop
    end do Sub_annual_loop !----------------------------------------------------

    !5.Calculate the average annual and summer water-table depths to be used for
    !the next layer litter calculation
    Water_table_depth_av: do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
           col_wt_depth_av(x, y) = col_wt_sum(x, y) / (real(sub_year_counter))
           if (net_rain_res /= 1) then
             col_wt_depth_av_su(x, y) = col_wt_sum_su(x, y) / &
                                      (real(sub_year_counter_su))
           end if
        end if
      end do
    end do Water_table_depth_av

    !Calculate the root depth for the next inputs
    !Root litter is distributed following the new layer calculation
    call rooting(x_extent, y_extent, no_layers, transmissivity, &
      max_root, root_depth, mineral, root_count, root_mask, bg_root_mass, &
      n_pft, new_layer_mass, pft_params, new_m_start, new_m_end, pft_prop, &
      activation_status, roots)

    !6.Calculate the mean timestep in minutes and reset counters
    mean_t_step = (1.0 / real(sub_year_counter)) * 525600
    write(130, '(20F20.8)') mean_t_step
    !Reset mean timestep for the next main loop
    mean_t_step = 0.0
    !Reset sub annual counter for new year
    sub_annual_counter = 0

    !7.Check for model write out
    output_counter = output_counter + 1
    if (output_counter == output_interval) then
      !Write results to file
      do x = 1, x_extent
        do y = 1, y_extent
          if(activation_status (x,y) == "off" &
            .or. activation_status (x,y) == "diri" &
            .or. activation_status (x,y) == "neu") then
            !Column height
            write(070, '(20F20.8)') -999.0
            !Water-table height
            write(080, '(20F20.8)') -999.0
            !Transmissivity
            write(090, '(20F20.8)') -999.0
            !Column mass per area
            write(100, '(20F20.8)') -999.0
            !Column mean water-table depth
            write(110, '(20F20.8)') -999.0
            !Summer mean column water-table depth
            write(160, '(20F20.8)') -999.0
          else
            write(070, '(20F20.8)') &
            transmissivity(x, y, (no_layers(x,y) - 1), 1)
            write(080, '(20F20.8)') water_table(x, y)
            write(090, '(20F20.8)') transmissivity(x, y, no_layers(x,y) - 1, 2)
            write(100, '(20F20.8)') col_mass_per_area(x, y)
            write(110, '(20F20.8)') col_wt_depth_av(x, y)
            write(160, '(20F20.8)') col_wt_depth_av_su(x, y)
          end if
        end do
      end do
    end if

    !8.Annual counting
    !Update counter for annual update
    year_counter = year_counter + 1 ! Counts annual timesteps

    !Exit main annual loop if year counter exceeds model timespan
    if (year_counter > t_extent + t_add) exit

    !If the simulation is a spin-up read, the simulation needs to begin here
    !Reset column water-table sum for next annual loop
    col_wt_sum = 0.0
    col_wt_sum_su = 0.0
    !Reset sub-year counting for the next annual loop
    sub_year_counter = 0
    sub_year_counter_su = 0

    !Add new a layer
    call add_layer(x_extent, y_extent, no_layers, layer_attributes, &
      pond_depth, activation_status)

    !Modify the mass of the new top layer according to the CO2 concentration
    call co2_conc(pft_params, beta_fact, co2_ppm_ref, new_m_start, &
        new_m_end, pft_params_ref, co2_ppm_curr)

    !Specify the above and belowground mass of new layer
    call new_layer(x_extent, y_extent, no_layers, layer_attributes, &
      activation_status, root_depth, root_count, layer_age, t_extent, t_add, &
      layer_mass, pft_layer_mass, n_pft, fract_recal, recal_start, &
      recal_end, new_m_start, new_m_end, pft_prop, leaves, wood, &
      roots, new_layer_mass, year_counter, layers_in_layer, col_wt_depth_av, &
      n_fractions, n_pools, density, bg_root_mass, transmissivity, porosity, &
      col_mass_per_area, k_param_a, k_param_b, pft_params)

    !3.Layer Lumping
    !3.1 Lumping routine
    !(1) a set time has elapsed;
    !(2) a layer is below an exclusion level (set from the top peat layer); and
    !(3) both the current layer and the indexed lumping layer (i.e. next active
    !layer) are less than a set thickness; e.g. if z1 and z2 have been combined
    !and z3 is available for combining, z3 can only be combined with z1.

    !(1) Don't start lumping until after n years.
    if (year_counter > layer_exclusion + 3) then
      lumping: do x = 1, x_extent !x extent do
        do y = 1, y_extent ! y extent do
          lump_index = 2
          if (activation_status(x, y) == "on") then
            do z = 3, (no_layers(x,y) - 1) !z extent do. Ignore the pond.

              !(2) If the current layer is within the excluded layers, and
              !therefore unavailable for lumping, move the lump index up 1 and
              !pass the layer properties of the current layer to the indexed
              !layer. For example, if z1-z3 have been combined into a single
              !layer (i.e. z1) and z4 is not available for lumping, then the
              !properties of z4 are passed to z2 so it becomes the new z2
              !(the original z2 is combined with z1).
              if (z > (no_layers(x,y) - layer_exclusion)) then
                lump_index = lump_index + 1
                !Layer thickness
                !Layer attributes should really come after the layer mass files
                !even though nothing is done in this loop because it is the
                !result of changes to layer mass
                layer_attributes(x,y,lump_index,1) = layer_attributes(x,y,z,1)

                !pft fraction mass and thicknesses
                pft_layer_mass(x, y, lump_index, :, :, :) = &
                  pft_layer_mass(x, y, z, :, :, :)

                !Layer total mass (current, initial and remaining)
                layer_mass(x, y, lump_index, :) = layer_mass(x, y, z, :)

                !Layer age
                layer_age(x, y, lump_index) = layer_age(x, y, z)

                !Pass the number of layers in the current layer to the active
                !Number of layers in layers
                layers_in_layer(x, y, lump_index) = layers_in_layer(x, y, z)

              !(3) If the layer is outside the layer exclusion limit then
              !layers can be lumped if they are thinner than the threshold.
              !If the current layer is thinner than the threshold and the
              !indexed (active) layer (lump index) is thinner than the threshold
              !then combine them.
              else if ((layer_attributes(x,y,z,1) < threshold) .and. &
                 (layer_attributes(x,y,lump_index,1) < threshold)) then

                !Combine layer PFTs mass and thickness
                do pft = 1, n_pft
                  do litter_frac = 1, n_fractions
                  !current mass
                  pft_layer_mass(x, y, lump_index, pft, litter_frac, 1:2) = &
                    pft_layer_mass(x, y, z, pft, litter_frac, 1:2) + &
                    pft_layer_mass(x, y, lump_index, pft, litter_frac, 1:2)

                  !Original mass
                  pft_layer_mass(x, y, lump_index, pft, litter_frac, 3:4) = &
                    pft_layer_mass(x, y, z, pft, litter_frac, 3:4) + &
                    pft_layer_mass(x, y, lump_index, pft, litter_frac, 3:4)

                  !Thickness
                  pft_layer_mass(x, y, lump_index, pft, litter_frac, 5) = &
                    pft_layer_mass(x, y, z, pft, litter_frac, 5) + &
                    pft_layer_mass(x, y, lump_index, pft, litter_frac, 5)

                  !Remaining mass
                  pft_layer_mass(x, y, lump_index, pft, litter_frac, 6) = &
            sum(pft_layer_mass(x, y, lump_index, pft, litter_frac,(/1,2/))) / &
                sum(pft_layer_mass(x, y, lump_index, pft, litter_frac,(/3,4/)))

                  end do
                end do

                !Combine layer aggregate properties ====
                !Combine layer total thickness
                layer_attributes(x,y,lump_index,1) = &
                   sum(pft_layer_mass(x, y, lump_index, :, :,5))

                !Combined layer total current mass
                layer_mass(x, y, lump_index, 1) = &
                  sum(pft_layer_mass(x, y, lump_index, :, :,(/1,2/)))

                !Combined layer initial mass
                layer_mass(x, y, lump_index, 2) = &
                  sum(pft_layer_mass(x, y, lump_index, :, :,(/3,4/)))

                !Combined layer remaining mass
                layer_mass(x, y, lump_index, 3) = &
                  layer_mass(x, y, lump_index, 1) / &
                  layer_mass(x, y, lump_index, 2)

                !Calculate the weighted avg. of the two layers
                weight_z = layer_mass(x,y,z,1) / (layer_mass(x,y,z,1) + &
                  layer_mass(x,y,lump_index,1))

                weight_z_minus = &
                  layer_mass(x,y,lump_index,1) / (layer_mass(x,y,z,1) + &
                  layer_mass(x,y,lump_index,1))

               !Calculate age of layer
               layer_age(x,y,lump_index) = ((weight_z * layer_age(x,y,z)) + &
                        (weight_z_minus * layer_age(x,y,lump_index)))

               !Add number of layers in both layers being combined
               layers_in_layer(x,y,lump_index) = layers_in_layer(x,y,z) + &
                     layers_in_layer(x,y,lump_index)

               write (*,*) "Layer was lumped"

               !If one of the layers is > the lumping threshold then no lumping
               !occur. Pass the properties of the current layer to the active
               !layer as (2) above.
              else
                lump_index = lump_index + 1
                !Layer thickness
                layer_attributes(x, y, lump_index, 1) = &
                  layer_attributes(x, y, z, 1)

                !PFTs
                pft_layer_mass(x, y, lump_index, :, :, :) = &
                  pft_layer_mass(x, y, z, :, :, :)

                !Layer  mass totals
                layer_mass(x, y, lump_index, :) = layer_mass(x, y, z, :)

                !Layer age
                layer_age(x, y, lump_index) = layer_age(x, y, z)

                !Pass the number of layers in the current layer to the active
                !Number of layers in layers
                layers_in_layer(x, y, lump_index) = layers_in_layer(x, y, z)
              end if
            end do ! end z extent do

            lump_index = lump_index + 1

            !Move ponding layer to be above the uppermost active peat layer
            layer_attributes(x, y, lump_index, :) = &
              layer_attributes(x, y, no_layers(x, y), :)

            !Blank other values that remain in the ponding layer
            wet_proportion(x, y, lump_index) = 0.0
            layer_age(x, y, lump_index) = 0.0
            layers_in_layer(x, y, lump_index) = 0
            pft_layer_mass(x, y, lump_index, :, :, :) = 0.0
            layer_mass(x, y, lump_index, 1:3) = 0.0
            layer_storage(x,y,lump_index, 1:2) = 0.0
            transmissivity(x,y,lump_index, 1:2) = 0.0

            !Blank layer properties that are now above the pond
            wet_proportion(x, y, no_layers(x, y) + 1) = 0.0
            layer_age(x, y, no_layers(x, y) + 1) = 0.0
            layers_in_layer(x, y, no_layers(x, y) + 1) = 0
            pft_layer_mass(x, y, no_layers(x, y) + 1, :, :, :) = 0.0
            layer_mass(x, y, no_layers(x, y) + 1, :) = 0.0
            layer_attributes(x, y, no_layers(x, y) + 1, :) = 0.0
            layer_storage(x, y, no_layers(x, y) + 1, :) = 0.0
            transmissivity(x, y, no_layers(x, y) + 1, :) = 0.0
            no_layers(x, y) = lump_index
          end if !end activation status check
        end do! end y extent do loop
      end do lumping !end x extent do loop
    end if

    !3.2 Re-calculate layer k, water storage, depth averaged k and
    !transmissivity. Calculate amount of water (expressed as a depth) which
    !may be stored in each layer within each column.

    !Recalculate layer k
    !For peat layers only
    do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
          do z = 2, (no_layers(x,y) - 1)
            layer_attributes(x, y, z, 2) = k_param_a * (exp (k_param_b &
                                           * layer_mass(x, y, z, 3)))
          end do
        end if
      end do
    end do

    !Water storage in a layer
    call layer_water_depth(x_extent, y_extent, no_layers, layer_attributes, &
                            layer_storage, activation_status)

    !Depth averaged k
    call wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)

    !Transmissivity profile
    call trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                      transmissivity, activation_status)

    !4.Write the number of layers in a layer
    do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
        write (180,*) no_layers(x,y)
        end if
      end do
    end do

    !12. Annual loop management
    !Re-set output counter
    output_counter = 0
    !Re-set the belowground root litter fraction
    root_count = 0
    bg_root_mass= 0
    root_mask = .false.

  !End of main annual loop
  end do Annual_loop
  !-----------------------------------------------------------------------------

  !13.Write array outputs
  !Other files could be added to this loop if required.
  do x = 1, x_extent
    do y = 1, y_extent
      if (activation_status(x, y) == "on") then
        do z = 2, (no_layers(x,y) -1)
          write(060, '(20F20.8)') layer_mass(x, y, z, 3)     !mass remaining
          write(200, '(20F20.8)') layer_mass(x, y, z, 2)     !initial mass
          write(210, '(20F20.8)') layer_attributes(x, y, z, 1) !layer thickness
          write(120, '(20F20.8)') wet_proportion(x, y, z)    !layer wet prop_n
          write(190,*) y, ",", z, ",", layer_mass(x, y, z, 3), ",", &
            layer_attributes(x,y,z,1), ",", transmissivity(x,y,z,1)
          write(220, 1100) x, y, z, layers_in_layer(x, y, z), &!number of layers
            layer_age(x, y, z), & !layer age
            layer_attributes(x, y, z, 1), & ! layer thickness
            layer_mass(x, y, z, 2),  layer_mass(x, y, z, 3) !layer mass
        end do
      end if
    end do
  end do

  if (spin_up == 1) then
    open(unit=510, file=data_file_name_510, status="replace", &
      form = "unformatted")
    open(unit=520, file=data_file_name_520, status="replace", &
      form = "unformatted")
    open(unit=530, file=data_file_name_530, status="replace", &
      form = "unformatted")
    open(unit=540, file=data_file_name_540, status="replace", &
      form = "unformatted")
    open(unit=550, file=data_file_name_550, status="replace", &
      form = "unformatted")
    open(unit=560, file=data_file_name_560, status="replace", &
      form = "unformatted")
    open(unit=570, file=data_file_name_570, status="replace", &
      form = "unformatted")
    open(unit=580, file=data_file_name_580, status="replace", &
      form = "unformatted")
    open(unit=590, file=data_file_name_590, status="replace", &
      form = "unformatted")
    open(unit=600, file=data_file_name_600, status="replace", &
      form = "unformatted")

    write (510) wet_proportion(2:(x_extent - 1),2:(y_extent - 1), 2:)
    write (520) layer_storage(2:(x_extent - 1),2:(y_extent - 1), 2:, :)
    write (530) layer_age(2:(x_extent - 1),2:(y_extent - 1), 2:)
    write (540) layers_in_layer(2:(x_extent - 1),2:(y_extent - 1), 2:)
    write (550) layer_attributes(2:(x_extent - 1),2:(y_extent - 1), 2:, :)
    write (560) layer_mass(2:(x_extent - 1),2:(y_extent - 1), 2:, :)
    write (570) pft_layer_mass(2:(x_extent - 1),2:(y_extent - 1), 2:, :, :, :)
    write (580) water_table(2:(x_extent - 1),2:(y_extent - 1))
    write (590) no_layers(2:(x_extent - 1),2:(y_extent - 1))
    write (600) transmissivity(2:(x_extent - 1),2:(y_extent - 1), 2:, :)
  end if

  !14.Write the model elapsed time
  call cpu_time(time_finish)
  time_elapsed = time_finish - time_start
  write(170, *) time_elapsed

end program DIGIBOG_CONGO
