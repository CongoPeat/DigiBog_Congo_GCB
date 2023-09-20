module hydro_procedures_transmax
use global_definition

  !   DigiBog_Congo. a model to simulate peat accumulation in interfluvial
  !   domed peatlands in the Central Congo Basin.
  !   Copyright (c) 2023  the authors (see the main program file)

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

!-------------------------------------------------------------------------------
! Section 1.0 Program header
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !Modification history of code
  !-----------------------------------------------------------------------------
  !Programmer           Date           Modifications
  !-----------------------------------------------------------------------------
  !Dylan m. Young       18/06/14  This version is for the fixed Dirichlet
  !                               condition model. The activation status in
  !                               trans_height (section 2.0) and wat_k_mean
  !                               (section 4.0) have been removed because the
  !                               'diri' condition is updated using only water-
  !                               table height.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       25/06/14  Transmax scalar added to wat_K_mean to enable
  !                               variable timesteps to be used. Transmax code
  !                               created by Andy Baird.
  !-----------------------------------------------------------------------------
  !Dylan m. Young       23/01/15  Water-table-update (6.0). Added condition to
  !                               set marker = no_layers where the water-table
  !                               exceeds the ponding layer. The layer marker
  !                               is initialised at 0.
  !-----------------------------------------------------------------------------

  !Declarations

  implicit none

  !Global type definition
  !Define a real kind type q with at least 8 decimal digits and an exponent
  !range from 10**30 to 10**(-30)
  !integer, parameter :: q = SELECTED_REAL_KIND(p = 8, r = 30)

  contains

!-------------------------------------------------------------------------------
! Section 2.0   Initial transmissivity calculations
!-------------------------------------------------------------------------------

  subroutine trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                          transmissivity, activation_status)

    !This subroutine calculates the peat column's transmissivity from the
    !impermeable base of each column to the top of each layer


    !Declarations

    implicit none

    !Subroutine arguments

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent

    !Array arguments with intent(in)
    integer, dimension(:,:), intent(in) :: no_layers
    real(kind=q), dimension(:,:,:,:), intent(in) :: layer_attributes
    character(8), dimension(:,:), intent(in) :: activation_status

    !Array arguments with intent(inout)
    real(kind=q), dimension(:,:,:,:), intent(inout) :: transmissivity

    !Local scalars
    integer :: x, y, z

    !-- End of subroutine header -----------------------------------------------

    !Calculations

    do x = 1, x_extent
      do y = 1, y_extent
        !Calculations only needed for active columns and Dirichlet boundary
        !columns
        if(activation_status(x, y)     == "on") then
          transmissivity(x, y ,1, 1) = layer_attributes(x, y, 1, 1)
          transmissivity(x, y, 1, 2) = layer_attributes(x, y, 1, 1) &
                                     * layer_attributes(x, y, 1, 2)
          do z = 2, (no_layers(x,y) - 1)!Ignore ponding layer
            transmissivity(x, y, z, 1) = transmissivity(x, y, z - 1, 1) &
                                       + layer_attributes(x, y, z, 1)
            transmissivity(x, y, z, 2) = transmissivity(x, y, z - 1, 2) &
                                       + (layer_attributes(x, y, z, 1) &
                                       * layer_attributes(x, y, z, 2))
          end do
        end if
      end do
    end do

  end subroutine trans_height


!-------------------------------------------------------------------------------
! Section 3.0   Water capacity initialisation
!-------------------------------------------------------------------------------

  subroutine layer_water_depth (x_extent, y_extent, no_layers, &
                                layer_attributes, layer_storage, &
                                activation_status)

    !This subroutine calculates the amount of water (expressed as a depth)
    !capable of being stored in each layer in each column.

    !Declarations

    implicit none

    !Subroutine arguments

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent

    !Array arguments with intent(in)
    integer, dimension(:,:), intent(in) :: no_layers
    real(kind=q), dimension(:,:,:,:), intent(in) :: layer_attributes
    character(8), dimension(:,:), intent(in) :: activation_status

    !Array arguments with intent(inout)
    real(kind=q), dimension(:,:,:,:), intent(inout) :: layer_storage

    !Local scalars
    integer :: x, y, z

    !-- End of subroutine header -----------------------------------------------

    !Calculations

    do x = 1, x_extent
      do y = 1, y_extent
        !Calculations only needed for active columns
        if (activation_status(x, y) == "on") then
          layer_storage(x, y, 1, 1) = layer_attributes(x, y, 1, 1)
          layer_storage(x, y, 1, 2) = layer_attributes(x, y, 1, 1) &
                                    * layer_attributes(x, y, 1, 3)
          do z = 2, no_layers(x,y) !Ponding layer included
            layer_storage(x, y, z, 1) = layer_storage(x, y, z - 1, 1) &
                                      + layer_attributes(x, y, z ,1)
            layer_storage(x, y, z, 2) = layer_attributes(x, y, z, 1) &
                                      * layer_attributes(x, y, z, 3)
          end do
        end if
      end do
    end do

  end subroutine layer_water_depth


!------------------------------------------------------------------------------
!Section 4.0 Depth-averaged hydraulic conductivity calculations
!------------------------------------------------------------------------------

  subroutine wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)

  !This subroutine calculates the depth-averaged k below the water-table
  !for each column. It also records the value of the transmissivity for
  !that model column with the highest transmissivity below the
  !water-table for use in the variable time-step calculation in the main
  !programme.

  !Declarations

  implicit none

  !Subroutine arguments

  !Scalar arguments with intent(in)
  integer, intent(in) :: x_extent, y_extent

  !Array arguments with intent(in)
  integer, dimension(:,:), intent(in) :: no_layers
  real(kind=q), dimension(:,:), intent(in) :: water_table
  real(kind=q), dimension(:,:,:,:), intent(in) :: layer_attributes, &
                                                  transmissivity
  character(8), dimension(:,:), intent(in) :: activation_status

  !Scalar arguments with intent(inout)
  real(kind=q), intent(inout) :: transmax

  !Array arguments with intent(inout)
  real(kind=q), dimension(:,:), intent(inout) :: wk_mean

  !Local scalars
  integer :: x, y, z

  !-- End of subroutine header-----------------------------------------------

  !Calculations
  do x = 1, x_extent
    do y= 1, y_extent
      !Calculations only needed for active columns
      if (activation_status(x, y) == "on") then
        !condition 1: is the water table at or above the peatland surface
        !(i.e. residing in the pond)?
        !If the water table is at or above the peatland surface then mean k is
        !that of the peat profile below the ponding layer
        if (water_table(x, y) &
            >= transmissivity(x, y, (no_layers(x,y) -1), 1)) then
          wk_mean(x, y) = transmissivity(x, y, (no_layers(x,y) -1), 2)
          if (wk_mean(x, y) > transmax) then
            transmax = wk_mean(x, y)
          end if
          wk_mean(x, y) = wk_mean(x, y) / &
                            transmissivity(x, y, (no_layers(x,y) -1), 1)
        !condition 2: is the water table below the top of the first layer?
        !**note: first layer is mineral soil**
        else if (water_table(x, y) < transmissivity(x, y, 1 , 1)) then
          !Depth avg k = k of first layer
          wk_mean(x, y) = layer_attributes(x, y, 1, 2)
          if (transmissivity(x, y, 1, 2) > transmax) then
            transmax = transmissivity(x, y, 1, 2)
          end if
        !condition 3: is the water table equal to or above the top of the
        !current layer and below the peatland surface?
        else
          do z = 1, (no_layers(x,y) - 2)
            !condition 3.1
            !If the water table is equal to the height of the current layer
            !Testing the equality of real numbers but this is ok here
            if (water_table(x, y) == transmissivity(x, y, z, 1)) then
              wk_mean(x, y) = transmissivity(x, y, z, 2)
              if (wk_mean(x, y) > transmax) then
                transmax = wk_mean(x, y)
              end if
              wk_mean(x, y) = wk_mean(x, y) / transmissivity(x, y, z, 1)
            !condition 3.2: is the water table between two layers
            !If the water table is above the top of layer z
            else if (water_table(x, y) > transmissivity(x, y, z, 1)) then
            !If the water table is below the top of z + 1
              if (water_table(x, y) < transmissivity(x, y, z + 1, 1)) then
              wk_mean(x, y) = (transmissivity(x, y, z, 2) &
                            + (layer_attributes(x, y, z + 1, 2) &
                            * (water_table(x, y) &
                            - transmissivity(x, y, z, 1))))
                if (wk_mean(x, y) > transmax) then
                transmax = wk_mean(x, y)
                end if
              wk_mean(x, y) = wk_mean(x, y) / water_table(x, y)
              end if
            end if
          end do
        end if
      end if
    end do
  end do

  end subroutine wat_k_mean


!-------------------------------------------------------------------------------
! Section 5.0   2-dimensional flux calculations
!-------------------------------------------------------------------------------

  subroutine move_water(x_extent, y_extent, timestep, spatial_step, rainfall, &
                        base_altitude, water_change, water_table, wk_mean, &
                        activation_status, &
                        pond_olf, height_olf, mineral, no_layers, layer_storage)

    !This subroutine calculates the net movement of water (expressed as a depth)
    !between columns. It is assumed that a harmonic mean operates between the
    !mean (depth-averaged) k of each column. The output is a change in depth
    !(positive or negative) of water (volume per unit area)to be added to the
    !previous water-table elevation.

    !Declarations

    implicit none

    !Subroutine arguments

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent
    real(kind=q), intent(in) :: spatial_step, timestep, rainfall, &
                                pond_olf, height_olf, mineral

    !Array arguments with intent(in)
    real(kind=q), dimension(:,:), intent(in) :: base_altitude, water_table, &
                                                wk_mean

    character(8), dimension(:,:), intent(in) :: activation_status

    !olf
    integer, dimension(:,:), intent(in) :: no_layers
    real(kind=q), dimension(:,:,:,:), intent(in) :: layer_storage

    !Array arguments with intent(inout)
    real(kind=q), dimension(:,:), intent(inout) :: water_change

    !Local scalars
    integer :: x, y
    real :: olf

    !Local arrays
    real(kind=q), dimension(x_extent, y_extent) :: x_flux, y_flux


    !-- End of subroutine header -----------------------------------------------


    !Calculations

    !Volume 'flux' (in x direction) between columns using Dupuit-Forchheimer
    !approximation
    !16.09.2013 Changed way in which flow to and from Diri cells works
    do x = 1, (x_extent - 1)
      do y = 2, (y_extent - 1)
        !Different rules apply to cells with different statuses
        !Case 1
        if (activation_status(x, y) == "neu" &
           .and. activation_status(x + 1, y) == "on") then
          x_flux(x, y) = 0.0
        !Case 2
        else if (activation_status(x, y) == "diri" &
                 .and. activation_status(x + 1, y) == "on") then
          x_flux(x, y) = wk_mean(x + 1, y) &
                       * ((water_table(x, y) + water_table(x + 1, y)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x + 1, y)  - water_table(x + 1, y)) &
                       * 2.0 * timestep
        !Case 3
        else if (activation_status(x, y) == "on" &
                 .and. activation_status(x + 1, y) == "neu") then
          x_flux(x, y) = 0.0
        !Case 4
        else if (activation_status(x, y) == "on" &
                 .and. activation_status(x + 1, y) == "diri") then
          x_flux(x, y) = wk_mean(x, y) &
                       * ((water_table(x, y) + water_table(x + 1, y)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x + 1, y)  - water_table(x + 1, y)) &
                       * 2.0 * timestep
        !Case 5
        else if (activation_status(x, y) == "on" &
                 .and. activation_status(x + 1, y) == "on") then
          x_flux(x, y) = (2 * wk_mean(x, y) * wk_mean(x + 1, y)) &
                       / (wk_mean(x, y) + wk_mean(x + 1, y)) &
                       * ((water_table(x, y) + water_table(x + 1, y)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x + 1, y)  - water_table(x + 1, y)) &
                       * timestep
        end if
      end do
    end do

    !Volume 'flux' (in y direction) between columns using Dupuit-Forchheimer
    !approximation
    do y = 1, (y_extent - 1)
      do x = 2, (x_extent - 1)
        !Different rules apply to cells with different statuses
        !Case 1
        if (activation_status(x, y) == "neu" &
           .and. activation_status(x, y + 1) == "on") then
          y_flux(x, y) = 0.0
        !Case 2
        else if (activation_status(x, y) == "diri" &
                 .and. activation_status(x, y + 1) == "on") then
          y_flux(x, y) = wk_mean(x, y + 1) &
                       * ((water_table(x, y) + water_table(x, y + 1)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x, y + 1)  - water_table(x, y + 1)) &
                       * 2.0 * timestep
        !Case 3
        else if (activation_status(x, y) == "on" &
                 .and. activation_status(x, y + 1) == "neu") then
          y_flux(x, y) = 0.0
        !Case 4
        else if (activation_status(x, y) == "on" &
                 .and. activation_status(x, y + 1) == "diri") then
          y_flux(x, y) = wk_mean(x, y) &
                       * ((water_table(x, y) + water_table(x, y + 1)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x, y + 1)  - water_table(x, y + 1)) &
                       * 2.0 * timestep
        !Case 5
        else if (activation_status(x, y) == "on" &
                 .and. activation_status(x, y + 1) == "on") then
          y_flux(x, y) = (2 * wk_mean(x, y) * wk_mean(x, y + 1)) &
                       / (wk_mean(x, y) + wk_mean(x, y + 1)) &
                       * ((water_table(x, y) + water_table(x, y + 1)) / 2.0) &
                       * (base_altitude(x, y) + water_table(x, y) &
                       - base_altitude(x, y + 1)  - water_table(x, y + 1)) &
                       * timestep
        end if
      end do
    end do

    !Convert volume into depth of water (in absence of soil matrix - i.e. the
    !depth of water if the drainable porosity were 1)
    !Ignore all edge cells (they have to be boundary conditions)
    do x = 2, (x_extent - 1)
      do y = 2, (y_extent - 1)
        if (activation_status(x, y) == "on") then
          !Is the water table at or below the peatland surface? If yes, then
          !there is no olf.
          if (water_table(x, y) <= layer_storage(x, y, no_layers(x, y) -1, 1)) &
            then
            water_change(x, y) = (x_flux(x - 1, y) - x_flux(x ,y) &
                                 +  y_flux(x, y - 1) - y_flux(x ,y)) &
                                 / (spatial_step ** 2)
            water_change(x, y) = water_change(x, y) + (rainfall * timestep)
          !If the water table resides in the pond, calculate olf.
          else
            olf = &
              !water-table depth x a
              ((water_table (x,y) - layer_storage(x,y,no_layers(x,y) -1, 1)) &
                * pond_olf) + &
              !peat thickness (excluding mineral) x b
              ((layer_storage(x,y,no_layers(x,y) -1, 1) - mineral) * height_olf)
            !Is olf greater than the water in pond?
            if (olf * timestep > &
             (water_table (x,y) - layer_storage(x,y,no_layers(x,y) -1, 1))) then
              olf = &
                water_table (x,y) - layer_storage(x,y,no_layers(x,y) -1, 1)
            else
              olf = olf * timestep
            end if
            !Add olf to flux
            water_change(x, y) = (x_flux(x - 1, y) - x_flux(x ,y) &
                                 +  y_flux(x, y - 1) - y_flux(x ,y)) &
                                 / (spatial_step ** 2)
            water_change(x, y) = water_change(x, y) + (rainfall * timestep) &
                                 !remove olf
                                 - olf
          end if
        end if
      end do
    end do

  end subroutine move_water


!-------------------------------------------------------------------------------
! Section 6.0   Water-table update
!-------------------------------------------------------------------------------

  subroutine water_table_update(x_extent, y_extent, no_layers, water_change, &
                                water_table, layer_attributes, layer_storage, &
                                activation_status)

    !This subroutine updates the position of the water table based on the
    !storage available in layers above or below the water table

    !Declarations

    implicit none

    !Subroutine arguments

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent

    !Array arguments with intent(in)
    integer, dimension(:,:), intent(in) :: no_layers
    real(kind=q), dimension(:,:,:,:), intent(in) :: layer_attributes, &
                                                    layer_storage
    character(8), dimension(:,:), intent(in) :: activation_status

    !Array arguments with intent(inout):
    real(kind=q), dimension(:,:), intent(inout) :: water_change, water_table

    !Local scalars
    integer :: x, y, z, marker
    real(kind=q) :: test_depth

    !-- End of subroutine header -----------------------------------------------

    !Calculations

    !Ignore all edge cells (they have to be boundary conditions)
    do x = 2, (x_extent - 1)
      do y = 2, (y_extent - 1)
        if (activation_status(x, y) == "on") then
          !Locate layer in which water table resides
          !Is the water table below the top of the first layer?
          if (water_table(x, y) <= layer_storage(x, y, 1, 1)) then
            marker = 1
          else
            do z = 1, (no_layers(x,y) - 1)
              !Is the water table above the top of the current layer?
              if (water_table(x, y) > layer_storage(x, y, z, 1)) then
                !Is the water table in the layer above (z+1)?
                if (water_table(x, y) <= layer_storage(x, y, z + 1, 1)) then
                  marker = z + 1
                !Is the water table above the surface of the pond?
                else if (water_table(x, y) > &
                  layer_storage(x, y, no_layers(x,y), 1)) then
                    marker = no_layers(x,y)
                  exit
                end if
              end if
            end do
          end if
          !Does water table rise?
          if (water_change(x, y) > 0.0) then
            !Ignore any water-table updates if water level already at top
            !of uppermost layer
            if (water_table(x, y) < layer_storage(x, y, no_layers(x,y), 1)) then
              z = marker
              !Is there any storage available in the current layer?
              !Note for debugging the following == is ok
              if (water_table(x, y) == layer_storage(x, y, z, 1)) then
                !No storage available - advance to next layer
                z = z + 1
              end if
              do
                test_depth = layer_storage(x, y, z, 2) &
                           * (layer_storage(x, y, z, 1) - water_table(x, y)) &
                           / layer_attributes(x, y, z, 1)
                if (water_change(x, y) <= test_depth) then
                  !Water table moves up within layer z
                  water_table(x, y) = water_table(x, y) + (water_change(x, y) &
                                    / layer_attributes(x, y, z, 3))
                  exit
                !Water table moves up to next layer or rises to the surface
                else
                  !Has top of ponding layer been reached?
                  if (z == no_layers(x,y)) then
                    water_table(x, y) = layer_storage(x, y, z, 1)
                    exit
                  else
                    water_table(x, y) = layer_storage(x, y, z, 1)
                    water_change(x, y) = water_change(x, y) - test_depth
                    z = z + 1
                  end if
                end if
              end do
            end if
          end if
          !Does the water table fall?
          if (water_change(x, y) < 0.0) then
            !Ignore any water-table updates if water table already at or
            !below base of column
            if (water_table(x, y) > 0.0 ) then
              z = marker
              do
                if (z == 1) then
                  test_depth = layer_storage(x, y, 1, 2) * water_table(x ,y) &
                             / layer_attributes(x, y, 1, 1)
                else
                  test_depth = layer_storage(x, y, z, 2) &
                             * (water_table(x, y) &
                             - layer_storage(x, y, z - 1, 1)) &
                             / layer_attributes(x, y, z, 1)
                end if
                if (abs(water_change(x, y)) <= test_depth) then
                  !Water table falls within layer z or to top of layer z - 1
                  water_table(x, y) = water_table(x, y) + (water_change(x, y) &
                                    / layer_attributes(x, y, z, 3) )
                  exit
                else !Water table falls to top of layer z - 1 (if it exists)
                     !and needs to fall further, provided a lower layer exits
                  !Has base of column been reached/exceeded?
                  if (z == 1) then
                    !If water table drops below the base, set it to half the
                    !height of the mineral layer.
                    water_table(x,y) = layer_attributes(x, y, 1, 1) / 2
                  else
                    water_table(x, y) = layer_storage(x, y, z - 1, 1)
                    water_change(x, y) = water_change(x, y) + test_depth
                    z = z - 1
                  end if
                end if
              end do
            end if
          end if
        end if
      end do
    end do

  end subroutine water_table_update

end module hydro_procedures_transmax
