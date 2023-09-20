module new_litter_procedures
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

  !For details of code ownership, model updates etc, see section 1.0 of
  !DigiBog_Hydro

  !-----------------------------------------------------------------------------
  !Modification history of code
  !-----------------------------------------------------------------------------
  !Programmer           Date           Modifications
  !-----------------------------------------------------------------------------
  !Dylan m. Young       24/01/23  Includes sub-routines needed to run the
  !                               model from the end of a spin-up.
  !-----------------------------------------------------------------------------

  !Declarations

  implicit none

  contains

!-------------------------------------------------------------------------------
! Section 2.0   Add a new layer
!------------------------------------------------------------------------------

  subroutine add_layer(x_extent, y_extent, no_layers, layer_attributes, &
    pond_depth, activation_status)

    !This subroutine adds a new layer at the top of a peat column and sets the
    !parameters of the ponding layer.

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent

    !Variable argumnets with intent(in)
    real(kind=q), intent(in) :: pond_depth

    !Array arguments with intent(in)
    character(8), dimension(:,:), intent(in) :: activation_status

    !Array arguments with intent(inout)
    integer, dimension(:,:), intent(inout) :: no_layers
    real(kind=q), dimension(:,:,:,:), intent(inout) :: layer_attributes

    !Local scalars
    integer :: x, y


    !Calculations
    !Add new layer for each additional year current layer counter +1 for pond
    do x = 1, x_extent
      do y = 1, y_extent
        no_layers(x,y) = no_layers(x,y) + 1
      end do
    end do

    !Initialise above-ground storage layer properties for new ponding layer
    !for active columns only
    Pond_props: do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
         !Initialise the properties of the pond
          layer_attributes(x, y, no_layers(x,y), 1) = pond_depth
          layer_attributes(x, y, no_layers(x,y), 2) = 0.0 !k
          layer_attributes(x, y, no_layers(x,y), 3) = 1.0 !s
        end if
      end do
    end do Pond_props

  end subroutine add_layer

!-------------------------------------------------------------------------------
! Section 3.0  Modify the litter fraction mass according to CO2 concentration
!------------------------------------------------------------------------------
  subroutine co2_conc(pft_params, beta_fact, co2_ppm_ref, new_m_start, &
      new_m_end, pft_params_ref, co2_ppm_curr)

  !This subroutine alters the mass of pft litter fractions according to
  !atmospheric CO2 concentrations in ppm.

    !Arguments with intent(in)
    integer, intent(in) :: new_m_start, new_m_end
    real(kind=q), intent(in) :: beta_fact, co2_ppm_ref
    real(kind=q), dimension(:,:), intent(in) :: pft_params_ref
    !Arguments with intent(inout)
    real(kind=q), dimension(:,:), intent(inout) :: pft_params
    !Local
    real(kind=q) :: co2_ppm_curr

    !Read in the current ppm values
    read(035, *) co2_ppm_curr
    !Resort to original pft parameters
    pft_params = pft_params_ref
    !Calculate new initial mass based on current ppm
    pft_params(:,new_m_start:new_m_end) = &
      pft_params(:,new_m_start:new_m_end) + &
        (pft_params(:,new_m_start:new_m_end) * &
          !beta is the proportion per 100 ppm. Divide by 100 for ppm-1
          ((beta_fact / 100) * (co2_ppm_curr - co2_ppm_ref)))

  end subroutine co2_conc

!-------------------------------------------------------------------------------
! Section 4.0   Set the values for the new layer
!------------------------------------------------------------------------------

  subroutine new_layer(x_extent, y_extent, no_layers, layer_attributes, &
    activation_status, root_depth, root_count, layer_age, t_extent, t_add, &
    layer_mass, pft_layer_mass, n_pft, fract_recal, recal_start, &
    recal_end, new_m_start, new_m_end, pft_prop, leaves, wood, &
    roots, new_layer_mass, year_counter, layers_in_layer, col_wt_depth_av, &
    n_fractions, n_pools, density, bg_root_mass, transmissivity, porosity, &
    col_mass_per_area, k_param_a, k_param_b, pft_params)

    !This subroutine adds a new layer at the top of a peat column and sets the
    !parameters of the ponding layer.

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent, t_add, t_extent, n_pft, &
                           new_m_start, new_m_end, recal_start, &
                           recal_end, leaves, wood, roots, pft_prop, &
                           n_pools, n_fractions, year_counter

    real(kind=q), intent(in) :: porosity, k_param_a, k_param_b

    !Variable argumnets with intent(inout)
    real(kind=q), intent(inout) :: root_depth

    !Array arguments with intent(in)
    real(kind=q), dimension(:,:), intent(in) :: density, pft_params, &
      col_wt_depth_av
    real(kind=q), dimension(:,:,:), intent(in) :: bg_root_mass
    integer, dimension(:,:), intent(in) :: no_layers, root_count
    character(8), dimension(:,:), intent(in) :: activation_status

    !Array arguments with intent(inout)
    real(kind=q), dimension(:), intent(inout) :: fract_recal, new_layer_mass
    real(kind=q), dimension(:,:), intent(inout) :: col_mass_per_area
    real(kind=q), dimension(:,:,:), intent(inout) :: layer_age
    integer, dimension(:,:,:), intent(inout) :: layers_in_layer
    real(kind=q), dimension(:,:,:,:), intent(inout) :: layer_attributes, &
     layer_mass, transmissivity
    real(kind=q), dimension(:,:,:,:,:,:), intent(inout) :: pft_layer_mass

    !Local scalars
    integer :: x = 1, y = 1, z = 1, peat_top = 1, pft = 1, litter_frac = 0
    real(kind=q) :: age_orig = 0, mass_orig = 0, weight_orig = 0, weight_new = 0

    !Calculations
    !Add the new litter fractions
    !11.Calculate mass of new annual litter layer using the gem plots data.
    New_layer_props: do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then

          !Set the peat surface for the column
          peat_top = no_layers(x, y) - 1

          !Update the root-depth for each column.
          !Calculate the thickness of all layers within the rooting depth.
          !These layers are between the layer below peat top (peat top - 1) and
          !the root count + 1 (because the root count was done before the new
          !later was added, e.g. what was z5 is now z6).

          root_depth = &
            sum(layer_attributes(x, y ,(peat_top - root_count(x, y)): &
                                      peat_top - 1, 1))

          !Layer age of most recent layer
          layer_age(x, y, no_layers(x, y) -1) = (t_extent + t_add) - &
            year_counter + 1
          !layers in layer
          layers_in_layer(x, y, no_layers(x, y) -1) = 1

          !1. If water table is deep, then all litter is added to the surface
          !Is average water table > x cm below surface? If so, assign v. low
          !productivity to new annual layer.
          !Note: this doesn't apply to this version of DB_C because litter input
          !does not vary with water-table depth. However, it remains in the
          !model code in case we use it to simulate other peatland types where
          !it is relevant. In which case, the hardcoded 200 cm water table
          !depth would need to be reviewed.

          !For all fractions (including roots).
          if (col_wt_depth_av(x, y) > 200.0) then
            !1E-7 (the original db value) divided over the number of
            !fractions x pools
            pft_layer_mass(x, y, peat_top, :, :, 1:4) = &
              1E-7 / (n_fractions * n_pools)

            !Total current and original mass
            layer_mass(x, y, peat_top, 1:2) = &
            !current and original mass are the same in a new layer so just
            !add the current mass for labile and recalcitrant pools and set
            !the sum as the total current and original mass for the layer.
            sum(pft_layer_mass(x, y, peat_top, :, :, (/1,2/)))

          !Or 2. If the water table greater than x cm
          !(i.e. within the 'normal' range).
          else
            !For just the leaves and wood litter fractions of the new surface
            !layer
            do pft = 1, n_pft
              !Set the recalcitrant proportions for the pft
              fract_recal = pft_params(pft,recal_start:recal_end)
              !Set the new litter mass for each fraction
              new_layer_mass = pft_params(pft,new_m_start:new_m_end) * &
                            pft_params(pft, pft_prop)
              do litter_frac = leaves, wood
                !Current mass and original mass of the two pools
                !Labile pool
                pft_layer_mass(x, y, peat_top, pft, litter_frac, (/1,3/)) = &
                  new_layer_mass(litter_frac) * (1 - fract_recal(litter_frac))
                !Recalcitrant pool
                pft_layer_mass(x, y, peat_top, pft, litter_frac, (/2,4/)) = &
                  new_layer_mass(litter_frac) * fract_recal(litter_frac)
                !Set thickness of both fractions
                pft_layer_mass(x, y, peat_top, pft, litter_frac, 5) = &
                  sum(pft_layer_mass(x,y,peat_top,pft,litter_frac, (/1,2/))) / &
                    density(pft,litter_frac)
                !Set remaining mass of both fractions (= 100%)
                pft_layer_mass(x, y, peat_top, pft, litter_frac, 6) = 1.0
              end do
            end do

            !Belowground roots for each pft
            !The rooting depth has been calculated before the new layer is
            !added. In this version of the model, roots are distriuted
            !according to the thickness of each layer withing the rooting
            !depth. Therefore, the total new root mass is multiplied by the
            !thickness of a layer divided by the root depth times the pool
            !fraction. This proportion of roots is added to the layer in
            !question.

            !bg_root_mass is calculated after the root_layers routine, which
            !counts the layers in the rooting zone.

            !Add a small amount of root mass to the surface layer for both
            !pools
            pft_layer_mass(x, y, peat_top, :, roots, 1:4) = (1E-8 / n_pft)

            !Distribute the remaining roots to the layers below the surface
            !layer (peat_top) and within the rooting depth.
            do z = peat_top - root_count(x, y), peat_top - 1
              !Save the age of the layer before the roots are added
              age_orig = layer_age(x, y, z)
              !Save the current mass of the layer
              mass_orig = layer_mass(x, y, z, 1)
              !Calculate the new mass for each litter fraction
              do pft = 1, n_pft
                !Set the recalcitrant proportions for the pft
                fract_recal = pft_params(pft, recal_start:recal_end)

                !Current mass and original mass of the two pools
                !Add the new root mass to the existing root mass
                !Labile pool ====
                !Add to current mass
                pft_layer_mass(x, y, z, pft, roots, 1) = &
                  pft_layer_mass(x, y, z, pft, roots, 1) + &
                    (bg_root_mass(x, y, pft) * ((layer_attributes(x, y, z, 1) /&
                    root_depth)) * (1 - fract_recal(roots)))
                !Add to original mass
                pft_layer_mass(x, y, z, pft, roots, 3) = &
                  pft_layer_mass(x, y, z, pft, roots, 3) + &
                    (bg_root_mass(x, y, pft) * ((layer_attributes(x, y, z, 1) /&
                    root_depth)) * (1 - fract_recal(roots)))

                !Recalcitrant pool ====
                !Add to current mass
                pft_layer_mass(x, y, z, pft, roots, 2) = &
                  pft_layer_mass(x, y, z, pft, roots, 2) + &
                    (bg_root_mass(x, y, pft) * ((layer_attributes(x, y, z, 1) /&
                    root_depth)) * fract_recal(roots))
                !Add to original mass
                pft_layer_mass(x, y, z, pft, roots, 4) = &
                  pft_layer_mass(x, y, z, pft, roots, 4) + &
                    (bg_root_mass(x, y, pft) * ((layer_attributes(x, y, z, 1) /&
                    root_depth)) * fract_recal(roots))

                !Recalculate the thickness of root fraction
                pft_layer_mass(x, y, z, pft, roots, 5) = &
                  sum(pft_layer_mass(x, y, z, pft, roots, (/1,2/))) / &
                    density(pft, roots)

                !Recalculate the remaining mass of root fraction
                pft_layer_mass(x, y, z, pft, roots, 6) = &
                  sum(pft_layer_mass(x, y, z, pft, roots, (/1,2/))) / &
                    sum(pft_layer_mass(x, y, z, pft, roots, (/3,4/)))

              end do !End pft loop
            end do !End of layers loop for new layer

            !Recalculation of thickness and elevation
            !rooting depth.
            do z = 2, peat_top ! z = 1 is the mineral soil
              !Recalculate total thickness of layers
              layer_attributes(x, y, z, 1) = &
                 sum(pft_layer_mass(x, y, z, :, :, 5))

              !Recalculate layer elevation
              transmissivity(x, y, z, 1) = &
                !layer elevation below layer + thickness of current layer
                transmissivity(x, y, z - 1, 1) + &
                  layer_attributes(x, y, z, 1)
            end do
          end if !End adding new surface layer and belowground roots

          !For whole new top layer ====
          !Layer total mass
          !Add the current mass for labile and recalcitrant pools for each
          !fraction, and set both current and original mass for total mass
          !array.
          layer_mass(x, y, peat_top, 1:2) = &
              sum(pft_layer_mass(x, y, peat_top, :, :, (/1,2/)))

          !Layer remaining mass of fractions in the new top layer = 100%
          !Remaining mass for whole layer
          layer_mass(x, y, peat_top, 3) = 1.0

          !Layer mass for the layers below the top because root mass has
          !been added (prevents remaining mass > 1). Overall thickness and
          !elevation have already been recalculated based on the pft mass.
          do z = 2, peat_top - 1
                  layer_mass(x, y, z, 1) = &
                    sum(pft_layer_mass(x, y, z, :, :,(/1,2/)))
                  layer_mass(x, y, z, 2) = &
                    sum(pft_layer_mass(x, y, z, :, :,(/3,4/)))
              layer_mass(x, y, z, 3) = &
                layer_mass(x, y, z, 1) / layer_mass(x, y, z, 2)
          end do

          !Calculate new age of layers. Weight the age depending on old and new
          !mass
          do z = peat_top - root_count(x, y), peat_top - 1
            weight_orig = mass_orig / layer_mass(x, y, z, 1)
            weight_new = (layer_mass(x, y, z, 1) - &
              mass_orig) / (layer_mass(x, y, z, 1))
            !The age of new material is 1. And the weights should add to 1
            layer_age(x, y, z) = (age_orig * weight_orig ) + (1 * weight_new)
          end do

          !Layer drainable porosity
          layer_attributes(x, y, peat_top, 3) = porosity
          !Layer k (recalculated after lumping so could be removed from here).
          layer_attributes(x, y, peat_top, 2) =  k_param_a * (exp (k_param_b))

          !Calculate transmissivity with new layer
          transmissivity(x, y, peat_top, 2) = &
              transmissivity(x, y, peat_top - 1, 2) &
              + (layer_attributes(x, y, peat_top, 1) &
              * layer_attributes(x, y, peat_top, 2))

          col_mass_per_area(x, y) = col_mass_per_area(x, y) + &
            layer_mass(x, y, peat_top, 1)

        end if !End of activation status check
      end do !End of y dim loop
    end do New_layer_props ! End of x dim loop

  end subroutine new_layer

!-------------------------------------------------------------------------------
! Section 4.0 Rooting
!------------------------------------------------------------------------------
!Belowground root litter ===================================================
    !This section determines the rooting depth and counts how many layers are
    !within it.

    !Set rooting depth for year.
    !In this version, roots are distributed downcore based on the
    !following condition.
    !1. Is the peat thickness < max_root? If it is, the rooting depth is equal
        !to the peat thickness.
    !2. Is the peat thickness > max_root? If so, the rooting depth is equal to
        !max_root

  subroutine rooting(x_extent, y_extent, no_layers, transmissivity, &
    max_root, root_depth, mineral, root_count, root_mask, bg_root_mass, &
     n_pft, new_layer_mass, pft_params, new_m_start, new_m_end, pft_prop, &
    activation_status, roots)

    !Scalar arguments with intent(in)
    integer, intent(in) :: x_extent, y_extent, n_pft, pft_prop, new_m_start, &
                           new_m_end, roots

    !Variable argumnets with intent(in)
    real(kind=q), intent(in) :: mineral, max_root
    real(kind=q), intent(inout) :: root_depth

    !Logical argument with intent(inout)
    logical, dimension(:, :, :) :: root_mask

    !Array arguments with intent(in)
    real(kind=q), dimension(:,:), intent(in) :: pft_params
    character(8), dimension(:,:), intent(in) :: activation_status
    real(kind=q), dimension(:,:,:,:), intent(in) :: transmissivity

    !Array arguments with intent(inout)
    real(kind=q), dimension(:), intent(inout) :: new_layer_mass
    integer, dimension(:,:), intent(inout) :: no_layers, root_count
    real(kind=q), dimension(:,:,:), intent(inout) :: bg_root_mass

    !Local scalars
    integer :: x, y, z, pft

    !Calculations
    do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
          if ((transmissivity(x,y, (no_layers(x,y)) - 1, 1) - mineral) &
                  < max_root) then
            root_depth = &
                    transmissivity(x, y, (no_layers(x, y)) - 1, 1) - mineral
          else
            root_depth = max_root
          end if
        end if
      end do
    end do

   !Count the layers in the root depth
    do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
          !Allows rooting into the first peat layer, but not the mineral soil.
          !The routine uses the top of the layer below the current to determine
          !if the current is within the rooting depth. i.e. the top of the
          !column minus the elevation of the layer below current.
          do z = 2, (no_layers(x,y) -1)
            root_mask(x, y, z) = &
                    (transmissivity(x, y, (no_layers(x,y) - 1), 1) - &
                    transmissivity(x, y, z - 1, 1)) <= root_depth
          end do
          !Count the layers
          root_count(x, y) = count(root_mask(x, y, :))
        end if
      end do
    end do

    !Save the fraction of root mass for downcore distribution (minus the
    !small mass added to the surface layer).
    !Root litter is distributed following the new layer calculation
    root_mass: do x = 1, x_extent
      do y = 1, y_extent
        if (activation_status(x, y) == "on") then
          do pft = 1, n_pft
          !Set the new root mass for each fraction
          new_layer_mass = pft_params(pft,new_m_start:new_m_end) * &
                             pft_params(pft, pft_prop)
            !Remove the very small amount of litter added to the surface layer
            !for each pft.
            bg_root_mass(x, y, pft) = (new_layer_mass(roots) - (1E-8 / n_pft))
          end do
        end if
      end do
    end do root_mass

  end subroutine rooting

end module new_litter_procedures
