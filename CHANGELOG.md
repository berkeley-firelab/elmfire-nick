# Changelog

We constantly update and upgrade ELMFIRE based on user feedback and new cases / requirements.

We want to expand the list of errors that can be caught early in the simulation.

Please create issues with your errors, requests or observations so they can be worked on for the next release.

## ELMFIRE-memopt.0331

### Bug Fixes

- Fixed a bug when the NODATA_VALUE of rasters was a positive value
- Reworked the BSQ metadata parser to find its target values more flexibly and not rely on fixed character ranges (caused issues depending on how the input tifs were converted to BSQ)
- Fixed a bug where the final ROS value in the Fire Potential mode (MODE = 2) was incorrectly calculated.
- Fixed a bug where user-specified firebrand inputs were overwritten in the code back to default.
- Fixed a bug for flame ellipse in building fire spread where nonuniform building input rasters is introduced and there are NODATA_VALUE in building size and structure separation distance rasters to have default values.
- Fixed a bug when structure separation distance inout raster has too high value (> 50m) to have default maximum value (50m) since the coarser input resolution is 30m so far.

### Removed

- Removed the COMPUTATIONAL_DOMAIN input namelist. The domain XLLCORNER, YLLCORNER, CELLSIZE and A_SRS values are set by ELMFIRE via GDAL commands.
- Removed the LATITUDE and LONGITUDE inputs required for diurnal fire adjustments. ELMFIRE calculates those automatically from the domain center dimensions.

### Added

- More examples: Additional examples have been included in the new examples/ folder.
- Verification scripts: In addition to the old verification scripts, a new set of verification cases has been added. These cases are described in detail in the guide. The cases can be ran automatically by copying the TEMPLATE folder, renaming it, and running each case within, with the results described here.
- Automated Validation: We created a brand new pipeline for generating historic validation cases automatically, and validating fire spread models on a large scale. The pipeline is WildfireAV (https://github.com/nick-cloudfire/WildfireAV
), and has been setup to be used by any wildfire spread model. The Validation section of this guide outlines the results of the Validation study. A research paper will be released soon comparing ELMFIRE and FARSITE in terms of their validation statistics, to give relative magnitude to the values in this section. ELMFIRE and FARSITE were shown to have near-matching overall similarity coefficients.
- Rolling Memory Buffer: Larger simulations (big domains, long simulation times) would require lots of RAM to run. ELMFIRE now only loads some of the data at a time in memory to reduce RAM usage by up to 80%, with a runtime penalty of 5% to 20%. Check WX_BANDS_KEPT_IN_MEM in the &SIMULATOR namelist for details.
- Barrier breaching: ELMFIRE can now simulate the effect of fire spread barriers like fuel breaks, roads and rivers. Specify USE_BARRIERS = .TRUE. in the &INPUTS namelist and then set BARRIER_FILENAME to the name of the barrier geotiff. The way this new file works is that, if a barrier exists in a cell, that cell contains the width value of the barrier. See the examples folder for more details.
- CFFDRS implementation: The Rothermel model is the standard fire rate of spread estimation model. To expand the potential use of ELMFIRE into Canada, the CFFDRS model (FBP) has been added, to calculate the head ROS values of the FBP fuels, along with the Build-up effect and the CFFDRS crown involvement algorithm. To run ELMFIRE with the CFFDRS model, set SURFACE_SPREAD_MODEL = "CFFDRS", START_DC = (starting Drought Code value), START_DMC = (starting Duff Moisture Code value), all in the &INPUTS namelist. To test the CFFDRS implementation, a set of 20 verification cases have been included, taken from the 2009 FBP model updates, along with automated scripts to run the cases and compare their outputs with the correct solutions. The Dogrib fire has also been included under Validation to compare the ELMFIRE implementation of the Dogrib fire to the Prometheus/WISE validation case (the standard CFFDRS fire spread model).
- Output options: Some options have been added to the &OUTPUTS namelist. DUMP_CRITICAL_FLIN outputs the critical fireline intensity for canopy fire ignition. DUMP_EMITIMES creates the emitimes.txt file required for a HYSPLIT smoke simulation. DUMP_SPREAD_DIRECTION outputs the direction of primary fire spread in degrees. SPREAD_RATE_IN_M converts the ROS output values from feet per minute to meters per minute.
- Input file checking: some input checks have been implemented at the start of a simulation. The simulation will now be stopped if not enough files are specified, if required input parameters are not set, if the input raster files have mismatching size / extent / band count / non-meter CRS, and others.
- Ellipse adjustment parameter: A new tunable parameter, HRR_ELLIPSE_ADJ, has been added to the &SIMULATOR namespace to scale the effective semi-major and semi-minor axes of the fire ellipse used in building fire spread calculations. This parameter multiplies the actual ellipse dimensions to account for non-uniform heat distribution within the ellipse. The default value is 0.5, but it can be calibrated for improved accuracy.

### Changed

- The smoke submodel of ELMFIRE has been reworked to make integration with HYSPLIT smoother. A simplified emission calculation model is used based on bulk flaming and smoldering wood emission values. The new parameters to be set are PM_EMISSION_FACTOR_FLAMING, PM_EMISSION_FACTOR_SMOLDERING, DRY_WOOD_CALORIFIC_VALUE, FLAMING_TIME, SMOLDERING_TIME. All have standard values that can be modified if needed. DUMP_EMITIMES has been added in the &OUTPUTS namelist to also create the emitimes.txt file required by HYSPLIT.
- Changed the default values of CURRENT_YEAR, HOUR_OF_YEAR, SUNRISE_HOUR, SUNSET_HOUR to -1, so an error occurs if these variables are needed for diurnal adjustment but are not set. SUNRISE_HOUR and SUNSET_HOUR are calculated automatically anyway but changing their values in the input files overwrites the automatic calculation.

## ELMFIRE-memopt.0427

### Bug Fixes

- Fixed WU-E accumulated heat-flux outputs so `DUMP_TOTAL_DFC_RECEIVED` and `DUMP_TOTAL_RAD_RECEIVED` write gridded totals for all receiving cells, rather than only values stored on tagged nodes.
- Fixed WU-E heat-flux time integration to use the active level-set timestep `DT`, instead of the user-specified `SIMULATION_DT`, so accumulated heat remains consistent with CFL-limited timesteps.
- Added safer WU-E building fuel model assignment for newly appended linked-list nodes, including fallback handling for missing building fuel model values.
- Fixed WU-E ignition initialization for random and scheduled ignitions so `PHIP`, `TIME_OF_ARRIVAL`, and `SURFACE_FIRE` are set before WUI tagging and structure HRR evaluation.
- Fixed WU-E transient HRR diagnostic flag usage to use `DUMP_HRR_TRANSIENT`.
- Guarded threshold interface-model state reads so `INTERFACE_MODEL_TYPE = 2` does not access WU-E arrays unless building spread model type 2 is active.
- Prevented WU-E structure nodes from being removed by the stale-node timeout before their design HRR curve has ended.
- Added single-node handling to linked-list deletion to keep list head/tail pointers valid after deleting the final node.

### Added

- Added WU-E transient output switches: `DUMP_TRANSIENT_DFC` for direct flame contact heat flux, `DUMP_TRANSIENT_RAD` for radiant heat flux, `DUMP_HRR_TRANSIENT` for transient HRRPUA, and `DUMP_FUEL_CONSUMPTION` for remaining fuel-load diagnostics.
- Added `BANDTHICKNESS_WUI` to the `&WUI` namelist to limit the neighborhood used for WU-E tagging and heat-flux updates. The default is 5 cells. This should be revised in the future with resolution-independent criterions (e.g. 100 meter).
- Added `CRITICL_HF_WUI` to the `&WUI` namelist as a critical heat-flux threshold for sustaining structure burning. The default is 0.
- Added WU-E state arrays for gridded transient and accumulated fields: `HRR_TRANSIENT_MAP`, `TOTAL_DFC_WUI`, `TOTAL_RADIATION_WUI`, `TRANSIENT_DFC_WUI`, `TRANSIENT_RADIATION_WUI`, `FUEL_LOAD_REMAIN`, and `ELLIPSE_PROPERTY_MAP`.
- Added gridded interface-model state arrays, `TEST_INTERFACE_WUI` and `WTU_SPREAD_WUI`, to carry threshold interface effects from WU-E heat-flux calculation back into linked-list node updates.
- Added `LIST_WUI_BURNING` plus `TAG_WUI` and `UNTAG_CELLS_WUI` to track WUI cells that need transient HRR and heat-flux updates.

### Changed

- Refactored the WU-E model so transient heat-flux calculation is handled by `CALC_WUI_HEATFLUX`, separate from `UMD_UCB_BLDG_SPREAD`.
- Moved WU-E heat-flux field updates out of `UX_AND_UY_ELLIPTICAL`; ROS calculation now reads the gridded transient heat-flux fields.
- Restored the threshold interface-model pathway for `INTERFACE_MODEL_TYPE = 2` by mapping vegetation-to-interface effects onto gridded WUI flags before node velocity and RK2 updates.
- Updated WU-E heat-flux calculation to use receiving-cell building fuel model properties for direct flame contact and radiation coefficients.
- Updated WU-E heat-flux normalization to use `HRR_ELLIPSE_ADJ` (in the &WUI namelist, previously in the &SIMULATOR namelist) when scaling the effective flame ellipse area.
- Simplified vegetative WU-E source preparation by using head-fire `FLIN_DMS_SURFACE` to estimate source HRRPUA, avoiding unnecessary normal-vector and ellipse-velocity updates for diagnostic/source-only vegetative cells. This is to substitute the constant 250 kW/m2 for vegetative cells in the previous versions. (Another option would be introducing a new wui vegetation class in the building fuel models).
- Replaced the previous WU-E dynamic node-pointer array with linked-list based WUI tracking and gridded state arrays.
- Updated WU-E transient HRR handling so structure and nearby vegetative WUI cells can update `HRR_TRANSIENT_MAP`.
- Limited WU-E transient and accumulated output rasters to building spread model type 2.
- Preserved isolated tagged pixels when the building spread model is enabled so isolated burning structures are not prematurely removed from the tagged list.
