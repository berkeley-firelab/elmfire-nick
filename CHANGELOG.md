# Changelog

We constantly update and upgrade ELMFIRE based on user feedback and new cases / requirements.

We want to expand the list of errors that can be caught early in the simulation.

Please create issues with your errors, requests or observations so they can be worked on for the next release.

## ELMFIRE 1.1

> Starting with this release, ELMFIRE uses semantic version numbers. Earlier builds were identified only by their compile date, which made it hard to reference or reproduce a specific version. This changelog covers all changes between ELMFIRE 1.1 and the original code this repository was forked from.

### Added
- **`DUMP_EVERY_STEP` (`&OUTPUTS`)**: writes transient outputs at every computation timestep.
- **Transient dump manifest**: dump-index filenames for transient rasters, plus a `dump_times_<case>.csv` manifest recording each dump index, its exact output time, and whether the dump is final.
- **`FEEDBACK_LEVEL` (`&SIMULATOR`)**: controls the verbosity of runtime feedback, with the options 0 to 3. This replaces the previous numeric `DEBUG_LEVEL` (see Removed).
- **Modular spotting model**: spotting is now broken into independent components through new switches — `GENERATION_MODEL` (`'RANDOM'`, `'PER-AREA'`, `'PER-MW'`), `SPOTTING_DISTANCE_MODEL` (`'UNIFORM'`, `'LOGNORMAL'`, `'EMPIRICAL'`), `ACCUMULATION_MODEL` (`'LAGRANGIAN'`, `'EULERIAN'`), and `IGNITION_MODEL` (`'DIRECT'`, `'SIMPLE'`, `'PHYSICAL'`). These replace the older boolean switches (see Removed and Superseded). The guide details which combinations are valid or invalid.
- **Line ignitions**: specify a line ignition through its start/end coordinates and a time using `X_LINE_IGN_START(:)`, `Y_LINE_IGN_START(:)`, `X_LINE_IGN_END(:)`, `Y_LINE_IGN_END(:)`, and `T_LINE_IGN(:)`. For example, `T_LINE_IGN(1)` applies to the line between (`X_LINE_IGN_START(1)`, `Y_LINE_IGN_START(1)`) and (`X_LINE_IGN_END(1)`, `Y_LINE_IGN_END(1)`). `NUM_IGNITIONS` does not need to be set when line ignitions are used.
- **WU-E transient output switches (`&OUTPUTS`)**: `DUMP_TRANSIENT_DFC` (direct flame contact heat flux), `DUMP_TRANSIENT_RAD` (radiant heat flux), and `DUMP_FUEL_CONSUMPTION` (remaining fuel-load diagnostics). (`DUMP_HRR_TRANSIENT` already existed and was corrected — see Fixed.)
- **`BANDTHICKNESS_WUI` (`&WUI`)**: limits the neighborhood used for WU-E tagging and heat-flux updates. Default is 5 cells. This should be revised in the future with a resolution-independent criterion (e.g. 100 m).
- **`CRITICL_HF_WUI` (`&WUI`)**: critical heat-flux threshold for sustaining structure burning. Default is 0.
- **`DUMP_EMBER_FLUX_TRANSIENT` (`&OUTPUTS`)**: transient ember-flux output in the Lagrangian accumulation path, consistent with the Eulerian path.
- **`HRR_ELLIPSE_ADJ` (`&WUI`)**: tunable parameter that scales the effective semi-major and semi-minor axes of the fire ellipse used in building fire spread, accounting for non-uniform heat distribution within the ellipse. Default is 0.5; it can be calibrated for improved accuracy.
- **More examples**: additional examples have been added to the new `examples/` folder.
- **Verification scripts**: in addition to the existing verification scripts, a new set of verification cases has been added (described in detail in the guide). Cases can be run automatically by copying the `TEMPLATE` folder, renaming it, and running each case within.
- **Automated validation**: a new pipeline, [WildfireAV](https://github.com/nick-cloudfire/WildfireAV), generates historic validation cases automatically and validates fire spread models at large scale. It is designed to work with any wildfire spread model. The Validation section of the guide outlines the results of the study; a forthcoming paper compares ELMFIRE and FARSITE validation statistics, which were shown to have near-matching overall similarity coefficients.
- **Rolling memory buffer**: large simulations (big domains, long run times) previously required large amounts of RAM. ELMFIRE now keeps only part of the data in memory at a time, reducing RAM usage by up to 80% with a runtime penalty of 5–20%. See `WX_BANDS_KEPT_IN_MEM` in the `&SIMULATOR` namelist.
- **Barrier breaching**: ELMFIRE can now model the effect of fire spread barriers such as fuel breaks, roads, and rivers. Set `USE_BARRIERS = .TRUE.` in `&INPUTS` and point `BARRIER_FILENAME` at a barrier GeoTIFF, where each barrier cell holds the barrier's width value. See the examples folder for details.
- **CFFDRS implementation**: to extend ELMFIRE into Canada, the CFFDRS (FBP) model has been added — head ROS for the FBP fuels, the build-up effect, and the CFFDRS crown involvement algorithm. Set `SURFACE_SPREAD_MODEL = "CFFDRS"`, `START_DC` (starting Drought Code), and `START_DMC` (starting Duff Moisture Code) in `&INPUTS`. Twenty verification cases from the 2009 FBP model updates are included, with automated scripts to run them and compare against the reference solutions. The Dogrib fire is included under Validation to compare the ELMFIRE CFFDRS implementation against the Prometheus/WISE reference case (the standard CFFDRS fire spread model).
- **Additional `&OUTPUTS` options**: `DUMP_CRITICAL_FLIN` (critical fireline intensity for canopy fire ignition), `DUMP_EMITIMES` (creates the `emitimes.txt` file required for a HYSPLIT smoke simulation), `DUMP_SPREAD_DIRECTION` (primary fire spread direction in degrees), `SPREAD_RATE_IN_M` (converts ROS output from feet per minute to meters per minute), and `DUMP_CFFDRS_DEBUG` (CFFDRS diagnostics).
- **Raster perturbation distributions (`&MONTE_CARLO`)**: each perturbed raster can now draw its Monte Carlo offset from a `PDF_TYPE` of `'UNIFORM'`, `'GAUSSIAN'`, or `'LOGNORMAL'`. `'UNIFORM'` uses `PDF_LOWER_LIMIT`/`PDF_UPPER_LIMIT` as before; `'GAUSSIAN'` and `'LOGNORMAL'` use the new `PDF_MEAN`/`PDF_SIGMA` parameters.
- **Input file checking**: input validation now runs at the start of a simulation. The run stops if too few files are specified, if required input parameters are unset, or if input rasters have mismatching size / extent / band count or a non-meter CRS, among other checks.

### Changed
- **Transient output scheduling**: fixed-`DTDUMP` output is no longer limited by the previous fixed-size dump-time array.
- **Final-output handling**: now uses an explicit final-dump flag, preserving final-only outputs (time of arrival, fireline intensity, spread rate, ember ignition) while allowing transient fields to dump repeatedly.
- **Smoke submodel**: reworked for smoother HYSPLIT integration. A simplified emission model based on bulk flaming and smoldering wood emission values is now used. The new parameters are `PM_EMISSION_FACTOR_FLAMING`, `PM_EMISSION_FACTOR_SMOLDERING`, `DRY_WOOD_CALORIFIC_VALUE`, `FLAMING_TIME`, and `SMOLDERING_TIME` (all have standard defaults that can be modified if needed). `DUMP_EMITIMES` (`&OUTPUTS`) creates the `emitimes.txt` file required by HYSPLIT.
- **`HRR_ELLIPSE_ADJ` moved** from `&SIMULATOR` to `&WUI`; WU-E heat-flux normalization now uses it when scaling the effective flame ellipse area.
- **Vegetative WU-E source preparation simplified**: head-fire `FLIN_DMS_SURFACE` is used to estimate source HRRPUA, avoiding unnecessary normal-vector and ellipse-velocity updates for diagnostic/source-only vegetative cells. This replaces the constant 250 kW/m² used for vegetative cells in earlier versions. (An alternative would be a dedicated WUI vegetation class in the building fuel models.)
- **WU-E heat-flux coefficients**: direct flame contact and radiation now use the receiving cell's building fuel model properties.
- **WU-E transient and accumulated output rasters** are now limited to building spread model type 2.
- **`NUM_METEOROLOGY_TIMES` default changed to -1**: leaving it unset makes it equal the number of bands in the wind speed file.
- **Diurnal-adjustment defaults changed to -1** for `CURRENT_YEAR`, `HOUR_OF_YEAR`, `SUNRISE_HOUR`, and `SUNSET_HOUR`, so an error is raised if these are needed for diurnal adjustment but are not set. `SUNRISE_HOUR` and `SUNSET_HOUR` are computed automatically; setting them in the input file overrides the automatic calculation.

### Removed
- Removed `DEBUG_LEVEL`, which has been replaced by `FEEDBACK_LEVEL` (see Added).
- Removed `NUM_NODES_OMP_THRESHOLD`, as it was never used in the code. 
- Removed various spotting input switches, which were replaced with new, clearer input settings (see Added and Superseded): `USE_EULERIAN_SPOTTING`, `USE_PHYSICAL_EMBER_NUMBER`, `USE_EMBER_IGNITION_MODEL`, `USE_SIMPLE_IGNITION_MODEL`, and the active use of `USE_UMD_SPOTTING_MODEL`.
- Removed the 7×7 burned-cell exclusion from the new Lagrangian direct-ignition spotting path. The legacy spotting path preserves its original neighborhood exclusion.
- Removed the `COMPUTATIONAL_DOMAIN` input namelist. The domain `XLLCORNER`, `YLLCORNER`, `CELLSIZE`, and `A_SRS` values are now set by ELMFIRE via GDAL commands.
- Removed the `LATITUDE` and `LONGITUDE` inputs required for diurnal fire adjustments. ELMFIRE now calculates them automatically from the domain center dimensions.

### Fixed
- Fixed `HRR_TRANSIENT` outputs so `DUMP_HRR_TRANSIENT = .TRUE.` writes HRRPUA for all actively burning fuels within residence time, not only WUI cells.
- Fixed simulation stop handling so runs are clipped to `SIMULATION_TSTOP` and exit after final cleanup instead of continuing to the end of the meteorology bands.
- Fixed WU-E accumulated heat-flux outputs so `DUMP_TOTAL_DFC_RECEIVED` and `DUMP_TOTAL_RAD_RECEIVED` write gridded totals for all receiving cells, rather than only values stored on tagged nodes.
- Fixed WU-E heat-flux time integration to use the active level-set timestep `DT` instead of the user-specified `SIMULATION_DT`, so accumulated heat remains consistent with CFL-limited timesteps.
- Added safer WU-E building fuel model assignment for newly appended linked-list nodes, including fallback handling for missing building fuel model values.
- Fixed WU-E ignition initialization for random and scheduled ignitions so `PHIP`, `TIME_OF_ARRIVAL`, and `SURFACE_FIRE` are set before WUI tagging and structure HRR evaluation.
- Fixed WU-E transient HRR diagnostic output to be controlled by `DUMP_HRR_TRANSIENT`.
- Guarded threshold interface-model state reads so `INTERFACE_MODEL_TYPE = 2` does not access WU-E arrays unless building spread model type 2 is active.
- Restored the threshold interface model (`INTERFACE_MODEL_TYPE = 2`) so vegetation-to-interface effects are applied before node velocity and RK2 updates.
- Preserved isolated tagged pixels when the building spread model is enabled, so isolated burning structures are not prematurely removed from the tagged list.
- Prevented WU-E structure nodes from being removed by the stale-node timeout before their design HRR curve has ended.
- Added single-node handling to linked-list deletion to keep list head/tail pointers valid after deleting the final node.
- Fixed a bug where a positive `NODATA_VALUE` in rasters was handled incorrectly.
- Reworked the BSQ metadata parser to find its target values more flexibly instead of relying on fixed character ranges (which caused issues depending on how the input TIFs were converted to BSQ).
- Fixed a bug where the final ROS value in Fire Potential mode (`MODE = 2`) was incorrectly calculated.
- Fixed a bug where user-specified firebrand inputs were overwritten back to their defaults.
- Fixed flame-ellipse building fire spread to fall back to default values when nonuniform building rasters contain `NODATA_VALUE` in the building-size or structure-separation-distance layers.
- Fixed handling of the structure-separation-distance input raster so values above 50 m fall back to the 50 m default (the coarsest input resolution is currently 30 m).

### Superseded
- Legacy spotting remains available with `USE_SUPERSEDED_SPOTTING = .TRUE.` and uses the legacy `SPOTTING_DISTRIBUTION_TYPE` setting.

### Developer / Internal
These changes are implementation details with no direct user-facing input or behavior change; they are listed for developers working on the source.

- Added a quick compile option (./make_gnu.sh -f) that only compiles the basic elmfire executable.
- Added WU-E gridded state arrays for transient and accumulated fields: `HRR_TRANSIENT_MAP`, `TOTAL_DFC_WUI`, `TOTAL_RADIATION_WUI`, `TRANSIENT_DFC_WUI`, `TRANSIENT_RADIATION_WUI`, `FUEL_LOAD_REMAIN`, and `ELLIPSE_PROPERTY_MAP`.
- Added gridded interface-model state arrays `TEST_INTERFACE_WUI` and `WTU_SPREAD_WUI` to carry threshold interface effects from the WU-E heat-flux calculation back into linked-list node updates.
- Added `LIST_WUI_BURNING`, `TAG_WUI`, and `UNTAG_CELLS_WUI` to track WUI cells that need transient HRR and heat-flux updates.
- Replaced the previous WU-E dynamic node-pointer array with linked-list based WUI tracking and gridded state arrays.
- Refactored the WU-E model so transient heat-flux calculation is handled by `CALC_WUI_HEATFLUX`, separate from `UMD_UCB_BLDG_SPREAD`.
- Moved WU-E heat-flux field updates out of `UX_AND_UY_ELLIPTICAL`; ROS calculation now reads the gridded transient heat-flux fields.
- Mapped vegetation-to-interface effects onto gridded WUI flags before node velocity and RK2 updates (the implementation behind the restored `INTERFACE_MODEL_TYPE = 2` pathway).
- Updated WU-E transient HRR handling so structure and nearby vegetative WUI cells can update `HRR_TRANSIENT_MAP`.
- Removed the active-model dependency on `USE_UMD_SPOTTING_MODEL` for allocating and updating `SPOTTING_STATS`.
- Renamed functions `SARDOY_PDF_PARAMETERS` → `EMPIRICAL_PDF_PARAMETERS`, `SARDOY_PDFINV` → `LOGNORM_CDF`, and `SARDOY_CDF` → `LOGNORM_CDF_DEFINITE` to clarify their roles.
- Renamed the variable `EMBER_TIGN` → `EMBER_TOA` (the time of arrival of the first ember).
