# HS_tool
Tool for plume quantification using data from high-resolution imaging spectrometers.

# HS_tool description

- It goes from L1 to L4.

- Considered missions: EnMAP, EMIT, Gaofen-5A, and PRISMA

- L1 readers for these missions are included.

- Radiance filtering is manually applied

- L2 concentration enhancement maps using the matched-filter method (Thompson et al., 2016) are obtained for ch4, co2, h2o, c2h4, and nh3 and are expressed in ppm·m units.

- The unit gas absorption spectrum required for the matched filter method is obtained from LookUp Tables, which are included in this repository. The LUTs were obtained by combining HITRAN2024 + ARTS + Libradtran under a default midlatituder-summer atmosphere from the RFM model (Dudhia et al., 2017). Gas enhancement within the first 500 m above surface are assumed with surface height at sea level and PWV = 3 cm.

- Plume quantification is carried out using the IME approach (Frankenberg et al., 2016) and wind speed values are extracted from GEOS-FP (Molod et al., 2012).

- Plume quantification is only allowed for ch4, nh3, and c2h4, although there might be potential extension to other gases in the future.

- L2 files are stored in ENVI format; L4-related values are stored into a csv file; GPC-based georreferenced retrievals and plume bands are saved as tif files.

- The tool involves user intervention.

# Required settings

You first need to install all the required dependencies (requirements.txt + gdal + georasters + spectral).

The tool needs user intervation through a GUI. If run directly through terminal with 'ipython', we sometimes found the following error: 'Segmentation default (core dumped)'. In order to avoid it, we used the Python IDE Spyder, where the 'Graphic backend' (tools/preference/ipython console/plotting) was set to automatic. Further commits might fix this issue in the future.

# Required files

· Included in the repository:
- Default atmosphere file: 'RFM_midlatitude_summer.dat'
- LUTs for ch4, co2, h2o, and nh3, and unit c2h4 absorption spectrum (.npy format). These are located in the 'lut_nc' folder.
Important: change the path of these files according to where you store them in 'variable_definitions.py'

· Not included in the repository:
- L1 files. You can get them from PRISMA (https://prisma.asi.it/missionselect/), EMIT (https://search.earthdata.nasa.gov/), and EnMAP (https://planning.enmap.org/) data portals. Regarding Gaofen-5A, data was obtained through personal communication with Chinese institutions.

# Tool steps (manual intervention):

1) If the L2 retrieval is already stored, answer 'There is already a retrieval from this scene. Do you want to reprocess it? y/n:'. If you want to reprocess, type 'y'. If not 'n', and then skip to the quantification step (Step 3)

2) Radiance filtering. You have the opportunity to set radiance thresholds (water, clouds, cloud shadows, etc.). First, a radiance band quicklook will be shown. Check whether there are bottom or upper radiance thresholds you want to set. Close the window. Then, 'Masking lower values? y/n' and 'Masking upper values? y/n'. If 'y', set the desired values. Once the masked is applied, a new figure will be shown for you to confirm/reject the masking: 'Do you accept the masking? y/n'. If 'y', L2 retrievals are run and stored and you will continue to the next step. If 'n', Step 2 will be repeated.

3) Plume delineation. A new figure is shown. It has 2 panels: the previous radiance band and the gas concentration enhancement retrieval. A new question arises in console: 'Check the map. Do you still want to proceed? y/n'. If no plume detection, type 'n' and the whole process ends; if there is plume detection, type 'y'. Now you can interact with the figure following the indications from the suptitle: a) zoom-in to the subset of the image that captures the plume - use the magnifying lens for that -, and once you are finished press 'Enter'; b) click in the pixel where you think the emmission source is - you only have one opportunity, do not miss it or you will have to re-run the whole process; c) delineate the plume by clicking with the left button of the mouse - if you want to delete points you can use the right button -, and once you are done just press 'Enter'. Now, the wind speed data is extracted and stored or uploaded if already stored. Finally, the quantification is made and the related values are stored in a csv file. GPC-based georreferenced retrievals and plume bands are also stored.
