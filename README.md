# COOLING CAPACITY INDEX FOR QGIS FOR LANDSAT IMAGERY

A QGIS tool based on the Stanford University [INVEST project](https://naturalcapitalproject.stanford.edu/software/invest) and their Urban Cooling Model.

### Where to Place the Plugin
Place the downloaded file into your QGIS installation into processing/scripts subfolders in specific QGIS profile.

## Inputs and Outputs

## Requird Minimal Amount of Input Data

- Landsat 8/9 Level 1 and Landsat 8/9 Level 2 Imagery with MTL file
- Digital Terrain Model in the same extent and spatial resolution as satellite imagery
- Weather Information (Air Temperature, Relative Humidity, Wind Speed, Atmospheric Pressure)

### Outputs

- CCI GeoTIFF file 
- Several GeoTIFF files as intermediate results (e.g. ET0, Albedo, Land Surface Temperature etc.)

## Usage

Access the tool from Processing Toolbar from Scripts Section. 
