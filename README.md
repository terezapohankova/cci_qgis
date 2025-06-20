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
- List of output Data:
  
| Filename          | Fullname      | Units
| -------------     | ------------- |-------------
| albedo        | Albedo | unitless
| cci           | Cooling Capacity Index | unitless
| et0           | Reference Evapotranspiration | mm/day
| eti           | Evapotranspiration Index | unitless
| g             | Soil Heat Flux | W/m²
| hillshade     | Hillshaded Terrain | 0-255
| kc            | Crop Coefficient | unitless
| longin        | Incoming Longwave Solar Radiation | W/m²
| longout       | Outgoing Longwave Solar Radiation | W/m²
| lse_b10       | Land Surface Emissivity from Band 10 | unitless
| lse_b11       | Land Surface Emissivity from Band 11 | unitless
| lsens_b10     | Sensor Radiance for Band 10 | W/(m²·sr·µm)
| lsens_b11     | Sensor Radiance for Band 11 | W/(m²·sr·µm)
| lst           | Land Surface Temeprature | °C
| ndvi          | Normalized Differential Vegetation Index | unitless
| pv            | Proportion of vegetation | unitless
| rnet          | Solar Net Radiation | W/m²
| shortout      | Outgoing Shortwabe Solar Radiation | W/m²
| tbright_b10   | Brightness Temperature for B10| K
| tbright_b11   | Brightness Temperature for B11| K

## Usage

Access the tool from Processing Toolbar from Scripts Section. 
