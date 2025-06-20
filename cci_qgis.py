"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""
import os, sys
import numpy as np
from osgeo import gdal
from typing import Optional, Dict, Any
import json
import math
from qgis.core import (
    QgsFeatureSink,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingFeedback,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterDefinition
    
)
from qgis import processing

from PyQt5.QtCore import QCoreApplication


class CCI(QgsProcessingAlgorithm):

    INPUT_FOLDER_L1 = "INPUT_L1"
    INPUT_FOLDER_L2 = "INPUT_L2"
    OUTPUT_FOLDER = "OUTPUT"
    DTM = "DTM"
    MTL_FILE = "MTL_FILE"
    
    VEG_COEFF_B10 = "VEG_COEFF_B10"
    SOIL_COEFF_B10 = "SOIL_COEFF_B10"
    VEG_COEFF_B11 = "VEG_COEFF_B11"
    SOIL_COEFF_B11 = "SOIL_COEFF_B11´"
    
    AIR_TEMP_C = "AIR_TEMP_C"
    REL_HUM = "REL_HUM"
    WIND_SP = "WIND_SP"
    ATM_PRESS = "ATM_PRESS"
    
    C0 = "C0"
    C1 = "C1"
    C2 = "C2"
    C3 = "C3"
    C4 = "C4"
    C5 = "C5"
    C6 = "C6"
    
    SOLAR_CONSTANT = "SOLAR_CONSTANT"
    ATM_TRANS = "ATM_TRANS"
    
    

    def name(self) -> str:
        return "cciqgis"

    def displayName(self) -> str:
        return "Cooling Capacity Index"

    def group(self) -> str:
        return "Cooling Capacity Index"

    def groupId(self) -> str:
        return "cci"

    def shortHelpString(self) -> str:
        return "Example algorithm short description"
    
    def tr(self, message):
        return QCoreApplication.translate("CCI", message)
        
        
    # === Utility Methods ===

    def read_and_scale_band_l2(self, path):
        """
        Open Level-2 band and convert DN to surface reflectance.

        Returns:
            Scaled reflectance array, and GDAL dataset object
        """
        dataset = gdal.Open(path)
        if dataset is None:
            raise FileNotFoundError(f"Could not open {path}")
        band_array = dataset.ReadAsArray().astype(np.float32)
        return (band_array * 0.0000275) - 0.2, dataset
    
    def read_l1(self, path):
        """
        Open Level-1 band and return raw digital numbers (DN).
        """
        dataset = gdal.Open(path)
        if dataset is None:
            raise FileNotFoundError(f"Could not open {path}")
        band_array = dataset.ReadAsArray().astype(np.float32)
        return band_array, dataset

    def create_output_path(self, parameters, context, filename):
        """
        Construct full output path for saving a raster.
        """
        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        return os.path.join(output_folder, f"{filename}.tif")

    def save_raster(self, output_path, array, reference_path, feedback):
        """
        Save a NumPy array as a GeoTIFF, using metadata from a reference raster.
        """
        ref_ds = gdal.Open(reference_path)
        if ref_ds is None:
            feedback.reportError(f"Cannot open reference file {reference_path} for saving")
            raise IOError(f"Cannot open reference {reference_path}")

        driver = ref_ds.GetDriver()
        out_ds = driver.Create(
            output_path,
            ref_ds.RasterXSize,
            ref_ds.RasterYSize,
            1,
            gdal.GDT_Float32
        )
        out_ds.SetGeoTransform(ref_ds.GetGeoTransform())
        out_ds.SetProjection(ref_ds.GetProjection())

        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(array)
        out_band.SetNoDataValue(-9999)
        out_band.FlushCache()

        feedback.pushInfo(f"Saved raster to {output_path}")
        out_ds = None  # Close file
        ref_ds = None

    def W_toMJday(self, var):
        """
        Convert W/m² to MJ/m²/day.
        """
        return var * 0.0864
    
    # === Remote Sensing Calculations ===

    def ndvi(self, b4, b5):
        """
        Compute NDVI from Red (B4) and NIR (B5) bands.
        """
        ndvi = (b5 - b4) / (b4 + b5)
        return np.clip(ndvi, -1, 1)

    def pv(self, ndvi):
        """
        Compute fractional vegetation cover (Pv) from NDVI.
        """
        return ((ndvi - np.nanmin(ndvi)) / (np.nanmax(ndvi) - np.nanmin(ndvi))) ** 2

    def lse(self, pv, emis_vege, emis_soil):
        """
        Compute Land Surface Emissivity from Pv and emissivity coefficients.
        """
        return emis_vege * pv + emis_soil * (1 - pv)

    def sensor_radiance(self, thermal_band, gain, offset):
        """
        Convert DN to radiance using radiometric rescaling factors.
        """
        return gain * thermal_band + offset

    def tbright(self, k1, k2, l_sens):
        """
        Convert radiance to brightness temperature using Planck's equation.
        """
        return (k2 / np.log(k1 / l_sens + 1))

    def lst_sw(self, tbright_b10, tbright_diff, lse_avg, lse_diff, water_vap, c0, c1, c2, c3, c4, c5, c6):
        """
        Compute Land Surface Temperature using the Split-Window algorithm.
        """
        return (tbright_b10 + c1 * tbright_diff + c2 * tbright_diff**2 +
                c0 + (c3 + c4 * water_vap) * (1 - lse_avg) +
                (c5 + c6 * water_vap) * lse_diff) - 273.15

    # === Atmospheric and Radiation Models ===

    def water_wap(self, e0, rh):
        """
        Estimate water vapor content.
        """
        return (0.0981 * (10 * e0 * rh) + 0.1697) / 10

    def e0(self, rh, es):
        """
        Compute actual vapor pressure from RH and saturation vapor pressure.
        """
        return (rh / 100) * es

    def es(self, ta):
        """
        Compute saturation vapor pressure (kPa) from air temperature.
        """
        return (610.78 * math.exp((17.27 * ta) / (ta + 237.3))) / 1000

    def albedo(self, b2, b4, b5, b6, b7): 
        """
        Calculate broadband surface albedo from multispectral bands.
        """
        return 0.356 * b2 + 0.130 * b4 + 0.373 * b5 + 0.085 * b6 + 0.072 * b7 - 0.0018

    def kc(self, b4, b5, ndvi):
        """
        Estimate crop coefficient from SAVI and LAI.
        """
        savi = ((1 + 0.5)*(b5 - b4)) / (b5 + b4 + 0.5)
        lai = (ndvi - np.nanmin(ndvi)) / 0.6
        return 1.1 * (1 - np.exp(-1.5 * lai))

    def emisatm(self, ta_k):
        """
        Estimate atmospheric emissivity from temperature (Kelvin).
        """
        return 0.0000092 * (ta_k ** 2)

    def zenithangle(self, sun_elevation):
        """
        Convert sun elevation to solar zenith angle.
        """
        return 90 - sun_elevation

    def invertSE(self, dES):
        """
        Compute inverse square of Earth-Sun distance (for irradiance correction).
        """
        return 1 / (dES ** 2)

    def long_rad(self, lst, emis):
        """
        Compute outgoing or incoming longwave radiation.
        """
        return emis * 5.6703e-8 * lst ** 4

    def shortin_rad(self, solar_cons, zenith_angle, invertSE, atm_trans):
        """
        Compute incoming shortwave radiation.
        """
        return solar_cons * math.cos(math.radians(zenith_angle)) * invertSE * atm_trans

    def shortout_rad(self, albedo, shortin):
        """
        Compute reflected shortwave radiation.
        """
        return albedo * shortin

    
    def net_radiation(self, albedo, lse, shortin, longin, longout, shortout):
        """
        Calculate net radiation at the land surface.
        """
        return (1 - albedo) * shortin + longin - longout - (1 - lse) * longin
        
    def soil_flux(self, lst, albedo, ndvi, rn):
        """
        Estimate soil heat flux.
        """
        g = lst / albedo * (0.0038 * albedo + 0.0074 * albedo ** 2) * (1 - 0.98 * ndvi ** 4) * rn
        g = np.where(ndvi < 0, rn * 0.5, g)  #assume water
        g = np.where((lst < 4) & (albedo > 0.45), rn * 0.5, g) #assume snow
        return g
        
    def psychrometric_cons(self, p):
        """
        Calculate psychrometric constant.
        """
        
        return 0.000665 * p
        
    def vapour_press_curv(self, ta):
        """
        Calculate slope of saturation vapor pressure curve.
        """
        return (4098 * 0.6108 * np.exp ((17.27 * ta)/(ta+237.3)) / (ta + 237.3) ** 2)

    
    def et0(self, e0, slopevappress, windsp, es, g_mj, psychrom, rnet_mj, ta):
        """
        Calculate reference evapotranspiration (ET0).
        """
        return 0.408 * slopevappress * (rnet_mj - g_mj) + (900 * psychrom * windsp * (es - e0)) /  (ta+273.15) / (slopevappress + psychrom * (1 + 0.34 * windsp))

    def eti(self, et0, kc):
        """
        Calculate evapotranspiration index.
        """
        eti = (kc * et0) /  np.nanmax(et0)
        eti = np.where(eti > 1, 1, eti)
        return eti
    
    def hillshade(self, dmr, output):
        
        """
        Generate hillshade raster.
        """
        # Run the hillshade algorithm
        result = processing.run("qgis:hillshade", {
            'INPUT': dmr,
            'Z_FACTOR': 1.0,
            'AZIMUTH': 315.0,
            'ALTITUDE': 45.0,
            'COMPUTE_EDGES': False,
            'ZEVENBERGEN': False,
            'MULTIDIRECTIONAL': False,
            'OUTPUT': output
        })
        return result['OUTPUT']
        
    def cci(self, albedo, eti, hillshade):
        """
        Compute Cooling Capacity Index (CCI).
        """
    
        dataset = gdal.Open(hillshade)
        if dataset is None:
            raise FileNotFoundError(f"Could not open {path}")
        band_array = dataset.ReadAsArray().astype(np.float32)
        normalization_hillshade = band_array / 255
        
        cci = (0.6 * normalization_hillshade) + (0.2 * albedo) + (0.2 * eti)
        return np.clip(cci,0,1)

            
    def initAlgorithm(self, config: Optional[Dict[str, Any]] = None):

        self.addParameter(
            QgsProcessingParameterFile(
            self.INPUT_FOLDER_L1,
            self.tr("Select Input Folder with Imagery -- Level 1 TP"),
            behavior=QgsProcessingParameterFile.Folder
            
            )
        )

        self.addParameter(
            QgsProcessingParameterFile(
            self.INPUT_FOLDER_L2,
            self.tr("Select Input Folder with Imagery -- Level 2 SP"),
            behavior=QgsProcessingParameterFile.Folder
            
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.MTL_FILE,
                self.tr("Level 2 Metadata File (MTL)")
                
                )
        )

        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.OUTPUT_FOLDER,
                self.tr("Select Output Folder")
                )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.DTM,
                self.tr("Digital Terrain Model [m]")
                )
        )
        
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.AIR_TEMP_C,
                description=self.tr("Air Temeprature [°C]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=22.3
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.REL_HUM,
                description=self.tr("Relative Air Humidity [%]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=69
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.WIND_SP,
                description=self.tr("Wind Speed [m/s]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=3.4
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.ATM_PRESS,
                description=self.tr("Atmospheric Pressure [kPa]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=101.5
            )
        )
                            
        param = QgsProcessingParameterNumber(
            self.VEG_COEFF_B10,
            description=self.tr("Vegetation Emissivity Coefficient for B10"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.987
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.SOIL_COEFF_B10,
            description=self.tr("Soil Emissivity Coefficient for B10"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.971
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.VEG_COEFF_B11,
            description=self.tr("Vegetation Emissivity Coefficient for B11"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.989
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.SOIL_COEFF_B11,
            description=self.tr("Soil Emissivity Coefficient for B11"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.977
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        
        
        param = QgsProcessingParameterNumber(
            self.C0,
            description=self.tr("Split Window Coefficient C0"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=-0.268
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.C1,
            description=self.tr("Split Window Coefficient C1"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=1.378
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.C2,
            description=self.tr("Split Window Coefficient C2"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.183
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.C3,
            description=self.tr("Split Window Coefficient C3"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=54.300
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.C4,
            description=self.tr("Split Window Coefficient C4"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=-2.238
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.C5,
            description=self.tr("Split Window Coefficient C5"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=-129.200
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)

        param = QgsProcessingParameterNumber(
            self.C6,
            description=self.tr("Split Window Coefficient C6"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=16.400
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        
        param = QgsProcessingParameterNumber(
            self.SOLAR_CONSTANT,
            description=self.tr("Solar Constant [W/m2]"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=1367
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)
        
        param = QgsProcessingParameterNumber(
            self.ATM_TRANS,
            description=self.tr("Transmissiity of the Atmosphere"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.98
        )
        param.setFlags(param.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(param)



        

    def processAlgorithm(self, parameters, context, feedback):
        # === Retrieve input folder and parameter values from user input ===
        input_folder_l1 = self.parameterAsString(parameters, self.INPUT_FOLDER_L1, context)
        input_folder_l2 = self.parameterAsString(parameters, self.INPUT_FOLDER_L2, context)
        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        dtm = self.parameterAsString(parameters, self.DTM, context)
        mtl = self.parameterAsString(parameters, self.MTL_FILE, context)

        # === Read user-defined emissivity coefficients for thermal bands B10 and B11 ===
        vege_coeff_b10 = float(self.parameterAsDouble(parameters, self.VEG_COEFF_B10, context))
        soil_coeff_b10 = float(self.parameterAsDouble(parameters, self.SOIL_COEFF_B10, context))
        vege_coeff_b11 = float(self.parameterAsDouble(parameters, self.VEG_COEFF_B11, context))
        soil_coeff_b11 = float(self.parameterAsDouble(parameters, self.SOIL_COEFF_B11, context))

        # === Read meteorological data parameters ===
        ta = self.parameterAsDouble(parameters, self.AIR_TEMP_C, context)       # Air temperature (°C)
        rh = self.parameterAsDouble(parameters, self.REL_HUM, context)          # Relative humidity (%)
        solar_cons = self.parameterAsDouble(parameters, self.SOLAR_CONSTANT, context)  # Solar constant (W/m²)
        atm_transmis = self.parameterAsDouble(parameters, self.ATM_TRANS, context)     # Atmospheric transmissivity
        wind_speed = self.parameterAsDouble(parameters, self.WIND_SP, context)         # Wind speed (m/s)
        atm_press = self.parameterAsDouble(parameters, self.ATM_PRESS, context)        # Atmospheric pressure (kPa)

        # === Read coefficients for the Split-Window LST algorithm ===
        c0 = self.parameterAsDouble(parameters, self.C0, context)
        c1 = self.parameterAsDouble(parameters, self.C1, context)
        c2 = self.parameterAsDouble(parameters, self.C2, context)
        c3 = self.parameterAsDouble(parameters, self.C3, context)
        c4 = self.parameterAsDouble(parameters, self.C4, context)
        c5 = self.parameterAsDouble(parameters, self.C5, context)
        c6 = self.parameterAsDouble(parameters, self.C6, context)

        # === Dictionaries for storing band file paths and metadata ===
        l2_dict = {}      # For Level 2 (reflectance) bands
        l1_dict = {}      # For Level 1 (thermal) bands
        l1_metadata = {}  # For metadata from the MTL JSON file

        # === Parse Level 2 files and store paths for required bands ===
        for filename in os.listdir(input_folder_l2):
            if not filename.lower().endswith('.tif'):
                continue
            lower_name = filename.lower()
            if 'b2' in lower_name and 'sr' in lower_name:
                l2_dict['b2'] = os.path.join(input_folder_l2, filename)
            elif 'b3' in lower_name and 'sr' in lower_name:
                l2_dict['b3'] = os.path.join(input_folder_l2, filename)
            elif 'b4' in lower_name and 'sr' in lower_name:
                l2_dict['b4'] = os.path.join(input_folder_l2, filename)
            elif 'b5' in lower_name and 'sr' in lower_name:
                l2_dict['b5'] = os.path.join(input_folder_l2, filename)
            elif 'b6' in lower_name and 'sr' in lower_name:
                l2_dict['b6'] = os.path.join(input_folder_l2, filename)
            elif 'b7' in lower_name and 'sr' in lower_name:
                l2_dict['b7'] = os.path.join(input_folder_l2, filename)

        # === Parse Level 1 files for thermal bands and metadata ===
        for filename in os.listdir(input_folder_l1):
            if not filename.lower().endswith('.tif'):
                continue
            lower_name = filename.lower()
            if 'b10' in lower_name:
                l1_dict['b10'] = os.path.join(input_folder_l1, filename)
            elif 'b11' in lower_name:
                l1_dict['b11'] = os.path.join(input_folder_l1, filename)

        # === Read thermal metadata from MTL JSON file ===
        with open(mtl, 'r') as f:
            mtl_data = json.load(f)

        # Radiometric and thermal constants from metadata
        l1_metadata['K1_B10'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_THERMAL_CONSTANTS']['K1_CONSTANT_BAND_10'])
        l1_metadata['K1_B11'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_THERMAL_CONSTANTS']['K1_CONSTANT_BAND_11'])
        l1_metadata['K2_B10'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_THERMAL_CONSTANTS']['K2_CONSTANT_BAND_10'])
        l1_metadata['K2_B11'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_THERMAL_CONSTANTS']['K2_CONSTANT_BAND_11'])

        l1_metadata['RAD_MULT_B10'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING']['RADIANCE_MULT_BAND_10'])
        l1_metadata['RAD_MULT_B11'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING']['RADIANCE_MULT_BAND_11'])
        l1_metadata['RAD_ADD_B10'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING']['RADIANCE_ADD_BAND_10'])
        l1_metadata['RAD_ADD_B11'] = float(mtl_data['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING']['RADIANCE_ADD_BAND_11'])

        l1_metadata['SUN_ELEV'] = float(mtl_data['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES']['SUN_ELEVATION'])
        l1_metadata['dES'] = float(mtl_data['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES']['EARTH_SUN_DISTANCE'])

        # === Log the paths for debugging ===
        for band, path in l2_dict.items():
            feedback.pushInfo(f"{band}: {path}")
        for band, path in l1_dict.items():
            feedback.pushInfo(f"{band}: {path}")
            
        feedback.pushInfo(f"{l1_dict}")

        # === Check that all required bands are present ===
        required_bands_l2 = ['b2','b3','b4','b5','b6','b7']
        for b in required_bands_l2:
            if b not in l2_dict:
                raise QgsProcessingException(f"Missing band {b} in Level 2 input folder.")
        required_bands_l1 = ['b10','b11']
        for b in required_bands_l1:
            if b not in l1_dict:
                raise QgsProcessingException(f"Missing band {b} in Level 1 input folder.")

        # === Read and scale all spectral bands ===
        b2, ds_ref = self.read_and_scale_band_l2(l2_dict['b2'])
        b3, _ = self.read_and_scale_band_l2(l2_dict['b3'])
        b4, _ = self.read_and_scale_band_l2(l2_dict['b4'])
        b5, _ = self.read_and_scale_band_l2(l2_dict['b5'])
        b6, _ = self.read_and_scale_band_l2(l2_dict['b6'])
        b7, _ = self.read_and_scale_band_l2(l2_dict['b7'])

        b10, _ = self.read_l1(l1_dict['b10'])
        b11, _ = self.read_l1(l1_dict['b11'])

        # === Compute NDVI, Pv, and emissivities ===
        ndvi = self.ndvi(b4, b5)
        self.save_raster(self.create_output_path(parameters, context, "ndvi"), ndvi, l2_dict['b2'], feedback)

        pv = self.pv(ndvi)
        self.save_raster(self.create_output_path(parameters, context, "pv"), pv, l2_dict['b2'], feedback)

        lse_b10 = self.lse(pv, vege_coeff_b10, soil_coeff_b10)
        self.save_raster(self.create_output_path(parameters, context, "lse_b10"), lse_b10, l2_dict['b2'], feedback)

        lse_b11 = self.lse(pv, vege_coeff_b11, soil_coeff_b11)
        self.save_raster(self.create_output_path(parameters, context, "lse_b11"), lse_b11, l2_dict['b2'], feedback)

        # === Convert thermal bands to radiance and brightness temperature ===
        l_sens_b10 = self.sensor_radiance(b10, l1_metadata['RAD_MULT_B10'], l1_metadata['RAD_ADD_B10'])
        self.save_raster(self.create_output_path(parameters, context, "lsens_b10"), l_sens_b10, l2_dict['b2'], feedback)

        l_sens_b11 = self.sensor_radiance(b11, l1_metadata['RAD_MULT_B11'], l1_metadata['RAD_ADD_B11'])
        self.save_raster(self.create_output_path(parameters, context, "lsens_b11"), l_sens_b11, l2_dict['b2'], feedback)

        tbright_b10 = self.tbright(l1_metadata['K1_B10'], l1_metadata['K2_B10'], l_sens_b10)
        self.save_raster(self.create_output_path(parameters, context, "tbright_b10"), tbright_b10, l2_dict['b2'], feedback)

        tbright_b11 = self.tbright(l1_metadata['K1_B11'], l1_metadata['K2_B11'], l_sens_b10)  # NOTE: Should use `l_sens_b11`?
        self.save_raster(self.create_output_path(parameters, context, "tbright_b11"), tbright_b11, l2_dict['b2'], feedback)

        # === Atmospheric corrections ===
        saturated = self.es(ta)
        actual = self.e0(rh, saturated)
        w = self.water_wap(actual, rh)
        atm_emis = self.emisatm(ta + 273.15)
        psychrometric_cons = self.psychrometric_cons(atm_press)
        vap_curve = self.vapour_press_curv(ta)

        # === Compute Land Surface Temperature using Split-Window algorithm ===
        tbright_diff = tbright_b10 - tbright_b11
        lse_avg = 0.5 * (lse_b10 + lse_b11)
        lse_diff = lse_b10 - lse_b11
        lst = self.lst_sw(tbright_b10, tbright_diff, lse_avg, lse_diff, w, c0, c1, c2, c3, c4, c5, c6)
        self.save_raster(self.create_output_path(parameters, context, "lst"), lst, l2_dict['b2'], feedback)

        # === Calculate surface albedo and crop coefficient ===
        alb = self.albedo(b2, b4, b5, b6, b7)
        self.save_raster(self.create_output_path(parameters, context, "albedo"), alb, l2_dict['b2'], feedback)

        crop_coeff = self.kc(b4, b5, ndvi)
        self.save_raster(self.create_output_path(parameters, context, "kc"), crop_coeff, l2_dict['b2'], feedback)

        # === Radiation budget components ===
        zenith_ang = self.zenithangle(l1_metadata['SUN_ELEV'])
        inver_des = self.invertSE(l1_metadata['dES'])
        longout = self.long_rad(lst+273.15, lse_b10)
        longin = self.long_rad(lst+273.15, atm_emis)
        shortin = self.shortin_rad(solar_cons, zenith_ang, inver_des, atm_transmis)
        shortout = self.shortout_rad(alb, shortin)
        rnet = self.net_radiation(alb, lse_b10, shortin, longin, longout, shortout)

        self.save_raster(self.create_output_path(parameters, context, "longout"), longout, l2_dict['b2'], feedback)
        self.save_raster(self.create_output_path(parameters, context, "longin"), longin, l2_dict['b2'], feedback)
        self.save_raster(self.create_output_path(parameters, context, "shortout"), shortout, l2_dict['b2'], feedback)
        self.save_raster(self.create_output_path(parameters, context, "rnet"), rnet, l2_dict['b2'], feedback)

        # === Soil heat flux ===
        g = self.soil_flux(lst, alb, ndvi, rnet)
        self.save_raster(self.create_output_path(parameters, context, "g"), g, l2_dict['b2'], feedback)

        # === Evapotranspiration estimation ===
        g_mj = self.W_toMJday(g)
        rnet_mj = self.W_toMJday(rnet)
        et0 = self.et0(actual, vap_curve, wind_speed, saturated, g_mj, psychrometric_cons, rnet_mj, ta)
        self.save_raster(self.create_output_path(parameters, context, "et0"), et0, l2_dict['b2'], feedback)

        evapotr_index = self.eti(et0, crop_coeff)
        self.save_raster(self.create_output_path(parameters, context, "eti"), evapotr_index, l2_dict['b2'], feedback)

        # === Terrain shading effect from DTM ===
        hs_terrain = self.hillshade(dtm, os.path.join(output_folder, 'hillshade.tif'))

        # === Cooling Capacity Index (CCI) ===
        cooling_capacity_index = self.cci(alb, evapotr_index, hs_terrain)
        self.save_raster(self.create_output_path(parameters, context, "cci"), cooling_capacity_index, l2_dict['b2'], feedback)

        
        return {}
        

    def createInstance(self):
        return self.__class__()
