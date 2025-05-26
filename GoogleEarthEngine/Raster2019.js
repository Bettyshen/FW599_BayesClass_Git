//===== The goal of the script is to export Landsat spectral band, temperature, and precipitation of 2019 as TIFF to GoogleDrive

// Load Enhanced vegetation index
var EVI = ee.Image("projects/shenf934044906/assets/EnvirPredict_2019/EVI_May_July_2019").select("EVI_median");
// Load PRISM data (precipitation & temperature)
var PRISM = ee.Image("projects/shenf934044906/assets/EnvirPredict_2019/PRISM_May_July_2019").select(['ppt_median','tmean_median']);



//======== Visualize on map ===========//

// PRISM
var precipitationVis = {
  min: 0.0,
  max: 300.0,
  palette: ['red', 'yellow', 'green', 'cyan', 'purple'],
};

Map.addLayer(PRISM.select('ppt_median'), precipitationVis, 'Precipitation_2019');
//====Import EVI and visualize it====//
  //A nice EVI palette
var palett = {'min': -0.1,
'max':1.0,
'palette':[
  'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718',
    '74A901', '66A000', '529400', '3E8601', '207401', '056201',
    '004C00', '023B01', '012E01', '011D01', '011301']};
    
Map.addLayer(EVI, palett, 'Enhanced Vegetation Index_MayJune_2019');

 
// Combine all bands into one image
//var combined = EVI.addBands(PRISM);
var combined = ee.Image.cat([EVI,PRISM]);

// Convert all bands to Float32 to ensure consistent data type
var combinedFloat32 = combined.toFloat();

// Export to drive
Export.image.toDrive({
  image: combinedFloat32,
  scale: 30,
  folder: 'GoogleEarthEngine',
  description: 'Composite_2019',
  crs: 'EPSG:4326',
  maxPixels: 1e12
});


