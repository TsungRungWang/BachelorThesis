/* 
################################################################################
AREA
################################################################################
*/ 
var Sub_AOI = subaoi.geometry();

// A backup for variable names, in case I have to assign boundaries manually
var AOI = Sub_AOI;
print('AOI', AOI);

var reduceRange = AOI;
var exportRange = AOI;


/* 
################################################################################
TIME
################################################################################
*/ 
var e = ee.Date('2018-06-28T00:00:00Z');

// Sentinel-1 time range
var preStart  = e.advance(-10000, 'day'),
    preEnd    = e.advance(0, 'day'),
    postStart = e.advance(10, 'day'), // July 8 is 10 days after June 28
    postEnd   = e.advance(10 + 4, 'day');

var PreEventPeriod = ee.Filter.date(preStart, preEnd),
    PostEventPeriod = ee.Filter.date(postStart, postEnd);

var exportFolder = 'S1Export';


/* 
################################################################################
NASADEM and masks
################################################################################
*/
var slope_threshold = 5; // Exclude areas with hillslope angle < slope_threshold, unit: degree
var curv_threshold = -0.005; // Exclude areas with hillslope curvature > curv_threshold, unit: m/m^2
// Negative value means it's a peak here
var valley_threshold = 0.003; // Re-include valley areas that had been excluded by flat slope

// extract elevation from NASA Digital Elevation Model (DEM)
var elevation = NASADEM_dataset.select('elevation');

var waterMask = NASADEM_dataset.select('swb').eq(0); // Create a binary water mask.
var slope = ee.Terrain.slope(elevation);
var slopeMask = slope.gte(slope_threshold); // slope mask with values 0 or 1, removes values less than or equal to threshold

// Calculate hillslope curvature
// Define a Gaussian kernel for smoothing. This step helps reduce noise in the curvature maps
var smooth_curv = ee.Kernel.gaussian({
  radius: 60,
  sigma: 30,
  units: 'meters',
  normalize: true,
});
// Smoothing the DEM with the gaussian kernel.
var elevation_smooth= elevation.convolve(smooth_curv).resample("bilinear");
var xyDemGrad = elevation_smooth.gradient().resample("bilinear");
var xGradient = xyDemGrad.select('x').gradient().resample("bilinear");
var yGradient = xyDemGrad.select('y').gradient().resample("bilinear");
var curvature = xGradient.select('x').add(yGradient.select('y'));

var curvatureMask = curvature.gte(curv_threshold);
var valleyMask = curvature.gte(valley_threshold);

var mask4in1 = waterMask.multiply(slopeMask.or(valleyMask)).multiply(curvatureMask);


/* 
################################################################################
Sentinel 1 radar image - Preprocessing
################################################################################
*/ 
var pmode = 'VH'; // Polarization Mode

var imgcol = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', pmode))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filterBounds(AOI)
  .filterDate(preStart,postEnd)
  .sort('system:time_start');
  
// Add Local Incidence Angle band to S1 images
var tclib = require('users/wcr/BS:getS1LIA');
imgcol = tclib.getS1LIA(imgcol);

imgcol = imgcol
  // Remove low edge values as suggested by GEE
  .map(function(image) {
    var edge = image.lt(-30.0);
    var maskedImage = image.mask().and(edge.not());
    return image.updateMask(maskedImage);
  })
  // Remove 4 in 1 masked areas
  .map(function(image) {
    return image.updateMask(image.mask().and(mask4in1));
  })
  // Removal no data areas due to terrain
  .map(function(image) {
    return image.updateMask(image.mask().and(image.select('no_data_mask'))
    )})
  /* Angle and no data mask bands no longer required, 
  remove unecessary bands to improve performance */
  .select([pmode]);

print(imgcol, 'All S1 images');
print(imgcol.filter(PreEventPeriod), 'All Pre S1 images');
print(imgcol.filter(PostEventPeriod), 'All Post S1 images');

// Just take a look at available paths. This list is no longer needed. 
var pathList = imgcol.aggregate_array('relativeOrbitNumber_start').distinct();
print(pathList, 'Path List');


/* 
################################################################################
Sentinel 1 radar image - Functions 
################################################################################
*/ 
// No median is taken in this step in the new method
var preimg = imgcol.filter(PreEventPeriod),
    postimg = imgcol.filter(PostEventPeriod);

function precolPostimg(precol, postimg) {
  var postimgPath = postimg.get('relativeOrbitNumber_start');
  var precolInPath = precol
    .filter(ee.Filter.eq('relativeOrbitNumber_start', postimgPath));
  var I_r = precolInPath.map(function(img) {
    return img.subtract(postimg);
  });
  var upperT = 90
  // Calculate pixels in an image where I ratio go beyond the upper threshold
  var Ibt = I_r.map(function(img) {
    var tr = img.reduceRegion({
      reducer: ee.Reducer.percentile([upperT, 10]),
      geometry: reduceRange,
      scale: 10, // Unit: meter
    });
    return img.gte(ee.Number(tr.get('VH_p' + upperT)));
  });
  var Ibt_m = Ibt.mean();
  return Ibt_m;
}

function precolPostcol(precol, postcol) {
  var SI = postcol.map(function(img) {
    return precolPostimg(precol, img);
  });
  return SI.mean();
}

var SI = precolPostcol(preimg, postimg);
print(SI, 'SI');
Map.addLayer(SI, {min:0, max:1}, 'SI');


/* 
################################################################################
Sentinel 1 radar image - Exporting results - SI
################################################################################
*/ 

Export.image.toDrive({
  image: SI,
  description: 'SI',
  folder: exportFolder,
  scale: 10, // Resolution in meters per pixel
  fileFormat: 'GeoTIFF',
  region: exportRange
});


/* 
################################################################################
Showing boundaries
################################################################################
*/ 

Map.centerObject(exportRange, 14);
Map.addLayer(AOI, [], 'AOI', false);
Map.addLayer(exportRange, [], 'Export Range', false);
Map.addLayer(GSI_AJG_Inventory, [], 'GSI_AJG_Inventory');