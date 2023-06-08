// correction function for radiometric slope correction on a
// Sentinel-1 image collection
var getS1LIA = function (collection,
                         options
                        ){
    // set defaults if undefined options
    options = options || {};
    var model = options.model || 'volume';
    var elevation = options.elevation || ee.Image('USGS/SRTMGL1_003');
    var buffer = options.buffer || 0;

    // we need a 90 degree in radians image for a couple of calculations
    var ninetyRad = ee.Image.constant(90).multiply(Math.PI/180);
    
    // buffer function (thanks Noel)
    function _erode(img, distance) {
      var d = (img.not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()));

      return img.updateMask(d.gt(distance));
    }

    // calculate masks
    function _masking(alpha_rRad, theta_iRad, proj, buffer){
        // layover, where slope > radar viewing angle
        var layover = alpha_rRad.lt(theta_iRad).rename('layover');

        // shadow
        var shadow = alpha_rRad.gt(ee.Image.constant(-1)
          .multiply(ninetyRad.subtract(theta_iRad))).rename('shadow');

        // combine layover and shadow
        var mask = layover.and(shadow);

        // add buffer to final mask
        if (buffer > 0)
            mask = _erode(mask, buffer);

        return mask.rename('no_data_mask');
   }

    function _getLIA(image){
        // get image geometry and projection
        var proj = image.select(1).projection();

        // get look direction angle
        var heading = (ee.Terrain.aspect(
            image.select('angle')).reduceRegion(ee.Reducer.mean()).get('aspect')
            );

        // Radar geometry
        var theta_iRad = image.select('angle').multiply(Math.PI/180);
        var phi_iRad = ee.Image.constant(heading).multiply(Math.PI/180);

        // Terrain geometry
        var alpha_sRad = ee.Terrain.slope(elevation).select('slope')
          .multiply(Math.PI/180).setDefaultProjection(proj);
        var phi_sRad = ee.Terrain.aspect(elevation).select('aspect')
          .multiply(Math.PI/180).setDefaultProjection(proj);

        //reduce to 3 angle
        var phi_rRad = phi_iRad.subtract(phi_sRad);

        // slope steepness in range
        var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();
        
        // Local incident angle, in degrees
        var LIA = (theta_iRad.subtract(alpha_rRad)).multiply(180/Math.PI).rename('LIA');
        var MSK = _masking(alpha_rRad, theta_iRad, proj, buffer);

        // return gamma_flat plus mask
        return image.addBands([LIA, MSK]).copyProperties(image);
    }
    // run correction function and return corrected collection
    return collection.map(_getLIA);
};

// export function
exports.getS1LIA = getS1LIA;