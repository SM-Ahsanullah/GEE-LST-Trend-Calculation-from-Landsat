// S M Ahsanullah
// Project LST for Cox's Bazar with Hossain Bin Idrish and Dr. Bhuiyan Alam

/****************  Cox’s Bazar ROI  ****************/
var coxsbazar = ee.FeatureCollection('projects/ee-smgge/assets/CoxsBazarStudyArea');

/****************  Years & sensor mapping  ***********/
var years = [1988, 1994, 1995, 2004, 2005, 2014, 2015, 2024, 2025];
var sensors = {
  1988: {id:'LANDSAT/LT05/C02/T1_L2', code:'LT05'},
  1994: {id:'LANDSAT/LT05/C02/T1_L2', code:'LT05'},
  1995: {id:'LANDSAT/LT05/C02/T1_L2', code:'LT05'},
  2004: {id:'LANDSAT/LT05/C02/T1_L2', code:'LT05'},
  2005: {id:'LANDSAT/LT05/C02/T1_L2', code:'LT05'},
  2014: {id:'LANDSAT/LC08/C02/T1_L2', code:'LC08'},
  2015: {id:'LANDSAT/LC08/C02/T1_L2', code:'LC08'},
  2024: {id:'LANDSAT/LC09/C02/T1_L2', code:'LC09'},
  2025: {id:'LANDSAT/LC09/C02/T1_L2', code:'LC09'}
};

/****************  Cloud & shadow mask  *******/
function maskL2(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)   // cloud
           .and(qa.bitwiseAnd(1 << 4).eq(0));  // shadow
  return image.updateMask(mask);
}

/****************  Metric & stats functions  *******/
function processImage(image) {
  var sensor = ee.String(image.get('SENSOR')),
      isLT05  = sensor.equals('LT05');

  // NDVI
  var nir = ee.Image(ee.Algorithms.If(isLT05, image.select('SR_B4'), image.select('SR_B5')))
              .multiply(0.0000275).add(-0.2),
      red = ee.Image(ee.Algorithms.If(isLT05, image.select('SR_B3'), image.select('SR_B4')))
              .multiply(0.0000275).add(-0.2),
      ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  // LST (°C)
  var thermalK = ee.Image(ee.Algorithms.If(isLT05, image.select('ST_B6'), image.select('ST_B10')))
                  .multiply(0.00341802).add(149.0),
      pv       = ndvi.subtract(0.2).divide(0.3).pow(2),
      emis     = pv.multiply(0.004).add(0.986),
      lst      = thermalK.expression(
                   '(T)/(1 + (0.00115 * (T/1.4388) * log(e))) - 273.15',
                   {T: thermalK, e: emis}
                 ).rename('LST');

  // NDBI
  var swir = ee.Image(ee.Algorithms.If(isLT05, image.select('SR_B5'), image.select('SR_B6')))
               .multiply(0.0000275).add(-0.2),
      ndbi = swir.subtract(red).divide(swir.add(red)).rename('NDBI');

  // UHI (null‐safe)
  var ruralMask    = ndbi.lt(0.1),
      ruralMeanRaw = lst.updateMask(ruralMask).reduceRegion({
                       reducer: ee.Reducer.mean(),
                       geometry: coxsbazar,
                       scale: 30,
                       maxPixels: 1e13
                     }).get('LST'),
      uhiImage    = ee.Algorithms.If(ruralMeanRaw, lst.subtract(ee.Number(ruralMeanRaw)), ee.Image(0)),
      uhi         = ee.Image(uhiImage).rename('UHI');

  // LULC
  var lulc = ee.Image(0)
             .where(ndvi.gt(0.3), 1)
             .where(ndbi.gt(0.2).and(ndvi.lt(0.2)), 2)
             .where(ndvi.lt(0), 3)
             .rename('LULC');

  return ee.Image.cat([ndvi, lst, ndbi, uhi, lulc])
           .set('system:time_start', image.get('system:time_start'))
           .set('SENSOR', sensor)
           .clip(coxsbazar);
}

function lstSceneStats(img) {
  var lstImg = processImage(img).select('LST'),
      stats  = lstImg.reduceRegion({
                 reducer: ee.Reducer.minMax(),
                 geometry: coxsbazar,
                 scale: 30,
                 maxPixels: 1e13
               });
  return ee.Feature(null, {
    date   : ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
    year   : ee.Date(img.get('system:time_start')).get('year'),
    id     : img.id(),
    minLST : stats.get('LST_min'),
    maxLST : stats.get('LST_max')
  });
}

/****************  Build stats & print image lists  *******/
var allStats = ee.FeatureCollection(
  years.map(function(year) {
    var cfg  = sensors[year],
        coll = ee.ImageCollection(cfg.id)
               .filterDate(year + '-01-01', year + '-12-31')
               .filter(ee.Filter.eq('WRS_PATH', 136))
               .filter(ee.Filter.eq('WRS_ROW', 45))
               .filterBounds(coxsbazar)
               .filterMetadata('CLOUD_COVER', 'less_than', 50)
               .map(maskL2)
               .map(function(img) { return img.set('SENSOR', cfg.code); });

    // ---- NEW: print out full scene names & IDs ----
    print('Year ' + year + ' — scene count:', coll.size());
    print('Year ' + year + ' — scene names:', coll.aggregate_array('system:index'));
    print('Year ' + year + ' — full asset IDs:', coll.aggregate_array('id'));

    return coll.map(lstSceneStats);
  })
).flatten();

print('Sample stats:', allStats.limit(10));

/****************  Export imagery  *******/
years.forEach(function(year) {
  var cfg = sensors[year];
  var coll = ee.ImageCollection(cfg.id)
    .filterDate(year + '-01-01', year + '-12-31')
    .filter(ee.Filter.eq('WRS_PATH', 136))
    .filter(ee.Filter.eq('WRS_ROW', 45))
    .filterBounds(coxsbazar)
    .filterMetadata('CLOUD_COVER', 'less_than', 50)
    .map(maskL2)
    .map(function(img) { return img.set('SENSOR', cfg.code); });

  var list = coll.toList(coll.size());

  var count = list.size().getInfo(); // get the number of images synchronously

  for (var i = 0; i < count; i++) {
  var img = ee.Image(list.get(i));
  var processed = processImage(img);
  var date = ee.Date(img.get('system:time_start')).format('YYYYMMdd');
  var baseName = year + '_' + date.getInfo();

  var bands = ['NDVI', 'LST', 'NDBI', 'UHI', 'LULC'];

  for (var j = 0; j < bands.length; j++) {
    var band = bands[j];
    var bandImage = processed.select(band);

    Export.image.toDrive({
      image: bandImage,
      description: band + '_' + baseName,
      fileNamePrefix: band + '_' + baseName,
      folder: 'GEE_Export_' + year,
      region: coxsbazar,
      scale: 30,
      maxPixels: 1e13
    });
  }
}
});


/****************  Annotated charts per year  *******/
/****************  Yearly average LST chart  ***************/
/****************  Yearly average LST chart  ***************/
var yearlyAvgStats = ee.FeatureCollection(
  years.map(function(year) {
    var stats = allStats.filter(ee.Filter.eq('year', year));
    
    var meanMin = ee.Number(stats.aggregate_stats('minLST').get('mean'));
    var meanMax = ee.Number(stats.aggregate_stats('maxLST').get('mean'));
    
    return ee.Feature(null, {
      year: year,
      meanMinLST: meanMin,
      meanMaxLST: meanMax
    });
  })
);

var avgChart = ui.Chart.feature.byFeature(yearlyAvgStats, 'year', ['meanMinLST', 'meanMaxLST'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Yearly Average LST (°C)',
    hAxis: {title: 'Year', format: '####'},
    vAxis: {title: 'LST (°C)'},
    pointSize: 6,
    lineWidth: 2,
    series: {
      0: {labelInLegend: 'Mean Min LST'},
      1: {labelInLegend: 'Mean Max LST'}
    },
    legend: {position: 'right'}
  });

print(avgChart);
