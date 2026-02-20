// =======================================================
// EXPORT: 30-yr Moving Trends (p05 DEF, p05 P, p95 PET, p95 VPD/TMMN/TMMX)
// TerraClimate 1958–2024
// Corrected: uses gte/lte for window selection and includes year band
// =======================================================

// --- CONUS mask + bbox (small payload)
var states = ee.FeatureCollection('TIGER/2018/States');
var lower48 = states.filter(ee.Filter.inList('STATEFP', [
  '01','04','05','06','08','09','10','11','12','13','16','17','18','19','20','21',
  '22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37',
  '38','39','40','41','42','44','45','46','47','48','49','50','51','53','54','55','56'
]));
var conusMask = ee.Image(0).byte().paint(lower48, 1).selfMask();
var conusBBox = ee.Geometry.Rectangle([-125, 24, -66.5, 50], null, false);

// --- TerraClimate & params
var ic = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE');
var startYear = 1958;
var endYear   = 2024;
var windowSize = 30;

// -----------------------------------------------------------------------------
// 1) Annual stats (full-year). Keep scalings consistent with TerraClimate table.
//    pr (p) has no scaling; pet/tmmn/tmmx/vpd scaled as requested.
// -----------------------------------------------------------------------------
var annualStats = ee.ImageCollection(
  ee.List.sequence(startYear, endYear).map(function(y) {
    var year = ee.Number(y);

    // full-year slice
    var filtered = ic.filter(ee.Filter.calendarRange(year, year, 'year'));

    // annual sums/means
    var p   = filtered.select('pr').sum().rename('p');                 // mm (no scaling)
    var pet = filtered.select('pet').sum().multiply(0.1).rename('pet'); // scaled *0.1 (mm)
    var def = p.subtract(pet).rename('def');                           // mm

    var vpd  = filtered.select('vpd').mean().multiply(0.01).rename('vpd');   // kPa
    var tmmn = filtered.select('tmmn').mean().multiply(0.1).rename('tmmn');  // °C
    var tmmx = filtered.select('tmmx').mean().multiply(0.1).rename('tmmx');  // °C

    // include year band (useful for debugging and ensures linearFit x-values are consistent)
    var yearBand = ee.Image.constant(year).rename('year').toFloat();

    return p.addBands([pet, def, vpd, tmmn, tmmx, yearBand])
      .set('year', year)
      .toFloat()
      .updateMask(conusMask); // mask early to reduce graph size
  })
);

// -----------------------------------------------------------------------------
// 2) 30-yr moving percentiles (end-year indexed).
//    p05 for def & p, p95 for pet, vpd, tmmn, tmmx
// -----------------------------------------------------------------------------
var movingPercentiles = ee.ImageCollection(
  ee.List.sequence(startYear + windowSize - 1, endYear).map(function(y) {
    var endW = ee.Number(y);
    var startW = endW.subtract(windowSize - 1);

    // CORRECT window selection using scalar comparisons
    var windowSlice = annualStats.filter(
      ee.Filter.and(
        ee.Filter.gte('year', startW),
        ee.Filter.lte('year', endW)
      )
    );

    // percentiles per band
    var p05_def = windowSlice.select('def').reduce(ee.Reducer.percentile([5])).rename('p05_def');
    var p05_p   = windowSlice.select('p')  .reduce(ee.Reducer.percentile([5])).rename('p05_p');
    var p95_pet = windowSlice.select('pet').reduce(ee.Reducer.percentile([95])).rename('p95_pet');

    var p95_vpd  = windowSlice.select('vpd') .reduce(ee.Reducer.percentile([95])).rename('p95_vpd');
    var p95_tmmn = windowSlice.select('tmmn').reduce(ee.Reducer.percentile([95])).rename('p95_tmmn');
    var p95_tmmx = windowSlice.select('tmmx').reduce(ee.Reducer.percentile([95])).rename('p95_tmmx');

    var yearImg = ee.Image.constant(endW).rename('year').toFloat();

    // return all percentile bands + year band (end-year)
    return ee.Image.cat([p05_def, p05_p, p95_pet, p95_vpd, p95_tmmn, p95_tmmx, yearImg])
      .set('end_year', endW)
      .updateMask(conusMask);
  })
);

// -----------------------------------------------------------------------------
// 3) Compute linear trends (slopes) for each moving-percentile metric
//    slope units will be: (units_of_metric) per year
// -----------------------------------------------------------------------------
var metrics = ['p05_def', 'p05_p', 'p95_pet', 'p95_vpd', 'p95_tmmn', 'p95_tmmx'];

metrics.forEach(function(metric) {
  var fit = movingPercentiles.select(['year', metric]).reduce(ee.Reducer.linearFit());

  // slope = 'scale' band from linearFit()
  var slope = fit.select('scale')
    .rename(metric + '_slope')
    .toFloat()
    .updateMask(conusMask);

  // export slope to an EE asset (one per metric)
  Export.image.toAsset({
    image: slope,
    description: 'TC_trend_30yr_moving_' + metric,
    assetId: 'projects/ee-zhoylman/assets/TC_trend_30yr_moving_' + metric,
    region: conusBBox,
    scale: 4638,
    maxPixels: 1e13
  });
});