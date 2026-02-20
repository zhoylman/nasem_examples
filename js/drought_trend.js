// =======================================================
// EXPORT: 30-yr Moving Trends (p05 DEF, p05 P, p95 PET, p95 VPD/TMMN/TMMX)
// TerraClimate 1958–2024
// =======================================================

var states = ee.FeatureCollection('TIGER/2018/States');
var lower48 = states.filter(ee.Filter.inList('STATEFP', [
  '01','04','05','06','08','09','10','11','12','13','16','17','18','19','20','21',
  '22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37',
  '38','39','40','41','42','44','45','46','47','48','49','50','51','53','54','55','56'
]));
var conusMask = ee.Image(0).byte().paint(lower48, 1).selfMask();
var conusBBox = ee.Geometry.Rectangle([-125, 24, -66.5, 50], null, false);

var ic = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE');
var startYear = 1958;
var endYear   = 2024;
var windowSize = 30;

// 1) Annual totals/means for P, PET, DEF, plus VPD/TMMN/TMMX
var annualStats = ee.ImageCollection(
  ee.List.sequence(startYear, endYear).map(function(y) {
    var year = ee.Number(y);
    var filtered = ic.filter(ee.Filter.calendarRange(year, year, 'year'));

    var p   = filtered.select('pr').sum();
    var pet = filtered.select('pet').sum().multiply(0.1);   // <-- KEEP PET scaling
    var def = p.subtract(pet);

    var vpd  = filtered.select('vpd').mean().multiply(0.01);
    var tmmn = filtered.select('tmmn').mean().multiply(0.1);
    var tmmx = filtered.select('tmmx').mean().multiply(0.1);

    return p.rename('p')
      .addBands(pet.rename('pet'))
      .addBands(def.rename('def'))
      .addBands(vpd.rename('vpd'))
      .addBands(tmmn.rename('tmmn'))
      .addBands(tmmx.rename('tmmx'))
      .set('year', year)
      .toFloat();
  })
);

// 2) 30-yr moving window percentiles (end-year indexed)
var movingPercentiles = ee.ImageCollection(
  ee.List.sequence(startYear + windowSize - 1, endYear).map(function(y) {
    var endW = ee.Number(y);
    var startW = endW.subtract(windowSize - 1);
    var windowSlice = annualStats.filter(ee.Filter.rangeContains('year', startW, endW));

    // Existing percentiles
    var p05_def = windowSlice.select('def').reduce(ee.Reducer.percentile([5])).rename('p05_def');
    var p05_p   = windowSlice.select('p')  .reduce(ee.Reducer.percentile([5])).rename('p05_p');
    var p95_pet = windowSlice.select('pet').reduce(ee.Reducer.percentile([95])).rename('p95_pet');

    // NEW: p95 for VPD/TMMN/TMMX
    var p95_vpd  = windowSlice.select('vpd') .reduce(ee.Reducer.percentile([95])).rename('p95_vpd');
    var p95_tmmn = windowSlice.select('tmmn').reduce(ee.Reducer.percentile([95])).rename('p95_tmmn');
    var p95_tmmx = windowSlice.select('tmmx').reduce(ee.Reducer.percentile([95])).rename('p95_tmmx');

    var yearImg = ee.Image.constant(endW).rename('year').toFloat();

    return ee.Image([p05_def, p05_p, p95_pet, p95_vpd, p95_tmmn, p95_tmmx, yearImg])
      .set('end_year', endW)
      .updateMask(conusMask);
  })
);

// 3) Compute linear trends (slopes) for each moving-percentile metric
var metrics = ['p05_def', 'p05_p', 'p95_pet', 'p95_vpd', 'p95_tmmn', 'p95_tmmx'];

metrics.forEach(function(metric) {
  var slope = movingPercentiles.select(['year', metric])
    .reduce(ee.Reducer.linearFit())
    .select('scale')
    .rename(metric + '_slope')
    .toFloat();

  // 4) Export to Assets
  Export.image.toAsset({
    image: slope,
    description: 'TC_trend_30yr_moving_' + metric,
    assetId: 'projects/ee-zhoylman/assets/TC_trend_30yr_moving_' + metric,
    region: conusBBox,
    scale: 4638,
    maxPixels: 1e13
  });
});