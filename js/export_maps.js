// =======================================================
// Export 30-yr moving trend assets to Google Drive
// (p05 DEF, p05 P, p95 PET, p95 VPD, p95 TMMN, p95 TMMX)
// =======================================================

// ---- Assets ----
var def  = ee.Image('projects/ee-zhoylman/assets/TC_trend_30yr_moving_p05_def');
var p    = ee.Image('projects/ee-zhoylman/assets/TC_trend_30yr_moving_p05_p');
var pet  = ee.Image('projects/ee-zhoylman/assets/TC_trend_30yr_moving_p95_pet');

var vpd  = ee.Image('projects/ee-zhoylman/assets/TC_trend_30yr_moving_p95_vpd');
var tmmn = ee.Image('projects/ee-zhoylman/assets/TC_trend_30yr_moving_p95_tmmn');
var tmmx = ee.Image('projects/ee-zhoylman/assets/TC_trend_30yr_moving_p95_tmmx');

// ---- CONUS mask (lower 48) ----
var states = ee.FeatureCollection('TIGER/2018/States');

var lower48 = states.filter(ee.Filter.inList('STATEFP', [
  '01','04','05','06','08','09','10','11','12','13','16','17','18','19','20','21',
  '22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37',
  '38','39','40','41','42','44','45','46','47','48','49','50','51','53','54','55','56'
]));

var conusMask = ee.Image(0).byte().paint(lower48, 1).selfMask();

// ---- Simple bounding box (small request payload) ----
var conusBBox = ee.Geometry.Rectangle([-125, 24, -66.5, 50], null, false);

// ---- Apply mask ----
var defMasked  = def.updateMask(conusMask);
var pMasked    = p.updateMask(conusMask);
var petMasked  = pet.updateMask(conusMask);
var vpdMasked  = vpd.updateMask(conusMask);
var tmmnMasked = tmmn.updateMask(conusMask);
var tmmxMasked = tmmx.updateMask(conusMask);

// =======================================================
// EXPORTS
// =======================================================

Export.image.toDrive({
  image: defMasked,
  description: 'TC_p05_def_trend',
  folder: 'earthengine',
  fileNamePrefix: 'TC_p05_def_trend_30yr_CONUS',
  region: conusBBox,
  scale: 4638,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: pMasked,
  description: 'TC_p05_p_trend',
  folder: 'earthengine',
  fileNamePrefix: 'TC_p05_p_trend_30yr_CONUS',
  region: conusBBox,
  scale: 4638,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: petMasked,
  description: 'TC_p95_pet_trend',
  folder: 'earthengine',
  fileNamePrefix: 'TC_p95_pet_trend_30yr_CONUS',
  region: conusBBox,
  scale: 4638,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: vpdMasked,
  description: 'TC_p95_vpd_trend',
  folder: 'earthengine',
  fileNamePrefix: 'TC_p95_vpd_trend_30yr_CONUS',
  region: conusBBox,
  scale: 4638,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: tmmnMasked,
  description: 'TC_p95_tmmn_trend',
  folder: 'earthengine',
  fileNamePrefix: 'TC_p95_tmmn_trend_30yr_CONUS',
  region: conusBBox,
  scale: 4638,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: tmmxMasked,
  description: 'TC_p95_tmmx_trend',
  folder: 'earthengine',
  fileNamePrefix: 'TC_p95_tmmx_trend_30yr_CONUS',
  region: conusBBox,
  scale: 4638,
  maxPixels: 1e13
});