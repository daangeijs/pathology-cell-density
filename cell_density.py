from wholeslidedata.annotation.wholeslideannotation import WholeSlideAnnotation
from wholeslidedata.annotation.callbacks import ScalingAnnotationCallback
from wholeslidedata import WholeSlideImage
import shapely
from rasterio import features
import rasterio
import geopandas as gpd
import numpy as np


def load_mask(path, spacing):
    mask_obj = WholeSlideImage(path, backend='asap')
    mask = mask_obj.get_slide(spacing)
    return mask.squeeze()


def mask_to_polygons_layer(mask):
    all_polygons = []
    for shape, value in features.shapes(mask.astype(np.int16), mask=(mask > 0),
                                        transform=rasterio.Affine(1.0, 0, 0, 0, 1.0, 0)):
        all_polygons.append(shapely.geometry.shape(shape).buffer(0))

    all_polygons = shapely.geometry.MultiPolygon(all_polygons)
    return all_polygons


def cell_density(tumor_poly, cells, spacing):
    cell_points = [cell._geometry for cell in cells]
    a = gpd.GeoDataFrame(cell_points).set_geometry(0)
    b = gpd.GeoDataFrame([tumor_poly]).set_geometry(0)
    check = gpd.tools.sjoin(a, b, how='left')
    counts = sum(check['index_right'].notna())
    return counts / (tumor_poly.area * spacing * 10 ** -6)


mask_path = 'files/tumor_map.tif'
cell_detections = 'files/cell_detections.xml'

# Load the tissue tumor mask. Its a pyramid tif and the lowest level contains pixel spacing 2.0. We will load this on spacing 2.0.
mask = load_mask(mask_path, 2.0)

# We add a callback to scale the coordinates of the annotations to the same level as the mask, so we go from 0.24 um/pixel to 2.0 um/pixel
scaler = ScalingAnnotationCallback(1 / 8)
cells = WholeSlideAnnotation(cell_detections, callbacks=(scaler,))
print('Cells loaded:', len(cells.annotations))

# Convert the tumor (pixel-value=2) mask to a polygon
tumor_poly = mask_to_polygons_layer(mask == 2)

# Now we can calculate the cell density
cell_dens = cell_density(tumor_poly, cells.annotations, 2.0)
print(f"Cells per mm2 tumor area: {int(cell_dens)}")
