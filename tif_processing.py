from osgeo import gdal
from osgeo import osr


# From a tif file and a list of latitudes and longitudes, return the data corresponding to each point in a text file.
def compute_tif_data(tif_file, output_file, lats, lons):
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    tif = gdal.Open(tif_file)
    printInfo(tif)
    rows_tif = tif.RasterYSize
    cols_tif = tif.RasterXSize
    # coordinates to get pixel values for
    res = open(output_file, 'w')
    for lat in lats:
        for lon in lons:
            pixel = latLonToPixel(tif, [lat, lon], True)
            if 0 <= pixel[0] < cols_tif and 0 <= pixel[1] < rows_tif:
                val = value(tif, pixel)
            else:
                val = "NA"
            res.write(str(lat))
            res.write("\t")
            res.write(str(lon))
            res.write("\t")
            res.write(str(val))
            res.write("\n")

    res.close()


# The following method translates given latitude/longitude pairs into pixel locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#      latLonPairs - The decimal lat/lon pairings to be translated in the form [[lat1,lon1],[lat2,lon2]]
# OUTPUT: The pixel translation of the lat/lon pairings in the form [[x1,y1],[x2,y2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough
#      image resolution for pixel size to be insignificant
def latLonToPixel(ds, latLon, lcl):
    # Load the image dataset
    # ds = gdal.Open(geotifAddr)
    # Get a geo-transform of the dataset
    gt = ds.GetGeoTransform()
    # Create a spatial reference object for the dataset
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    # Set up the coordinate transformation object
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srsLatLong, srs)
    # Change the point locations into the GeoTransform space
    # It seems there is an error with that method in case the DS is already defined in degrees we do not need this !!
    if (lcl): (latLon[1], latLon[0], holder) = ct.TransformPoint(latLon[1], latLon[0])
    # Translate the x and y coordinates into pixel values
    x = (latLon[1] - gt[0]) / gt[1]
    y = (latLon[0] - gt[3]) / gt[5]
    return [int(x), int(y)]


# The following method translates given pixel locations into latitude/longitude locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#      pixelPairs - The pixel pairings to be translated in the form [[x1,y1],[x2,y2]]
# OUTPUT: The lat/lon translation of the pixel pairings in the form [[lat1,lon1],[lat2,lon2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough
#      image resolution for pixel size to be insignificant
def pixelToLatLon(ds, pixel):
    # Get a geo-transform of the dataset
    gt = ds.GetGeoTransform()
    # Create a spatial reference object for the dataset
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    # Set up the coordinate transformation object
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srs, srsLatLong)
    # Translate the pixel pairs into untranslated points
    ulon = pixel[0] * gt[1] + gt[0]
    ulat = pixel[1] * gt[5] + gt[3]
    # Transform the points to the space
    (lon, lat, holder) = ct.TransformPoint(ulon, ulat)
    return [lat, lon]


def value(ds, pixel):
    band = ds.GetRasterBand(1)  # 1-based index
    # read data and add the value to the string
    data = band.ReadAsArray(pixel[0], pixel[1], 1, 1)
    return data[0, 0]


def printInfo(ds):
    # get georeference info
    transform = ds.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    print
    'Transform =', transform
    print
    'Origin = (', transform[0], ',', transform[3], ')'
    print
    'Pixel Size = (', transform[1], ',', transform[5], ')'