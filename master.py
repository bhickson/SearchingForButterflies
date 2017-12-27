import urllib.request, os, zipfile, ogr, numpy, random, shutil, struct, sys, osr
from osgeo.gdalnumeric import *
from gdalconst import *
from random import randint

from datetime import datetime, timedelta

weeks = 6

date = (datetime.today() - timedelta(days=6))  # prism data less than 5 days old is considered very preliiminary

# DIRECTORIES
prismDir = r"../DataLayers/PRISM"
npnDir = r"../DataLayers/NPN_Yearly"
staticDir = r"../DataLayers/SA/StaticRasters"

# SHAPEFILES
fc_merge = r"../DataLayers/SA/iNat_GBIF_merge_rural.shp"
backgrPoints = r"../DataLayers/SA/backgroundPoints_rural.shp"
fishnet = r"../DataLayers/SA/Fishnet4k_Selection.shp"

# VARIABLES
#  average first leaf out and first bloom dates (National Phrenology Network data)
npn_variables = ["avleaf", "avbloom"]
# ppt, tmin, tmax, tmean, tdmean (mean dewpoint temp), vpdmin (minimum vapor pressure deficit), vpdmax
prism_variables = ["ppt", "tmin", "tmax", "tmean", "tdmean", "vpdmin", "vpdmax"]
# elevation, % tree cover, distance from watercourse (NHD derived), GAP (2011) Formation Class, 5m TPI, Slope (SRTM 30m)
static_variables = {"Elevation_SRTM1_SA.tif": "Int16",
                    "TreeCover_NLCD_SA.tif": "UInt16",
                    "EucDist_NHD_SA.tif": "Float",
                    "FRM_CLASS_GAPv221_SA.tif": "Byte",
                    "TPI_5m.tif": "Float",
                    "Slope_SRTM1_SA.tif": "Float"}

point_sr = 102003  # EPSG Code for Point Feature Class

def downloadData(fname, down_addr):
    #Download the file via wcs for the relevant data

    with urllib.request.urlopen(down_addr) as response, open(fname, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)

    #u = urllib.reurlopen(down_addr)
    """f = open(fname, 'wb')
    meta = u.info()

    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)

    f.close()"""

def unzipFiles(in_dir):
    extension = ".zip"
    count=0
    print(in_dir)
    #os.chdir(in_dir) # change directory from working dir to dir with files

    for item in os.listdir(in_dir): # loop through items in dir
        if item.endswith(extension): # check for ".zip" extension
            count+=1
            #print("Starting on file", count,":",item)
            file_name = os.path.join(in_dir, item) # get full path of files
            zip_ref = zipfile.ZipFile(file_name) # create zipfile object
            zip_ref.extractall(in_dir) # extract file to dir
            zip_ref.close() # close file
            #os.remove(file_name) # delete zipped file

    print("Finished Unzipping File in", in_dir, datetime.now())

def doyFormat(doy):
    if len(doy) == 1:
        doy ="00" + doy
    elif len(doy) == 2:
        doy = "0" + doy
    return doy

def moveRasters(directory):
    print("starting move")
    if not os.path.exists(directory):
        os.makedirs(directory)

    for root, dirs, files in os.walk(directory):
        for file in files:
            fpath = os.path.join(root, file)

            if file.endswith(".bil") or file.endswith(".hdr") or file.endswith("prj"):
                extension = file[-4:]

                # Get year and day of year from file name
                file_parse = file.split("_")
                try:
                    file_date = file_parse[4]
                except:
                    continue
                dateobj = datetime.strptime(file_date, '%Y%m%d')
                year = dateobj.year
                doy = str(dateobj.timetuple().tm_yday)
                # set doy format to consistent length, e.g. 001 instead of 1
                doy = doyFormat(doy)

                # Create output directory
                year_dir = directory + "/" + str(year)
                if not os.path.exists(year_dir):
                    os.makedirs(year_dir)

                new_name = doy + extension
                opath = year_dir + "/" + new_name

                shutil.move(fpath, opath)
                    # move the raster and all dependent files with same name
                #parentdir = os.path.abspath(fpath + "/../")

                #for r, ds, fs in os.walk(parentdir):
                #    for f in fs:
                #        if f.split(".")[0] == file.split(".")[0]:
                #            fp = os.path.join(r, f)

            else:
                if not file.endswith(".zip"):
                    os.remove(fpath)    # Get rid of xml, txt, csv, and stx files

    print("Finished Raster Organization -", datetime.now())

def getPRISMData_ForDay(element): # ppt, tmin, tmax, tmea, tdmean (mean dewpoint temp), vpdmin (minimum vapor pressure deficit), vpdmax
    global prismDir, weeks
    daysprior = weeks * 7

    date = datetime.today() - timedelta(days=5) # prism data less than 5 days old is considered very preliiminary

    prismaddr = r'http://services.nacse.org/prism/data/public/4km/'

    outdir = prismDir + "/" + element + "/" #""PRISM_6w_" + today.strftime('%Y%m%d')
    print(outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Created output directory - ", outdir)

    prior_count = 0
    while daysprior>0:
        prevday = date - timedelta(days=daysprior)
        prevday = prevday.strftime('%Y%m%d')

        newFileName = "PRISM_" + element + "_" + prevday + ".zip"

        downloadaddr = prismaddr + element + "/" + prevday
        print("\t", downloadaddr)

        file_name = outdir + "/" + newFileName

        downloadData(file_name, downloadaddr)
        #Download the file via wcs for the relevant data

        zip_ref = zipfile.ZipFile(file_name)  # create zipfile object
        zip_ref.extractall(outdir)  # extract file to dir
        zip_ref.close()  # close file
        os.remove(file_name) # delete zipped file

        daysprior-=1

    print("FINISHED DOWNLOADS")

# USE TO GET BULK DATA FOR YEARS
def getPRISMData_BULK(element,start,end):
    global prismDir

    start = datetime.strptime(start, '%Y%m%d')
    end = datetime.strptime(end, '%Y%m%d')

    prismaddr = r'http://services.nacse.org/prism/data/public/4km/'

    outdir = prismDir + "/" + element
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Created output directory - ", outdir)

    day = end
    while day >= start:
        daystring = day.strftime("%Y%m%d")

        newFileName = element + "_" + daystring + ".zip"

        downloadaddr = prismaddr + element + "/" + daystring
        year = day.strftime("%Y")

        filepath = outdir + "/" + newFileName
        if not os.path.exists(filepath):
            # Download the file via wcs for the relevant data
            downloadData(filepath, downloadaddr)
            print("Finished Download:", day)

        day -= timedelta(days=1)

    print("FINISHED ALL YEAR DOWNLOADS FOR", element)

    #if element != "tmin":
    unzipFiles(outdir)
    moveRasters(outdir)

def getNPNData(year):
    datatypes = {"bloom":"si-x:average_bloom_prism",
                 "leaf":"si-x:average_leaf_prism"}

    outdir = npnDir + "/" + str(year)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for k,v in datatypes.items():
        fname = outdir + "/av" + k +"_" + str(year) + ".tif"
        print(k,v)
        geoserver_base_addr = 'https://geoserver.usanpn.org/geoserver/wcs?service=WCS&version=2.0.1&request=GetCoverage&coverageId='
        timestring ='&SUBSET=time("' + str(year) + '-01-01T00:00:00.000Z")'
        formstring = '&format=geotiff'
        address = geoserver_base_addr + v + timestring + formstring
        print(address)
        downloadData(fname,address)

def setRandomDates(fc):
    feat_count = 0
    # open vector layer
    drv = ogr.GetDriverByName('ESRI Shapefile')  # assuming shapefile?
    ds = drv.Open(fc, 1)  # open for editing
    lyr = ds.GetLayer(0)
    date_field = ogr.FieldDefn("DATE", ogr.OFTDate)

    lyr.CreateField(date_field)

    format = "%Y/%m/%d"

    beg_date = datetime.strptime("2007/01/01", format).date()
    end_date = datetime.strptime("2016/09/30", format).date()
    elapse = (end_date - beg_date).days

    for feat in lyr:
        random_day = randint(0, elapse)

        date = beg_date + timedelta(days=random_day)
        date_str = datetime.strftime(date, "%Y/%m/%d")

        feat.SetField("DATE", date_str)

        lyr.SetFeature(feat)
        feat_count += 1
        if feat_count % 100 == 0:
            print(feat_count)

def getTransform(raster_spatialref):
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(point_sr)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(raster_spatialref)
    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    return coordTransform

def computeFieldValues_NPN(fc, raster_dir, element, df):
    raster_sr = 4269  # EPSG CODE FOR RASTER PROJECTION
    print("Starting Compute of", element)
    feat_count = 0

    data_type = "Int16"
    # Open Vector Layer
    drv=ogr.GetDriverByName('ESRI Shapefile')   # Assuming Shapefile?
    ds=drv.Open(fc,1)   # Open For Editing
    lyr=ds.GetLayer(0)
    layerDef = lyr.GetLayerDefn()

    field_name = "DA_" + element

    fields = []
    for i in range(layerDef.GetFieldCount()):
        fname = layerDef.GetFieldDefn(i).GetName()
        fields.append(fname)
    if field_name not in fields:
        new_field = ogr.FieldDefn(field_name, ogr.OFTInteger64)
        lyr.CreateField(new_field)

    transformation = getTransform(raster_sr)

    for feat in lyr:
        geom = feat.GetGeometryRef()
        mx = geom.Centroid().GetX()
        my = geom.Centroid().GetY()

        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(mx, my)

        point.Transform(transformation)
        # Get feature Geometry
        mx = point.GetX()
        my = point.GetY()

        try:
            try:
                obs_date = datetime.strptime(feat.GetField(df), '%Y/%m/%d')
            except:
                obs_date = datetime.strptime(feat.GetField(df), '%Y-%m-%d')
        except:
            print("unable to find date field value")

        year = str(obs_date.year)
        doy = obs_date.timetuple().tm_yday

        raster = raster_dir + "/" + element + "_" + year + ".tif"

        if not os.path.exists(raster):
            print("Unable to find raster for", year, "day", doy)
            exit()

        vap = getPointValue(mx, my, raster, data_type)  # value at point of raster
        if vap > 365 or vap < 0:
            print("bad value for", ", year:", year, ", doy:", doy, vap)
            val = -9999
        else:
            val = int(doy - vap)    # Difference between the day of observation and estimated first leaf or bloom day of that year

        feat.SetField(field_name, val)

        lyr.SetFeature(feat)
        feat_count += 1
        if feat_count%1000 == 0:
            print(feat_count)
    lyr = None
    ds = None
    drv = None

def computeFieldValues_PRISM(fc, raster_dir, element, df):
    print("Starting Compute of", element, "for", fc)
    raster_sr = 4269
    data_type = "Float"
    feat_count = 0
    daysprior = weeks * 7
    # Open Vector Layer
    drv = ogr.GetDriverByName('ESRI Shapefile') #assuming shapefile?
    ds = drv.Open(fc,1) # open for editing
    lyr = ds.GetLayer()
    layerDef = lyr.GetLayerDefn()

    field_name = element.upper() + str(weeks) + 'W'

    if element != "ppt":
        field_name += "AV"

    fields = []
    for i in range(layerDef.GetFieldCount()):
        fname = layerDef.GetFieldDefn(i).GetName()
        fields.append(fname)
    if field_name not in fields:
        new_field = ogr.FieldDefn(field_name, ogr.OFTReal)
        lyr.CreateField(new_field)

    transformation = getTransform(raster_sr)

    for feat in lyr:
        #obs_date_1 = feat.GetField("obs_date_1")

        #if obs_date_1 == "0":
            geom = feat.GetGeometryRef()
            mx = geom.Centroid().GetX()
            my = geom.Centroid().GetY()

            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(mx, my)

            point.Transform(transformation)
            # Get feature Geometry
            mx = point.GetX()
            my = point.GetY()

            try:
                try:
                    obs_date = datetime.strptime(feat.GetField(df), '%Y/%m/%d')
                except:
                    obs_date = datetime.strptime(feat.GetField(df), '%Y-%m-%d')
            except:
                print("unable to find date field value")

            sum_value = 0
            prior_count = 0

            while prior_count < daysprior:
                prior_date = obs_date - timedelta(days=prior_count)
                year = str(prior_date.year)
                doy = str(prior_date.timetuple().tm_yday)
                doy = doyFormat(doy)
                raster = raster_dir + "/" + element + "/" + year + "/" + doy + ".bil"
                if not os.path.exists(raster):
                    print("Unable to find raster for", year, "day", doy)
                    exit()
                vap = getPointValue(mx, my, raster, data_type) # value at point of raster'
                if vap != -9999 and vap != None:
                    #print(element,year,doy)
                    sum_value += vap
                else:
                    print("bad value for", feat.GetField("OBJECTID"), ", year:", year, ", doy:", doy)

                prior_count+=1
            if element != "ppt":
                avgval = sum_value / daysprior
                feat.SetField(field_name, avgval)
            else:
                feat.SetField(field_name, sum_value)

            lyr.SetFeature(feat)

            feat_count+=1

            if feat_count%100 == 0:
                print(element, feat_count)

    lyr = None
    ds = None
    drv = None

def computeFieldValues_Static(fc, raster_dir, raster, data_type):
    print("Starting Compute of", raster)
    raster_sr = 102003
    field_name = raster.split("_")[0]
    raster_path = raster_dir + "/" + raster
    feat_count = 0
    # Open Vector Layer
    drv = ogr.GetDriverByName('ESRI Shapefile') #assuming shapefile?
    print(fc)
    ds = drv.Open(fc,1) #open for editing
    lyr = ds.GetLayer()
    layerDef = lyr.GetLayerDefn()

    fields = []
    for i in range(layerDef.GetFieldCount()):
        fname = layerDef.GetFieldDefn(i).GetName()
        fields.append(fname)
    if field_name not in fields:
        new_field = ogr.FieldDefn(field_name, ogr.OFTReal)
        lyr.CreateField(new_field)

    transformation = getTransform(raster_sr)

    for feat in lyr:
        geom = feat.GetGeometryRef()
        mx = geom.Centroid().GetX()
        my = geom.Centroid().GetY()

        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(mx, my)

        point.Transform(transformation)
        # Get feature Geometry
        mx = point.GetX()
        my = point.GetY()

        vap = getPointValue(mx, my, raster_path, data_type) # value at point of raster'

        feat.SetField(field_name, vap)
        lyr.SetFeature(feat)

        feat_count+=1
        if feat_count%100 == 0:
            print(raster, feat_count)
    lyr = None
    ds = None
    drv = None

def getPointValue(mx, my, rast, dt):
    system_byte_order = sys.byteorder
    if system_byte_order == "little":
        fmt = "<"
    elif system_byte_order == "big":
        fmt = ">"
    else:
        print("Unknown system byte order. Exiting")

    # open raster layer
    src_ds = gdal.Open(rast)
    gt = src_ds.GetGeoTransform()

    px = int((mx - gt[0]) / gt[1])  #x pixel (longitude - left coordinate) / raster resolution
    py = int((my - gt[3]) / gt[5])  #y pixel (latitdue - top coordinate) / raster resolution
    rb = src_ds.GetRasterBand(1)
    gdal.UseExceptions()    # so it doesn't print(to screen everytime point is outside grid

    if dt == "Int16":
        structval = rb.ReadRaster(px, py, 1, 1, buf_type=gdal.GDT_Int16)
        fmt += "h"
    elif dt == "Float":
        structval = rb.ReadRaster(px, py, 1, 1, buf_type=gdal.GDT_Float32)
        fmt += "f"
    elif dt == "UInt16":
        structval = rb.ReadRaster(px, py, 1, 1, buf_type=gdal.GDT_UInt16)
        fmt += "H"
    elif dt == "Byte":
        structval = rb.ReadRaster(px, py, 1, 1, buf_type=gdal.GDT_Byte)
        fmt += "b"
    #print(mx,my,structval)
    intval = struct.unpack(fmt, structval)
    #print(intval)
    val = intval[0] # stuct returns list of values even if only one
    src_ds = None

    return val

def computePredicionDayValues(date, prismVariables, npnVariables): # 20170615, prism_variables

    year = str(date.year)
    doy = date.timetuple().tm_yday

    prior_count = weeks * 7
    prior_date = date- timedelta(days=prior_count)

    startday = date.strftime("%Y%m%d")
    endday = prior_date.strftime("%Y%m%d")

    for variable in prismVariables:
        print("Starting", variable)
        getPRISMData_BULK(variable, endday, startday)
        prismVariablesDir = prismDir + "/" + variable
        computePRISMData(startday, variable, prismVariablesDir)
    for variable in npnVariables:
        getNPNData(year)
        computeNPNData(year, doy, variable, npnDir)

def computeNPNData(year, doy, variable, variablesDir):
    yearDir = variablesDir + "/" + year
    day = datetime(int(year), 1, 1) + timedelta(int(doy) -1)
    day = datetime.strftime(day, "%Y%m%d")

    for root,dirs,files in os.walk(yearDir):
        for file in files:
            fvariable = file.split("_")[0]  # <--- NPN Downloads follow avbloom_2017.tif format
            if variable == fvariable and file.endswith(".tif"):  # <--- default for PRISM DATA
                print(file)
                fpath = os.path.join(yearDir, file)
                ds = gdal.Open(fpath, GA_ReadOnly)
                band = ds.GetRasterBand(1)
                bandArray = BandReadAsArray(band)
                dataOut = doy - bandArray
                dataOut[dataOut >= 9999] = -9999

                newOutDir = variablesDir + "/" + day
                if not os.path.exists(newOutDir):
                    os.mkdir(newOutDir)

                outfile = newOutDir + "/" + "DA_" + variable.upper() + "_" + str(day) + ".tif"
                driver = gdal.GetDriverByName("GTiff")
                dsOut = driver.Create(outfile, ds.RasterXSize, ds.RasterYSize, 1, band.DataType)
                CopyDatasetInfo(ds, dsOut)
                bandOut = dsOut.GetRasterBand(1)
                BandWriteArray(bandOut, dataOut)

                # Close the variables/datasets
                ds = None
                band = None
                bandArray = None
                dataOut = None

def computePRISMData(day, variable, variablesDir):
    prior_count = weeks * 7
    for root, dirs, files in os.walk(variablesDir):
        for file in files:
            if file.endswith(".bil"):  # <--- default for PRISM DATA
                print(file)
                fpath = os.path.join(root, file)
                sample_file = fpath
                ds = gdal.Open(fpath, GA_ReadOnly)
                band = ds.GetRasterBand(1)
                bandArray = BandReadAsArray(band)
                if 'dataOut' not in locals():  # or dataOut == None:
                    dataOut = bandArray
                else:
                    dataOut += bandArray

                band = None
                ds = None
                bandArray = None

    # values without data (e.g. pacific ocean) are set to -9999 in each raster. Summming gives large negative number. Set all negatives.

    if variable != "ppt":
        dataOut = dataOut / (prior_count + 1)  # <--- Getting average for tmin/tmax/tmean/tdmean
    else:
        sumbadvalues = -9999 * ((prior_count) + 1)  # <-- 42 days plus start day
        dataOut[dataOut == sumbadvalues] = -9999

    # Get generic info from sample file
    sample_ds = gdal.Open(sample_file, GA_ReadOnly)
    sample_band = sample_ds.GetRasterBand(1)
    # Write the out file
    newOutDir = prismDir + "/" + day
    if not os.path.exists(newOutDir):
        os.mkdir(newOutDir)

    outfile = newOutDir + "/" + variable.upper() + "_6W_" + day + ".tif"
    driver = gdal.GetDriverByName("GTiff")
    dsOut = driver.Create(outfile, sample_ds.RasterXSize, sample_ds.RasterYSize, 1, sample_band.DataType)
    CopyDatasetInfo(sample_ds, dsOut)
    bandOut = dsOut.GetRasterBand(1)
    BandWriteArray(bandOut, dataOut)

    # Close the datasets
    bandOut = None
    dsOut = None
    sample_ds = None
    sample_band = None
    del dataOut


setRandomDates(backgrPoints)


startyear = 2007
while startyear < 2017:
    getNPNData(startyear)
    startyear += 1

for variable in npn_variables:
    print("Starting", variable)
    computeFieldValues_NPN(fc_merge, npnDir, variable, "obs_day")
    computeFieldValues_NPN(backgrPoints, npnDir, variable, "DATE")
for k, v in static_variables.items():
    print("Starting", k)
    computeFieldValues_Static(fc_merge, staticDir, k, v,)
    computeFieldValues_Static(backgrPoints, staticDir, k, v)
    computeFieldValues_Static(fishnet, staticDir, k,v)
for variable in prism_variables:
    print("Starting", variable)
    getPRISMData_BULK(variable, "20070101", "20161231")
    computeFieldValues_PRISM(fc_merge, prismDir, variable, "obs_day")
    computeFieldValues_PRISM(backgrPoints, prismDir, variable, "DATE")


# COMPUTE VALUES FOR A DAY OF YEAR
#  GET PRISM DATA FIRST
#  NPN DATA layers will be a leaf and bloom raster calculated as the doy (say june 15, 2017 with avleaf and avbloom for 2017 subtracted from it
#  Static variables are static

computePredicionDayValues(date, prism_variables, npn_variables)
