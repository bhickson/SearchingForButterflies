from geopy.distance import vincenty
from osgeo import ogr
import os

az = r"D:\eButterfly\DataLayers\SA\Extracts_Aug4th\tucson_extracts_FS.shp"
wa = r"C:\Users\hicksonb\Google Drive\eButterfly\DataLayers\SA\Extracts_July15\seattle_extracts.shp"
mt = r"C:\Users\hicksonb\Google Drive\eButterfly\DataLayers\SA\Extracts_July8\boseman_extracts.shp"
ca = r"C:\Users\hicksonb\Google Drive\eButterfly\DataLayers\SA\Extracts_July15\sacramento_extracts.shp"
ut = r"C:\Users\hicksonb\Google Drive\eButterfly\DataLayers\SA\Extracts_July15\slc_extracts.shp"

shapefiles = [az]#, nv, ut]

for shapefile in shapefiles:
    print("Starting",shapefile)
    drv = ogr.GetDriverByName('ESRI Shapefile')  # assuming shapefile?
    ds = drv.Open(shapefile, 1)  # open for editing
    layer = ds.GetLayer(0)
    layerDef = layer.GetLayerDefn()

    fields = []
    for i in range(layerDef.GetFieldCount()):
        fname = layerDef.GetFieldDefn(i).GetName()
        fields.append(fname)
    if "Set" not in fields:
        set_field = ogr.FieldDefn("Set", ogr.OFTInteger)
        layer.CreateField(set_field)
        print("Created Field")


    def getHighestValue(FIDs):
        values = {}
        for fid in FIDs:
            subFeat = layer.GetFeature(int(fid))
            values[fid] = subFeat.GetField("GRID_CODE")
        try:
            maxFID = max(values.keys(), key=(lambda key: values[key]))
        except:
            print("no max for fids", FIDs)
            return(0)
        return maxFID


    def getXY(feat):
        geom = feat.GetGeometryRef()
        tpy = geom.Centroid().GetY()
        tpx = geom.Centroid().GetX()
        xy2 = (tpy, tpx)
        return(xy2)

    FIDlist = []
    for feat in layer:
        featID = feat.GetFID()
        FIDlist.append(str(feat.GetFID()))
    tempList = FIDlist[:]
    print("TOTAL FEATURES:",len(FIDlist))
    setCount = 1
    while setCount <= 15:
            three2five = []
            highPointFID = getHighestValue(FIDlist)
            tempList.remove(highPointFID)
            FIDlist = tempList[:]

            try:
                feat1 = layer.GetFeature(int(highPointFID))
            except:
                print("FID not in list", highPointFID)

            xy1 = getXY(feat1)

            for n in FIDlist:
                feat2 = layer.GetFeature(int(n))
                xy2 = getXY(feat2)
                dist = vincenty(xy1, xy2).miles
                if dist >= 2 and dist <= 4:
                    three2five.append(n)
                if dist <= 1: # Remove all points within 1 mile of first point in set
                    tempList.remove(n)
                    layer.DeleteFeature(int(n))

                feat = None

            FIDlist = tempList[:]

            highPointFID_3to5 = getHighestValue(three2five)
            try:
                tempList.remove(highPointFID_3to5)
                FIDlist = tempList[:]
            except:
                print("NO HIGH VALUE", highPointFID_3to5, len(tempList), " remaining")
                continue
            if highPointFID_3to5 == 0:
                print("NO CONTINUE")

            #Iterate over all points and remove those within 1 miles of the second point in set
            feat2 = layer.GetFeature(int(highPointFID_3to5))
            xy1_2 = getXY(feat2)
            for fid in FIDlist:
                feat2_2 = layer.GetFeature(int(fid))
                xy2_2 = getXY(feat2_2)
                dist_2 = vincenty(xy1_2, xy2_2).miles
                if dist_2 <= 1:
                    tempList.remove(fid)
                    layer.DeleteFeature(int(fid))
                    layer.DeleteFeature

            FIDlist = tempList[:]

            feat1.SetField("Set", setCount)
            layer.SetFeature(feat1)
            feat2.SetField("Set", setCount)
            layer.SetFeature(feat2)

            print("SetCount = ", setCount, "-", len(FIDlist))

            setCount += 1

    print("FINISHED",shapefile)
    drv = None
    layer = None
    ds = None

print("FINISHED")