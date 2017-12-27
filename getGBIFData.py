import requests, json, csv, os

years = "2007 TO 2016"

outpath = r"../bisonGBIF_datatemp.csv"

#Families dictionary of butterfly families and their hierarchy homonym string code
families = {#"Hedylidae":"None",
            "Hesperiidae": "117278",
            "Lycaenidae": "117236",
            "Nymphalidae": "117276",
            "Papilionidae": "117268",
            "Pieridae": "117267",
            "Riodinidae": "117266"}

# set number of rows to something very high.  Will be max query return
numRows = str(10000000)

# construct RESTful query
url = "https://bison.usgs.gov/solr/occurrences/select?"
url += "q=hierarchy_homonym_string:(*-h_h_s-*)" \
       " AND year:[" + years + "]" \
       " AND decimalLongitude:[-124.762795 TO -98.210199]" \
       " AND decimalLatitude:[26.391708 TO 49.003595]" \
       " AND coordinatePrecision:[* TO 10]" \
       " AND coordinateUncertaintyInMeters:[* TO 10]"
url += "&wt=json&stored=true&rows=" + numRows
url += "&fl=" \
        "decimalLatitude," \
        "decimalLongitude," \
        "eventDate," \
        "ownerInstitutionCollectionCode"\
        "ITISscientificName," \
        "coordinateUncertaintyInMeters," \
        "coordinatePrecision," \
        "collectionID"

for k,v in families.items():
    url = url.replace("h_h_s",v)
    #print(url)

    response = requests.get(url, stream=True)
    if response.status_code == 200:
        dict = json.loads(response.content.decode())
        return_dict = dict["response"]["docs"]

        columns = [x for row in return_dict for x in row.keys()]
        columns = list(set(columns))
        #print(columns)
        def writevalues(csvf, dict):
            for i_r in dict:
                #print("i_r",i_r)
                csvf.writerow(map(lambda x: i_r.get(x, ""), columns))

        if os.path.exists(outpath):
            with open(outpath, 'a+') as ofile:
                csv_w = csv.writer(ofile, lineterminator='\n')
                writevalues(csv_w, return_dict)
        else:
            with open(outpath, 'w') as ofile:
                csv_w = csv.writer(ofile, lineterminator='\n')
                csv_w.writerow(columns)
                writevalues(csv_w, return_dict)

    print("Done with", k)

    # Loading the response data into a dict variable
    # json.loads takes in only binary or string variables so using content to fetch binary content
    # Loads (Load String) takes a Json file and converts into python data structure (dict or list, depending on JSON)

    #print("\n")
    #print(jData)
    #for key in jData:
    #    print(key + " : ", jData[key])
else:
    # If response code is not ok (200), print the resulting http error code with description
    response.raise_for_status()

del response

print("FINISHED")