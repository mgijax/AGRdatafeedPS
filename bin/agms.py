
import db
import json
import re

def getAGMs () :
    q = '''
        SELECT
            aa.accid
        FROM
            GXD_Genotype g,
            ACC_Accession aa
        WHERE
            g._genotype_key = aa._object_key
            and aa._mgitype_key = 12
            and aa._logicaldb_key = 1
            and aa.preferred = 1
        '''
    return db.sql(q, 'auto')

def getJsonObject (r) :
    obj = {
        "curie" : r["accid"],
        "taxon": "NCBITaxon:10090",
        "internal": False,
        "subtype" : "genotype",
    }
    return obj

def main () :
    print('{')
    print('"agm_ingest_set": [')
    for j,r in enumerate(getAGMs()):
        if j: print(',', end='')
        o = getJsonObject(r)
        print(json.dumps(o))
    print(']')
    print('}')

main()
