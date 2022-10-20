
import db
import json
import re
from adfLib import symbolToHtml

APPROVED_ALLELE_STATUS = 847114
AUTOLOAD_ALLELE_STATUS = 3983021

def getAlleles () :
    q = '''
        SELECT
            aa.accid, a.symbol, a.name
        FROM
            ALL_Allele a,
            ACC_Accession aa
        WHERE a._allele_status_key in (%d,%d)
            and a._allele_key = aa._object_key
            and aa._mgitype_key = 11
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and aa.private = 0
        ''' % (APPROVED_ALLELE_STATUS, AUTOLOAD_ALLELE_STATUS)
    return db.sql(q, 'auto')

def getJsonObject (r) :
    obj = {
        "curie" : r["accid"],
        "taxon": "NCBITaxon:10090",
        "internal": False,
        "symbol" : symbolToHtml(r["symbol"]),
        "name" : r["name"],
    }
    return obj

def main () :
    print('{')
    print('"allele_ingest_set": [')
    for j,r in enumerate(getAlleles()):
        if j: print(',', end='')
        o = getJsonObject(r)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
