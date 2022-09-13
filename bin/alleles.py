
import db
import json
import re

def getAlleles () :
    q = '''
        SELECT
            aa.accid, a.symbol, a.name
        FROM
            ALL_Allele a,
            ACC_Accession aa
        WHERE
            a._marker_key in (
                select _marker_key
                from mrk_marker
                where _organism_key = 1
                and _marker_type_key = 1
                and _marker_status_key = 1
                )
            and a._allele_key = aa._object_key
            and aa._mgitype_key = 11
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and aa.private = 0
        '''
    return db.sql(q, 'auto')

def getJsonObject (r) :
    obj = {
        "curie" : r["accid"],
        "taxon": "NCBITaxon:10090",
        "internal": False,
        "symbol" : r["symbol"],
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

main()
