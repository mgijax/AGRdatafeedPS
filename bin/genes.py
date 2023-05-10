
import db
import json
import re
from adfLib import LINKML_VERSION, symbolToHtml

def getGenes () :
    q = '''
        SELECT
            aa.accid, mm.symbol, mm.name
        FROM
            MRK_Marker mm,
            ACC_Accession aa
        WHERE
            mm._organism_key = 1
            and mm._marker_type_key = 1
            and mm._marker_status_key = 1
            and mm._marker_key = aa._object_key
            and aa._mgitype_key = 2
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and aa.private = 0
        '''
    return db.sql(q, 'auto')

def getJsonObject (r) :
    obj = {
        "curie" : r["accid"],
        "taxon_curie": "NCBITaxon:10090",
        "internal": False,
        "gene_symbol_dto" : {
	    "name_type_name" : "nomenclature_symbol",
	    "format_text" : r["symbol"],
	    "display_text" : r["symbol"],
	    "internal" : False,
	    },
        "gene_full_name_dto" : {
	    "name_type_name" : "full_name",
	    "format_text" : r["name"],
	    "display_text" : r["name"],
	    "internal" : False
	}
    }
    return obj

def main () :
    print('{')
    print(LINKML_VERSION)
    print('"gene_ingest_set": [')
    for j,r in enumerate(getGenes()):
        if j: print(',', end='')
        o = getJsonObject(r)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
