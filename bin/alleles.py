
import db
import json
import re
from adfLib import LINKML_VERSION, symbolToHtml, indexResults

APPROVED_ALLELE_STATUS = 847114
AUTOLOAD_ALLELE_STATUS = 3983021

def getAlleleRefs () :
    q = '''
        SELECT ra._object_key as _allele_key, a.accid as mgiid
        FROM MGI_RefAssocType rat,
          MGI_Reference_Assoc ra,
          ACC_Accession a
        WHERE ra._refassoctype_key = rat._refassoctype_key
        AND rat._mgitype_key = 11
        AND ra._refs_key = a._object_key
        AND a._mgitype_key = 1
        AND a._logicaldb_key = 1
        AND a.preferred = 1
        AND a.prefixpart = 'MGI:'
        '''
    return indexResults(db.sql(q), '_allele_key', 'mgiid', multi=True)

def getAlleles () :
    q = '''
        SELECT
            a._allele_key,
            aa.accid,
            a.symbol,
            a.name,
            CASE
              WHEN m.term = 'Not Specified' THEN ''
              WHEN m.term = 'Not Applicable' THEN ''
              WHEN m.term = 'Not Curated' THEN ''
              ELSE m.term
            END as mode,
            CASE
              WHEN c.term = 'Not Specified' THEN ''
              ELSE c.term
            END as collection,
            a.isextinct 
        FROM
            ALL_Allele a,
            ACC_Accession aa,
            VOC_Term m,
            VOC_Term c
        WHERE a._allele_status_key in (%d,%d)
            and a._allele_key = aa._object_key
            and aa._mgitype_key = 11
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and aa.private = 0
            and a._mode_key = m._term_key
            and a._collection_key = c._term_key
        ''' % (APPROVED_ALLELE_STATUS, AUTOLOAD_ALLELE_STATUS)
    return db.sql(q, 'auto')

def getJsonObject (r, ak2refs) :
    refs = list(set(ak2refs.get(r["_allele_key"], [])))
    refs.sort()
    obj = {
        "curie" : r["accid"],
        "taxon_curie": "NCBITaxon:10090",
        "internal": False,
        "allele_symbol_dto" : {
	    "name_type_name" : "nomenclature_symbol",
	    "format_text" : symbolToHtml(r["symbol"]),
	    "display_text" : symbolToHtml(r["symbol"]),
	    "internal" : False,
	},
        "allele_full_name_dto" : {
	    "name_type_name" : "full_name",
	    "format_text" : r["name"],
	    "display_text" : r["name"],
	    "internal" : False
	},
	"allele_inheritance_mode_dtos": [{
          "inheritance_mode_name": r["mode"].lower(),
          "internal": False
        }],
        "is_extinct" : (r["isextinct"] == 1),
        "reference_curies" : refs
    }
    collName = r["collection"].replace(" ","_")
    if collName:
    	obj["in_collection_name"] = collName
    return obj

def main () :
    print('{')
    print(LINKML_VERSION)
    print('"allele_ingest_set": [')
    ak2refs = getAlleleRefs()
    for j,r in enumerate(getAlleles()):
        if j: print(',', end='')
        o = getJsonObject(r, ak2refs)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
