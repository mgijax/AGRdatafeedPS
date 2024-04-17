
import db
import json
import re

from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, setCommonFields
#-----------------------------------
# Mapping from MCV term key to SO id.
# Initialize with hard-coded mappings then load what's in the db (which is incomplete).
#
MCV2SO = { 
    #"complex/cluster/region"
    6238175 : "SO:0000110",
    #"cytogenetic marker"
    6238176 : "SO:0000110",
    #"BAC/YAC end"
    6238177 : "SO:0000110",
    #"other genome feature"
    6238178 : "SO:0000110",
    #"DNA segment"
    6238179 : "SO:0000110",
    #"unclassified gene"
    6238184 : "SO:0000704",
    #"other feature type"
    6238185 : "SO:0000110",
    #"unclassified non-coding RNA gene"
    6238186 : "SO:0000704",
    #"unclassified cytogenetic marker"
    7222413 : "SO:0000110",
    #"unclassified other genome feature"
    7648969 : "SO:0000110",
    #"mutation defined region"
    11928467 : "SO:0000110",
}

def initMCV2SO () :
    so_re = re.compile(r'SO:[0-9]+')
    q = '''
            SELECT _term_key, term, note
            FROM VOC_Term
            WHERE _vocab_key = 79 /* MCV */
            ''' 
    for r in db.sql(q):
        m = so_re.search(r['note'])
        if m:
            MCV2SO[r['_term_key']] = m.group(0)

initMCV2SO()

def getGenes () :
    q = '''
        SELECT
            aa.accid, mm.*, mc._mcvterm_key
        FROM
            MRK_Marker mm,
            ACC_Accession aa,
            MRK_MCV_Cache mc
        WHERE
            mm._organism_key = 1
            and mm._marker_type_key = 1
            and mm._marker_status_key = 1
            and mm._marker_key = mc._marker_key
            and mc.qualifier = 'D'
            and mm._marker_key = aa._object_key
            and aa._mgitype_key = 2
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and aa.private = 0
        '''
    return db.sql(q, 'auto')

def getJsonObject (r) :
    obj = {
        "mod_entity_id" : r["accid"],
        "gene_type_curie" : MCV2SO[r['_mcvterm_key']],
        "data_provider_dto": getDataProviderDto(r["accid"], "gene"),
        "taxon_curie": "NCBITaxon:10090",
        "internal": False,
        "gene_symbol_dto" : {
            "name_type_name" : "nomenclature_symbol",
            "format_text" : symbolToHtml(r["symbol"]),
            "display_text" : symbolToHtml(r["symbol"]),
            "internal" : False,
	    },
        "gene_full_name_dto" : {
            "name_type_name" : "full_name",
            "format_text" : symbolToHtml(r["name"]),
            "display_text" : symbolToHtml(r["name"]),
            "internal" : False
        }
    }
    setCommonFields(r, obj)
    return obj

def main () :
    print('{')
    print(getHeaderAttributes())
    print('"gene_ingest_set": [')
    for j,r in mainQuery(getGenes()):
        if j: print(',', end='')
        o = getJsonObject(r)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
