
import db
import json
import re
from adfLib import symbolToHtml

def getAGMnames () :
    q = '''
        SELECT g._genotype_key, n.note as alleles
        FROM GXD_Genotype g,
            MGI_Note n
        WHERE
            g._genotype_key = n._object_key
        AND n._notetype_key = 1016
        '''
    d = {}
    for r in db.sql(q, 'auto'):
        d[r['_genotype_key']] = symbolToHtml(r['alleles']).replace('\n', ' ')
    return d

def getAGMs () :
    q = '''
        SELECT
            aa.accid,
            g._genotype_key,
            s.strain
        FROM
            GXD_Genotype g,
            ACC_Accession aa,
            PRB_Strain s
        WHERE
            g._genotype_key > 0  /* skip the not applicable and not specified genotypes */
            and g._genotype_key = aa._object_key
            and aa._mgitype_key = 12
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and g._strain_key = s._strain_key
        '''
    return db.sql(q, 'auto')

def getJsonObject (r, agmKey2name) :
    
    obj = {
        "curie" : r["accid"],
        "name" : (agmKey2name.get(r["_genotype_key"], "") + " [background:] " + r["strain"]).strip(),
        "taxon": "NCBITaxon:10090",
        "internal": False,
        "subtype" : "genotype",
    }
    return obj

def main () :
    agmKey2name = getAGMnames()
    print('{')
    print('"agm_ingest_set": [')
    for j,r in enumerate(getAGMs()):
        if j: print(',', end='')
        o = getJsonObject(r, agmKey2name)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
