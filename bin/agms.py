
import db
import json
import re
import argparse
from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, setCommonFields

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
            g.*,
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

#
# valid GENO term ids for Alliance submissions
# See: http://www.ontobee.org/ontology/GENO
#
genoZygosityTerms = [ 
("GENO:0000602", "homoplasmic"),
("GENO:0000603", "heteroplasmic"),
("GENO:0000604", "hemizygous X-linked"),
("GENO:0000605", "hemizygous Y-linked"),
("GENO:0000606", "hemizygous insertion-linked"),
("GENO:0000135", "heterozygous"),
("GENO:0000136", "homozygous"),
("GENO:0000137", "unspecified zygosity"),
]

mgi2geno = { 
 "Homozygous" : "GENO:0000136", 
 "Heterozygous" : "GENO:0000135",
 "Hemizygous Insertion" : "GENO:0000606",
 "Indeterminate" : "GENO:0000137",
 "Hemizygous X-linked" : "GENO:0000604",
 "Hemizygous Y-linked" : "GENO:0000605",
 "Hemizygous Deletion" : "GENO:0000134", # no term for hemizygous deletion. Use plain hemizygous.
 "Homoplasmic" : "GENO:0000602",
 "Heteroplasmic" : "GENO:0000603"
}

def getAGMComponents () :
    q = '''
        SELECT
          ap._genotype_key,
          aa.accid,
          ps.term
        FROM 
          GXD_AllelePair ap,
          ALL_Allele a1,
          ACC_Accession aa,
          VOC_Term ps
        WHERE ap._allele_key_1 = a1._allele_key
        AND a1._allele_key = aa._object_key
        AND aa._mgitype_key = 11
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        AND ap._pairstate_key = ps._term_key
        '''
    gk2comps = {}
    for r in db.sql(q, 'auto'):
        gk2comps.setdefault(r["_genotype_key"],[]).append({
            "allele_identifier" : r["accid"],
            "zygosity_curie" : mgi2geno[r["term"]],
            "internal" : False
        })
    return gk2comps

def getJsonObject (r, agmKey2name) :
    
    obj = {
        "primary_external_id" : r["accid"],
        "name" : (agmKey2name.get(r["_genotype_key"], "") + " [background:] " + symbolToHtml(r["strain"])).strip(),
        "taxon_curie": "NCBITaxon:10090",
        "data_provider_dto": getDataProviderDto(r["accid"], "genotype"),
        "internal": False,
        "subtype_name" : "genotype"
    }
    setCommonFields(r, obj)
    return obj

def getOpts () :
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--type',choices=['genotypes','associations'],help="What to output.")
    return parser.parse_args()

def main () :
    opts = getOpts()
    agmKey2name = getAGMnames()
    genoKey2comps = getAGMComponents()
    print('{')
    print(getHeaderAttributes())
    if opts.type == "genotypes":
        print('"agm_ingest_set": [')
    else:
        print('"agm_allele_association_ingest_set": [')  
    first=True          
    for j,r in mainQuery(getAGMs()):
        if opts.type == "associations":
            obj = {
                "agm_subject_identifier": r["accid"],
                "relation_name": "has_component",
            }
            objs = genoKey2comps.get(r["_genotype_key"], [])
            if objs:
                for obj in objs:
                    obj["agm_subject_identifier"] = r["accid"]
                    obj["relation_name"] = "has_component"
                    if not first: print(",", end=' ')
                    first = False
                    print(json.dumps(obj))                                        
            continue

        # else opts.type == "genotypes"...
        if j: print(',', end='')
        o = getJsonObject(r, agmKey2name)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
