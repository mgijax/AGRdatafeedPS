
import db
import json
import re
import argparse
from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, setCommonFields

clonalityMap = {
    "Polyclonal" : "polyclonal",
    "Monoclonal" : "monoclonal",
    "Not Specified" : "not_specified",
}

def getAntibodies () :
    q = '''
        SELECT
            a.*,
            aa.accid,
            c.term as clonality,
            cls.term as heavy_chain
        FROM
            GXD_Antibody a,
            ACC_Accession aa,
            VOC_Term c,
            VOC_Term cls
        WHERE a._antibody_key = aa._object_key
        AND aa._mgitype_key = 6
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        AND a._antibodytype_key = c._term_key
        AND a._antibodyclass_key = cls._term_key
        '''
    return db.sql(q, 'auto')

def getAntibodyGenes () :
    abk2geneIds = {}
    q = '''
        SELECT am._antibody_key, am._marker_key, aa.accid
        FROM gxd_antibodymarker am, acc_accession aa
        WHERE am._marker_key = aa._object_key
        AND aa._mgitype_key = 2
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        '''
    for r in db.sql(q):
        abk2geneIds.setdefault(r['_antibody_key'],[]).append(r['accid'])
    return abk2geneIds

def getAntibodyRefs () :
    abk2refs = {}
    q = '''
        SELECT ra._object_key as _antibody_key, cc.mgiid, cc.pubmedid
        FROM mgi_reference_assoc ra, bib_citation_cache cc
        WHERE ra._refassoctype_key in (1026,1027)
        AND ra._refs_key = cc._refs_key
        '''
    for r in db.sql(q):
        rid = r['mgiid']
        if r['pubmedid'] != None:
            rid = 'PMID:' + str(r['pubmedid'])
        abk2refs.setdefault(r['_antibody_key'],[]).append(rid)
    return abk2refs

def addNoteDTOs (r, obj) :
    mapping = [
        ("regioncovered", "antigen_sequence_note"),
        ("antibodynote", "antibody_note"),
        ("antigennote", "antigen_note"),
    ]
    note_dtos = []
    for (n, n2) in mapping:
        if r[n]:
            note_dtos.append({
                "free_text" : r[n],
                "note_type_name" : n2,
                "internal": False,
            })
    if len(note_dtos) > 0 :
        obj["note_dtos"] = note_dtos

    return obj

def getJsonObject (r, abk2geneIds, abk2refs) :
    obj = {
        "primary_external_id" : r["accid"],
        "name" : r["antibodyname"],
        "clonality_name" : clonalityMap[r["clonality"]],
        "data_provider_dto": getDataProviderDto(r["accid"], "antibody"),
        "internal": False,
        "antibody_target_gene_identifiers" : abk2geneIds.get(r['_antibody_key'],[]),
    }
    refs = abk2refs.get(r['_antibody_key'],[])
    if len(refs) > 0:
        obj["original_reference_curie"] = refs[0]
    if len(refs) > 1:
        obj["reference_curies"] = refs[1:]

    hci = r["heavy_chain"]
    if hci != "Not Applicable" :
        if hci == "Not Specified" :
            hci = "not_specified"
        obj["heavy_chain_isotype_name"] = hci

    addNoteDTOs(r, obj)
    setCommonFields(r, obj)
    return obj

def main () :
    print('{')
    print(getHeaderAttributes())
    print('"antibody_ingest_set": [')  
    first=True          
    abk2geneIds = getAntibodyGenes()
    abk2refs = getAntibodyRefs()
    for j,r in mainQuery(getAntibodies()):
        if j: print(',', end='')
        o = getJsonObject(r, abk2geneIds, abk2refs)
        print(json.dumps(o))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
