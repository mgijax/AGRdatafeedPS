
import db
import json
import re
import argparse
from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, setCommonFields, getPreferredRefId
from allianceTaxonCV import getTaxonByOrganismKey

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
            cls.term as heavy_chain,
            p._organism_key as _antigen_organism_key
        FROM
            GXD_Antibody a,
            PRB_Source p,
            ACC_Accession aa,
            VOC_Term c,
            VOC_Term cls
        WHERE a._antibody_key = aa._object_key
        AND aa._mgitype_key = 6
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        AND a._antibodytype_key = c._term_key
        AND a._antibodyclass_key = cls._term_key
        AND a._source_key = p._source_key
        '''
    return db.sql(q, 'auto')

def getAntibodySynonyms () :
    abk2synonyms = {}
    q = '''
        SELECT _antibody_key, alias as synonym, _refs_key
        FROM gxd_antibodyalias
        '''
    for r in db.sql(q) :
        abk2synonyms.setdefault(r['_antibody_key'],[]).append(r)
    return abk2synonyms

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
        FROM mgi_reference_assoc ra, bib_citation_cache cc, bib_refs br
        WHERE ra._refassoctype_key in (1026,1027)
        AND ra._refs_key = cc._refs_key
        AND ra._refs_key = br._refs_key
        ORDER BY br.year, br.creation_date
        '''
    for r in db.sql(q):
        rid = r['mgiid']
        if r['pubmedid'] != None:
            rid = 'PMID:' + str(r['pubmedid'])
        abk2refs.setdefault(r['_antibody_key'],[]).append(rid)
    return abk2refs

def getOrganismTaxons () :
    ok2taxonid = {}
    q = '''
        SELECT o._organism_key, a.accid
        FROM MGI_Organism o, ACC_Accession a
        WHERE a._object_key = o._organism_key
        AND a._mgitype_key = 20
        AND a._logicaldb_key = 32
        '''
    for r in db.sql(q) :
        ok2taxonid[r['_organism_key']] = 'NCBITaxon:' + r['accid']
    #
    return ok2taxonid

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

def addSynonymDTOs(r, obj, abk2synonyms) :
    syns = abk2synonyms.get(r['_antibody_key'], None)
    if syns:
        obj["antibody_synonym_dtos"] = []
        for s in syns:
            dto = {
                "name_type_name": "synonym",
                "format_text": s["synonym"],
                "display_text": s["synonym"],
                "synonym_scope_name": "exact",
                "internal": False
            }
            if s["_refs_key"] :
                dto["evidence_curies"] = [
                    getPreferredRefId(s["_refs_key"])
                ]
            obj["antibody_synonym_dtos"].append(dto)
    return obj

def getJsonObject (r, abk2geneIds, abk2refs, abk2synonyms, ok2taxonid) :
    obj = {
        "primary_external_id" : r["accid"],
        "name" : r["antibodyname"],
        "clonality_name" : clonalityMap[r["clonality"]],
        "data_provider_dto": getDataProviderDto(),
        "internal": False,
        "antibody_target_gene_identifiers" : abk2geneIds.get(r['_antibody_key'],[]),
    }
    #
    taxid = getTaxonByOrganismKey(r['_organism_key'])
    if taxid:
        obj['host_taxon_term_name'] = taxid
    #
    agtaxid = getTaxonByOrganismKey(r['_antigen_organism_key'])
    if agtaxid:
        obj['antigen_taxon_term_name'] = agtaxid
    #
    #addSynonymDTOs(r, obj, abk2synonyms)

    #
    refs = abk2refs.get(r['_antibody_key'],[])
    if len(refs) > 0:
        obj["original_reference_curie"] = refs[0]
    if len(refs) > 1:
        obj["reference_curies"] = refs[1:]
    #
    hci = r["heavy_chain"]
    if hci != "Not Applicable" :
        if hci == "Not Specified" :
            hci = "not_specified"
        obj["heavy_chain_isotype_name"] = hci
    #
    addNoteDTOs(r, obj)
    #
    setCommonFields(r, obj)
    return obj

def main () :
    print('{')
    print(getHeaderAttributes())
    print('"antibody_ingest_set": [')  
    first=True          
    abk2geneIds = getAntibodyGenes()
    abk2refs = getAntibodyRefs()
    abk2synonyms = getAntibodySynonyms()
    ok2taxonid = getOrganismTaxons()
    for j,r in mainQuery(getAntibodies()):
        if j: print(',', end='')
        o = getJsonObject(r, abk2geneIds, abk2refs, abk2synonyms, ok2taxonid)
        print(json.dumps(o, indent=2))
    print(']')
    print('}')

if __name__ == "__main__":
    main()
