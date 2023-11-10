#
# contructs.py
#
# Script to generate construct objects for alleles in MGI for submission to the Alliance curation site.
# Copied and adapted from the original AGRdatafeed product.
#
# Alleles that have expressed components and/or driver genes are considered to
# have constructs. One construct is created per allele, which  may contain multiple components.
# (e.g., multiple expressed components, or an expressed component and a driver, or ...)
#
# The alliance curation schema wants an identifier for each construct. MGI does not have
# actual construct objects, much less identifiers for them. So this script creates ersatz
# IDs for the constructs, based on the owning allele's id. If an allele's id is MGI:123456,
# its contruct has the id MGI:123456_con
#

import sys
import db
import json
import re
import argparse
from adfLib import getHeaderAttributes, symbolToHtml, indexResults, getDataProviderDto, mainQuery, log

EXPRESSES_cat_key = 1004
DRIVER_cat_key = 1006

mk2nmdId = {} # marker key -> non-mouse ID
def loadNonMouseGeneIds () :
    for r in db.sql(qConstructNonMouseComponents):
        mk2nmdId[r['_marker_key']] = r['accid']

rk2id = {}
def loadRefIds () :
    for r in db.sql(qRefs):
        rk2id[r['_refs_key']] = ('PMID:' + r['pmid']) if r['pmid'] else r['mgiid']
        
def loadRelationship (key) :
    rels = []
    # read the relationships 
    for r in db.sql(tConstructRelationships % key):
        rk = r['_relationship_key']
        rels.append(r)
    return rels

# Given an allele-gene relationship record, returns a JSON object
# If the gene has a primary ID, returns a ConstructGenomicEntityAssociationDTO
# Otherwise, returns a ConstructComponentSlotAnnotationDTO
# (types defined in the Alliance schema)
def rel2constrComp (r, construct_id) :
    symbol = r["genesymbol"]
    gid = r["mgiid"]
    if gid is None:
        gid = mk2nmdId.get(r["_marker_key"],None)
        if gid :
            if gid.startswith('ZDB'):
                gid = 'ZFIN:' + gid
            elif gid.startswith('XB-'):
                gid = 'Xenbase:' + gid

    if r["_category_key"] == EXPRESSES_cat_key:
        reln = "expresses"
    elif r["_category_key"] == DRIVER_cat_key:
        reln = "is_regulated_by"
    else:
        raise RuntimeError("Internal error: unknown _category_key: " + str(r))

    if gid:
        # we have an id for the gene. Return a ConstructGenomicEntityAssociationDTO
        rval = {
          "construct_identifier": construct_id,
          "genomic_entity_relation_name": reln,
          "genomic_entity_curie": gid,
          "evidence_curies": [
            rk2id[r['_refs_key']]
          ],  
          #"note_dtos": [
          #  {   
          #    "free_text": "note about association",
          #    "note_type_name": "comment",
          #    "internal": false
          #      }   
          #],  
          "internal": False,
          "obsolete": False,
          #"date_created": "2015-04-09T10:15:30+00:00",
          #"date_updated": "2022-07-11T13:12:51+00:00",
          #"created_by_curie": "WB:WBPerson002314",
          #"updated_by_curie": "WB:WBPerson010002"
        }
        return ("ConstructGenomicEntityAssociationDTO", rval)
    else:
        # don't have a curie for the gene. Return a ConstructComponentSlotAnnotationDTO
        rval = {
            "relation_name" : reln,
            "component_symbol" : symbol,
            "internal" : False,
            "evidence_curies": [
              rk2id[r['_refs_key']]
            ],  
        }
        if r["taxonid"]:
            rval["taxon_curie"] = "NCBITaxon:" + r["taxonid"]
        if r["commonname"]:
            rval["taxon_text"] = r["commonname"]
        return ("ConstructComponentSlotAnnotationDTO", rval)
    
def getOpts () :
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--type',choices=['constructs','associations'],help="What to output.")
    return parser.parse_args()

def main () :
    opts = getOpts()
    loadNonMouseGeneIds()
    loadRefIds()
    #
    aid2rels = {}
    for r in loadRelationship(EXPRESSES_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    for r in loadRelationship(DRIVER_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    aids = list(aid2rels.keys())
    #
    print('{')
    print(getHeaderAttributes())
    if opts.type == "constructs":
        print('"construct_ingest_set": [')
    else:
        print('"construct_genomic_entity_association_ingest_set": [')
    first=True
    for aid in aids:
        construct_id = aid + '_con'
        arels = aid2rels[aid]
        ccomps = []
        cgassocs = []
        for arel in arels:
           tp, obj = rel2constrComp(arel, construct_id)
           if tp == "ConstructComponentSlotAnnotationDTO":
               ccomps.append(obj)
           elif tp == "ConstructGenomicEntityAssociationDTO":
               cgassocs.append(obj)
           else:
               raise RuntimeError("Unknown type: " + str(tp))

        if opts.type == "associations":
            for a in cgassocs:
                if not first: print(",", end=' ')
                first = False
                print(json.dumps(a, indent=2))
            continue

        # if opts.type == "constructs"...
        symbol = arels[0]["allelesymbol"] + ' construct'
        obj = {
          "internal" : False,
          "mod_internal_id" : construct_id,
          "construct_symbol_dto" : {
              "name_type_name": "nomenclature_symbol",
              "format_text": symbol,
              "display_text": symbol,
              # "evidence_curies": [], #PMIDs
              "internal": False,
          },
          "construct_component_dtos": ccomps,
          "data_provider_dto": getDataProviderDto(aid, "allele"),
        }
        if not first: print(",", end=' ')
        print(json.dumps(obj, indent=2))
        first=False
    print("]}")

#

# ------------------------------------------------
# QUERIES.
# ------------------------------------------------
tConstructRelationships = ''' 
    SELECT 
        r._relationship_key,
        r._category_key,
        al._allele_key,
        aa.accid as allele,
        al.symbol as alleleSymbol,
        mm._organism_key,
        mo.commonname,
        moa.accid as taxonid,
        mm._marker_key,
        am.accid as mgiid,
        mm.symbol as genesymbol,
        rr.term as relationship,
        q.term as qualifier,
        e.abbreviation as evidencecode,
        r._refs_key
    FROM
        MGI_Relationship r
        LEFT JOIN ACC_Accession am 
            /* Have to outer join to get the MGI id since not all of these will be mouse genes. */
            ON r._object_key_2 = am._object_key
            AND am._mgitype_key = 2
            AND am._logicaldb_key = 1
            AND am.preferred = 1
        JOIN VOC_Term q
            ON r._qualifier_key = q._term_key
        JOIN VOC_Term e
            ON r._evidence_key = e._term_key
        JOIN VOC_Term rr
            ON r._relationshipterm_key = rr._term_key
        JOIN ACC_Accession aa
            ON r._object_key_1 = aa._object_key
            AND aa._mgitype_key = 11
            AND aa._logicaldb_key = 1
            AND aa.preferred = 1
        JOIN ALL_Allele al
            ON r._object_key_1 = al._allele_key
        JOIN MRK_Marker mm
            ON r._object_key_2 = mm._marker_key
        LEFT JOIN MGI_Organism mo
            ON mm._organism_key = mo._organism_key
        LEFT JOIN ACC_Accession moa
            ON mo._organism_key = moa._object_key
            AND moa._mgitype_key = 20
            AND moa._logicaldb_key = 32
    WHERE r._category_key = %s
    ORDER BY r._relationship_key
    '''
# query for relationship properties. Will get attached as list to relationship.
# Arg: category key
tConstructProperties = '''
    SELECT r._relationship_key, t.term as property, p.value
    FROM MGI_Relationship r
      JOIN MGI_Relationship_Property p
        ON r._relationship_key = p._relationship_key
      JOIN VOC_Term t
        ON p._propertyname_key = t._term_key
    WHERE r._category_key = %s
    ORDER BY r._relationship_key, p.sequenceNum
    '''

# query to get accids for non-mouse drivers 
qConstructNonMouseComponents = '''
    SELECT distinct a.accid, m._marker_key
    FROM
        MGI_Relationship r,
        MRK_Marker m,
        ACC_Accession a
    WHERE
        r._category_key in (1004,1006)
    AND r._object_key_2 = m._marker_key
    AND m._organism_key != 1
    AND a._object_key = m._marker_key
    AND a._mgitype_key = 2
    AND a._logicaldb_key in (64,47,172,225) /* HGNC, RGD, ZFIN, Xenbase */
    AND a.preferred = 1
    '''
qRefs = '''
    SELECT DISTINCT r._refs_key, a1.accid as mgiid, a2.accid as pmid
    FROM
        MGI_Relationship r
        JOIN ACC_Accession a1
            ON r._refs_key = a1._object_key
            AND a1._mgitype_key = 1
            AND a1._logicaldb_key = 1
            AND a1.preferred = 1
            AND a1.prefixPart = 'MGI:'
        LEFT JOIN ACC_Accession a2
            ON r._refs_key = a2._object_key
            AND a2._mgitype_key = 1
            AND a2._logicaldb_key = 29
            AND a2.preferred = 1
    WHERE r._category_key in (1004,1006)
'''

main()
