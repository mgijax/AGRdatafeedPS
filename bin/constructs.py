#
# contructs.py
#
# Script to generate construct objects for alleles in MGI for submission to the Alliance curation site.
# Copied and adapted from the original AGRdatafeed product.
#
# A Construct object connects an allele with its engineered components, specifically, expressed and/or driver genes.
# The Alliance model divides components for an allele into those having curie IDs in the db (e.g., human genes) and
# those without (e.g. bacterial genes). A Construct object contains a list of the latter kind
# (ConstructComponentSlotAnnotation), and has a set of associations to the former kind (ConstructGenomicEntityAssociation).
#
# This script produces one of two outputs, depending on the -t command line arg:
#  -t constructs       ConstructDTO object, one per allele having constructs. Contains list of constructs lacking curies.
#  -t associations     ConstructGenomicEntityAssociation, associations between constructs and components with curies.
# (Yes, you have to run the script twice to get the complete output.)
#
# The Alliance model wants an identifier for each construct. MGI does not have
# actual construct objects, much less identifiers for them. So this script creates ersatz
# IDs for the constructs, based on the owning allele's id. If an allele's id is MGI:123456,
# its contruct has the id MGI:123456_con. This works because an allele only has at most a single construct 
# (the way we're doing things).
#

import sys
import db
import json
import re
import argparse
from adfLib import getHeaderAttributes, symbolToHtml, indexResults, getDataProviderDto, mainQuery, log, setCommonFields

EXPRESSES_cat_key = 1004
DRIVER_cat_key = 1006

mk2nmdId = {} # marker key -> non-mouse ID
def loadNonMouseGeneIds () :
    for r in db.sql(qConstructNonMouseComponents):
        mk2nmdId[r['_marker_key']] = r['accid']

rk2id = {} # _refs_key -> either PMID or MGI id
def loadRefIds () :
    for r in db.sql(qConstructRefs):
        rk2id[r['_refs_key']] = ('PMID:' + r['pmid']) if r['pmid'] else r['mgiid']
        
rk2note = {} # _relationship_key -> note obj
def loadConstructNotes () :
    for n in db.sql(qConstructNotes):
        rk2note[n['_relationship_key']] = n
        # some notes mistakenly surrounded by double quote characters
        # remove them here
        if n['note'].startswith('"') and n['note'].endswith('"'):
            n['note'] = n['note'][1:-1]
    log("Loaded notes for %d constructs." % len(rk2note))

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

    note = rk2note.get(r["_relationship_key"], None)
    note_dto = None
    if note:
        note_dto = {
          "free_text": note["note"],
          "note_type_name": "comment",
        }
        setCommonFields(note, note_dto)


    # NOTE: additional constraint that the marker type is gene. Only needed while 
    # Alliance persistent store contains only genes. When Allianec contains all markers,
    # can remove the added constraint.
    if gid and r["_marker_type_key"] == 1:
        # we have an id for the gene. Return a ConstructGenomicEntityAssociationDTO
        rval = {
          "construct_identifier": construct_id,
          "genomic_entity_relation_name": reln,
          "genomic_entity_identifier": gid,
          "evidence_curies": [
            rk2id[r['_refs_key']]
          ],  
        }
        if note_dto: rval["note_dtos"] = [ note_dto ]
        setCommonFields(r, rval)
        return ("ConstructGenomicEntityAssociationDTO", rval)
    else:
        # don't have a curie for the gene. Return a ConstructComponentSlotAnnotationDTO
        rval = {
            "relation_name" : reln,
            "component_symbol" : symbol,
            "evidence_curies": [
              rk2id[r['_refs_key']]
            ],  
        }
        if note_dto: rval["note_dtos"] = [ note_dto ]
        setCommonFields(r, rval)
        if r["taxonid"]:
            rval["taxon_curie"] = "NCBITaxon:" + r["taxonid"]
        if r["commonname"]:
            rval["taxon_text"] = r["commonname"]
        return ("ConstructComponentSlotAnnotationDTO", rval)
    
def getOpts () :
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--type',choices=['constructs','associations'],help="What to output.")
    return parser.parse_args()

# Returns mapping from allele MGI id to list of relationship records for components.
def getAlleleConstructRelationships () :
    aid2rels = {}
    for r in loadRelationship(EXPRESSES_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    for r in loadRelationship(DRIVER_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    return aid2rels
    
def main () :
    opts = getOpts()
    loadNonMouseGeneIds()
    loadRefIds()
    loadConstructNotes()
    # Get all relationship records for alleles (expresses-component and driven-by).
    # Then aggregate them into a single list of components per allele.
    aid2rels = getAlleleConstructRelationships()
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

        minCreatedDate = None
        minCreatedBy = None
        maxUpdatedDate = None
        maxUpdatedBy = None
       
        for arel in arels:
           tp, obj = rel2constrComp(arel, construct_id)
           if tp == "ConstructComponentSlotAnnotationDTO":
               ccomps.append(obj)
           elif tp == "ConstructGenomicEntityAssociationDTO":
               cgassocs.append(obj)
           else:
               raise RuntimeError("Unknown type: " + str(tp))

           if not minCreatedDate or (obj["date_created"] < minCreatedDate):
               minCreatedDate = obj["date_created"]
               minCreatedBy = obj["created_by_curie"]

           if not maxUpdatedDate or (obj["date_updated"] > maxUpdatedDate):
               maxUpdatedDate = obj["date_updated"]
               maxUpdatedBy = obj["created_by_curie"]

        if opts.type == "associations":
            for a in cgassocs:
                if not first: print(",", end=' ')
                first = False
                print(json.dumps(a, indent=2))
            continue

        # else opts.type == "constructs"...
        symbol = arels[0]["allelesymbol"] + ' construct'
        obj = {
          "internal" : False,
          "obsolete" : False,
          "date_created" : minCreatedDate,
          "created_by_curie" : minCreatedBy,
          "date_updated" : maxUpdatedDate,
          "updated_by_curie" : maxUpdatedBy,
          "mod_internal_id" : construct_id,
          "construct_symbol_dto" : {
              "name_type_name": "nomenclature_symbol",
              "format_text": symbol,
              "display_text": symbol,
              "internal": False,
          },
          "data_provider_dto": getDataProviderDto(aid, "allele"),
        }
        if len(ccomps): obj["construct_component_dtos"] = ccomps
        #
        if not first: print(",", end=' ')
        print(json.dumps(obj, indent=2))
        first=False
        #
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
        mm._marker_type_key,
        am.accid as mgiid,
        mm.symbol as genesymbol,
        rr.term as relationship,
        q.term as qualifier,
        e.abbreviation as evidencecode,
        r._refs_key,

        r.creation_date,
        r.modification_date,
        r._createdby_key,
        r._modifiedby_key
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

# query to get accids for non-mouse component genes 
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

# query to return all construct association references
qConstructRefs = '''
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

# query to returns notes attached to construct associations
qConstructNotes = '''
    SELECT r._relationship_key, n.*
    FROM mgi_note n, mgi_relationship r
    WHERE n._notetype_key = 1042
    AND n._object_key = r._relationship_key
    AND r._category_key IN (1004,1006)
'''

if __name__ == "__main__":
    main()
