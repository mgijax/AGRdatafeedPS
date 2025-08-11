
import sys
import db
import json
import re
import argparse
from adfLib import getHeaderAttributes, symbolToHtml, indexResults, getDataProviderDto, mainQuery, log, setCommonFields, getPreferredRefId, getNotesOfType, getNoteDTO
from genes import getSubmittedGeneIds
from constructs import getAlleleConstructRelationships

# Currently only uploading alleles where status is pproved and autoload.
APPROVED_ALLELE_STATUS = 847114
AUTOLOAD_ALLELE_STATUS = 3983021
IN_PROGRESS_ALLELE_STATUS = 847111
DELETED_ALLELE_STATUS = 847112
RESERVED_ALLELE_STATUS = 847113

def getAlleleRefs () :
    q = '''
        SELECT distinct ra._object_key as _allele_key, ra._refs_key, ra._refassoctype_key
        FROM MGI_RefAssocType rat,
          MGI_Reference_Assoc ra
        WHERE ra._refassoctype_key = rat._refassoctype_key
        AND rat._mgitype_key = 11
        AND rat._refassoctype_key != 1014
        '''

    def mapper (r) :
        r['preferredRefId'] = getPreferredRefId(r["_refs_key"])
        return r

    return indexResults(db.sql(q), '_allele_key', None, multi=True, mapper=mapper)

# Returns index from allele MGI id to the _refs_key of its original reference
def getOriginalRefs () :
    q = '''
        SELECT aa.accid as alleleId, aa._object_key as _allele_key, ra._refs_key, ra._refassoctype_key
        FROM 
          MGI_Reference_Assoc ra,
          ACC_Accession aa
        WHERE 1 = 1
        AND ra._refassoctype_key = 1011
        AND aa._object_key = ra._object_key
        AND aa._mgitype_key = 11
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        '''
    return indexResults(db.sql(q), 'alleleId', None, multi=False, mapper=lambda x:getPreferredRefId(x["_refs_key"]) )


def getAlleleTransmission () :
    q = '''
        SELECT a._allele_key, t.term
        FROM all_allele a, voc_term t
        WHERE a._transmission_key = t._term_key
        '''
    return indexResults(db.sql(q), '_allele_key', 'term', multi=False, mapper=lambda t: t.lower())

def getAlleleSynonyms () :
    q = '''
        SELECT s._object_key as _allele_key, s.synonym, s._refs_key
        FROM MGI_Synonym s
        WHERE s._synonymtype_key = 1016
        '''
    mapper = lambda r : (r['synonym'], getPreferredRefId(r['_refs_key']))
    return indexResults(db.sql(q), '_allele_key', None, multi=True, mapper=mapper)

def getAlleleMolecularNotes () :
    return getNotesOfType(1021)

def getAlleleAttributes () :
    q = '''
        SELECT va._object_key as _allele_key, vt.term
        FROM VOC_Annot va, VOC_Term vt
        WHERE va._annottype_key = 1014
        AND va._term_key = vt._term_key
        '''
    return indexResults(db.sql(q), '_allele_key', 'term', multi=True)

def getAlleleMutations () :
    q = '''
        SELECT m._allele_key, t.term
        FROM all_allele_mutation m, voc_term t
        WHERE m._mutation_key = t._term_key
        '''
    return indexResults(db.sql(q), '_allele_key', 'term', multi=True, mapper = lambda s: MUTATION_2_SOID.get(s, None))

def getAlleleSecondaryIds () :
    q = '''
        SELECT _object_key as _allele_key, accid
        FROM acc_accession
        WHERE _mgitype_key = 11
        AND _logicaldb_key = 1
        AND preferred = 0
        '''
    return indexResults(db.sql(q), '_allele_key', 'accid', multi=True)

def getAlleles () :
    q = '''
        SELECT
            a.*,
            aa.accid,
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
            CASE
              WHEN st.term = 'Autoload' THEN 'autoloaded'
              ELSE st.term
            END as status
        FROM
            ALL_Allele a,
            ACC_Accession aa,
            VOC_Term m,
            VOC_Term c,
            VOC_Term st
        WHERE a._allele_status_key in (%d,%d)
            and a._allele_key = aa._object_key
            and aa._mgitype_key = 11
            and aa._logicaldb_key = 1
            and aa.preferred = 1
            and aa.private = 0
            and a._mode_key = m._term_key
            and a._collection_key = c._term_key
            and a._allele_status_key = st._term_key
        ''' % (APPROVED_ALLELE_STATUS, AUTOLOAD_ALLELE_STATUS)
    return db.sql(q, 'auto')

def getAlleleJsonObject (r, ak2refs, ak2trans, ak2syns, ak2attrs, ak2muts, ak2secids, ak2mnotes) :
    refs = ak2refs.get(r["_allele_key"], [])
    origRefs  = list(map(lambda r: r["preferredRefId"], filter(lambda r: r["_refassoctype_key"] == 1011, refs)))
    molecRefs = list(map(lambda r: r["preferredRefId"], filter(lambda r: r["_refassoctype_key"] == 1012, refs)))
    allrefids = list(set(map(lambda r: r["preferredRefId"], refs)))
    allrefids.sort()
    obj = {
        "primary_external_id" : r["accid"],
        "data_provider_dto": getDataProviderDto(r["accid"], "allele"),
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
            "format_text" : symbolToHtml(r["name"]),
            "display_text" : symbolToHtml(r["name"]),
            "internal" : False
        },
        "is_extinct" : (r["isextinct"] == 1),
        "reference_curies" : allrefids,
        "note_dtos" : [],
    }
    setCommonFields(r, obj)
    # inheritance mode
    mode = INHERITANCE_MODE.get(r["mode"], None)
    if mode:
        dto = {
          "inheritance_mode_name": mode,
          "internal": False
        }
        obj["allele_inheritance_mode_dtos"] = [ dto ]
    # collection
    collName = r["collection"].replace(" ","_")
    if collName:
    	obj["in_collection_name"] = collName
    # transmission status
    trans = ak2trans.get(r["_allele_key"], None)
    if trans not in [None, 'not applicable', 'not specified']:
        dto = {
            "germline_transmission_status_name": GERMLINE_TRANS[trans],
            "internal": False
        }
        transRefs = list(filter(lambda r: r["_refassoctype_key"] == 1023, refs))
        if len(transRefs) == 1:
           dto["evidence_curies"] = [ transRefs[0]["preferredRefId"] ]
        obj["allele_germline_transmission_status_dto"] = dto
    # synonyms
    syns = ak2syns.get(r["_allele_key"], None)
    if syns:
        obj["allele_synonym_dtos"] = []
        for s in syns:
            stext = symbolToHtml(s[0])
            pubid = s[1]
            dto = {
                "name_type_name" : "unspecified",
                "format_text" : stext,
                "display_text" : stext,
                "internal" : False,
            }
            if pubid:
                dto["evidence_curies"] = [ pubid ]
            obj["allele_synonym_dtos"].append(dto)
    # status
    st = r['status'].lower()
    obj["allele_database_status_dto"] = {
        "database_status_name": st,
        "internal": False
    }
    # allele attributes (annotation type 1014)
    attrs = ak2attrs.get(r["_allele_key"], [])
    attrs = list(filter(lambda x:x, map(lambda a: FUNC_IMPACT[a], attrs)))
    if len(attrs):
        dto = {
            "functional_impact_names" : attrs,
            "internal" : False
        }
        if len(molecRefs) == 1:
            dto["evidence_curies"] = molecRefs
        obj["allele_functional_impact_dtos"] = [ dto ]
    # allele mutations
    muts = list(filter(lambda x:x, ak2muts.get(r["_allele_key"], [])))
    if len(muts):
        dto = {
            "mutation_type_curies" : muts,
            "internal" : False,    
        }
        if len(molecRefs) == 1:
            dto["evidence_curies"] = molecRefs
        obj["allele_mutation_type_dtos"] = [ dto ]
    # secondary ids
    secids = ak2secids.get(r["_allele_key"], None)
    if secids:
        obj["allele_secondary_id_dtos"] = [{
            "secondary_id" : s,
            "internal" : False,
        } for s in secids]
    # molecular note
    for n in ak2mnotes.get(r["_allele_key"],[]):
        n["evidence_curies"] = molecRefs
        dto = getNoteDTO(n, "mutation_description")
        obj["note_dtos"].append(dto)
    #
    return obj

def getOpts () :
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--type',choices=['alleles','associations'],help="What to output.")
    return parser.parse_args()

def outputAlleles () :
    print('{')
    print(getHeaderAttributes())
    print('"allele_ingest_set": [')
    ak2refs = getAlleleRefs()
    ak2trans = getAlleleTransmission()
    ak2syns = getAlleleSynonyms()
    ak2attrs = getAlleleAttributes()
    ak2muts = getAlleleMutations()
    ak2mnotes = getAlleleMolecularNotes()
    ak2secids = getAlleleSecondaryIds ()
    for j,r in mainQuery(getAlleles()):
        if j: print(',', end='')
        o = getAlleleJsonObject(r, ak2refs, ak2trans, ak2syns, ak2attrs, ak2muts, ak2secids, ak2mnotes)
        print(json.dumps(o))
    print(']')
    print('}')

def getAlleleOfAssociations () :
    q = '''
        SELECT aa.accid as alleleId, ma.accid as markerId, 'is_allele_of' as relationship, null as _refs_key
        FROM ALL_Allele a, MRK_Marker m, ACC_Accession aa, ACC_Accession ma
        WHERE a._marker_key = m._marker_key
        AND m._marker_status_key = 1
        AND m._organism_key = 1
        AND a._allele_key = aa._object_key
        AND aa._mgitype_key = 11
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        AND aa.private = 0
        AND m._marker_key = ma._object_key
        AND ma._mgitype_key = 2
        AND ma._logicaldb_key = 1
        AND ma.preferred = 1
        AND ma.private = 0
        AND a._allele_status_key in (%d,%d)
        ''' % (APPROVED_ALLELE_STATUS, AUTOLOAD_ALLELE_STATUS)
    return db.sql(q)

def getMutationInvolvesAssociations () :
    q = '''
        SELECT aa.accid as alleleId, ma.accid as markerId, rt.term as relationship, ec.abbreviation as evidence, r.*
        FROM MGI_Relationship r, ALL_Allele a, MRK_Marker m, ACC_Accession aa, ACC_Accession ma, VOC_Term rt, VOC_Term ec
        WHERE r._category_key = 1003
        AND r._relationshipterm_key = rt._term_key
        AND r._evidence_key = ec._term_key
        AND r._object_key_1 = a._allele_key
        AND r._object_key_2 = m._marker_key
        AND m._marker_status_key = 1
        AND m._organism_key = 1
        AND a._allele_key = aa._object_key
        AND aa._mgitype_key = 11
        AND aa._logicaldb_key = 1
        AND aa.preferred = 1
        AND aa.private = 0
        AND m._marker_key = ma._object_key
        AND ma._mgitype_key = 2
        AND ma._logicaldb_key = 1
        AND ma.preferred = 1
        AND ma.private = 0
        AND a._allele_status_key in (%d,%d)
        ''' % (APPROVED_ALLELE_STATUS, AUTOLOAD_ALLELE_STATUS)
    return db.sql(q)

def getAlleleGeneAssociations () :
    return getAlleleOfAssociations() + getMutationInvolvesAssociations()

def getAlleleConstructAssociations () :
    aid2rels = getAlleleConstructRelationships()
    aid2oref = getOriginalRefs()
    
    aids = list(aid2rels.keys())
    aids.sort()
    jobjs = []
    for a in aids:
        jobj = {
            "allele_identifier" : a,
            "construct_identifier" : a + "_con",
            "relation_name": "contains",
            "evidence_curies" : [],
            "internal" : False,
        }
        oref = aid2oref.get(a, None)
        if oref:
            jobj["evidence_curies"] = [ oref ]
        jobjs.append(jobj)
    return jobjs

def getGeneAssociationJsonObject (r, geneIds) :
    if not r["markerId"] in geneIds:
        return None

    jobj = {
        "allele_identifier" : r["alleleId"],
        "gene_identifier" : r["markerId"],
        "evidence_curies" : [],
        "internal" : False,
        "obsolete" : False,
        "relation_name" : r["relationship"].replace(" ", "_"),
    }
    if r.has_key("_refs_key"):
        rid = getPreferredRefId(r["_refs_key"])
        if rid:
            jobj["evidence_curies"] = [ rid ]
    if r.has_key("evidence"):
        eco = EVIDENCE_2_ECO.get(r["evidence"], None)
        if eco:
            jobj["evidence_code_curie"] = eco
    setCommonFields(r, jobj)
    return jobj

def outputAssociations () :
    geneIds = getSubmittedGeneIds()
    print('{')
    print(getHeaderAttributes())
    print('"allele_gene_association_ingest_set": [')

    # allele-gene associations
    sep = ''
    for j,r in mainQuery(getAlleleGeneAssociations()):
        o = getGeneAssociationJsonObject(r, geneIds)
        if o:
            print(sep, end='')
            print(json.dumps(o))
            sep = ','
    print(']')

    # allele-construct associations
    sep = ''
    print(',')
    print('"allele_construct_association_ingest_set": [')
    for jobj in getAlleleConstructAssociations():
        print(sep, end='')
        print(json.dumps(jobj))
        sep = ','
    print(']')

    print('}')

def main () :
    opts = getOpts()
    if opts.type == "alleles" :
        outputAlleles()
    elif opts.type == "associations" :
        outputAssociations()

######
# Translation tables 
######

FUNC_IMPACT = dict([
        ("Dominant negative", "dominant_negative_(antimorphic)"),
        ("Constitutively active", "constitutively_active"),
        ("Hypomorph", "hypomorphic_(reduction_of_function)"),
        ("Knockdown", "knockdown"),
        ("Null/knockout", "amorphic_(null/knockout)"),
        ("Not Specified", ""),
        ("Inducible", "inducible"),
        ("Inserted expressed sequence", "inserted_expressed_sequence"),
        ("Modified isoform(s)", "modified_isoform(s)"),
        ("No functional change", "no_functional_change"),
        ("Recombinase", "recombinase"),
        ("Reporter", "reporter"),
        ("RMCE-ready", "RMCE-ready"),
        ("Transactivator", "transactivator"),
        ("Transposase", "transposase"),
        ("Transposon concatemer", "transposon_concatemer"),
        ("Epitope tag", "epitope_tag"),
        ("Endonuclease", "endonuclease"),
        ("Modified regulatory region", "modified_regulatory_region"),
        ("Not Applicable", ""),
        ("Conditional ready", "conditional_ready"),
        ("Humanized sequence", "humanized_sequence"),
        ("Lineage barcode", "lineage_barcode"),
        ("Altered localization", "altered_localization"),
        ("Inducible degradation", "inducible_degradation"),
    ])

GERMLINE_TRANS = dict([
        ("germline", "germline"),
        ("chimeric", "chimeric"),
        ("cell line", "cell_line"),
    ])

INHERITANCE_MODE = {
    "Codominant" : "codominant",
    "Semidominant" : "semi-dominant",
    "Recessive" : "recessive",
    "Dominant" : "dominant",
}

MUTATION_2_SOID = dict([
    ("Intergenic deletion",     "SO:0000159"),
    ("Intragenic deletion",     "SO:0000159"),
    ("Duplication",     "SO:1000035"),
    ("Insertion",       "SO:0000667"),
    ("Insertion of gene trap vector",   "SO:0001218"),
    ("Viral insertion", "SO:0000667"),
    ("Transposon insertion",    "SO:0001837"),
    ("Inversion",       "SO:1000036"),
    ("Nucleotide substitutions",        "SO:0002007"),
    ("Single point mutation",   "SO:1000008"),
    ("Translocation",   "SO:0000199"),
    ("Nucleotide repeat expansion",     "SO:0002162"),
    ("Not Applicable",  None),
    ("Not Specified",   None),
    ("Other",   None),
    ("Undefined",       None),
    ])

EVIDENCE_2_ECO = dict([
    ("IGC", "ECO:0000317"),
    ("IDA", "ECO:0000314"),
    ("ISO", "ECO:0000266"),
    ("IMP", "ECO:0000315"),
    ("EXP", "ECO:0000269"),
    ("IEA", "ECO:0000501"),
    ("IGI", "ECO:0000316"),
    ])

if __name__ == "__main__":
    main()
