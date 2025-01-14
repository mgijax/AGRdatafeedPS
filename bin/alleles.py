
import sys
import db
import json
import re
from adfLib import getHeaderAttributes, symbolToHtml, indexResults, getDataProviderDto, mainQuery, log, setCommonFields

APPROVED_ALLELE_STATUS = 847114
AUTOLOAD_ALLELE_STATUS = 3983021

def getReferenceIds () :
    q = '''
        SELECT a1._object_key as _refs_key, a1.accid as mgiid, a2.accid as pubmedid
        FROM acc_accession a1
          LEFT JOIN acc_accession a2
          ON a1._object_key = a2._object_key
          AND a2._mgitype_key = 1
          AND a2._logicaldb_key = 29
          AND a2.preferred = 1
        WHERE a1._mgitype_key = 1
          AND a1.prefixpart = 'MGI:'
          AND a1._logicaldb_key = 1
          AND a1.preferred = 1
        '''
    return indexResults(db.sql(q), '_refs_key', None, multi=False)

rk2ids = getReferenceIds()
def getPreferredRefId (rk):
    if rk is None: return None
    ids = rk2ids[rk] # every reference must have an entry, else error
    if ids["pubmedid"]:
        return "PMID:" + ids["pubmedid"]
    else:
        return ids["mgiid"]

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

def getJsonObject (r, ak2refs, ak2trans, ak2syns, ak2attrs, ak2muts, ak2secids) :
    refs = ak2refs.get(r["_allele_key"], [])
    molecRefs = list(map(lambda r: r["preferredRefId"], filter(lambda r: r["_refassoctype_key"] == 1012, ak2refs.get(r['_allele_key'],[]))))
    allrefids = list(set(map(lambda r: r["preferredRefId"], refs)))
    allrefids.sort()
    obj = {
        "mod_entity_id" : r["accid"],
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
        "reference_curies" : allrefids
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
    #
    return obj

def main () :
    print('{')
    print(getHeaderAttributes())
    print('"allele_ingest_set": [')
    rk2ids = getReferenceIds()
    ak2refs = getAlleleRefs()
    ak2trans = getAlleleTransmission()
    ak2syns = getAlleleSynonyms()
    ak2attrs = getAlleleAttributes()
    ak2muts = getAlleleMutations()
    ak2secids = getAlleleSecondaryIds ()
    for j,r in mainQuery(getAlleles()):
        if j: print(',', end='')
        o = getJsonObject(r, ak2refs, ak2trans, ak2syns, ak2attrs, ak2muts, ak2secids)
        print(json.dumps(o))
    print(']')
    print('}')

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

if __name__ == "__main__":
    main()
