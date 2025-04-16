
import sys
import db
import json
import re
from genes import qGenes
from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, getTimeStamp, setCommonFields

def getDiseaseAnnotations (cfg) :
    q = '''
        SELECT
            va._annot_key,
            ve._annotevidence_key,
            av.accid AS doid,
            vt.term AS doterm,
            qt.term AS qualifier,
            ag.accid AS subjectid,
            ra.accid AS mgipubid,
            pma.accid as pmid,
            ve.creation_date,
            ve.modification_date
        FROM
            VOC_Annot va,
            VOC_Term vt,
            ACC_Accession ag,
            ACC_Accession av,
            VOC_Term qt,
            ACC_Accession ra,
            VOC_Evidence ve LEFT JOIN ACC_Accession pma
                ON pma._object_key = ve._refs_key
                AND pma._mgitype_key = 1
                AND pma._logicaldb_key = 29
        WHERE va._annottype_key                  = %(_annottype_key)d
        AND va._qualifier_key = qt._term_key
        AND va._term_key = vt._term_key
        AND ag._object_key = va._object_key
        AND ag._mgitype_key                      = %(_mgitype_key)d
        AND ag._logicaldb_key = 1
        AND ag.preferred = 1
        AND av._object_key = vt._term_key
        AND av._mgitype_key = 13
        AND av._logicaldb_key = 191
        AND av.preferred = 1
        AND va._annot_key = ve._annot_key
        AND ve._refs_key = ra._object_key
        AND ra._mgitype_key = 1
        AND ra._logicaldb_key = 1
        AND ra.accid like 'MGI:%%'
        ''' % cfg
    return db.sql(q, 'auto')

# When we set inferred_gene, want to only include ids for things we actually submitted.
# This function imports/uses the query from genes.py to get this set.
def getSubmittedGeneIds ():
    ids = set()
    for g in db.sql(qGenes, 'auto'):
        ids.add(g['accid'])
    return ids

# Returns a mapping from _annot_key to inferred_allele/inferred_gene.
# The query works by seeing if an _annot_key matches a back-reference for
# a derived annotation. If so, the derived annotation's gene or allele is
# the inferred_gene/inferred_allele for the _annot_key.
def getRollups ():
    # builds a map from annotation key to inferred gene and/or allele ids.

    #
    q = '''
        select _marker_key
        from mrk_marker
        where _marker_type_key = 12
        '''
    mouseTgKeys = set()
    for r in db.sql(q):
        mouseTgKeys.add(r['_marker_key'])

    annotKey2inferred = {}
    #
    # The query runs once for derived gene annotations and once for
    # derived allele annotations.
    cfgs = [{
        '_annottype_key' : 1023,
        '_mgitype_key' : 2,
        '_logicaldb_key' : '1',
        'fieldname' : 'inferred_gene',
    },{
        '_annottype_key' : 1029,
        '_mgitype_key' : 11,
        '_logicaldb_key' : '1',
        'fieldname' : 'inferred_allele',
    },{
        '_annottype_key' : 1032,
        '_mgitype_key' : 2,
        '_logicaldb_key' : '47,64,172,225',
        'fieldname' : 'inferred_gene',
    }]
    #
    q = '''
        select distinct
          va._object_key,
          a1.accid,
          cast(vep.value as integer) as _annot_key
        from
          voc_annot va,
          voc_evidence ve,
          voc_evidence_property vep,
          voc_term vept,
          acc_accession a1
        where va._annot_key = ve._annot_key
        and ve._annotevidence_key = vep._annotevidence_key
        and vep._propertyterm_key = vept._term_key
        and vept.term = '_SourceAnnot_key'
        and va._annottype_key = %(_annottype_key)d
        and a1._object_key = va._object_key
        and a1._mgitype_key = %(_mgitype_key)d
        and a1._logicaldb_key in (%(_logicaldb_key)s)
        and a1.preferred = 1
        '''
    #
    for cfg in cfgs:
        fieldname = cfg['fieldname']
        for r in db.sql(q % cfg, 'auto'):
            ak = r['_annot_key']
            mgiid = r['accid']
            inferreds = annotKey2inferred.setdefault(ak, {})
            if fieldname in inferreds:
                if r['_object_key'] in mouseTgKeys:
                    #sys.stderr.write("Skipping replacement because inferred gene is a Tg.\n")
                    continue
                #sys.stderr.write("Replacing " + fieldname + " " + inferreds[fieldname] + " with " + mgiid + "\n")
            inferreds[fieldname] = mgiid
    #
    return annotKey2inferred
            

def getPrivateCuratorNotes (cfg) :
    q = '''
        SELECT
            n._object_key,
            n.note
        FROM
            MGI_Note n,
            VOC_Evidence ve,
            VOC_Annot va
        WHERE
            n._notetype_key = 1008
        AND n._object_key = ve._annotevidence_key
        AND ve._annot_key = va._annot_key
        AND va._annottype_key                   = %(_annottype_key)d
        ''' % cfg
    ek2note = {}
    for r in db.sql(q, 'auto'):
        ek2note[r['_object_key']] = r['note']
    return ek2note

def getJsonObject (cfg, r, ek2note, annotKey2inferred, submittedGeneIds) :
    unique_id = "MGI:diseaseannotation_%s_%s" % (r['_annot_key'], r['_annotevidence_key'])
    obj = {
      "mod_internal_id" : unique_id,
      "evidence_code_curies": [ "ECO:0000033" ],  # all disease annots use TAS
      "annotation_type_name" : "manually_curated",
      "evidence_curie": "PMID:"+r["pmid"] if r["pmid"] else r["mgipubid"],
      "data_provider_dto": getDataProviderDto(r["subjectid"], cfg["subjecttype"]),
      "do_term_curie": r["doid"],
      "disease_relation_name": cfg["predicate"],
      "negated" : r["qualifier"] == "NOT",
    }
    setCommonFields(r, obj)
    #
    obj[cfg["curie_field"]] = r["subjectid"]
    #
    if r['_annotevidence_key'] in ek2note:
        obj['note_dtos'] = [{
                'free_text' : ek2note[r['_annotevidence_key']],
                'internal'  : True,
                'note_type_name' : 'disease_note'
            }]
    #
    inferred = annotKey2inferred.get(r['_annot_key'], {})
    igene = inferred.get('inferred_gene', None)
    iallele = inferred.get('inferred_allele', None)
    if igene and (igene in submittedGeneIds or not igene.startswith('MGI:')):
        obj['inferred_gene_identifier'] = igene
    if iallele:
        obj['inferred_allele_identifier'] = iallele
    #
    return obj

def main () :
    submittedGeneIds = getSubmittedGeneIds()
    annotKey2inferred = getRollups()
    cfg = {
        "disease_agm_ingest_set": {
            "_annottype_key" : 1020,
            "_mgitype_key"   : 12,
            "predicate"      : "is_model_of",
            "curie_field"    : "agm_identifier",
            "subjecttype"    : "genotype",
        },
        "disease_allele_ingest_set": {
            "_annottype_key" : 1021,
            "_mgitype_key"   : 11,
            "predicate"      : "is_implicated_in",
            "curie_field"    : "allele_identifier",
            "subjecttype"    : "allele",
        }
    }

    print('{')
    print(getHeaderAttributes())
    for i, (section, scfg) in enumerate(cfg.items()):
        ek2note = getPrivateCuratorNotes(scfg)
        if i: print(',', end='')
        print('"%s": [' % section)
        for j,r in mainQuery(getDiseaseAnnotations(scfg)):
            if j: print(',', end='')
            o = getJsonObject(scfg, r, ek2note, annotKey2inferred, submittedGeneIds)
            print(json.dumps(o))
        print(']')
    print('}')

if __name__ == "__main__":
    main()
