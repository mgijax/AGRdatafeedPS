
import db
import json
import re

def getDiseaseAnnotations (cfg) :
    q = '''
        SELECT
            va._annot_key,
            ve._annotevidence_key,
            av.accid AS doid,
            vt.term AS doterm,
            qt.term AS qualifier,
            ag.accid AS genoid,
            ra.accid AS mgipubid,
            pma.accid as pmid,
            ve.creation_date,
            ve.modification_date
        FROM
            VOC_Annot va,
            GXD_Genotype gg,
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
        AND va._object_key = gg._genotype_key
        AND va._term_key = vt._term_key
        AND ag._object_key = gg._genotype_key
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

def formatDate (s) :
    s = re.sub('\.[0-9]*$', '', s)
    s = s.replace(" ", "T")
    s = s + "Z"
    return s

def getJsonObject (cfg, r, ek2note) :
    unique_id = "MGI:diseaseannotation_%s_%s" % (r['_annot_key'], r['_annotevidence_key'])
    obj = {
      "mod_entity_id" : unique_id,
      "internal": False,
      "evidence_codes": [ "ECO:0000033" ],  # all disease annots use traceable author statement (TAS) codes
      "annotation_type" : "manually_curated",
      "single_reference": "PMID:"+r["pmid"] if r["pmid"] else r["mgipubid"],
      "data_provider": "MGI",  
      "object": r["doid"],
      "created_by": "MGI:curation_staff",
      "updated_by": "MGI:curation_staff",
      "subject": r["genoid"],
      "predicate": cfg["predicate"],
      "negated" : r["qualifier"] == "NOT",
      "date_created" : formatDate(r["creation_date"]),
      "date_updated" : formatDate(r["modification_date"])
    }
    if r['_annotevidence_key'] in ek2note:
        obj['related_notes'] = [{
                'free_text' : ek2note[r['_annotevidence_key']],
                'internal'  : True,
                'note_type' : 'disease_note'
            }]
    return obj

def main () :
    cfg = {
        "disease_agm_ingest_set": {
            "_annottype_key" : 1020,
            "_mgitype_key"   : 12,
            "predicate"      : "is_model_of"
        },
        "disease_allele_ingest_set": {
            "_annottype_key" : 1021,
            "_mgitype_key"   : 11,
            "predicate"      : "is_implicated_in"
        }
    }

    print('{')
    for i, (section, scfg) in enumerate(cfg.items()):
        ek2note = getPrivateCuratorNotes(scfg)
        if i: print(',', end='')
        print('"%s": [' % section)
        for j,r in enumerate(getDiseaseAnnotations(scfg)):
            if j: print(',', end='')
            o = getJsonObject(scfg, r, ek2note)
            print(json.dumps(o))
        print(']')
    print('}')

main()
