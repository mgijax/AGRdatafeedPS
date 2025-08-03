
import db
import json
import re

from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, setCommonFields

# ----------------------------------------------------------
# Mapping from MCV term key to SO id.
# Initialize with hard-coded mappings then load what's in the db (which is incomplete).
# ----------------------------------------------------------

MCV2SO = { 
    #"complex/cluster/region" -> "biological_region"
    6238175 : "SO:0001411",
    #"cytogenetic marker"
    #6238176 : "SO:0000110",
    #"BAC/YAC end" -> "chromosome_part"
    6238177 : "SO:0000830",
    #"other genome feature" -> "biological_region"
    6238178 : "SO:0001411",
    #"DNA segment" -> "biological_region"
    6238179 : "SO:0001411",
    #"unclassified gene" -> "gene"
    6238184 : "SO:0000704",
    #"other feature type" -> "biological_region"
    6238185 : "SO:0001411",
    #"unclassified non-coding RNA gene" -> "ncRNA_gene"
    6238186 : "SO:0001263",
    #"unclassified cytogenetic marker"
    #7222413 : "SO:0000110",
    #"unclassified other genome feature" -> "biological_region"
    7648969 : "SO:0001411",
    #"mutation defined region" -> "heritable_phenotypic_marker"
    11928467 : "SO:0001500",
}

def initMCV2SO () :
    so_re = re.compile(r'SO:[0-9]+')
    q = '''
            SELECT _term_key, term, note
            FROM VOC_Term
            WHERE _vocab_key = 79 /* MCV */
            ''' 
    for r in db.sql(q):
        m = so_re.search(r['note'])
        if m and r['_term_key'] not in MCV2SO:
             MCV2SO[r['_term_key']] = m.group(0)

# ----------------------------------------------------------
# ----------------------------------------------------------

# Returns a dictionary from name to set of marker keys.
def getGeneSets () :
    sets = [
        (qGeneHasPhenotype,"hasPhenotype"),
        (qGeneHasImpc,"hasImpc"),
        (qGeneHasExpression,"hasExpression"),
        (qGeneHasExpressionImage, "hasExpressionImage"),
        ]
    rval = {}
    for (q,n) in sets:
        rval[n] = qset = set()
        for r in db.sql(q):
            qset.add(r['_marker_key'])
    return rval

# ----------------------------------------------------------
# Cross References
# ----------------------------------------------------------

LDB2PREFIX = {
    "Ensembl Gene Model" : "ENSEMBL:",
    "Entrez Gene" : "NCBI_Gene:",
    "TrEMBL" : "UniProtKB:",
    "SWISS-PROT" : "UniProtKB:",
}

# build map from marker key to the list of xrefs for that marker
def getXrefs () :
    mk2xrefs = {}
    for r in db.sql(qXrefs, 'auto'):
        mk2xrefs.setdefault(r['_marker_key'], []).append(r)
    return mk2xrefs

def getFormattedXrefs (mkey, xrefs, gsets) :
    xrs = xrefs.get(mkey, None)
    xrs2 = []
    if not xrs:
        return xrs2
    for xr in xrs:
        if xr["dbname"] == "MGI" :
            pgs = ["gene","gene/references"]
            if mkey in gsets["hasExpression"]:
                pgs.append("gene/expression")
            if mkey in gsets["hasExpressionImage"]:
                pgs.append("gene/expression_images")
            if mkey in gsets["hasPhenotype"]:
                pgs.append('gene/phenotypes')
            if mkey in gsets["hasImpc"]:
                pgs.append('gene/phenotypes_impc')
            for pg in pgs:
                xrDto = {
                    "display_name" : xr["accid"],
                    "referenced_curie" : xr["accid"],
                    "page_area" : pg,
                    "prefix" : "MGI",
                    "internal" : False,
                } 
                xrs2.append(xrDto)
        else :
            prefix = LDB2PREFIX[xr["dbname"]]
            xrDto = {
                "display_name" : prefix + xr["accid"],
                "referenced_curie" : prefix + xr["accid"],
                "page_area" : "default",
                "prefix" : prefix[:-1],
                "internal" : False,
            }
            xrs2.append(xrDto)
    return xrs2

# build map from marker key to the list of notes for that marker
def getGeneNotes () :
    gene_notes = {}
    for r in db.sql(qGeneNotes, 'auto'):
        gene_notes.setdefault(r['_marker_key'], []).append(r)
    return gene_notes

# get gene notes in dto format
def getGeneNoteDtos (mkey, gnotes) :
    note_dtos = []
    if not mkey in gnotes:  # no notes, just return the empty array
        return note_dtos
    
    notes = gnotes[mkey]
    if len(notes) < 1:   # no notes, just return the empty array
        return note_dtos
    
    for note in notes:
        note_dto = {
            "note_type_name" : "MOD_provided_gene_description",
            "free_text" : note["note"]
        }
        setCommonFields(note, note_dto)
        note_dtos.append(note_dto)
    return note_dtos

# build map from marker key to list of synonyms
def getGeneSynonyms () :
    gene_syns = {}
    for r in db.sql(qGeneSynonyms, 'auto'):
        gene_syns.setdefault(r['_marker_key'], []).append(r)
    return gene_syns

# get synonyms for a gene in dto format
def getGeneSynonymDtos(mkey, gsyns):
    syn_dtos = []
    syns = gsyns.get(mkey, [])
    for s in syns:
        dto = {   
          "name_type_name": "unspecified",
          "format_text": s["synonym"],
          "display_text": s["synonym"],
          "synonym_scope_name": "exact",
          "internal": False
        }   
        syn_dtos.append(dto)
    return syn_dtos

# ----------------------------------------------------------
# ----------------------------------------------------------

def getJsonObject (r, xrefs, gsets, gnotes, gsynonyms) :
    obj = {
        "primary_external_id" : r["accid"],
        "gene_type_curie" : MCV2SO[r['_mcvterm_key']],
        "data_provider_dto": getDataProviderDto(r["accid"], "gene"),
        "taxon_curie": "NCBITaxon:10090",
        "internal": False,
        "gene_symbol_dto" : {
            "name_type_name" : "nomenclature_symbol",
            "format_text" : symbolToHtml(r["symbol"]),
            "display_text" : symbolToHtml(r["symbol"]),
            "internal" : False,
	    },
        "gene_full_name_dto" : {
            "name_type_name" : "full_name",
            "format_text" : symbolToHtml(r["name"]),
            "display_text" : symbolToHtml(r["name"]),
            "internal" : False
        }
    }
    setCommonFields(r, obj)
    obj["cross_reference_dtos"] = getFormattedXrefs(r["_marker_key"], xrefs, gsets)

    # add gene notes
    note_dtos = getGeneNoteDtos(r["_marker_key"], gnotes)
    if len(note_dtos) > 0:  # add gene dtos if exist
         obj["note_dtos"] = note_dtos

    # add gene synonyms
    syn_dtos = getGeneSynonymDtos(r["_marker_key"], gsynonyms)
    if len(syn_dtos) > 0:
        obj["gene_synonym_dtos"] = syn_dtos

    return obj

def main () :
    initMCV2SO()
    xrefs = getXrefs()
    gsets = getGeneSets()
    gnotes = getGeneNotes()
    gsynonyms = getGeneSynonyms()

    print('{')
    print(getHeaderAttributes())
    print('"gene_ingest_set": [')
    for j,r in mainQuery(db.sql(qGenes, 'auto')):
        if j: print(',', end='')
        o = getJsonObject(r, xrefs, gsets, gnotes, gsynonyms)
        print(json.dumps(o, indent=2))
    print(']')
    print('}')

# Returns the set of MGI ids for submitted genes
def getSubmittedGeneIds () :
    ids = set()
    for r in db.sql(qGenes):
        ids.add(r['accid'])
    return ids
# ----------------------------------------------------------
# Queries
# ----------------------------------------------------------

# Do not generate gene records for the following MCV types.
excludeMcvKeys = ",".join([str(x) for x in [
    # this group is being excluded because it was an MGI practice to always have a marker
    # for an allele to refer to. We will not carry these markers over to the Alliance.
    6238174, # Transgenes
    7196768, # chromosomal deletion
    7196774, # chromosomal duplication
    7196770, # chromosomal inversion
    7196773, # chromosomal translocation
    7196775, # chromosomal transposition
    7196769, # insertion
    7196772, # reciprocal chromosomal translocation
    7222413, # unclassified cytogenetic marker
    7196771, # Robertsonian fusion
    # This group is being excluded for now because we and the Alliance have not 
    # decided how to deal with them yet.
    6238173, # QTL
    #97015609, # CTCF binding site
    #36700088, # TSS cluster
    #97015607, # enhancer
    #15406207, # promoter
    #15406205, # CpG island
]])

# Basic info for each gene
qGenes = '''
    SELECT
        aa.accid, mm.*, mc._mcvterm_key
    FROM
        MRK_Marker mm,
        ACC_Accession aa,
        MRK_MCV_Cache mc
    WHERE
        mm._organism_key = 1
        and mm._marker_status_key = 1
        and mm._marker_key = mc._marker_key
        and mc.qualifier = 'D'
        and mm._marker_key = aa._object_key
        and aa._mgitype_key = 2
        and aa._logicaldb_key = 1
        and aa.preferred = 1
        and aa.private = 0
        and mc._mcvterm_key not in (%s)
    ''' % excludeMcvKeys

# Gene synonyms
qGeneSynonyms = '''
    SELECT DISTINCT ml._marker_key, ml.label as synonym
    FROM mrk_label ml, mrk_marker m
    WHERE ml._marker_key = m._marker_key
    AND ml.labeltype in ('MS','MN','MY')
    AND ml.labeltypename != 'human synonym'
    AND ml.label != m.symbol
    AND ml.label != m.name
    '''

# Cross references
qXrefs = '''
    SELECT m._marker_key, a.accid, a._logicaldb_key, d.name as dbname, a.creation_date, a.modification_date
    FROM ACC_Accession a, MRK_Marker m, ACC_LogicalDB d
    WHERE a._object_key = m._marker_key
    AND a._mgitype_key = 2
    AND a._logicaldb_key in (1, 13, 41, 55, 60) /* MGI, SwissProt, Trembl, EntrezGene, Ensembl gene model*/
    and a._logicaldb_key = d._logicaldb_key
    and a.preferred = 1
    '''
# genes with phenotype annots
qGeneHasPhenotype = ''' 
    SELECT distinct _object_key as _marker_key
    FROM VOC_Annot
    WHERE _annottype_key = 1015 /* MP-Gene */
    '''

# genes for alleles in the IMPC collection
qGeneHasImpc = ''' 
    SELECT distinct _marker_key
    FROM ALL_Allele
    WHERE _collection_key = 24755824 /* IMPC */
    '''

# genes that have expression data
qGeneHasExpression = '''
    SELECT distinct _marker_key
    FROM GXD_Expression
    '''

# genes that have expression data
qGeneHasExpressionImage = '''
    SELECT distinct _marker_key
    FROM GXD_Expression
    WHERE hasimage = 1
    '''

# genes with notes
qGeneNotes = ''' 
    SELECT m._marker_key, n.note, n.creation_date, n.modification_date
    FROM MRK_marker m, MRK_notes n
    WHERE m._marker_key = n._marker_key
    AND m._organism_key = 1 /* mouse, laboratory */
    '''

#########################################################

if __name__ == "__main__":
    main()
