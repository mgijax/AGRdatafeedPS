
import db
import json
import re

from adfLib import getHeaderAttributes, symbolToHtml, getDataProviderDto, mainQuery, setCommonFields

# ----------------------------------------------------------
# Mapping from MCV term key to SO id.
# Initialize with hard-coded mappings then load what's in the db (which is incomplete).
# ----------------------------------------------------------

MCV2SO = { 
    #"complex/cluster/region"
    6238175 : "SO:0000110",
    #"cytogenetic marker"
    6238176 : "SO:0000110",
    #"BAC/YAC end"
    6238177 : "SO:0000110",
    #"other genome feature"
    6238178 : "SO:0000110",
    #"DNA segment"
    6238179 : "SO:0000110",
    #"unclassified gene"
    6238184 : "SO:0000704",
    #"other feature type"
    6238185 : "SO:0000110",
    #"unclassified non-coding RNA gene"
    6238186 : "SO:0000704",
    #"unclassified cytogenetic marker"
    7222413 : "SO:0000110",
    #"unclassified other genome feature"
    7648969 : "SO:0000110",
    #"mutation defined region"
    11928467 : "SO:0000110",
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
        if m:
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

# ----------------------------------------------------------
# ----------------------------------------------------------


def getJsonObject (r, xrefs, gsets) :
    obj = {
        "mod_entity_id" : r["accid"],
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
    return obj

def main () :
    initMCV2SO()
    xrefs = getXrefs()
    gsets = getGeneSets()

    print('{')
    print(getHeaderAttributes())
    print('"gene_ingest_set": [')
    for j,r in mainQuery(db.sql(qGenes, 'auto')):
        if j: print(',', end='')
        o = getJsonObject(r, xrefs, gsets)
        print(json.dumps(o, indent=2))
    print(']')
    print('}')

# ----------------------------------------------------------
# Queries
# ----------------------------------------------------------

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
        and mm._marker_type_key = 1
        and mm._marker_status_key = 1
        and mm._marker_key = mc._marker_key
        and mc.qualifier = 'D'
        and mm._marker_key = aa._object_key
        and aa._mgitype_key = 2
        and aa._logicaldb_key = 1
        and aa.preferred = 1
        and aa.private = 0
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
#########################################################

if __name__ == "__main__":
    main()
