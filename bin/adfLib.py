
import re
import os
import sys
import time
import datetime
import db

#----------------------------------
# See: http://henry.precheur.org/projects/rfc3339 
from rfc3339 import rfc3339
date_re = re.compile(r'(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)')
#----------------------------------
# RFC 3339 timestamps
#
# With no arguments, returns the current date-time in RFC-3339 format
# With a string argument of the form 'yyyy-mm-dd', converts to rfc3339 format and returns.
# Note that simply concatenating a string such as 'T00:00:00-05:00' to the date is not
# sufficient because of DST differences (i.e., for part of the year, the offset is -04:00
# rather that "-05:00"
# Examples:
#   getTimeStamp() --> "2017-01-26T15:00:42-05:00"
#   getTimeStamp("2007-01-15") --> "2007-01-15T00:00:00-05:00"
#   getTimeStamp("2014-05-01") --> "2014-05-01T00:00:00-04:00"
#   getTimeStamp("2023-11-15 09:21:20") --> 
#
def getTimeStamp(s = None):
    if s:
        m = date_re.match(s)
        d = datetime.datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)),  int(m.group(4)), int(m.group(5)), int(m.group(6)))
        return rfc3339(d)
    else:
        return rfc3339(time.time())
#-----------------------------------

#
def log (msg, addTimestamp=True, lineTerminator="\n") :
    if addTimestamp:
        timestamp = time.asctime(time.localtime(time.time()))
        sys.stderr.write(timestamp + ' ')
    sys.stderr.write(str(msg))
    if lineTerminator:
        sys.stderr.write(lineTerminator)

# Returns a string to be used as the "alliance_member_release_version" field for data updates
# (Added in curation_schema v1.11.0)
# Concatendates public_version and lastdump_date from mgi_dbinfo.
#
def getReleaseVersion () :
    dbi = db.sql('select * from mgi_dbinfo')[0]
    return '%s %s' % (dbi['public_version'],dbi['lastdump_date'])
#
def getHeaderAttributes () :
    linkml_version = f'''"linkml_version": "{os.environ.get('AGR_CURATION_SCHEMA_VERSION','default')}",'''
    release_version = f'''"alliance_member_release_version" : "{getReleaseVersion()}",'''
    return '%s\n%s\n' % (linkml_version, release_version)



# Wraps the execution of the main query so that we can implement a "sample" option, that outputs
# a small sample of records. Useful for development. At the moment,a sample is just the first 20 records.
# Args:
#   results (iterable) Somethat can be iterated over to get the result records
# Yields:
#   Tuples (n,r) where n is a 0-based count and r is a record
#
DO_SAMPLE = os.environ.get('DO_SAMPLE','')
SAMPLE_COUNT=200
def mainQuery(results) :
    n = 0
    for r in results:
        if DO_SAMPLE and n >= SAMPLE_COUNT:
            break;
        yield n, r
        n += 1

# Converts "Foo<Bar>" to "Foo<sup>Bar</sup>"
def symbolToHtml (s) :
    return re.sub(r'<([^>]*)>', r'<sup>\1</sup>', s)

# Builds a key-value mapping from a list of objects.
# Args:
#  results - query result objects
#  keyfield - name of the field holding the key
#  valuefield - name of the field holding the value. If None, the whole record is the value
#  multi - if true, returns multivalued index; if false, single valued
#  mapper - optional. function to map result values to index entries.
#       Default maps a term to itself.
# Returns:
#  If multi is false, a mapping (dict) from keys to single values
#  If multi is true, a mapping from keys to lists of values
def indexResults (results, keyfield, valuefield, multi=False, mapper=None) :
    index = {}
    mapper = mapper if mapper else lambda x:x
    for r in results:
        v = r if valuefield is None else r[valuefield]
        v = mapper(v)
        if multi:
            index.setdefault(r[keyfield],[]).append(v)
        else:
            index[r[keyfield]] = v
    return index

# Returns a "data_provider_dto" object for the given curie and page area.
# These were introduced in schema version 1.6.0 for diseaseAnnotations, and 
# for genes, alleles, and agms in 1.7.0.
# Examples:
#       getDataProviderDto("DOID:123456", "disease/mgi")
#       getDataProviderDto("MGI:567890", "gene")
def getDataProviderDto (curie, pageArea):
      i = curie.find(":")
      prefix = curie[:i]
      dto = {
      	"source_organization_abbreviation": "MGI",  
	"internal" : False,
        "cross_reference_dto": {
          "referenced_curie": curie,
          "prefix": prefix,
          "page_area": pageArea,
          "display_name": curie,
          "internal": False
        }  
      }
      return dto 

# Sets common fields in the object, based on values in the db record
def setCommonFields (rec, obj, internal=False, obsolete=False) :
    obj["internal"] = internal
    obj["obsolete"] = obsolete
    if rec.has_key("creation_date"): obj["date_created"] = getTimeStamp(rec["creation_date"])
    if rec.has_key("modification_date"): obj["date_updated"] = getTimeStamp(rec["modification_date"])
    obj["created_by_curie"] = "MGI:curation_staff"
    obj["updated_by_curie"] = "MGI:curation_staff"

# Return a mapping from reference key to the MGI id and PMID (if available).
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

# Return the preferred ID for a reference, given it ref key.
# The preferred ID is the PMID, if there is one, otherswise the MGI id.
rk2ids = None
def getPreferredRefId (rk):
    global rk2ids
    if rk2ids is None:
        rk2ids = getReferenceIds()
    if rk is None: return None
    ids = rk2ids[rk] # every reference must have an entry, else error
    if ids["pubmedid"]:
        return "PMID:" + ids["pubmedid"]
    else:
        return ids["mgiid"]


# Get the MGI notes of the given note type.
# Return an index of _object_key to the note record(s) for that object
def getNotesOfType (noteTypeKey) :
    q = '''SELECT *
        FROM MGI_Note
        WHERE _noteType_key = %s
        ''' % noteTypeKey
    k2n = {}
    for r in db.sql(q) :
        k2n.setdefault(r['_object_key'], []).append(r)
    return k2n

# Turns a note record (from the MGI_Note table) into a NoteDTO object for submission.
# 
def getNoteDTO ( r , noteTypeName="comment") :
    rr = {
    "note_type_name": noteTypeName,
    "free_text": r["note"],
    }
    if r.has_key("evidence_curies"):
        rr["evidence_curies"] = r["evidence_curies"]
    setCommonFields(r, rr)
    return rr
