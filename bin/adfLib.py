
import re
import os
import sys
import time

#
def log (msg, addTimestamp=True, lineTerminator="\n") :
    if addTimestamp:
        timestamp = time.asctime(time.localtime(time.time()))
        sys.stderr.write(timestamp + ' ')
    sys.stderr.write(str(msg))
    if lineTerminator:
        sys.stderr.write(lineTerminator)

#
LINKML_VERSION = f'''"linkml_version": "{os.environ.get('AGR_CURATION_SCHEMA_VERSION','default')}",''' 

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

