
import re
import os

#
LINKML_VERSION = f'''"linkml_version": "{os.environ.get('AGR_CURATION_SCHEMA_VERSION','default')}",''' 

# Converts "Foo<Bar>" to "Foo<sup>Bar</sup>"
def symbolToHtml (s) :
    return re.sub(r'<([^>]*)>', r'<sup>\1</sup>', s)

# Builds a key-value mapping from a list of objects.
# Args:
#  results - query result objects
#  keyfield - name of the field holding the key
#  valuefield - name of the field holding the value
#  multi - if true, returns multivalued index; if false, single valued
# Returns:
#  If multi is false, a mapping (dict) from keys to single values
#  If multi is true, a mapping from keys to lists of values
def indexResults (results, keyfield, valuefield, multi=False) :
    index = {}
    for r in results:
        if multi:
            index.setdefault(r[keyfield],[]).append(r[valuefield])
        else:
            index[r[keyfield]] = r[valuefield]
    return index
