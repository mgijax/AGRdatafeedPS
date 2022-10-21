
import re

def symbolToHtml (s) :
    return re.sub(r'<([^>]*)>', r'<sup>\1</sup>', s)

def indexResults (results, keyfield, valuefield, multi=False) :
    index = {}
    for r in results:
        if multi:
            index.setdefault(r[keyfield],[]).append(r[valuefield])
        else:
            index[r[keyfield]] = r[valuefield]
    return index
