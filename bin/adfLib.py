
import re

def symbolToHtml (s) :
    return re.sub(r'<([^>]*)>', r'<sup>\1</sup>', s)
