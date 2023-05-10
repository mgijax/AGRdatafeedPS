import jsonschema
import json
import sys
import argparse

def getArgs () :
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="schemaFile", help="The JSON schema file.")
    parser.add_argument("-d", dest="dataFile", help="The JSON data file to test.")
    return parser.parse_args()

def main ():
    args = getArgs()
    with open(args.dataFile, 'r') as fd:
        dataObj = json.load(fd)
    with open(args.schemaFile, 'r') as sd:
        schemaObj = json.load(sd)

    try:
        jsonschema.validate(instance=dataObj, schema=schemaObj)
        sys.exit(0)
    except jsonschema.exceptions.ValidationError as err:
        print(err.message)
        print(err)
        sys.exit(1)

main()
