
import db

# 
NS = "not_specified"
TNS = "taxon not specified"

# MGI organism "names" that map to not_specified at the Alliance
NS2 = "Not Specified"
NS3 = "Other (see notes)"

# indexes
species2taxonTerm = None
okey2taxonTerm = None

# This function reads/parses the taxon term file.
# Returns a mapping from latinName to NCBITaxon ID.
# For terms without NCBITaxon ID, has a mapping from
# common name to common name
#
# Alliance taxon term list has this format:
# Column 1 = the term, either an NCBITaxon ID, or a species common name
# Column 2 = the "definition" field for the alliance term.
#   If col 1 is an NCBITaxon ID, Has the form "common name (Latin name)"
#   Otherwise, has the form "common name (taxon not specified)"
# There is one special line for the "not_specified" taxon. Has the "not_specified"
# in both columns.
# Example lines:
#   NCBITaxon:10029    hamster, Chinese (Cricetulus griseus)
#   monkey             monkey (taxon not specified)
#   not_specified      not_specified
def loadTaxonTerms () :
    name2term = {}
    for line in ALLIANCE_ANTIBODY_TAXON_LIST.strip().split('\n'):
        pieces = line[:-1].split('\t')
        term = pieces[0]
        if term == NS:
            commonName = term
            latinName = term
        else:
            pieces2 = pieces[1].split(' (')
            commonName = pieces2[0]
            latinName = pieces2[1].replace(')','')
            if latinName == TNS: latinName = commonName
        name2term[latinName] = term
    return name2term

def getTaxonTermForOrganismName(org):
    global species2taxonTerm
    if species2taxonTerm is None:
        species2taxonTerm = loadTaxonTerms()
    term = species2taxonTerm.get(org['latinname'], '')
    if not term:
        term = species2taxonTerm.get(org['commonname'], None)
    if not term:
        if org['commonname'] in [NS2, NS3]:
            return NS
    return term

def loadMgiOrganisms():
    global okey2taxonTerm
    okey2taxonTerm = {}
    #preload special case
    okey2taxonTerm[1] = "NCBITaxon:10090"
    #
    for org in db.sql("SELECT * from mgi_organism"):
        term = getTaxonTermForOrganismName(org)
        if term:
            okey2taxonTerm[org['_organism_key']] = term
    return okey2taxonTerm

def getTaxonByOrganismKey(okey) :
    global okey2taxonTerm
    if okey2taxonTerm is None:
        okey2taxonTerm = loadMgiOrganisms()
    return okey2taxonTerm.get(okey, None)

ALLIANCE_ANTIBODY_TAXON_LIST = '''
bacteria	bacteria (taxon not specified)
bacteriophage	bacteriophage (taxon not specified)
butterfly	butterfly (taxon not specified)
camel	camel (taxon not specified)
carp	carp (taxon not specified)
chicken	chicken (taxon not specified)
donkey	donkey (taxon not specified)
frog	frog (taxon not specified)
goat	goat (taxon not specified)
guinea_pig	guinea pig (taxon not specified)
hamster	hamster (taxon not specified)
horse	horse (taxon not specified)
monkey	monkey (taxon not specified)
mouse	mouse (taxon not specified)
NCBITaxon:10029	hamster, Chinese (Cricetulus griseus)
NCBITaxon:10036	hamster, golden (Mesocricetus auratus)
NCBITaxon:10090	mouse (Mus musculus)
NCBITaxon:10116	rat (Rattus norvegicus)
NCBITaxon:10141	guinea pig (Cavia porcellus)
NCBITaxon:105401	Zoanthus (Zoanthus)
NCBITaxon:1870830	electric ray, ocellated (Diplobatis ommata)
NCBITaxon:3055	green algae (Chlamydomonas reinhardtii)
NCBITaxon:31033	fugu (Takifugu rubripes)
NCBITaxon:3122392	hamster, Armenian (Cricetulus migratorius)
NCBITaxon:34765	larvacean (Oikopleura dioica)
NCBITaxon:36177	Atlantic sturgeon (Acipenser oxyrinchus)
NCBITaxon:4564	wheat (Triticum)
NCBITaxon:4932	yeast (Saccharomyces cerevisiae)
NCBITaxon:6100	jellyfish (Aequorea victoria)
NCBITaxon:6239	nematode (Caenorhabditis elegans)
NCBITaxon:69293	three-spined stickleback (Gasterosteus aculeatus)
NCBITaxon:7227	fruit fly (Drosophila melanogaster)
NCBITaxon:7668	purple sea urchin (Strongylocentrotus purpuratus)
NCBITaxon:7719	sea squirt (Ciona intestinalis)
NCBITaxon:7740	amphioxus (Branchiostoma lanceolatum)
NCBITaxon:7757	sea lamprey (Petromyzon marinus)
NCBITaxon:7769	Atlantic hagfish (Myxine glutinosa)
NCBITaxon:7787	Pacific electric ray (Torpedo californica)
NCBITaxon:7797	spiny dogfish (Squalus acanthias)
NCBITaxon:7924	bowfin (Amia calva)
NCBITaxon:7955	zebrafish (Danio rerio)
NCBITaxon:7957	goldfish (Carassius auratus)
NCBITaxon:7962	carp (Cyprinus carpio)
NCBITaxon:8022	rainbow trout (Oncorhynchus mykiss)
NCBITaxon:8030	salmon (Salmo salar)
NCBITaxon:8078	mummichog killifish (Fundulus heteroclitus)
NCBITaxon:8090	medaka (Oryzias latipes)
NCBITaxon:8353	Xenopus (Xenopus)
NCBITaxon:8355	African clawed frog (Xenopus laevis)
NCBITaxon:8364	Western clawed frog (Xenopus tropicalis)
NCBITaxon:8400	bullfrog (Rana catesbeiana)
NCBITaxon:86599	coral anemone (Discosoma)
NCBITaxon:9031	chicken (Gallus gallus)
NCBITaxon:9091	common quail (Coturnix coturnix)
NCBITaxon:9135	common canary (Serinus canaria)
NCBITaxon:9606	human (Homo sapiens)
NCBITaxon:9615	dog, domestic (Canis lupus familiaris)
NCBITaxon:9685	cat, domestic (Felis catus)
NCBITaxon:9796	horse (Equus caballus)
NCBITaxon:9823	pig (Sus scrofa)
NCBITaxon:9825	pig, domestic (Sus scrofa domestica L.)
NCBITaxon:9838	Arabian camel (Camelus dromedarius)
NCBITaxon:9843	vicugna (Vicugna vicugna)
NCBITaxon:9844	llama (Lama glama)
NCBITaxon:9913	cow (Bos taurus)
NCBITaxon:9925	goat (Capra hircus)
NCBITaxon:9940	sheep (Ovis aries)
NCBITaxon:9986	rabbit, European (Oryctolagus cuniculus)
NCBITaxon:99883	spotted green pufferfish (Tetraodon nigroviridis)
newt	newt (taxon not specified)
not_specified	not_specified
pigeon	pigeon (taxon not specified)
quail	quail (taxon not specified)
rabbit	rabbit (taxon not specified)
rat	rat (taxon not specified)
salamander	salamander (taxon not specified)
sheep	sheep (taxon not specified)
virus	virus (taxon not specified)
'''

if __name__ == "__main__":
    print(1, getTaxonByOrganismKey(1))
    print(40, getTaxonByOrganismKey(40))
    print(84, getTaxonByOrganismKey(84))
    print(76, getTaxonByOrganismKey(76))
    print(77, getTaxonByOrganismKey(77))
    print(60, getTaxonByOrganismKey(60))
    print(126, getTaxonByOrganismKey(126))

