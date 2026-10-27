
import db

# mapping from MGI organism key to Alliance term
okey2allianceTaxonTerm = None

def parseTaxonMapping () :
    global okey2allianceTaxonTerm
    okey2allianceTaxonTerm = {}
    lines = TAXON_MAPPING.strip().split('\n')
    for l in lines:
        okey, commonName, latinName, allianceTerm = l.split('\t')
        okey2allianceTaxonTerm[int(okey)] = allianceTerm

def getTaxonByOrganismKey(okey) :
    if okey2allianceTaxonTerm is None:
        parseTaxonMapping()
    return okey2allianceTaxonTerm[okey]

'''
This SQL query returns the distinct mgi_organisms used for antibodies and antigens

select a.accid, o._organism_key, o.commonname, o.latinname
  from mgi_organism o
  left join ACC_Accession a on o._organism_key = a._object_key
  and a._mgitype_key = 20
  and a._logicaldb_key = 32
where o._organism_key in
  (select distinct _organism_key from gxd_antibody
  UNION
  select distinct s._organism_key
  from gxd_antibody ga join prb_source s on ga._source_key = s._source_key
  )
order by commonname
'''
# Format: 
# _organism_key commonName  latinName   allianceTerm
#
TAXON_MAPPING = '''
60	bacteria	Not Specified	bacteria
90	bacteriophage	Not Specified	bacteriophage
61	butterfly	Not Specified	butterfly
62	carp	Not Specified	carp
9	cat, domestic	Felis catus	NCBITaxon:9685
11	cattle	Bos taurus	NCBITaxon:9913
88	C. elegans	Caenorhabditis elegans	NCBITaxon:6239
63	chicken	Gallus gallus	NCBITaxon:9031
93	coral anemone	Discosoma	NCBITaxon:86599
13	dog, domestic	Canis lupus familiaris	NCBITaxon:9615
64	Drosophila melanogaster	Drosophila melanogaster	NCBITaxon:7227
65	electric ray, ocellated	Diplobatis ommata	NCBITaxon:1870830
96	electric ray, Pacific	Torpedo californica	NCBITaxon:7787
92	frog	Xenopus	NCBITaxon:8353
86	frog, African clawed	Xenopus laevis	NCBITaxon:8355
18	goat	Capra hircus	NCBITaxon:9925
58	guinea pig	Cavia porcellus	NCBITaxon:10141
69	hamster	Cricetidae	hamster
89	hamster, Armenian	Nothocricetulus migratorius	NCBITaxon:3122392
20	hamster, Chinese	Cricetulus griseus	NCBITaxon:10029
47	hamster, golden	Mesocricetus auratus	NCBITaxon:10036
21	horse	Equus caballus	NCBITaxon:9796
2	human	Homo sapiens	NCBITaxon:9606
87	jellyfish	Aequorea victoria	NCBITaxon:6100
72	monkey	Not Specified	monkey
1	mouse, laboratory	Mus musculus/domesticus	NCBITaxon:10090
76	Not Specified	Not Specified	not_specified
77	Other (see notes)	Not Specified	not_specified
35	pig, domestic	Sus scrofa domestica L.	NCBITaxon:9823
140	pigeon	Not Specified	pigeon
91	quail	Not Specified	quail
79	rabbit	Not Specified	rabbit
40	rat	Rattus norvegicus	rat
44	sheep	Ovis aries	sheep
82	virus	Not Specified	virus
137	wheat	Triticum	NCBITaxon:4564
84	zebrafish	Danio rerio	NCBITaxon:7955
136	Zoanthus	Zoanthus	NCBITaxon:105401
74	Not Applicable	Not Applicable	not_specified
'''

if __name__ == "__main__":
    print(1, getTaxonByOrganismKey(1))
    print(40, getTaxonByOrganismKey(40))
    print(84, getTaxonByOrganismKey(84))
    print(76, getTaxonByOrganismKey(76))
    print(77, getTaxonByOrganismKey(77))
    print(60, getTaxonByOrganismKey(60))

