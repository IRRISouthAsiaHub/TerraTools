"""
From a trait name to the words that count as trait-relevant text.

  direct   words of the trait itself (and its plain synonyms): "drought", "water deficit" for Drought tolerance
  related  vocabulary known to belong to the trait ("abscisic acid", "osmotic" …): weaker evidence than a direct hit

The curated list below is a starting point and is meant to be extended: add a dict to CONCEPTS. A trait with no curated
entry still works from the words of its own name, and the run reports that no synonyms were used.
"""

from __future__ import annotations

import re

GENERIC = {
    "tolerance", "tolerant", "resistance", "resistant", "response", "responses", "trait", "traits", "stress", "of", "to",
    "and", "the", "in", "for", "use", "related", "level", "levels", "a", "an", "under", "against", "type", "index",
    "efficiency", "content",
    # words that name the crop, an organ or a stage rather than the trait: they would match nearly every record
    "rice", "wheat", "maize", "corn", "plant", "plants", "grain", "grains", "seed", "seeds", "stage", "stages", "date", "time",
    "early", "late", "low", "high", "rate", "days", "number", "per", "area", "total", "quality", "male", "female", "brown",
    "white", "red", "black", "yellow", "depth", "cooking", "bacterial", "leaf", "leaves", "gene", "genes",
}

CONCEPTS = [
    # ------------------------------------------------------------------ abiotic stress
    dict(key="drought", aliases=["drought", "water deficit", "water stress", "dehydration", "water use", "water-limited"],
         direct=["drought", "water deficit", "water stress", "dehydration", "desiccation"],
         related=["abscisic acid", "ABA", "osmotic", "osmoprotect", "stomata", "stomatal", "water use", "aquaporin",
                  "DREB", "late embryogenesis", "proline", "relative water content", "leaf rolling", "wax", "cuticle"]),
    dict(key="salt", aliases=["salt", "salinity", "saline", "sodium"],
         direct=["salt", "salinity", "saline", "NaCl", "sodium"],
         related=["ion homeostasis", "potassium", "HKT", "SOS1", "Na+", "osmotic", "vacuolar", "ion transport", "NHX"]),
    dict(key="heat", aliases=["heat", "high temperature", "thermotolerance", "thermal"],
         direct=["heat", "high temperature", "thermotolerance", "thermal"],
         related=["heat shock", "HSF", "HSP", "chaperone", "temperature", "unfolded protein", "spikelet fertility"]),
    dict(key="cold", aliases=["cold", "chilling", "freezing", "low temperature", "frost", "winter hardiness", "winter survival"],
         direct=["cold", "chilling", "freezing", "low temperature", "frost"],
         related=["cold acclimation", "CBF", "COLD1", "temperature", "cryoprotect", "antifreeze", "vernalization"]),
    dict(key="submergence", aliases=["submergence", "flooding", "waterlogging", "anaerobic", "anoxia", "flood", "anaerobic germination"],
         direct=["submergence", "flooding", "waterlogging", "anaerobic", "anoxia", "hypoxia"],
         related=["ethylene", "SUB1", "fermentation", "aerenchyma", "alcohol dehydrogenase", "coleoptile elongation", "trehalose-6-phosphate"]),
    dict(key="oxidative", aliases=["oxidative", "reactive oxygen", "antioxidant", "ozone"],
         direct=["oxidative stress", "reactive oxygen", "antioxidant", "ozone", "ultraviolet"],
         related=["superoxide dismutase", "catalase", "peroxidase", "ascorbate", "glutathione", "H2O2", "hydrogen peroxide"]),
    dict(key="metal", aliases=["cadmium", "arsenic", "heavy metal", "aluminium", "aluminum", "toxicity", "iron toxicity", "zinc toxicity"],
         direct=["cadmium", "arsenic", "heavy metal", "aluminium", "aluminum", "metal tolerance", "metal toxicity"],
         related=["metal transporter", "phytochelatin", "metallothionein", "vacuolar sequestration", "ALMT", "MATE", "citrate efflux", "Nramp"]),
    dict(key="alkalinity", aliases=["alkalinity", "alkaline", "sodicity", "acid soil", "soil acidity"],
         direct=["alkaline", "alkalinity", "sodic", "soil acidity", "acid soil"],
         related=["proton", "pH", "H+-ATPase", "organic acid", "aluminium"]),
    dict(key="nitrogen", aliases=["nitrogen", "nue", "nitrate", "ammonium", "n use", "n-use"],
         direct=["nitrogen", "nitrate", "ammonium", "nitrate transporter", "NUE"],
         related=["glutamine synthetase", "glutamate", "amino acid transport", "nitrate reductase", "NRT", "assimilation", "protein content", "NLP"]),
    dict(key="phosphorus", aliases=["phosphorus", "phosphate", "phosphorous", "pup1", "p use", "p-use"],
         direct=["phosphate", "phosphorus", "phosphorous"],
         related=["PHT", "PHO", "phosphate starvation", "SPX", "purple acid phosphatase", "root", "PSTOL1"]),
    dict(key="potassium", aliases=["potassium", "k use", "k deficiency"],
         direct=["potassium", "K+ transporter"],
         related=["HAK", "AKT", "K+ uptake", "ion transport", "CBL", "CIPK"]),
    dict(key="micronutrient", aliases=["zinc", "iron", "biofortification", "micronutrient", "mineral", "iron deficiency", "zinc deficiency"],
         direct=["zinc", "iron", "micronutrient", "mineral", "biofortification"],
         related=["ZIP", "YSL", "nicotianamine", "ferritin", "phytosiderophore", "IRT", "heavy metal transporter", "phytate", "NAS"]),
    # ------------------------------------------------------------------ biotic stress
    # The first entries of `direct` also go into the PubMed and Plant Reactome searches (disease name first, then the pathogen or pest).
    # --- fungal and oomycete diseases of leaves, panicles and ears
    dict(key="blast", aliases=["blast", "magnaporthe", "pyricularia", "neck blast", "leaf blast", "panicle blast", "wheat blast", "pot2"],
         direct=["blast", "Magnaporthe", "Pyricularia", "neck blast", "leaf blast", "panicle blast"],
         related=["Pi9", "Pib", "Pita", "Pi-ta", "Pi54", "Pik", "Pi5", "avirulence", "AVR-Pi", "appressorium", "disease resistance", "defense response", "defence response",
                  "NBS-LRR", "NLR", "pathogen", "hypersensitive", "salicylic acid", "jasmonic acid", "PR protein", "resistance gene", "innate immun", "fungal"]),
    dict(key="brown_spot", aliases=["brown spot", "bipolaris", "cochliobolus", "helminthosporium", "drechslera", "narrow brown leaf spot", "narrow brown spot"],
         direct=["brown spot", "Bipolaris oryzae", "Cochliobolus miyabeanus", "Helminthosporium", "Bipolaris", "Cochliobolus", "narrow brown leaf spot", "leaf spot"],
         related=["necrotroph", "necrotrophic", "lesion", "cell death", "reactive oxygen", "silicon", "OsWRKY", "OsMPK", "MAP kinase", "ethylene", "abscisic acid",
                  "fungal", "pathogen", "disease resistance", "defense response", "defence response", "NBS-LRR", "NLR", "PR protein", "resistance gene", "chitinase"]),
    dict(key="leaf_spot_blotch", aliases=["leaf spot", "spot blotch", "tan spot", "septoria", "stagonospora", "blotch", "net blotch", "cercospora", "gray leaf spot",
                                          "grey leaf spot", "northern corn leaf blight", "southern corn leaf blight", "turcicum", "exserohilum", "sorokiniana",
                                          "pyrenophora", "zymoseptoria", "mycosphaerella", "ncb", "nclb", "sclb", "stb", "gls", "phaeosphaeria", "eyespot",
                                          "leaf scald", "rhynchosporium"],
         direct=["leaf spot", "leaf blight", "spot blotch", "tan spot", "septoria", "blotch", "gray leaf spot", "northern corn leaf blight", "southern corn leaf blight",
                 "Exserohilum turcicum", "Bipolaris sorokiniana", "Cercospora", "Zymoseptoria", "Pyrenophora", "Stagonospora", "Rhynchosporium"],
         related=["Ht1", "Htn1", "Stb6", "Stb16q", "Tsn1", "Lr34", "Lr67", "quantitative resistance", "necrotrophic effector", "toxin sensitivity", "lesion", "cell death",
                  "disease resistance", "defense response", "defence response", "NLR", "NBS-LRR", "kinase", "ABC transporter", "hexose transporter", "pathogen", "fungal"]),
    dict(key="rust_wheat", aliases=["rust", "stripe rust", "yellow rust", "leaf rust", "stem rust", "puccinia", "common rust", "southern rust", "polysora", "uromyces", "ug99",
                                    "brown rust", "black rust"],
         direct=["rust", "stripe rust", "yellow rust", "leaf rust", "stem rust", "Puccinia", "common rust", "southern rust", "Ug99"],
         related=["adult plant resistance", "Lr34", "Lr67", "Yr36", "Yr15", "Sr35", "Sr50", "Rp1", "Rpp", "slow rusting", "NLR", "NBS-LRR", "disease resistance",
                  "defense response", "ABC transporter", "kinase", "START domain", "hexose transporter"]),
    dict(key="mildew", aliases=["powdery mildew", "mildew", "mlo", "downy mildew", "blumeria", "peronosclerospora", "sclerophthora", "crazy top", "sorghum downy"],
         direct=["mildew", "powdery mildew", "downy mildew", "Blumeria", "Peronosclerospora", "Sclerophthora", "crazy top"],
         related=["MLO", "Pm", "Pm3", "Pm8", "Pm21", "Mla", "NLR", "disease resistance", "defense response", "hypersensitive", "callose", "papilla"]),
    dict(key="sheath_blight", aliases=["sheath blight", "rhizoctonia", "sheath rot", "banded leaf", "banded leaf and sheath blight", "sheath spot", "thanatephorus", "sarocladium",
                                       "sheath brown rot"],
         direct=["sheath blight", "Rhizoctonia solani", "Rhizoctonia", "banded leaf and sheath blight", "sheath rot", "Sarocladium oryzae", "Thanatephorus"],
         related=["quantitative resistance", "disease resistance", "defense response", "pathogen", "fungal", "chitinase", "plant height", "lodging", "necrotroph", "lesion"]),
    dict(key="false_smut", aliases=["false smut", "ustilaginoidea", "green smut", "rice false smut", "kernel smut of rice", "udbatta"],
         direct=["false smut", "Ustilaginoidea virens", "Ustilaginoidea", "green smut", "rice false smut"],
         related=["ustiloxin", "mycotoxin", "effector", "biotroph", "floral infection", "spikelet", "disease resistance", "defense response", "pathogen", "fungal", "NLR"]),
    dict(key="smut_bunt", aliases=["smut", "bunt", "ustilago", "tilletia", "karnal bunt", "loose smut", "head smut", "sporisorium", "kernel smut", "tilletia indica", "covered smut",
                                   "common smut", "ergot", "claviceps", "neovossia"],
         direct=["smut", "bunt", "Ustilago", "Tilletia", "Karnal bunt", "loose smut", "head smut", "Sporisorium", "common smut", "ergot", "Claviceps"],
         related=["Ustilago maydis", "effector", "biotroph", "gall", "teliospore", "disease resistance", "defense response", "pathogen", "fungal", "NLR"]),
    dict(key="fusarium", aliases=["fusarium", "head blight", "scab", "fhb", "ear rot", "gibberella", "stalk rot", "bakanae", "fujikuroi", "verticillioides", "graminearum", "culmorum",
                                  "fusarium wilt", "crown rot", "seedling blight", "pokkah boeng", "foot rot"],
         direct=["Fusarium", "head blight", "Fusarium head blight", "scab", "ear rot", "stalk rot", "Gibberella", "bakanae", "Fusarium fujikuroi", "Fusarium graminearum",
                 "Fusarium verticillioides", "crown rot", "Fusarium wilt"],
         related=["Fhb1", "Fhb7", "Qfhb1", "deoxynivalenol", "mycotoxin", "DON", "fumonisin", "zearalenone", "UDP-glucosyltransferase", "type II resistance", "kernel infection",
                  "gibberellin", "disease resistance", "defense response", "fungal"]),
    dict(key="mycotoxin_aspergillus", aliases=["aflatoxin", "aspergillus", "mycotoxin", "fumonisin", "ochratoxin", "grain mold", "grain mould", "kernel rot", "ear mold",
                                               "discoloured grain", "discolored grain", "penicillium", "storage fungi"],
         direct=["aflatoxin", "Aspergillus", "Aspergillus flavus", "mycotoxin", "fumonisin", "ochratoxin", "grain mold", "kernel rot", "discolored grain"],
         related=["kernel pericarp", "wax", "phytoalexin", "resistance to infection", "lipoxygenase", "pathogenesis-related", "chitinase", "fungal", "defense response",
                  "disease resistance", "toxin"]),
    dict(key="anthracnose", aliases=["anthracnose", "colletotrichum", "red rot", "glomerella", "stalk anthracnose"],
         direct=["anthracnose", "Colletotrichum", "Colletotrichum graminicola", "red rot", "Glomerella"],
         related=["Rcg1", "hemibiotroph", "lesion", "disease resistance", "defense response", "defence response", "NLR", "fungal", "pathogen", "stalk"]),
    dict(key="soilborne_rot", aliases=["root rot", "pythium", "damping off", "damping-off", "phytophthora", "take-all", "take all", "gaeumannomyces", "verticillium",
                                       "soilborne", "soil-borne", "seedling rot", "stem rot", "sclerotium", "sclerotinia", "magnaporthiopsis", "oomycete"],
         direct=["root rot", "Pythium", "damping off", "Phytophthora", "take-all", "Gaeumannomyces", "Verticillium", "stem rot", "Sclerotium", "Sclerotinia", "oomycete"],
         related=["root", "seedling", "rhizosphere", "root exudate", "soil", "hypersensitive", "disease resistance", "defense response", "pathogen", "fungal", "NLR"]),
    dict(key="wilt", aliases=["bacterial wilt", "ralstonia", "stewart", "goss", "vascular wilt", "verticillium wilt", "fusarium wilt", "wilt disease", "clavibacter", "pantoea"],
         direct=["wilt", "bacterial wilt", "Ralstonia", "Fusarium wilt", "Verticillium wilt", "vascular wilt", "Stewart's wilt", "Goss's wilt", "Clavibacter", "Pantoea stewartii"],
         related=["xylem", "vascular", "tylose", "disease resistance", "defense response", "pathogen", "receptor kinase", "NLR", "effector"]),
    # --- bacterial diseases
    dict(key="bacterial_blight", aliases=["bacterial blight", "xanthomonas", "xoo", "bacterial leaf streak", "leaf streak", "xoc", "blb", "bacterial leaf blight", "kresek",
                                          "bls", "pseudomonas", "acidovorax", "bacterial stripe", "bacterial brown stripe", "burkholderia", "bacterial panicle blight",
                                          "panicle blight", "bacterial sheath brown rot", "bacterial grain rot", "bacterial leaf"],
         direct=["bacterial blight", "Xanthomonas oryzae", "Xanthomonas", "bacterial leaf streak", "leaf streak", "bacterial panicle blight", "Burkholderia glumae", "Pseudomonas",
                 "Acidovorax", "bacterial stripe"],
         related=["Xa21", "Xa4", "Xa5", "Xa7", "Xa13", "Xa23", "Xa26", "Xa27", "xa25", "SWEET", "SWEET11", "SWEET13", "SWEET14", "TAL effector", "TALE", "EBE", "executor gene",
                  "receptor kinase", "PAMP", "pattern-triggered", "disease resistance", "defense response", "NBS-LRR", "NLR", "pathogen", "type III effector", "quorum sensing"]),
    dict(key="phytoplasma", aliases=["phytoplasma", "yellows", "witches' broom", "witches broom", "stunt disease", "grassy stunt phytoplasma", "bushy stunt", "spiroplasma",
                                     "white leaf", "sugarcane white leaf"],
         direct=["phytoplasma", "yellows", "witches' broom", "bushy stunt", "Spiroplasma", "white leaf"],
         related=["leafhopper", "planthopper", "vector", "phloem", "sieve element", "effector", "SAP11", "disease resistance", "defense response", "pathogen"]),
    # --- viruses
    dict(key="virus", aliases=["virus", "viral", "tungro", "rice stripe", "yellow mottle", "streak virus", "mosaic", "rice dwarf", "grassy stunt", "ragged stunt",
                               "black-streaked dwarf", "black streaked dwarf", "stripe virus", "maize streak", "maize dwarf mosaic", "mdmv", "barley yellow dwarf", "bydv",
                               "wheat streak mosaic", "wsmv", "sugarcane mosaic", "scmv", "potyvirus", "geminivirus", "begomovirus", "luteovirus", "tenuivirus", "reovirus",
                               "yellow dwarf", "chlorotic mottle", "lethal necrosis", "mcmv", "rsv", "rdv", "rbsdv", "rtbv", "rymv", "rhdv", "rice hoja blanca", "hoja blanca"],
         direct=["virus", "viral", "tungro", "rice stripe virus", "rice dwarf virus", "rice grassy stunt", "rice ragged stunt", "rice black-streaked dwarf", "Rice yellow mottle virus",
                 "maize streak virus", "maize dwarf mosaic", "barley yellow dwarf", "wheat streak mosaic", "sugarcane mosaic", "mosaic", "yellow dwarf", "lethal necrosis",
                 "potyvirus", "geminivirus", "stripe virus"],
         related=["RNA silencing", "RNA interference", "siRNA", "Dicer", "Argonaute", "AGO", "RDR", "eIF4G", "eIF4E", "rymv1", "translation initiation", "movement protein",
                  "coat protein", "suppressor of silencing", "vector", "leafhopper", "planthopper", "aphid", "mite vector", "resistance gene", "disease resistance",
                  "defense response", "pathogen", "NLR"]),
    # --- nematodes
    dict(key="nematode", aliases=["nematode", "root-knot", "root knot", "cyst nematode", "meloidogyne", "root-lesion", "root lesion", "pratylenchus", "heterodera",
                                  "ditylenchus", "hirschmanniella", "aphelenchoides", "white tip", "stem nematode", "ufra", "globodera", "tylenchulus"],
         direct=["nematode", "root-knot", "root-knot nematode", "Meloidogyne graminicola", "Meloidogyne", "cyst nematode", "Heterodera", "root-lesion", "Pratylenchus",
                 "Ditylenchus", "Hirschmanniella", "Aphelenchoides", "white tip", "ufra"],
         related=["gall", "giant cell", "syncytium", "feeding site", "effector", "root", "Mi-1", "Hs1", "resistance gene", "defense response", "jasmonic acid", "ethylene",
                  "hypersensitive", "disease resistance"]),
    # --- insects and mites
    dict(key="planthopper_leafhopper", aliases=["planthopper", "brown planthopper", "bph", "whitebacked", "white-backed", "wbph", "nilaparvata", "sogatella", "laodelphax",
                                                "small brown planthopper", "sbph", "leafhopper", "nephotettix", "green leafhopper", "glh", "zigzag", "hopperburn", "cicadellid",
                                                "delphacid", "hopper"],
         direct=["planthopper", "brown planthopper", "Nilaparvata lugens", "whitebacked planthopper", "Sogatella furcifera", "Laodelphax striatellus", "leafhopper",
                 "Nephotettix", "green leafhopper", "hopperburn"],
         related=["Bph1", "Bph3", "Bph6", "Bph9", "Bph14", "Bph15", "Bph18", "Bph26", "Bph29", "Bph30", "Bph32", "antixenosis", "antibiosis", "callose",
                  "sieve tube", "phloem", "honeydew", "feeding", "salivary", "jasmonic acid", "salicylic acid", "volatile", "lectin", "trypsin inhibitor", "defense response",
                  "NBS-LRR", "lectin receptor kinase", "insect", "herbivore"]),
    dict(key="stem_borer", aliases=["stem borer", "stemborer", "yellow stem borer", "scirpophaga", "chilo", "striped stem borer", "pink stem borer", "sesamia", "corn borer",
                                    "european corn borer", "ostrinia", "borer", "dead heart", "whitehead", "diatraea", "busseola", "eldana", "chilo partellus", "maruca"],
         direct=["stem borer", "borer", "yellow stem borer", "Scirpophaga incertulas", "striped stem borer", "Chilo suppressalis", "corn borer", "European corn borer",
                 "Ostrinia nubilalis", "Sesamia", "dead heart", "whitehead", "Diatraea", "Busseola", "Chilo partellus"],
         related=["Bt", "Cry1Ab", "Cry1Ac", "benzoxazinoid", "DIMBOA", "maysin", "cell wall", "silica", "lignin", "trypsin inhibitor", "protease inhibitor", "volatile",
                  "jasmonic acid", "defense response", "herbivore", "insect", "antibiosis"]),
    dict(key="lepidoptera_defoliator", aliases=["armyworm", "fall armyworm", "spodoptera", "faw", "helicoverpa", "earworm", "cutworm", "leaf folder", "leaffolder", "cnaphalocrocis",
                                                "defoliator", "caterpillar", "lepidoptera", "leaf feeder", "rice skipper", "case worm", "caseworm", "nymphula", "mythimna",
                                                "heliothis", "agrotis", "defoliat", "chewing insect"],
         direct=["armyworm", "fall armyworm", "Spodoptera frugiperda", "Spodoptera", "Helicoverpa", "earworm", "cutworm", "leaf folder", "Cnaphalocrocis medinalis", "defoliator",
                 "caterpillar", "Lepidoptera", "Mythimna", "case worm"],
         related=["maize proteinase inhibitor", "Mir1", "Mir1-CP", "cysteine protease", "benzoxazinoid", "DIMBOA", "flavonoid", "maysin", "jasmonic acid", "wound", "volatile",
                  "terpene", "trypsin inhibitor", "defense response", "herbivore", "insect", "antibiosis", "Bt"]),
    dict(key="sucking_pest", aliases=["aphid", "russian wheat aphid", "diuraphis", "greenbug", "schizaphis", "bird cherry", "rhopalosiphum", "sitobion", "whitefly", "bemisia",
                                      "thrips", "mite", "spider mite", "tetranychus", "stink bug", "rice bug", "leptocorisa", "gundhi bug", "bug", "jassid", "mealybug",
                                      "scale insect", "sugarcane aphid", "melanaphis", "sap sucking", "sap-sucking", "piercing", "frankliniella", "stenotarsonemus", "panicle mite"],
         direct=["aphid", "Russian wheat aphid", "Diuraphis noxia", "greenbug", "Schizaphis graminum", "Rhopalosiphum", "Sitobion avenae", "whitefly", "thrips", "mite",
                 "spider mite", "stink bug", "rice bug", "Leptocorisa", "sugarcane aphid", "Melanaphis sacchari"],
         related=["Dn", "Dn4", "Gb3", "Rag1", "Sg1", "antixenosis", "antibiosis", "phloem", "sieve element", "honeydew", "callose", "feeding", "salivary",
                  "effector", "jasmonic acid", "salicylic acid", "defense response", "NBS-LRR", "insect", "herbivore"]),
    dict(key="midge_fly", aliases=["gall midge", "orseolia", "hessian fly", "mayetiola", "wheat midge", "sitodiplosis", "whorl maggot", "hydrellia", "shoot fly", "atherigona",
                                   "stem fly", "frit fly", "oscinella", "cecidomyiid", "midge", "sawfly", "cephus", "wheat stem sawfly", "stem maggot", "fly"],
         direct=["gall midge", "Asian rice gall midge", "Orseolia oryzae", "Hessian fly", "Mayetiola destructor", "wheat midge", "Sitodiplosis mosellana", "whorl maggot", "shoot fly",
                 "Atherigona", "sawfly", "Cephus"],
         related=["Gm1", "Gm2", "Gm4", "Gm8", "H13", "Sm1", "H9", "H25", "hypersensitive", "gall", "feeding site", "effector", "salivary", "antibiosis", "NBS-LRR",
                  "defense response", "insect", "herbivore"]),
    dict(key="storage_pest", aliases=["weevil", "sitophilus", "rice weevil", "storage pest", "stored grain", "stored-grain", "bruchid", "callosobruchus", "storage insect",
                                      "grain borer", "prostephanus", "larger grain borer", "rhyzopertha", "tribolium", "flour beetle", "sitotroga", "angoumois", "khapra",
                                      "post-harvest pest", "postharvest pest", "post harvest insect"],
         direct=["weevil", "rice weevil", "Sitophilus oryzae", "Sitophilus", "storage pest", "stored grain", "bruchid", "Callosobruchus", "grain borer", "Prostephanus truncatus",
                 "Rhyzopertha dominica", "Tribolium", "Sitotroga"],
         related=["seed coat", "seed hardness", "kernel hardness", "amylase inhibitor", "protease inhibitor", "lectin", "phenolic", "tannin", "antibiosis", "insect", "defense response"]),
    dict(key="root_pest", aliases=["rootworm", "diabrotica", "wireworm", "white grub", "grub", "root aphid", "termite", "root feeder", "root herbivore", "root-feeding", "chafer",
                                   "root weevil", "root maggot"],
         direct=["rootworm", "western corn rootworm", "Diabrotica", "wireworm", "white grub", "termite", "root herbivore", "root feeder", "root aphid"],
         related=["Bt", "Cry3Bb", "root", "volatile", "(E)-beta-caryophyllene", "terpene synthase", "benzoxazinoid", "jasmonic acid", "herbivore", "insect", "defense response"]),
    dict(key="insect", aliases=["insect", "pest", "herbivore", "herbivory", "pest resistance", "insect resistance", "gall midge", "antixenosis", "antibiosis", "entomo"],
         direct=["insect", "pest resistance", "insect resistance", "herbivore", "herbivory", "pest", "antixenosis", "antibiosis"],
         related=["jasmonic acid", "defense response", "volatile", "trypsin inhibitor", "protease inhibitor", "callose", "phloem", "wound", "glandular trichome", "trichome",
                  "terpene", "benzoxazinoid", "Bt", "lectin", "secondary metabolite", "Bph", "antixenosis", "antibiosis", "herbivore-induced"]),
    # --- parasitic plants
    dict(key="parasitic_weed", aliases=["striga", "witchweed", "orobanche", "broomrape", "parasitic weed", "parasitic plant", "cuscuta", "dodder", "rhamphicarpa", "alectra",
                                        "phelipanche"],
         direct=["Striga", "witchweed", "Striga hermonthica", "Striga asiatica", "Orobanche", "broomrape", "parasitic weed", "parasitic plant", "Cuscuta", "dodder", "Rhamphicarpa"],
         related=["strigolactone", "germination stimulant", "LGS1", "SL biosynthesis", "D27", "CCD7", "CCD8", "MAX1", "haustorium", "low germination stimulant", "root exudate",
                  "resistance", "SL", "karrikin", "MAX2", "D14"]),
    # --- general defence and immunity
    dict(key="disease", aliases=["disease", "pathogen", "blight", "mildew", "bacterial", "fungal", "immun", "biotic", "defense", "defence", "pathogenesis", "phytopathogen",
                                 "oomycete", "necrotroph", "biotroph", "lesion mimic", "spotted leaf", "systemic acquired resistance", "induced resistance", "quantitative resistance",
                                 "durable resistance", "broad-spectrum resistance", "broad spectrum resistance", "race-specific", "r gene", "resistance gene", "host resistance",
                                 "avirulence", "effector"],
         direct=["disease resistance", "biotic stress", "pathogen", "blight", "immunity", "defense response", "susceptibility", "pathogenesis", "broad-spectrum resistance"],
         related=["defence response", "NBS-LRR", "NLR", "hypersensitive", "hypersensitive response", "salicylic acid", "jasmonic acid", "ethylene", "PR protein", "pathogenesis-related",
                  "resistance gene", "innate immun", "receptor kinase", "receptor-like kinase", "RLK", "RLCK", "BIK1", "FLS2", "CERK1", "chitin", "PAMP", "pattern-triggered",
                  "effector-triggered", "effector", "systemic acquired resistance", "NPR1", "EDS1", "PAD4", "WRKY", "MAPK", "MAP kinase", "reactive oxygen", "respiratory burst",
                  "oxidative burst", "callose", "phytoalexin", "momilactone", "phytocassane", "chitinase", "glucanase", "thaumatin", "defensin", "thionin", "lectin", "lesion mimic",
                  "cell death", "RBOH", "lignin", "cell wall reinforcement", "silicon", "lipid transfer protein"]),
    dict(key="weed", aliases=["weed", "herbicide", "allelopathy", "weed competitiveness", "glyphosate"],
         direct=["weed", "herbicide", "allelopathy", "weed competitive"],
         related=["early vigor", "seedling vigor", "canopy", "detoxification", "cytochrome P450", "ALS", "acetolactate synthase"]),
    # ------------------------------------------------------------------ yield and its components
    dict(key="yield", aliases=["yield", "productivity"],
         direct=["yield", "grain yield", "productivity"],
         related=["grain filling", "tillering", "tiller", "biomass", "harvest index", "sink", "source", "panicle", "seed set"]),
    dict(key="grain_size", aliases=["grain size", "grain weight", "seed size", "seed weight", "kernel weight", "1000-grain", "thousand grain", "tgw", "grain length", "grain width", "kernel size"],
         direct=["grain size", "grain weight", "seed size", "seed weight", "kernel weight", "grain length", "grain width", "kernel size"],
         related=["GS3", "GW2", "GW5", "TGW6", "GL3", "cell cycle", "cell expansion", "brassinosteroid", "ubiquitin", "grain filling", "endosperm", "BG1"]),
    dict(key="grain_number", aliases=["grain number", "spikelet number", "panicle", "spike", "kernel number", "grains per", "branching", "panicle architecture", "spikelets", "kernel row"],
         direct=["grain number", "spikelet number", "panicle", "spike", "kernel number", "branching", "inflorescence"],
         related=["Gn1a", "cytokinin oxidase", "DEP1", "SPL14", "IPA1", "inflorescence meristem", "meristem", "floral organ", "rachis", "tassel"]),
    dict(key="tillering", aliases=["tiller", "tillering", "branching number", "shoot branching", "tiller number", "productive tiller"],
         direct=["tiller", "tillering", "shoot branching", "branching"],
         related=["strigolactone", "MOC1", "D3", "D10", "TB1", "axillary bud", "auxin", "cytokinin", "bud outgrowth"]),
    dict(key="harvest_index", aliases=["harvest index", "biomass", "dry matter", "straw", "fodder"],
         direct=["harvest index", "biomass", "dry matter", "straw"],
         related=["source", "sink", "partitioning", "carbon assimilation", "sucrose transport", "SWEET", "starch remobilization"]),
    dict(key="filling", aliases=["grain filling", "seed filling", "ripening", "filling rate"],
         direct=["grain filling", "seed filling", "ripening", "filling rate"],
         related=["starch", "sucrose", "SUT", "sugar transporter", "endosperm", "cell wall invertase", "GIF1", "amylopectin", "senescence"]),
    dict(key="heterosis", aliases=["heterosis", "hybrid vigor", "hybrid vigour", "combining ability"],
         direct=["heterosis", "hybrid vigor", "hybrid vigour"],
         related=["hybrid", "heterozygous", "overdominance", "dominance", "epistasis", "complementation"]),
    # ------------------------------------------------------------------ phenology and plant architecture
    dict(key="flowering", aliases=["flowering", "heading date", "heading", "photoperiod", "maturity", "days to", "earliness", "flowering time", "anthesis", "phenology", "vernalization"],
         direct=["flowering", "heading date", "photoperiod", "florigen", "vernalization", "flowering time", "anthesis"],
         related=["Hd1", "Hd3a", "FT", "circadian", "inflorescence", "floral", "MADS", "Ghd7", "VRN1", "PPD", "Ehd1", "CONSTANS", "phytochrome"]),
    dict(key="height", aliases=["plant height", "height", "dwarf", "semi-dwarf", "stature", "short stature", "internode", "plant stature"],
         direct=["plant height", "dwarf", "semi-dwarf", "stature", "internode"],
         related=["gibberellin", "GA20", "GA2ox", "brassinosteroid", "cell elongation", "culm", "stem elongation", "strigolactone", "SD1", "Rht", "reduced height"]),
    dict(key="lodging", aliases=["lodging", "culm strength", "stem strength", "stalk strength", "stalk lodging", "stem breaking", "root lodging"],
         direct=["lodging", "culm strength", "stem strength", "stalk strength"],
         related=["cell wall", "cellulose", "lignin", "internode", "dwarf", "brittle", "culm diameter", "sclerenchyma"]),
    dict(key="leaf_angle", aliases=["leaf angle", "leaf erectness", "erect leaf", "plant type", "plant architecture", "canopy", "leaf inclination", "lamina joint", "leaf size", "leaf width"],
         direct=["leaf angle", "leaf inclination", "lamina joint", "erect leaf", "plant architecture", "leaf width", "leaf size", "canopy architecture"],
         related=["brassinosteroid", "auxin", "LAZY1", "lamina", "ligule", "narrow leaf", "NAL1", "cell division", "cell expansion", "tiller angle"]),
    dict(key="root", aliases=["root", "rooting", "rhizosphere", "root length", "root depth", "deep root", "root system", "root hair"],
         direct=["root", "root architecture", "root hair", "lateral root", "root system", "root length", "root depth"],
         related=["auxin", "gravitropism", "root meristem", "adventitious", "crown root", "xylem", "nutrient uptake", "DRO1", "cytokinin", "crown root"]),
    dict(key="seedling", aliases=["seedling", "early vigor", "early vigour", "emergence", "establishment", "coleoptile", "mesocotyl", "seedling vigor"],
         direct=["seedling vigor", "seedling vigour", "early vigor", "seedling", "emergence", "coleoptile", "mesocotyl"],
         related=["germination", "seed reserve", "alpha-amylase", "gibberellin", "elongation", "etiolation", "photomorphogenesis"]),
    dict(key="senescence", aliases=["senescence", "stay-green", "stay green", "leaf senescence", "chlorophyll retention", "remobilization"],
         direct=["senescence", "stay-green", "chlorophyll degradation", "leaf senescence"],
         related=["NAC", "chlorophyll", "SGR", "nutrient remobilization", "cytokinin", "autophagy", "ethylene", "ageing"]),
    dict(key="photosynthesis", aliases=["photosynthesis", "photosynthetic", "chlorophyll", "rubisco", "carbon fixation", "c4", "leaf photosynthesis", "radiation use", "light use", "stomatal conductance"],
         direct=["photosynthesis", "photosynthetic", "chlorophyll", "Rubisco", "carbon fixation", "stomatal conductance"],
         related=["photosystem", "chloroplast", "light harvesting", "Calvin cycle", "carbonic anhydrase", "PEP carboxylase", "electron transport", "CO2 assimilation", "stomatal"]),
    # ------------------------------------------------------------------ seed biology and reproduction
    dict(key="dormancy", aliases=["dormancy", "pre-harvest sprouting", "preharvest sprouting", "sprouting", "seed longevity", "seed ageing", "storability", "germination", "viviparous"],
         direct=["dormancy", "pre-harvest sprouting", "preharvest sprouting", "seed longevity", "germination", "viviparous", "seed storability"],
         related=["abscisic acid", "ABA", "gibberellin", "DOG1", "MFT", "Sdr4", "seed coat", "alpha-amylase", "after-ripening", "VP1"]),
    dict(key="shattering", aliases=["shattering", "seed shattering", "pod shattering", "threshability", "free threshing", "domestication", "abscission"],
         direct=["shattering", "abscission", "threshability", "free threshing", "domestication"],
         related=["abscission layer", "sh4", "qSH1", "Q gene", "AP2", "dehiscence", "lignin", "cell wall degradation", "pod"]),
    dict(key="sterility", aliases=["male sterility", "sterility", "fertility", "tgms", "pollen", "restorer", "anther", "pollen sterility", "spikelet fertility", "seed set", "hybrid breeding", "wide compatibility"],
         direct=["male sterility", "sterility", "fertility restorer", "pollen", "anther", "cytoplasmic male", "pollen fertility", "spikelet fertility"],
         related=["tapetum", "pollen wall", "exine", "mitochondrial", "PPR", "restorer of fertility", "temperature-sensitive", "photoperiod-sensitive", "meiosis", "S5", "hybrid incompatibility"]),
    dict(key="apomixis", aliases=["apomixis", "haploid induction", "parthenogenesis", "clonal seed"],
         direct=["apomixis", "haploid induction", "parthenogenesis", "clonal"],
         related=["BBM", "BABY BOOM", "MTL", "matrilineal", "egg cell", "embryo sac", "meiosis", "PAIR1"]),
    # ------------------------------------------------------------------ grain and product quality
    dict(key="amylose", aliases=["amylose", "cooking quality", "eating quality", "eating", "texture", "gel consistency", "gelatinization", "gelatinisation", "stickiness", "waxy", "glutinous", "starch"],
         direct=["amylose", "cooking quality", "eating quality", "gel consistency", "gelatinization", "gelatinisation", "waxy", "glutinous", "starch"],
         related=["Wx", "GBSS", "granule-bound starch synthase", "SSIIa", "starch synthase", "amylopectin", "starch branching", "endosperm", "ADP-glucose", "pullulanase", "texture"]),
    dict(key="glycemic", aliases=["glycemic", "glycaemic", "resistant starch", "starch digestib", "digestibility"],
         direct=["glycemic", "glycaemic", "resistant starch", "starch digestibility", "amylose", "granule-bound starch synthase"],
         related=["starch", "starch synthase", "waxy", "GBSS", "amylopectin", "starch branching", "ADP-glucose", "endosperm",
                  "gelatinization", "sucrose", "sugar", "carbohydrate"]),
    dict(key="protein", aliases=["protein content", "protein", "gluten", "storage protein", "grain protein", "prolamin", "glutelin", "zein", "gliadin", "glutenin", "lysine", "essential amino"],
         direct=["protein content", "grain protein", "gluten", "storage protein", "glutelin", "prolamin", "zein", "gliadin", "glutenin", "lysine"],
         related=["amino acid", "nitrogen remobilization", "endosperm", "NAC", "opaque", "seed storage", "GPC", "amino acid transporter"]),
    dict(key="oil", aliases=["oil content", "oil", "fatty acid", "oleic", "linoleic", "lipid", "seed oil", "fat content", "oil quality"],
         direct=["oil content", "fatty acid", "oleic", "linoleic", "lipid", "seed oil", "triacylglycerol"],
         related=["FAD2", "FAD3", "DGAT", "WRI1", "oleosin", "desaturase", "acyl-CoA", "lipid biosynthesis", "stearoyl"]),
    dict(key="aroma", aliases=["aroma", "fragrance", "scent", "2-acetyl", "basmati", "badh2"],
         direct=["aroma", "fragrance", "scent", "2-acetyl-1-pyrroline", "betaine aldehyde dehydrogenase", "BADH2"],
         related=["pyrroline", "volatile", "proline", "aldehyde dehydrogenase", "badh2"]),
    dict(key="grain_appearance", aliases=["chalkiness", "chalk", "grain appearance", "translucency", "head rice", "milling", "grain shape", "grain color", "pericarp", "red rice", "pigment", "kernel hardness", "test weight"],
         direct=["chalkiness", "chalk", "grain appearance", "translucent", "head rice", "milling quality", "pericarp", "grain color", "kernel hardness"],
         related=["endosperm", "starch granule", "amyloplast", "storage protein", "Rc", "proanthocyanidin", "flavonoid", "PPDK", "hardness", "puroindoline"]),
    dict(key="bran", aliases=["anthocyanin", "polyphenol", "flavonoid", "carotenoid", "vitamin", "provitamin", "antioxidants", "phenolic", "tocopherol", "folate", "lutein", "nutritional quality", "nutrient content"],
         direct=["anthocyanin", "polyphenol", "flavonoid", "carotenoid", "vitamin", "provitamin A", "phenolic", "tocopherol", "folate", "lutein"],
         related=["phytoene synthase", "PSY", "lycopene", "beta-carotene", "chalcone synthase", "MYB", "bHLH", "pigment", "secondary metabolite", "biosynthesis"]),
    dict(key="antinutrient", aliases=["phytate", "phytic acid", "antinutrient", "low phytic", "lpa", "raffinose", "tannin", "glucosinolate", "gossypol", "lectin", "trypsin inhibitor", "mycotoxin", "aflatoxin"],
         direct=["phytate", "phytic acid", "antinutrient", "raffinose", "tannin", "glucosinolate", "aflatoxin"],
         related=["myo-inositol", "inositol", "MRP", "ABC transporter", "inositol phosphate", "detoxification", "oxalate"]),
    dict(key="sugar", aliases=["sugar", "sucrose", "sweetness", "brix", "sugar content", "sweet corn", "sugarcane", "starch content", "sugar accumulation"],
         direct=["sugar", "sucrose", "sweetness", "Brix", "sweet corn", "starch content"],
         related=["invertase", "sucrose synthase", "SUT", "SWEET", "vacuolar", "sugar transporter", "sucrose phosphate synthase", "hexose"]),
    dict(key="fibre", aliases=["fibre", "fiber", "lignin", "cellulose", "cell wall", "digestibility", "forage quality", "bioenergy", "biomass quality"],
         direct=["fibre", "fiber", "lignin", "cellulose", "cell wall", "forage quality"],
         related=["brown midrib", "bm", "COMT", "CAD", "cellulose synthase", "xylan", "hemicellulose", "secondary cell wall", "NAC", "MYB"]),
    # ------------------------------------------------------------------ regulation and signalling (useful as stand-alone searches)
    dict(key="hormone", aliases=["auxin", "gibberellin", "cytokinin", "brassinosteroid", "strigolactone", "ethylene", "jasmonate", "jasmonic", "abscisic", "salicylic", "hormone"],
         direct=["auxin", "gibberellin", "cytokinin", "brassinosteroid", "strigolactone", "ethylene", "jasmonic acid", "abscisic acid", "salicylic acid"],
         related=["hormone signaling", "hormone biosynthesis", "receptor", "transcription factor", "response regulator", "AUX/IAA", "ARF", "DELLA", "GID1"]),
    dict(key="transcription", aliases=["transcription factor", "regulatory gene"],
         direct=["transcription factor", "NAC", "MYB", "WRKY", "bZIP", "AP2", "ERF", "bHLH"],
         related=["DNA-binding", "regulation of transcription", "transcriptional regulator", "zinc finger", "homeobox", "MADS", "GRAS"]),
    dict(key="domestication", aliases=["domestication", "landrace", "wild relative", "introgression", "selective sweep"],
         direct=["domestication", "adaptation", "landrace", "wild relative", "introgression"],
         related=["shattering", "dormancy", "tillering", "plant architecture", "grain size", "flowering time", "pigment", "seed dispersal"]),
]


def _words(trait):
    return [w for w in re.findall(r"[A-Za-z][A-Za-z\-]+", trait) if w.lower() not in GENERIC and len(w) > 2]


def expand(trait):
    """{'trait', 'direct', 'related', 'concepts', 'note'} for a trait name."""
    low = trait.lower()
    direct, related, concepts = [], [], []
    for c in CONCEPTS:
        if any(re.search(r"(?<![a-z])" + re.escape(a), low) for a in c["aliases"]):
            concepts.append(c["key"])
            direct += c["direct"]
            related += c["related"]
    own = _words(trait)
    if not concepts:                      # a curated concept already holds the trait's words; otherwise use the words typed
        direct = list(dict.fromkeys(own))
    related = [r for r in dict.fromkeys(related) if r.lower() not in {d.lower() for d in direct}]
    note = ""
    if not concepts:
        note = (f"The trait '{trait}' has no built-in list of related terms, so only its own words ({', '.join(own) or trait}) were looked for. "
                "Genes whose records do not use these words get no trait evidence; try the wording the databases use (a Trait Ontology name such as 'amylose content' "
                "or a function such as 'starch synthase').")
    if not direct:
        direct = [trait.strip()]
    return {"trait": trait.strip(), "direct": direct, "related": related, "concepts": concepts, "note": note}


def _pattern(term):
    t = re.escape(term)
    if len(term) <= 4 and term.upper() == term:          # abbreviations (ABA, HSP, NUE): whole word, upper case
        return re.compile(r"(?<![A-Za-z0-9])" + t + r"(?![A-Za-z0-9])")
    return re.compile(r"(?<![a-z0-9])" + t, re.I)


def find_terms(text, terms):
    """Terms of `terms` found in `text` (a list of the terms that matched)."""
    if not text:
        return []
    return [t for t in terms if _pattern(t).search(text)]


def match(text, tt):
    """[(term, kind)] with kind 'direct' or 'related' for the terms of `tt` found in `text`."""
    return [(t, "direct") for t in find_terms(text, tt["direct"])] + [(t, "related") for t in find_terms(text, tt["related"])]


def search_words(tt, limit=6):
    """Words for the PubMed and Plant Reactome searches."""
    return tt["direct"][:limit]
