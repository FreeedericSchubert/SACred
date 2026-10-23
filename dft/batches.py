"""Motif lists for the DFT batches (motif string as in materials/graphene/sites.py, enzyme it comes from).

Batch 1 was fixed from the textbook first shells of the HER/ORR/CO2RR/NRR enzymes so the Curta runs could
start the first night; the PDB sweep (bio/) then checks it and supplies batch 2 (motifs it finds that are
missing here).
"""

GRAPHENE_BATCH1 = [
    # HER: hydrogenases
    ("Ni-S4:DV", "[NiFe] hydrogenase Ni (4 Cys)"),
    ("Fe-C2S2:DV:cis5", "[NiFe]/[FeFe] hydrogenase Fe (CN/CO carbons + 2 bridging S)"),
    # CO2RR: acetyl-CoA synthase, CO dehydrogenase, formate dehydrogenase, methyl-CoM reductase, carbonic anhydrase
    ("Ni-N2S2:DV:cis5", "acetyl-CoA synthase Ni_d (2 amide N, 2 Cys, cis)"),
    ("Ni-N2S2:DV:cis6", "same, other cis isomer (isomer sensitivity)"),
    ("Ni-S3:SV", "Ni-Fe CO dehydrogenase C-cluster Ni"),
    ("Mo-S4:DV", "formate dehydrogenase Mo-bisMGD (4 dithiolene S)"),
    ("W-S4:DV", "W-formate dehydrogenase / formylmethanofuran dehydrogenase"),
    ("Ni-N4:DV", "methyl-coenzyme M reductase F430"),
    ("Zn-N3:SV", "carbonic anhydrase (3 His)"),
    # NRR: nitrogenase
    ("Mo-NS3:DV", "nitrogenase FeMo-co Mo (3 S + His N)"),
    ("Fe-S3:SV", "nitrogenase FeMo-co belt Fe (3 S)"),
    ("Fe-S4:DV", "Fe-S clusters, P-cluster (4 S)"),
    # ORR: heme-copper oxidases, multicopper oxidases
    ("Fe-N4:DV", "heme a3 / heme b (porphyrin N4)"),
    ("Fe-N4:DV:ax=NH3", "heme with proximal His (NH3 as the axial N)"),
    ("Cu-N3:SV", "CuB of cytochrome c oxidase, T3 Cu of laccase (3 His)"),
    ("Cu-N2S:SV", "T1 Cu of laccase (2 His + Cys)"),
    # benchmark with many literature values
    ("Co-N4:DV", "benchmark: Co-N4 is the most studied HER single-atom site"),
]

# Batch 2: motifs of data/motifs_enzyme.csv (PDB sweep, 2026-10-20) among the 25 most frequent that batch 1 lacks;
# the count is the number of enzyme entries. Structural Zn (Zn-S3, Zn-N2O2) is left out. O donors at a vacancy
# edge are unusual in graphene; the O-rich motifs are run anyway, and also belong on silica.
GRAPHENE_BATCH2 = [
    ("Fe-NS2:SV", "His-ligated Fe of Fe-S clusters, 34 entries (4U9H, 3U7Q)"),
    ("Zn-NS2:SV", "beta-carbonic anhydrase (Cys2His), 25 entries (6LAC)"),
    ("Cu-N2O:SV", "T2 Cu of multicopper oxidases, 16 entries (6RI6)"),
    ("Cu-NS2:SV", "CuA of cytochrome c oxidase, 16 entries (5B1A)"),
    ("Fe-O3:SV", "non-heme Fe, carboxylate/water shell, 12 entries (6B8F)"),
    ("Fe-S2C:SV", "FeMo-co Fe bound to the carbide, 11 entries / 83 sites (3U7Q)"),
    ("Fe-O4:DV", "ferritin ferroxidase centre, 10 entries (5J93)"),
    ("Fe-NO3:DV", "non-heme diiron Fe, 7 entries (5JSK)"),
    ("Fe-NO2:SV", "non-heme Fe, 7 entries (8PP3)"),
    ("Mo-NS2O:DV:SNSO", "nitrogenase Mo as the sweep maps it (cis S2 + His N + homocitrate O), 7 entries (3U7Q)"),
    ("Ni-S3O:DV", "[NiFe]/[NiFeSe] hydrogenase Ni with an O ligand, 6 entries (4UQL)"),
]
