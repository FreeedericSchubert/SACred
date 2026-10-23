# SACred: turning enzyme active sites into usable single-atom catalysts

![From enzyme to graphene: [NiFe] hydrogenase (PDB 4U9H), the first shell of its Ni, and the graphene sites built from it](sacred_enzyme_to_graphene.png)

## Approach

Enzymes already run hydrogen evolution, O₂ and CO₂ reduction and N₂ fixation on single metal atoms. SACred uses them as templates for single-atom catalysts (SACs): an agentic AI pipeline and open dataset that turns every metal site in the Protein Data Bank into SAC candidates and checks whether they can actually be made.

- **PDB sweep:** 291,977 metal sites, each reduced to its first coordination shell and rebuilt as a site in graphene, MoS₂ or silica.
- **Descriptor funnel:** a fast semi-empirical screen (GFN2-xTB) with descriptors for activity and stability narrows 1,155 sites down to the best candidates for DFT.
- **DFT and machine learning:** the DFT data calibrate the xTB screen and train a MACE force field, which simulates the high-temperature synthesis.
- **Multi-agent literature screen (Valency):** agents searched for simulated and experimental SAC data to benchmark our models, check novelty, and find the biggest problems in SAC synthesis.
- **Visualisation (marimo):** interactive apps for the datasets and results, with a 3D view from enzyme to catalyst and a chat agent over the database.

## Results

Preliminary results:

- **24 usable motifs** bind hydrogen near the optimum in DFT, out of 71 computed; several have no earlier report.
- **Stability is the bottleneck:** most enzyme-like motifs would lose their metal during synthesis; N and C donors hold the metal, S and O let it go. The pipeline ranks motifs by activity and stability at once.
- **Validated:** our DFT reproduces published O₂-reduction data, and the screen independently finds a site (Ni–N₃) that has been made and tested for H₂ evolution.
- **Simulated synthesis:** MD at 1500 K with our force field shows a weak site losing its metal within a fraction of a picosecond while a promising motif keeps it, as the screen predicted.
- **Open dataset:** every PDB metal site with its first shell, ~11,000 DFT frames and ~5,000 xTB structures (energies and forces), released on Zenodo.

## Next steps

- **Scale:** all transition metals, more supports and reactions; the pipeline reruns with every PDB release, and the force field takes over more of the DFT.
- **Closer to reality:** the enzyme's second coordination sphere, pH and electrode potential, and competing particles and sulfides.
- **Make it:** synthesis of the most promising motifs with experimental partners, and a machine-learned stability model trained on the open dataset.
