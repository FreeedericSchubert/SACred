# SACred – pitch script (3:30 min, ~470 words)

Cues in [brackets] follow the sections of the pitch notebook (notebooks/sac_pitch.py).

**[0:00 · Title]**
Nature solved single-atom catalysis long before we did. Enzymes make hydrogen, reduce oxygen and CO₂, even fix nitrogen, each on a single metal atom held by a few well-chosen neighbours. We are SACred, and we turn these enzyme active sites into usable single-atom catalysts.

**[0:20 · 1 · The problem]**
Single-atom catalysts use every metal atom, so you need far less of scarce and expensive metals. But the field has two problems. First, almost everyone makes the same few sites, mostly a metal held by four nitrogens in graphene. The design space is huge and barely explored. Second, many sites don't survive: during synthesis or in the reactor, the metal atoms escape and clump into particles. So we asked: can we use nature's library of designs, and check right away which ones can actually be made?

**[0:50 · 2 · What we built – start the zoom]**
We built an agentic AI pipeline. We screened the entire Protein Data Bank, almost 300,000 metal sites. For each one, we take the atoms that directly hold the metal and rebuild them in graphene, MoS₂ or silica. Here you see it: from the enzyme, to its active site, to the same site in graphene.
Then a funnel: descriptors and fast semi-empirical calculations narrow over a thousand candidates down to the best ones, which we treat with DFT. And the DFT data works twice: it calibrates the fast screen, and it trains a machine-learned force field with MACE. So the pipeline gets cheaper and better with every round.

**[1:35 · 3 · Agents and tools]**
AI runs through the whole project. Claude wrote and ran the pipeline and managed more than 230 DFT calculations on our cluster. With Valency we ran a multi-agent literature search: about a hundred agent runs collected simulated and experimental data on known single-atom catalysts to benchmark our models, and found the biggest problems in making them. marimo brings it all together as interactive apps, including this pitch.

**[2:00 · 4 · What it finds – overview and gallery]**
What do we find? Still preliminary, but already 24 usable motifs that bind hydrogen close to the optimum, several of them never reported before, plus first candidates on MoS₂ and silica and for oxygen and CO₂ reduction. And without being told, our screen found a nickel site that experimentalists have already made and tested for hydrogen evolution. That tells us the approach works.

**[2:25 · Will it survive synthesis? – play the MD video]**
But can you make them? This is our force field at 1500 Kelvin, simulating synthesis. The weak site loses its metal in a fraction of a picosecond. Our candidate keeps it, exactly as the screen predicted. Activity and stability, ranked together.

**[2:45 · 5 · What comes next]**
And this scales. The pipeline reruns with every PDB release, for more metals, supports and reactions, and the force field takes over more of the expensive DFT. Next we add the enzyme's second coordination sphere and real conditions like pH and potential, and then we make the best motifs with experimental partners. Everything is open: the dataset, with over 11,000 DFT structures for force-field training, is on Zenodo.

**[3:20 · Close]**
SACred: from enzymes to single-atom catalysts that can actually be made. Thank you.

---

# Short version for the video (2:00, ~240 words)

**[0:00 · Title]**
Enzymes run hydrogen evolution, oxygen and CO₂ reduction on single metal atoms. We are SACred, and we turn these enzyme active sites into usable single-atom catalysts.

**[0:10 · The problem]**
Single-atom catalysts save scarce metals, but almost everyone makes the same few sites, and many don't survive: the metal escapes and clumps during synthesis. Nature has thousands of other designs. Which of them can we actually make?

**[0:25 · What we built – zoom]**
We built an agentic AI pipeline. It takes all 300,000 metal sites of the Protein Data Bank and rebuilds each one in graphene, MoS₂ or silica. Descriptors and fast xTB calculations funnel them down to the best candidates for DFT. The DFT data then calibrates the screen and trains a MACE force field.

**[0:50 · Agents and tools]**
Claude wrote and ran the pipeline. Valency agents searched the literature for benchmark data and the biggest synthesis problems. marimo turns it all into interactive apps.

**[1:05 · What it finds – MD video]**
Preliminary, but already 24 usable motifs, several never reported, and our screen independently found a nickel site that has already been made for hydrogen evolution. Our force field at 1500 Kelvin shows a weak site losing its metal within a fraction of a picosecond, while our candidate holds on, as predicted.

**[1:35 · What comes next]**
The pipeline scales with every PDB release, to more metals, supports and reactions. Next: real operating conditions, and synthesis with experimental partners. The full dataset is open on Zenodo.

**[1:50 · Close]**
SACred: from enzymes to single-atom catalysts that can actually be made.
