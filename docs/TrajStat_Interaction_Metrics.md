# TrajStat Interaction Metrics

This document records the hydrogen-bond and ionic-interaction definitions used by
`amia/trajstat.py`.

## Hydrogen Bonds

Hydrogen bonds are calculated with
[MDAnalysis `HydrogenBondAnalysis`](https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html).
The analysis reports one row for every observed bond in every analyzed frame.

### Geometric criteria

| Metric | Value in TrajStat | Meaning |
|---|---:|---|
| Donor-hydrogen distance | 1.2 A | MDAnalysis `d_h_cutoff` default; used when donor-hydrogen pairs must be identified from coordinates |
| Donor-acceptor distance | 4.5 A or 3.5 A | MDAnalysis `d_a_cutoff`; 4.5 A for protein-other-molecule and molecule-molecule analyses, 3.5 A for protein-nucleic analyses |
| D-H-A angle | 150 degrees | MDAnalysis `d_h_a_angle_cutoff` default |
| Selection updates | `False` or `True` | Fixed selections are used for ligand-like analyses; nucleic and molecule-pair analyses update selections each frame |

The donor-hydrogen distance is not a hydrogen-bond cutoff when the topology
contains bond information. In that case MDAnalysis uses the topology to identify
donor-hydrogen pairs. A PSF, TPR, or PRMTOP topology is therefore preferred.

### Donors, hydrogens, and acceptors

TrajStat does not hard-code individual donor and acceptor atom names for these
analyses. It uses MDAnalysis selection guessers:

- **Hydrogens:** `guess_hydrogens("protein")` or a molecule-specific selection.
  MDAnalysis identifies candidate hydrogens using mass and charge.
- **Acceptors:** `guess_acceptors(...)`. MDAnalysis identifies candidate
  acceptors using atomic charge and the supplied selection group.
- **Donor atoms:** donor-hydrogen relationships are obtained from topology bond
  information when available. If bond information is absent, donor assignment
  depends on the donor selection and the 1.2 A donor-hydrogen cutoff.

The current analysis paths use these selections:

| Function | Hydrogen selection | Acceptor selection | D-A cutoff |
|---|---|---|---:|
| `hbond_calc` | `guess_hydrogens("protein")` | `guess_acceptors("segid <other molecule>")` | 4.5 A |
| `nucleic_prot_hbonds` | `guess_hydrogens("protein")` | `guess_acceptors("nucleic")` | 3.5 A |
| `hetatm_calc` | `guess_hydrogens("segid <molecule>")` | `guess_acceptors("segid <other molecule>")` | 4.5 A |

All three paths use the 150-degree D-H-A angle criterion unless explicitly
changed in the `HydrogenBondAnalysis` constructor.

### Hydrogen-bond output metrics

For each analysis, TrajStat writes:

- `<name>_hbonds.csv`: time in ns and `count_by_time()`, the number of observed
  hydrogen bonds at each analyzed frame.
- `<name>_hbonds_observations.csv`: frame, donor index, hydrogen index, acceptor
  index, D-A distance in A, and D-H-A angle in degrees from
  `results.hbonds`.
- `<name>_hbonds_by_id.csv`: `count_by_ids()`, the number of observations for
  each donor-hydrogen-acceptor atom combination.
- `<name>_hbonds_by_type.csv`: `count_by_type()`, the number of observations for
  each donor/acceptor type combination.
- `<name>_hbonds_lifetime.csv`: the hydrogen-bond time autocorrelation from
  `lifetime()`, reported as tau in analyzed frames and autocorrelation value.

The lifetime calculation is bounded by the number of analyzed frames. Runs with
fewer than two analyzed frames produce an empty lifetime table because an
autocorrelation window cannot be calculated.

## Ionic Interactions

Ionic interactions are treated as contact ion pairs based on the article
["Dynamics of Ionic Interactions at Protein-Nucleic Acid Interfaces"](https://pmc.ncbi.nlm.nih.gov/articles/PMC7497705/).
The article describes contact ion pairs using direct heavy-atom contact and
reports protein-phosphate O...N distances below 6 A. TrajStat follows the
MDAnalysis contact-analysis example and uses its 4.5 A atom-contact cutoff for
the primary per-frame count. The 6 A article distance remains a useful broader
contact-ion-pair definition, but is not the primary TrajStat count.

Ionic interactions do **not** use a donor-hydrogen angle. The 150-degree angle
criterion applies only to hydrogen bonds.

### Protein-protein salt bridges

The selections are:

- **Cationic group:** `(resname ARG LYS) and (name NH* NZ)`
- **Anionic group:** `(resname ASP GLU) and (name OE* OD*)`
- **Contact criterion:** cation-anion atom pairs are counted when closer than
  4.5 A
- **Counting unit:** unique cation-residue/anion-residue pair per frame, so
  multiple qualifying atom pairs within one residue pair count once

### Protein-nucleic ionic interactions

The selections are:

- **Cationic group:** `(resname ARG LYS) and (name NH* NZ)`
- **Anionic group:** `nucleic and (name OP1 OP2 O1P O2P)`
- **Contact criterion:** selected cation-phosphate-oxygen atom pairs are counted
  when closer than 4.5 A
- **Counting unit:** unique cation-residue/anion-residue pair per frame

Using phosphate oxygen atoms rather than phosphorus measures the relevant
protein-side-chain-to-phosphate heavy-atom distance.

### Ionic-interaction output metrics

For each ionic analysis, TrajStat writes:

- `<name>_ionic_interactions.csv`: frame, time in ns, and the number of atom
  contacts within 4.5 A for that frame, matching the MDAnalysis
  `contact_matrix()` method.
- `<name>_ionic_pair_occupancy.csv`: cation residue ID/name, anion residue
  ID/name, frames in contact, occupancy, minimum distance in A, and mean distance
  in A.
- The legacy salt-bridge CSV (`*_salt_bridges.csv` or
  `*_Nucleic_Acid_saltbridges.csv`) retains the per-frame count used by the
  existing plotting functions.

The additional residue-pair occupancy table uses the same 4.5 A cutoff and is
calculated as:

```text
occupancy = frames in contact / total analyzed frames
```

## Units and Interpretation

- Distances are in Angstroms (A).
- Hydrogen-bond angles are in degrees.
- Trajectory times in the generated interaction CSV files are in ns.
- Hydrogen-bond lifetime tau is expressed in analyzed frames, not ns. Convert
  tau to time using the trajectory frame spacing when needed.
- A hydrogen-bond observation is an atom-level geometric event. An ionic contact
  count is a residue-pair-level event. These counts should not be compared as if
  they represented the same interaction definition.
