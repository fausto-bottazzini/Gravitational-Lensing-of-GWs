# Bibliography

This folder is where the source papers live on a local checkout. **The PDFs
themselves are not in the repository**: they are published articles, and their
copyright belongs to the journals (APS, IOP/AAS, OUP, Elsevier, Annual
Reviews), not to this project. Redistributing them here would be republishing
someone else's copyrighted work, so `bibliography/*.pdf` is in `.gitignore`
and the folder arrives empty.

Download them yourself from the links below. Every entry has a DOI; where a
free preprint exists, the arXiv link is given and is the one to use, since it
needs no subscription. Authoritative bibliographic data — the same DOIs, in
BibTeX — is in `theory/refs.bib`, which is what `theory/theory.tex` actually
compiles against; this file is the human-readable version of it plus the
download links.

Nothing in this repository reads these PDFs. They are here for a reader who
wants to check a citation against its source, and every equation cited in
`theory/theory.pdf` names the paper, the equation number and the section, so
a reader can go straight to the right page.

If you refill the folder, the naming convention already used is

```
<Title> - <Author surnames, comma-separated> - <Year>.pdf
```

## Free full text (arXiv)

Eight of the twenty are on arXiv. The preprint is the same physics as the
published version; equation numbering can differ slightly between the two, and
`theory/theory.pdf` cites the **published** numbering throughout.

| Reference | DOI | arXiv |
|---|---|---|
| Takahashi & Nakamura (2003), *Wave effects in the gravitational lensing of gravitational waves from chirping binaries*, ApJ **595**, 1039 | [10.1086/377430](https://doi.org/10.1086/377430) | [astro-ph/0305055](https://arxiv.org/abs/astro-ph/0305055) |
| D'Orazio & Loeb (2020), *Repeated gravitational lensing of gravitational waves in hierarchical black hole triples*, PRD **101**, 083031 | [10.1103/PhysRevD.101.083031](https://doi.org/10.1103/PhysRevD.101.083031) | [1910.02966](https://arxiv.org/abs/1910.02966) |
| Ulmer & Goodman (1995), *Femtolensing: beyond the semiclassical approximation*, ApJ **442**, 67 | [10.1086/175422](https://doi.org/10.1086/175422) | [astro-ph/9406042](https://arxiv.org/abs/astro-ph/9406042) |
| Buonanno, Iyer, Ochsner, Pan & Sathyaprakash (2009), *Comparison of post-Newtonian templates for compact binary inspiral signals in gravitational-wave detectors*, PRD **80**, 084043 | [10.1103/PhysRevD.80.084043](https://doi.org/10.1103/PhysRevD.80.084043) | [0907.0700](https://arxiv.org/abs/0907.0700) |
| Khan, Husa, Hannam, Ohme, Pürrer, Jiménez Forteza & Bohé (2016), *Frequency-domain gravitational waves from nonprecessing black-hole binaries. II*, PRD **93**, 044007 | [10.1103/PhysRevD.93.044007](https://doi.org/10.1103/PhysRevD.93.044007) | [1508.07253](https://arxiv.org/abs/1508.07253) |
| Berti, Cardoso & Will (2006), *Gravitational-wave spectroscopy of massive black holes with the space interferometer LISA*, PRD **73**, 064030 | [10.1103/PhysRevD.73.064030](https://doi.org/10.1103/PhysRevD.73.064030) | [gr-qc/0512160](https://arxiv.org/abs/gr-qc/0512160) |
| Naoz (2016), *The eccentric Kozai–Lidov effect and its applications*, ARA&A **54**, 441 | [10.1146/annurev-astro-081915-023315](https://doi.org/10.1146/annurev-astro-081915-023315) | [1601.07175](https://arxiv.org/abs/1601.07175) |
| Blaes, Lee & Socrates (2002), *The Kozai mechanism and the evolution of binary supermassive black hole systems*, ApJ **578**, 775 | [10.1086/342655](https://doi.org/10.1086/342655) | [astro-ph/0203370](https://arxiv.org/abs/astro-ph/0203370) |

## DOI only

These predate arXiv (or, for Mardling & Aarseth and Nakamura & Deguchi, were
never posted there — checked, not assumed), so the DOI is the only pointer.
The AAS journals (ApJ, AJ) and MNRAS have scanned full text, free to read, on
[NASA ADS](https://ui.adsabs.harvard.edu/); search the DOI there. The APS
articles (Phys. Rev., PRD, PRL) and the Elsevier one are paywalled and need an
institutional subscription.

| Reference | DOI |
|---|---|
| Deguchi & Watson (1986), *Diffraction in gravitational lensing for compact objects of low mass*, ApJ **307**, 30 | [10.1086/164389](https://doi.org/10.1086/164389) |
| Nakamura & Deguchi (1999), *Wave optics in gravitational lensing*, Prog. Theor. Phys. Suppl. **133**, 137 | [10.1143/PTPS.133.137](https://doi.org/10.1143/PTPS.133.137) |
| Paczyński (1986), *Gravitational microlensing by the galactic halo*, ApJ **304**, 1 | [10.1086/164140](https://doi.org/10.1086/164140) |
| Peters (1974), *Index of refraction for scalar, electromagnetic, and gravitational waves in weak gravitational fields*, PRD **9**, 2207 | [10.1103/PhysRevD.9.2207](https://doi.org/10.1103/PhysRevD.9.2207) |
| Peters (1964), *Gravitational radiation and the motion of two point masses*, Phys. Rev. **136**, B1224 | [10.1103/PhysRev.136.B1224](https://doi.org/10.1103/PhysRev.136.B1224) |
| Peters & Mathews (1963), *Gravitational radiation from point masses in a Keplerian orbit*, Phys. Rev. **131**, 435 | [10.1103/PhysRev.131.435](https://doi.org/10.1103/PhysRev.131.435) |
| Isaacson (1968), *Gravitational radiation in the limit of high frequency. I*, Phys. Rev. **166**, 1263 | [10.1103/PhysRev.166.1263](https://doi.org/10.1103/PhysRev.166.1263) |
| Regge & Wheeler (1957), *Stability of a Schwarzschild singularity*, Phys. Rev. **108**, 1063 | [10.1103/PhysRev.108.1063](https://doi.org/10.1103/PhysRev.108.1063) |
| Zerilli (1970), *Effective potential for even-parity Regge–Wheeler gravitational perturbation equations*, PRL **24**, 737 | [10.1103/PhysRevLett.24.737](https://doi.org/10.1103/PhysRevLett.24.737) |
| Teukolsky (1973), *Perturbations of a rotating black hole. I*, ApJ **185**, 635 | [10.1086/152444](https://doi.org/10.1086/152444) |
| Mardling & Aarseth (2001), *Tidal interactions in star cluster simulations*, MNRAS **321**, 398 | [10.1046/j.1365-8711.2001.03974.x](https://doi.org/10.1046/j.1365-8711.2001.03974.x) |
| Kozai (1962), *Secular perturbations of asteroids with high inclination and eccentricity*, AJ **67**, 591 | [10.1086/108790](https://doi.org/10.1086/108790) |
| Lidov (1962), *The evolution of orbits of artificial satellites of planets under the action of gravitational perturbations of external bodies*, Planet. Space Sci. **9**, 719 | [10.1016/0032-0633(62)90129-0](https://doi.org/10.1016/0032-0633(62)90129-0) |

## Books

Cited by chapter and section rather than by equation number, and not
distributed in any form here.

| Reference |
|---|
| Maggiore (2008), *Gravitational Waves, Volume 1: Theory and Experiments*, Oxford University Press |
| Schneider, Ehlers & Falco (1992), *Gravitational Lenses*, Springer-Verlag |
| Chandrasekhar (1983), *The Mathematical Theory of Black Holes*, Oxford University Press |
| Murray & Dermott (1999), *Solar System Dynamics*, Cambridge University Press |

## Software

| Reference |
|---|
| PyCBC — Nitz, Harry, Brown, Biwer *et al.*, <https://github.com/gwastro/pycbc> (open source, GPL) |
