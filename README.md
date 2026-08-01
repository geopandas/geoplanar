
# geoplanar
Planar enforcement for polygon geoseries

![mexico-us](https://i.imgur.com/CFgnecL.png)

[![unittests](https://github.com/sjsrey/geoplanar/workflows/.github/workflows/unittests.yml/badge.svg)](https://github.com/sjsrey/geoplanar/actions?query=workflow%3A.github%2Fworkflows%2Funittests.yml)
|[![Documentation Status](https://readthedocs.org/projects/geoplanar/badge/?version=latest)](https://geoplanar.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/382492314.svg)](https://zenodo.org/badge/latestdoi/382492314)




`geoplanar` supports the detection and correction of violations of
[planar enforcement](https://ibis.geog.ubc.ca/courses/klink/gis.notes/ncgia/u12.html#SEC12.6)
for polygon geoseries including:


- [gaps](https://github.com/sjsrey/geoplanar/blob/main/notebooks/gaps.ipynb)
- [nonplanar coincident edges](https://github.com/sjsrey/geoplanar/blob/main/notebooks/nonplanaredges.ipynb)
- [nonplanar touches](https://github.com/sjsrey/geoplanar/blob/main/notebooks/nonplanartouches.ipynb)
- [overlaps](https://github.com/sjsrey/geoplanar/blob/main/notebooks/overlaps.ipynb)
- [holes](https://github.com/sjsrey/geoplanar/blob/main/notebooks/holes.ipynb)


## Status

`geoplanar` is currently in beta and is open to contributions. It is being developed in
support of
[research ](https://nsf.gov/awardsearch/showAward?AWD_ID=1759746&HistoricalAwards=false)
and is likely to be undergoing changes as the project evolves.

## Contributing

`geoplanar` development uses a
[git-flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow)
model. Contributions following this model are welcomed.


## Funding

`geoplanar` is partially supported by
[NSF Award #1759746, Comparative Regional Inequality Dynamics: Multiscalar and Multinational Perspectives](https://nsf.gov/awardsearch/showAward?AWD_ID=1759746&HistoricalAwards=false)
and by the Charles University’s Primus programme through the project _“Influence of
Socioeconomic and Cultural Factors on Urban Structure in Central Europe,”_ project
reference `PRIMUS/24/SCI/023`.


## How to cite

To cite GeoPlanar please use the following
[software paper](https://doi.org/10.1177/23998083261470518) published in the Environment
and Planning B: Urban Analytics and City Science.

> Rey, S., Fleischmann, M., & Winkler, L. (2026). GeoPlanar: Planar enforcement and
> coverage topology repairing for Python. _Environment and Planning B: Urban Analytics
> and City Science_, 23998083261470518. https://doi.org/10.1177/23998083261470518

BibTeX:

```bibtex
@article{rey2026GeoPlanar,
  title = {{{GeoPlanar}}: {{Planar}} Enforcement and Coverage Topology Repairing for {{Python}}},
  shorttitle = {{{GeoPlanar}}},
  author = {Rey, Sergio and Fleischmann, Martin and Winkler, Lisa},
  year = 2026,
  month = jul,
  journal = {Environment and Planning B: Urban Analytics and City Science},
  pages = {23998083261470518},
  issn = {2399-8083, 2399-8091},
  doi = {10.1177/23998083261470518},
}
```