### QFitsView2 (Ott 2012)

It has all the features (and possible more), but it is **pre-compiled**.

### PINGSOF (Rosales-Ortega 2011)

It is written in **IDL** requiring an expensive license

### CubeViz (Ferguson 2019)

It is a visualisation and analysis toolbox for data cubes. It is designed to work with data cubes from the NIRSpec and MIRI instruments on JWST: https://cubeviz.readthedocs.io/en/0.3.0/index.html#
**If it can work with any IFU, it is a problem** -> It can work with JWST and MaNGA

The available formats are:

|      Format      |Read |Write     | Auto-identify|
|-----------------:|----:|---------:|--------:|
|      6dFGS-split | Yes |   No     |      Yes|
|    6dFGS-tabular | Yes |   No     |      Yes|
|    APOGEE apStar | Yes |   No     |      Yes|
|   APOGEE apVisit | Yes |   No     |      Yes|
|APOGEE aspcapStar | Yes |   No     |      Yes|
|            ASCII | Yes |   No     |      Yes|
|             ECSV | Yes |   No     |      Yes|
|          HST/COS | Yes |   No     |      Yes|
|         HST/STIS | Yes |   No     |      Yes|
|             IPAC | Yes |   No     |      Yes|
|         JWST c1d | Yes |   No     |      Yes|
|         JWST s2d | Yes |   No     |      Yes|
|         JWST s3d | Yes |   No     |      Yes|
|         JWST x1d | Yes |   No     |      Yes|
|      MUSCLES SED | Yes |   No     |      Yes|
|       MaNGA cube | Yes |   No     |      Yes|
|        MaNGA rss | Yes |   No     |      Yes|
|     SDSS spPlate | Yes |   No     |      Yes|
| SDSS-I/II spSpec | Yes |   No     |      Yes|
| SDSS-III/IV spec | Yes |   No     |      Yes|
| Subaru-pfsObject | Yes |   No     |      Yes|
|             iraf | Yes |   No     |      Yes|
|     tabular-fits | Yes |  Yes     |      Yes|
|       wcs1d-fits | Yes |  Yes     |      Yes|

From the CubeViz webpage: *Some of Jdaviz’s dependencies require non-Python packages to work (particularly the front-end stack that is part of the Jupyter ecosystem).*

### MPDAF (Bacon et al. 2016)

The MUSE Python Data Analysis Framework. It has been developed and used in the MUSE Consortium for several years, and is available freely for the community: https://mpdaf.readthedocs.io/en/latest/
**It is not interactive**

### MARVIN (Cherinka et al. 2019, https://www.sdss.org/dr15/manga/marvin/)

It is a complete ecosystem designed for overcoming the challenge of searching, accessing, and visualizing the MaNGA data.  It is designed **only for MaNGA data**

### Mapviewer  from GIST (Bittner et al. 2019)

It doesn't allow to visualize only the input data. **You have to run GIST** and analyse kinematics/stellar populations



# Advantage of our program:

1. It is simple  
2. It can visualize only the input data
3. It is written in Python 
4. It is interactive 
5. It is not specific 



|              | QFitsView2 |   PINGSOF  |   CubeViz  |    MPDAF   |   MARVIN   |  Mapviewer |
|--------------|:----------:|-----------:|------------|:----------:|-----------:|-----------:|
| Simple       |            |            |            |            |            |            |
| Input data   |     X      |     X      |     X      |      X     |     X      |            |
| Python       |            |            |            |      X     |            |     X      |
| interactive  |     X      |     X      |     X      |            |     X      |     X      |
| not specific |     X      |     X      |            |            |            |     X      |



