share
=====

share estimates sources that are SHared Across a REgion.  

Install using `devtools':

```S
devtools::install_github("share", "kralljr")
```

This package depends on the handles package, which also must be installed via github:

```S
devtools::install_github("handles", "kralljr")
```

share includes the function "share" that estimates unknown, latent sources across several datasets (ambient monitors) and the function "sharehealth" that takes in observed chemical data and health data and estimates health effects using either APCA (Thurston and Spengler 1985) or mAPCA (Thurston et al. 2011).  

share also includes a function "mAPCA" to perform mAPCA.  
