#Manuscript Title: Relatively large wings facilitate life at higher elevations in Nearctic dragonflies

#Year: 2023

# https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13946

#File Name: elevation.dat.drags.csv
#Description: Species’ traits of Nearctic Libelluloidea dragonflies that were used for the analyses to test the relationship between a species’ relative wing size and features of its elevational distribution. Species in this analysis are those that are present in a pruned, maximum credibility phylogeny that was produced by Rocha-Ortega et al. 2020 (Proceedings of the Royal Society of London B). The pruned version of the phylogeny used in the analyses for these data has the file name ‘phylo.drags.elev.tre’. The R scripts that use this dataset are named “code for dragonfly elevational distribution.R”. 
#Rows: 302, excluding the header. Each row corresponds to a species 
#Columns: 17. binom, body.size rel.wing, resid.wing, range.size, bio1_mean, bio6_mean, range.elev, max.mean.elev, min.mean.elev, max.elev, min.elev, max.min.elev, min.max.elev, max.mean.elev.finer, min.mean.elev.finer, range.elev.finer

#binom: Name of the dragonfly species

#body.size: adult body length. Taken as the mid-point of the total length (mm) listed in Paulson’s comprehensive field guides of North American odonates (2009, 2012) 

#rel.wing: relative wing size. Calculated by dividing the mid-point of the ln-transformed wing length (mm) by the ln-transformed total length (mm) listed in Paulson’s comprehensive field guides of North American odonates (2009, 2012) 

#resid.wing: residual wing size. Calculated as the residuals of a regression of ln-transformed wing size on ln-transformed body length. 

#range.size: size of a species’ geographic range (km2). Calculated from the occupied 1º x 1º grid cells from a presence/absence matrix across the Nearctic

#bio1_mean: Average Mean Annual Temperature (bioclim 1) from across the species’ range. 

#bio6_mean: Average Minimum Temperature of Coldest month (Bioclim 6) from across the species’ range

#range.elev: difference between the mean elevation of a species’ highest occupied 1º x 1º grid cell and the mean elevation of a species lowest occupied 1º x 1º grid cell (meters). 

#max.mean.elev: the mean elevation of a species’ highest occupied 1º x 1º grid cell (meters above sea level)

#min.mean.elev: the mean elevation of a species lowest occupied 1º x 1º grid cell (meters above sea level)

#max.elev: the highest elevation of a species’ highest occupied 1º x 1º grid cell
min.elev: the lowest elevation of a species’ lowest occupied 1º x 1º grid cell (meters above sea level)


#max.min.elev: the minimum elevation of a species’ highest occupied 1º x 1º grid cell (meters above sea level)


#min.max.elev: the maximum elevation of a species’ lowest occupied 1º x 1º grid cell (meters above sea level)


#max.mean.elev.finer: the mean elevation of a species highest occupied 0.25º x 0.25º grid cell (meters above sea level)


#min.mean.elev.finer: the mean elevation of a species lowest occupied 0.25º x 0.25º grid cell (meters above sea level)


#range.elev.finer: the difference between the mean elevation of a species’ highest occupied grid cell and the mean elevation of a species lowest occupied 0.25º x 0.25º grid cell (meters)
