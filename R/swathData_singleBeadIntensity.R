singleBeadIntensity <- function(tiffFile, coordinates) {
    
    ## function to calculate intensity for a singe bead in an image
    ## probably not perfect (as usual with Illumia's intensity algorithm)
    ## might need more testing but seems to do a reasonable job
    
    tiffSection <- readTIFF(tiffFile, xlim = c(floor(coordinates[1]) - 10, floor(coordinates[1]) + 9), ylim = c(floor(coordinates[2]) - 10, floor(coordinates[2]) + 9) )
    xfrac <- coordinates[1] - floor(coordinates[1]);
    yfrac <- coordinates[2] - floor(coordinates[2]);
    bg <- illuminaBackground(tiffSection, c(9+xfrac, 9+yfrac));
    tiffSection <- illuminaSharpen(tiffSection);
    fg <- illuminaForeground(tiffSection, c(9+xfrac, 9+yfrac));
    
    return(fg - bg);
}


singleBeadIntensity_6x6 <- function(tiffFile, coordinates) {
    
    ## function to calculate intensity for a singe bead in an image
    ## probably not perfect (as usual with Illumia's intensity algorithm)
    ## might need more testing but seems to do a reasonable job
    
    tiffSection <- readTIFF(tiffFile, xlim = c(floor(coordinates[1]) - 10, floor(coordinates[1]) + 9), ylim = c(floor(coordinates[2]) - 10, floor(coordinates[2]) + 9) )
    xfrac <- coordinates[1] - floor(coordinates[1]);
    yfrac <- coordinates[2] - floor(coordinates[2]);
    bg <- illuminaBackground(tiffSection, c(9+xfrac, 9+yfrac));
    fg <- illuminaForeground_6x6(tiffSection, c(9+xfrac, 9+yfrac));
    
    return(fg - bg);
}
