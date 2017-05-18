# Approximate circle in lat long (returns coordinates of n points)
# Code from http://mathforum.org/library/drmath/view/51816.html

# Take longitude and latitude in degrees, distance d in kilometers.
# Must convert lon/lat to radians, and distance to arc in radians (divide by Earth's radius in km)
# At the end convert lon/lat back to degrees and return a matrix
# Also had to switch the signs because west longitudes were treated as positive in the original.

latLongCircle <- function(lon1, lat1, d, n = 1000) {
  
  lat1 <- lat1 * pi/180
  lon1 <- lon1 * pi/180
  d <- d/6371
  
  thetas <- seq(0, 2*pi, length.out = n)
  lat <- asin(sin(lat1)*cos(d)+cos(lat1)*sin(d)*cos(thetas))
  dlon <- atan2(sin(thetas)*sin(d)*cos(lat1), cos(d)-sin(lat1)*sin(lat))
  lon = ((lon1+dlon+pi) %% (2*pi)) - pi
  
  cbind(lon, lat) * 180/pi
  
}
