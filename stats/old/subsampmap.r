focalpoly_fia <- make_map_polygons(states = c('California','Oregon','Washington'))
focalpoly_fia <- gBuffer(focalpoly_fia, width = 0)
focalpoly_fia <- spTransform(focalpoly_fia, CRSobj = CRS(aea_crs))

sample_idx <- SRS_iterative(focal_points = fia_aea_noedge, dist_mat = fia_pnw_dist, radius = 10 * 1000, n = 2000, show_progress = TRUE)

plot(focalpoly_fia, main = 'FIA Pacific Northwest')
points(fia_aea_noedge[sample_idx], pch=19, cex=0.5, col='blue')
