# Script for 2nd round of GEB revision
# 1. Find the number of neighbors used to calculate 50-km beta-diversity for each BBS route (based on whether midpoints are in the 50-km radius)
# 2. Find the percentage of stops within the 50-km radius for each of those neighbors.
# QDR/NASABIOXGEO/11 Sep 2019

# 1. get numbers of neighbors for BBS
# -----------------------------------

# Load bbs route coordinates
load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')
bbscov_oneyear_byroute <- bbscov_oneyear

library(dplyr)
library(sp)
library(rgdal)
library(purrr)


# for each BBS route, find the IDs of the routes within 50 km.
get_neighbors <- function(x, r = 50) {
	focalpointindex <- which(bbscov_oneyear_byroute$rteNo == x$rteNo)
	neighbordists <- spDistsN1(pts = cbind(bbscov_oneyear_byroute$lon, bbscov_oneyear_byroute$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	neighbor_ids <- bbscov_oneyear_byroute$rteNo[neighbordists > 0 & neighbordists <= r]
	neighbor_ids
}


bbs_neighbor_list <- bbscov_oneyear_byroute %>%
	group_by(rteNo) %>%
	do(nhb = get_neighbors(.))
	
# Get numbers of neighbors
n_neighbs <- map_int(bbs_neighbor_list$nhb, length)
table(n_neighbs)

# Result pasted in:
#   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
# 237 581 675 526 398 287 144 101  61  36  23   8   3   6   2   1

# We had to exclude 237 routes for not having beta-diversity neighbors, 581 routes have only 1 neighbor, but the remaining 2271 routes have 2 or more neighbors up to a max of 15.

# 2. find percentage of routes enclosed in circles for BBS
# --------------------------------------------------------

# For each of the neighbors identified, where we have mapped stop coordinates, find the distance of each of the STOPS from the midpoint of the single focal route.

# Load the stop coordinates.
stopcoords <- read.csv('/mnt/home/qdr/bbs_stop_coords.csv')

# For each route, find the coordinate of its midpoint, the coordinates of all the STOPS only in neighbor routes who have midpoints at least 100 km from the route. 
# That will account for the ones that are "sticking into" the circle even if their midpoint is not <= 50km.

# Get the 100 km neighbors
bbs_neighbor_list_100 <- bbscov_oneyear_byroute %>%
	rowwise %>%
	do(nhb = get_neighbors(., r = 100))

# Define function to get distances for all the neighboring stops
get_stop_distances <- function(x) {
	focalpointindex <- which(bbscov_oneyear_byroute$rteNo == x$rteNo)
	neighborstops <- stopcoords[stopcoords$rteNo %in% bbs_neighbor_list_100$nhb[[focalpointindex]], ] # Find the coordinates of all stops of routes with midpoints within 100 km
	stopdists <- spDistsN1(pts = cbind(neighborstops$lon, neighborstops$lat), pt = c(x$lon, x$lat), longlat = TRUE) # Find the distances to all those stops from focal route midpoint
	
	# Tally the neighbor stops by route and whether they are within 50 km
	stoptable <- neighborstops %>% 
		mutate(dist = stopdists) %>%
		group_by(rteNo) %>%
		summarize(n_within = sum(dist <= 50)) %>%
		rename(neighbor_rteNo = rteNo)
	
	stoptable
}

bbs_neighbor_stop_distances <- bbscov_oneyear_byroute %>%
	group_by(rteNo) %>%
	do(get_stop_distances(.)) %>%
	filter(n_within > 0)

# Convert the 50 km list into a data frame so we can join it with the stop distances data frame
# We can use this to tell us how many of the incomplete routes are "sticking in" to the circle and how many are not.

bbs_neighbor_df_50 <- data.frame(rteNo = rep(bbs_neighbor_list$rteNo, times = map_int(bbs_neighbor_list$nhb, length)), nhb = unlist(bbs_neighbor_list$nhb)) 

# Join the two data frames to identify which are within 50.
bbs_neighbor_stop_distances <- bbs_neighbor_stop_distances %>%
	group_by_all %>%
	do(data.frame(midpoint_within_50 = .$neighbor_rteNo %in% bbs_neighbor_df_50$nhb[bbs_neighbor_df_50$rteNo == .$rteNo]))
	
# Let's tabulate the stops
# For each route, we want to know (1) number of stops in routes where all 50 are in the 50 km radius, (2) number of stops in routes with the midpoint in and not all 50 in the radius, and (3) number of stops in routes with the midpoint out but at least some in the radius

stop_tallies <- bbs_neighbor_stop_distances %>%
	mutate(completely_within = n_within == 50,
		   n_outside = 50 - n_within) %>%
	group_by(rteNo, completely_within, midpoint_within_50) %>%
	summarize(n_unique_routes = length(unique(neighbor_rteNo)),
			  total_n_stops_within = sum(n_within),
			  total_n_stops_excluded = sum(n_outside))
			  
# Name the 3 possible cases
stop_tallies <- stop_tallies %>%
	ungroup %>%
	mutate(status = case_when(
		completely_within ~ 'entire neighbor route inside circle',
		!completely_within & midpoint_within_50 ~ 'midpoint inside circle, some stops outside circle',
		!midpoint_within_50 ~ 'midpoint outside circle, some stops inside circle'))
		
# Get grand totals by status
stop_grandtotals <- stop_tallies %>%
	group_by(status) %>%
	summarize(total_n_stops_within = sum(total_n_stops_within),
			  total_n_stops_excluded = sum(total_n_stops_excluded))
write.table(stop_grandtotals, row.names = FALSE, quote = FALSE, sep = '\t')			  
write.csv(stop_grandtotals, '~/stop_grandtotals.csv', row.names = FALSE)
# Result pasted in
# status total_n_stops_within total_n_stops_excluded
# entire neighbor route inside circle 258650 0
# midpoint inside circle, some stops outside circle 126412 48788
# midpoint outside circle, some stops inside circle 48882 140668


# Try another criterion, say that each route has to be at least 50% within the circle
twocriteria <- bbs_neighbor_stop_distances %>%
	ungroup %>%
	mutate(half_within_50 = n_within >= 25) %>%
	select(midpoint_within_50, half_within_50) %>%
	filter(midpoint_within_50 | half_within_50)
	
addmargins(table(twocriteria))

# Result pasted in
#                  half_within_50
#midpoint_within_50 FALSE TRUE  Sum
#             FALSE     0  181  181
#             TRUE    157 8520 8677
#             Sum     157 8701 8858

# Use a more stringent criterion of 75%

twocriteria_threequarters <- bbs_neighbor_stop_distances %>%
	ungroup %>%
	mutate(threequarters_within_50 = n_within > 37) %>%
	select(midpoint_within_50, threequarters_within_50) %>%
	filter(midpoint_within_50 | threequarters_within_50)
	
addmargins(table(twocriteria_threequarters))
# Result pasted in
#                   threequarters_within_50
#midpoint_within_50 FALSE TRUE  Sum
#             FALSE     0   21   21
#             TRUE   1948 6729 8677
#             Sum    1948 6750 8698
