#!/bin/bash
# QSUB commands to do geodiversity variable extraction
# 30 November 2017

# taxon = one of the following: bbs, fia
# geovar = one of the following: bioclim1k, bioclim5k, biocloud1k, biocloud5k, dhi, elevation, aspect, slope, tpi, hf, gea, night, soil

### FIA

# Elevation

qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 1-250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 251-500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 501-750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 751-1000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 1001-1250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 1251-1500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 1501-1750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 1751-2000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 2001-2250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 2251-2500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 2501-2750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 2751-3000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 3001-3250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 3251-3500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 3501-3750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 3751-4000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 4001-4250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 4251-4500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 4501-4750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 4751-5000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 5001-5250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 5251-5500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 5501-5750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 5751-6000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 6001-6250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 6251-6500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 6501-6750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 6751-7000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 7001-7250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 7251-7500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 7501-7750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 7751-8000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 8001-8250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 8251-8500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 8501-8750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 8751-9000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 9001-9250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 9251-9500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 9501-9750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 9751-10000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 10001-10250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 10251-10500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 10501-10750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 10751-11000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 11001-11250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 11251-11500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 11501-11750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 11751-12000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 12001-12250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 12251-12500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 12501-12750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 12751-13000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 13001-13250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 13251-13500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 13501-13750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 13751-14000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 14001-14250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 14251-14500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 14501-14750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 14751-15000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 15001-15250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 15251-15500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 15501-15750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 15751-16000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 16001-16250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 16251-16500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 16501-16750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 16751-17000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 17001-17250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 17251-17500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 17501-17750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 17751-18000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 18001-18250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 18251-18500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 18501-18750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 18751-19000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 19001-19250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 19251-19500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 19501-19750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 19751-20000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 20001-20250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 20251-20500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 20501-20750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 20751-21000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 21001-21250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 21251-21500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 21501-21750
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 21751-22000
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 22001-22250
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 22251-22500
qsub geoextract.sh -N elev_fia -v taxon=fia,geovar=elevation -t 22501-22531

# Bioclim 1k

qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 1-250
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 251-500
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 501-750
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 751-1000
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 1001-1250
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 1251-1500
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 1501-1750
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 1751-2000 
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 2001-2250 
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 2251-2500 
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 2501-2750 
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 2751-3000
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 3001-3250
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 3251-3500
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 3501-3750
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 3751-4000
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 4001-4250
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 4251-4500
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 4501-4750
qsub geoextract.sh -N clim1k_fia -v taxon=fia,geovar=bioclim1k -t 4751-5000

# Bioclim 5k

qsub geoextract.sh -N clim5k_fia -v taxon=fia,geovar=bioclim5k -t 1-250
qsub geoextract.sh -N clim5k_fia -v taxon=fia,geovar=bioclim5k -t 251-500
qsub geoextract.sh -N clim5k_fia -v taxon=fia,geovar=bioclim5k -t 501-750
qsub geoextract.sh -N clim5k_fia -v taxon=fia,geovar=bioclim5k -t 751-1000

# Biocloud 1k

qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 1-250
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 251-500
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 501-750
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 751-1000
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 1001-1250
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 1251-1500
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 1501-1750
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 1751-2000
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 2001-2250
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 2251-2500
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 2501-2750
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 2751-3000
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 3001-3250
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 3251-3500
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 3501-3750
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 3751-4000
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 4001-4250
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 4251-4500
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 4501-4750
qsub geoextract.sh -N cloud1k_fia -v taxon=fia,geovar=biocloud1k -t 4751-5000

# Biocloud 5k

qsub geoextract.sh -N cloud5k_fia -v taxon=fia,geovar=biocloud5k -t 1-250
qsub geoextract.sh -N cloud5k_fia -v taxon=fia,geovar=biocloud5k -t 251-500
qsub geoextract.sh -N cloud5k_fia -v taxon=fia,geovar=biocloud5k -t 501-750
qsub geoextract.sh -N cloud5k_fia -v taxon=fia,geovar=biocloud5k -t 751-1000

# DHI

qsub geoextract.sh -N dhi_fia -v taxon=fia,geovar=dhi -t 1-250

# Aspect

qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 1-250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 251-500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 501-750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 751-1000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 1001-1250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 1251-1500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 1501-1750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 1751-2000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 2001-2250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 2251-2500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 2501-2750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 2751-3000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 3001-3250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 3251-3500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 3501-3750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 3751-4000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 4001-4250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 4251-4500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 4501-4750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 4751-5000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 5001-5250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 5251-5500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 5501-5750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 5751-6000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 6001-6250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 6251-6500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 6501-6750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 6751-7000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 7001-7250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 7251-7500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 7501-7750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 7751-8000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 8001-8250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 8251-8500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 8501-8750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 8751-9000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 9001-9250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 9251-9500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 9501-9750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 9751-10000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 10001-10250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 10251-10500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 10501-10750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 10751-11000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 11001-11250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 11251-11500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 11501-11750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 11751-12000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 12001-12250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 12251-12500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 12501-12750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 12751-13000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 13001-13250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 13251-13500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 13501-13750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 13751-14000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 14001-14250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 14251-14500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 14501-14750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 14751-15000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 15001-15250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 15251-15500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 15501-15750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 15751-16000
#4h
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 16001-16250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 16251-16500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 16501-16750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 16751-17000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 17001-17250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 17251-17500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 17501-17750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 17751-18000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 18001-18250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 18251-18500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 18501-18750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 18751-19000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 19001-19250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 19251-19500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 19501-19750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 19751-20000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 20001-20250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 20251-20500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 20501-20750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 20751-21000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 21001-21250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 21251-21500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 21501-21750
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 21751-22000
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 22001-22250
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 22251-22500
qsub geoextract.sh -N aspect_fia -v taxon=fia,geovar=aspect -t 22501-22531

# Slope

qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 1-250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 251-500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 501-750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 751-1000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 1001-1250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 1251-1500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 1501-1750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 1751-2000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 2001-2250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 2251-2500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 2501-2750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 2751-3000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 3001-3250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 3251-3500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 3501-3750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 3751-4000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 4001-4250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 4251-4500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 4501-4750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 4751-5000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 5001-5250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 5251-5500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 5501-5750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 5751-6000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 6001-6250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 6251-6500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 6501-6750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 6751-7000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 7001-7250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 7251-7500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 7501-7750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 7751-8000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 8001-8250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 8251-8500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 8501-8750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 8751-9000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 9001-9250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 9251-9500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 9501-9750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 9751-10000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 10001-10250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 10251-10500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 10501-10750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 10751-11000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 11001-11250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 11251-11500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 11501-11750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 11751-12000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 12001-12250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 12251-12500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 12501-12750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 12751-13000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 13001-13250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 13251-13500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 13501-13750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 13751-14000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 14001-14250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 14251-14500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 14501-14750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 14751-15000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 15001-15250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 15251-15500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 15501-15750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 15751-16000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 16001-16250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 16251-16500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 16501-16750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 16751-17000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 17001-17250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 17251-17500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 17501-17750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 17751-18000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 18001-18250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 18251-18500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 18501-18750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 18751-19000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 19001-19250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 19251-19500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 19501-19750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 19751-20000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 20001-20250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 20251-20500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 20501-20750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 20751-21000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 21001-21250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 21251-21500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 21501-21750
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 21751-22000
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 22001-22250
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 22251-22500
qsub geoextract.sh -N slope_fia -v taxon=fia,geovar=slope -t 22501-22531

# TPI

qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 1-250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 251-500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 501-750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 751-1000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 1001-1250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 1251-1500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 1501-1750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 1751-2000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 2001-2250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 2251-2500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 2501-2750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 2751-3000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 3001-3250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 3251-3500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 3501-3750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 3751-4000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 4001-4250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 4251-4500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 4501-4750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 4751-5000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 5001-5250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 5251-5500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 5501-5750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 5751-6000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 6001-6250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 6251-6500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 6501-6750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 6751-7000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 7001-7250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 7251-7500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 7501-7750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 7751-8000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 8001-8250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 8251-8500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 8501-8750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 8751-9000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 9001-9250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 9251-9500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 9501-9750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 9751-10000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 10001-10250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 10251-10500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 10501-10750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 10751-11000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 11001-11250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 11251-11500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 11501-11750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 11751-12000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 12001-12250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 12251-12500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 12501-12750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 12751-13000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 13001-13250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 13251-13500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 13501-13750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 13751-14000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 14001-14250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 14251-14500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 14501-14750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 14751-15000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 15001-15250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 15251-15500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 15501-15750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 15751-16000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 16001-16250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 16251-16500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 16501-16750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 16751-17000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 17001-17250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 17251-17500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 17501-17750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 17751-18000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 18001-18250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 18251-18500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 18501-18750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 18751-19000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 19001-19250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 19251-19500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 19501-19750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 19751-20000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 20001-20250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 20251-20500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 20501-20750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 20751-21000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 21001-21250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 21251-21500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 21501-21750
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 21751-22000
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 22001-22250
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 22251-22500
qsub geoextract.sh -N tpi_fia -v taxon=fia,geovar=tpi -t 22501-22531

# Human footprint

qsub geoextract.sh -N dhi_fia -v taxon=fia,geovar=hf -t 1-250

# Geological age

qsub geoextract.sh -N gea_fia -v taxon=fia,geovar=gea -t 1-250

# Nightlights

qsub geoextract.sh -N night_fia -v taxon=fia,geovar=night -t 1-250

# Soil type

qsub geoextract.sh -N soil_fia -v taxon=fia,geovar=soil -t 1-250

# Check whether jobs are done

./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv bioclim1k_ .r 1 5000
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv bioclim5k_ .r 1 1000
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv biocloud1k_ .r 1 5000
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv biocloud5k_ .r 1 1000
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv dhi_ .r 1 250
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv elevation_ .r 1 22531
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv slope_ .r 1 22531
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv aspect_ .r 1 22531
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv tpi_ .r 1 22531
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv hf_ .r 1 250
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv gea_ .r 1 250
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv night_ .r 1 250
./didjob.sh /mnt/research/nasabio/data/fia/allgeodiv soil_ .r 1 250


###################################################

### BBS

# note: only a few bbs jobs were run. Everything will have to be restarted later.

# Bioclim 1k

qsub geoextract.sh -N clim1k_bbs -v taxon=bbs,geovar=bioclim1k -t 1-250
qsub geoextract.sh -N clim1k_bbs -v taxon=bbs,geovar=bioclim1k -t 251-500
qsub geoextract.sh -N clim1k_bbs -v taxon=bbs,geovar=bioclim1k -t 501-750
qsub geoextract.sh -N clim1k_bbs -v taxon=bbs,geovar=bioclim1k -t 751-1000

# Bioclim 5k

qsub geoextract.sh -N clim5k_bbs -v taxon=bbs,geovar=bioclim5k -t 1-250
qsub geoextract.sh -N clim5k_bbs -v taxon=bbs,geovar=bioclim5k -t 251-500

# Biocloud 1k

qsub geoextract.sh -N cloud1k_bbs -v taxon=bbs,geovar=biocloud1k -t 1-250
qsub geoextract.sh -N cloud1k_bbs -v taxon=bbs,geovar=biocloud1k -t 251-500
qsub geoextract.sh -N cloud1k_bbs -v taxon=bbs,geovar=biocloud1k -t 501-750
qsub geoextract.sh -N cloud1k_bbs -v taxon=bbs,geovar=biocloud1k -t 751-1000

# Biocloud 5k

qsub geoextract.sh -N cloud5k_bbs -v taxon=bbs,geovar=biocloud5k -t 1-250
qsub geoextract.sh -N cloud5k_bbs -v taxon=bbs,geovar=biocloud5k -t 251-500

# DHI

qsub geoextract.sh -N dhi_bbs -v taxon=bbs,geovar=dhi -t 1-100

# Elevation dem

qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 1-250
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 251-500
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 501-750
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 751-1000
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 1001-1250
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 1251-1500
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 1501-1750
qsub geoextract.sh -N elev_bbs -v taxon=bbs,geovar=elevation -t 1751-2000

# Aspect

qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 1-250
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 251-500
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 501-750
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 751-1000
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 1001-1250
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 1251-1500
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 1501-1750
qsub geoextract.sh -N aspect_bbs -v taxon=bbs,geovar=aspect -t 1751-2000

# Slope

qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 1-250
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 251-500
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 501-750
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 751-1000
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 1001-1250
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 1251-1500
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 1501-1750
qsub geoextract.sh -N slope_bbs -v taxon=bbs,geovar=slope -t 1751-2000

# TPI

qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 1-250
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 251-500
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 501-750
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 751-1000
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 1001-1250
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 1251-1500
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 1501-1750
qsub geoextract.sh -N tpi_bbs -v taxon=bbs,geovar=tpi -t 1751-2000

# Human footprint

qsub geoextract.sh -N dhi_bbs -v taxon=bbs,geovar=hf -t 1-100

# Geological age

qsub geoextract.sh -N gea_bbs -v taxon=bbs,geovar=gea -t 1-100

# Nightlights

qsub geoextract.sh -N night_bbs -v taxon=bbs,geovar=night -t 1-100

# Soil type

qsub geoextract.sh -N soil_bbs -v taxon=bbs,geovar=soil -t 1-100

# Check whether jobs are done

./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv bioclim1k_ .r 1 1000
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv bioclim5k_ .r 1 500
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv biocloud1k_ .r 1 1000
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv biocloud5k_ .r 1 500
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv dhi_ .r 1 100
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv elevation_ .r 1 2000
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv slope_ .r 1 2000
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv aspect_ .r 1 2000
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv tpi_ .r 1 2000
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv hf_ .r 1 100
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv gea_ .r 1 100
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv night_ .r 1 100
./didjob.sh /mnt/research/nasabio/data/bbs/allgeodiv soil_ .r 1 100

################################################
# FIA old-school beta-diversity
qsub fiabd.sh -t 3-250
qsub fiabd.sh -t 251-500
qsub fiabd.sh -t 501-750
qsub fiabd.sh -t 751-1000
qsub fiabd.sh -t 1001-1250
qsub fiabd.sh -t 1251-1500
qsub fiabd.sh -t 1501-1750
qsub fiabd.sh -t 1751-2000
qsub fiabd.sh -t 2001-2250
qsub fiabd.sh -t 2251-2500
qsub fiabd.sh -t 2501-2750
qsub fiabd.sh -t 2751-3000
qsub fiabd.sh -t 3001-3250
qsub fiabd.sh -t 3251-3500
qsub fiabd.sh -t 3501-3750
qsub fiabd.sh -t 3751-4000
qsub fiabd.sh -t 4001-4250
qsub fiabd.sh -t 4251-4500
qsub fiabd.sh -t 4501-4750
qsub fiabd.sh -t 4751-5000
qsub fiabd.sh -t 5001-5250
qsub fiabd.sh -t 5251-5500
qsub fiabd.sh -t 5501-5750
qsub fiabd.sh -t 5751-6000
qsub fiabd.sh -t 6001-6250
qsub fiabd.sh -t 6251-6500
qsub fiabd.sh -t 6501-6750
qsub fiabd.sh -t 6751-7000
qsub fiabd.sh -t 7001-7250
qsub fiabd.sh -t 7251-7500
qsub fiabd.sh -t 7501-7750
qsub fiabd.sh -t 7751-8000
qsub fiabd.sh -t 8001-8250
qsub fiabd.sh -t 8251-8500
qsub fiabd.sh -t 8501-8750
qsub fiabd.sh -t 8751-9000
qsub fiabd.sh -t 9001-9250
qsub fiabd.sh -t 9251-9500
qsub fiabd.sh -t 9501-9750
qsub fiabd.sh -t 9751-10000

./didjob.sh /mnt/research/nasabio/data/fia/diversity/unfuzzed beta_ .r 1 10000