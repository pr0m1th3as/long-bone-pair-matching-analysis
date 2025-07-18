## Load required packages
pkg load csg-dataset
pkg load csg-toolkit

## Aggregate CSG data for each bone
assemble_data

## Process each bone and compute pair-matching statistics
Femur_stats
Humerus_stats
Tibia_stats
Ulna_stats

## Create randomized assemblages as testing datasets
create_datasets
clear

## Evaluate sorting algorithm on randomized assemblages
evaluate_datasets
clear -x Femur Humerus Tibia Ulna

## Display summary results
Femur
Tibia
Humerus
Ulna
