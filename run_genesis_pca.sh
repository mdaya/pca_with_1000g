#!/bin/bash

#Set parameters
plink_bed_file=$1
plink_bim_file=`echo $plink_bed_file | sed 's/.bed/.bim/' | sed 's/.BED/.BIM/'`
plink_fam_file=`echo $plink_bed_file | sed 's/.bed/.fam/' | sed 's/.BED/.FAM/'`
gds_file=`echo $plink_bed_file | sed 's/.bed/.gds/' | sed 's/.BED/.gds/'`
pcvec_file=`echo $plink_bed_file | sed 's/.bed/.pcvec/' | sed 's/.BED/.pcvec/'`
pcvec_file=`echo $plink_bed_file | sed 's/.bed/.eigenval/' | sed 's/.BED/.eigenval/'`
slide_max_bp=$2
r2_ld_threshold=$3

#Call R script

cat /home/analyst/GENESIS_PCA.R | R --vanilla --args \
   $plink_bed_file \
   $plink_bim_file \
   $plink_fam_file \
   $gds_file \
   $pcvec_file \
   $pcvec_file \
   $slide_max_bp \
   $r2_ld_threshold

