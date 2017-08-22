#!/usr/bin/env bash

source bioscripts/vcfsh.sh

###############################################################################
 # Define our data

females=(BjF Daf DmajF KHF LuF MuF RoF SLF)
males=(M S)
all=("${females[@]}" "${males[@]}")

###############################################################################

# Add sample labels for samples
for sample_name in ${all[@]}; do
  zcat $sample_name.variant.vcf.gz \
  | vcf_add_info SN $sample_name \
  > $sample_name.variant.labelled.vcf
done

###############################################################################
# REMOVE SNP clusters!

# Make a file with all SNPs
cat `for sample in ${all[@]}; do echo $sample.variant.labelled.vcf; done` \
 | vcf_merge_same \
 > all_variants.vcf

# Remove SNPs that are within 10bp of another SNP in any another sample (but not on exactly the same position)
for sample in ${all[@]}; do
  vcf_filter_nearby all_variants.vcf $sample.variant.labelled.vcf 10 \
  > $sample.variant.cleaned.vcf
done

# How many were removed?
cat `for sample in ${all[@]}; do echo $sample.variant.cleaned.vcf; done` \
 | vcf_merge_same \
 > all_variants_no_snp_clusters.vcf

###############################################################################
# FEMALES

# All merged variants
cat `for female in ${females[@]}; do echo $female.variant.cleaned.vcf; done` \
 | vcf_merge_same \
 > females.merged_variants.vcf

# All common variants
cat females.merged_variants.vcf \
 | vcf_filter_info_length SN 7 7 \
 | vcf_add_info SN female \
 > females.common_variants.vcf

###############################################################################
#MALES

# All merged variants
cat `for male in ${males[@]}; do echo $male.variant.cleaned.vcf; done` \
 | vcf_merge_same \
 > males.merged_variants.vcf

# All common variants
cat males.merged_variants.vcf \
 | vcf_filter_info_length SN 2 2 \
 | vcf_add_info SN male \
 > males.common_variants.vcf

###############################################################################
# ALL variants in three or more samples

cat `for sample in ${all[@]}; do echo $sample.variant.cleaned.vcf; done` \
 | vcf_merge_same \
 | vcf_filter_info_length SN 3 7 \
 > common_variants.vcf

###############################################################################
# Overlaps

# Present in ALL males and NO females
vcf_diff females.merged_variants.vcf males.common_variants.vcf \
 | vcf_filter_info_value SN male \
 > males.unique_variants.vcf

# Present in ALL females and NO males
vcf_diff males.merged_variants.vcf females.common_variants.vcf \
 > females.unique_variants.vcf

###############################################################################
###############################################################################
# Summary

echo "Total number of SNPs observed: "  `grep -v '^#' all_variants.vcf | wc -l`
echo "Total remaining without SNP clusters: " `grep -v '^#' all_variants_no_snp_clusters.vcf | wc -l`
echo "Common to only ALL males: " `grep -v '^#' males.unique_variants.vcf | wc -l`
echo "Common to only ALL females: " `grep -v '^#' females.unique_variants.vcf | wc -l`
echo "Present in at least 3: " `grep -v '^#' common_variants.vcf | wc -l`

