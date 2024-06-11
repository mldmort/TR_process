#!/bin/bash

VCF_IN=../PB_ONT/longtr_merged_PB_ONT_addEND.vep.trimSamples.addInfo.vcf.gz
ZS_THR=3
OUT=zs_large_$ZS_THR.txt
N_SAM=$(bcftools query -l $VCF_IN | wc -l)
echo "number of samples: $N_SAM"

bcftools query -f '%ID[\t%ZS:%GB]' $VCF_IN | 
awk -v n_sam=$N_SAM \
	-v zs_thr=$ZS_THR \
'BEGIN{FS="\t";OFS="\t"}\
{\
	for (col=2;col<=NF;col++) {\
		split($col, p_col, ":");
		zs = p_col[1];\
		gb = p_col[2];\
		#if ($col != "." && $col != "nan|nan") {\
		if (zs != "." && zs != "nan|nan") {\
			#split($col, p, "|");\
			split(zs, p, "|");\
			split(gb, pp, "|");\
			for (i=1;i<=2;i++) {\
				if (p[i] >= zs_thr) {\
					print $1, p[i], pp[i];\
				}\
			}\
		}\
	}\
}' > $OUT
