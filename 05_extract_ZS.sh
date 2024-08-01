#!/bin/bash

VCF_IN=../../PB_ONT/longtr_merged_PB_ONT_addEND.vep.trimSamples.addInfo.vcf.gz
ZS_THR=3
OUT=zs_large_$ZS_THR.txt
N_SAM=$(bcftools query -l $VCF_IN | wc -l)
echo "number of samples: $N_SAM"

bcftools query -f '%ID[\t%ZS:%GB]' $VCF_IN | 
awk -v n_sam=$N_SAM \
	-v zs_thr=$ZS_THR \
'BEGIN{FS="\t";OFS="\t"}\
{\
	if (NR==FNR) {\
		sample_dict[NR] = $1;\
	}\
	else {\
		for (col=2;col<=NF;col++) {\
			split($col, p_col, ":");\
			zs = p_col[1];\
			gb = p_col[2];\
			if (zs != "." && zs != "nan|nan") {\
				n = split(zs, p_z, "|");\
				split(gb, p_g, "|");\
				for (i=1;i<=n;i++) {\
					abs_zs = p_z[i]; if (abs_zs<0) {abs_zs=-1*abs_zs};\
					if (abs_zs >= zs_thr) {\
						print $1, p_z[i], p_g[i], sample_dict[col-1];\
					}\
				}\
			}\
		}\
	}\
}' <(bcftools query -l $VCF_IN) - > $OUT
