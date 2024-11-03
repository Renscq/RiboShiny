rpf_CDT -l /mnt/t64/database/sce/fungidb/norm/sce.norm.txt \
 --rna /mnt/t64/rensc/serp/GSE67387/mrna/mrna/step05-merge/sce_rpf_merged.txt \
 --rpf /mnt/t64/rensc/serp/GSE67387/ribo/ribo/step05-merge/sce_rpf_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s A \
 --tis 10 \
 --tts 5 \
 -o A_ncr

rpf_CDT -l /mnt/t64/database/sce/fungidb/norm/sce.norm.txt \
 --rna /mnt/t64/rensc/serp/GSE67387/mrna/mrna/step05-merge/sce_rpf_merged.txt \
 --rpf /mnt/t64/rensc/serp/GSE67387/ribo/ribo/step05-merge/sce_rpf_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s P \
 --tis 10 \
 --tts 5 \
 -o P_ncr

rpf_CDT -l /mnt/t64/database/sce/fungidb/norm/sce.norm.txt \
 --rna /mnt/t64/rensc/serp/GSE67387/mrna/mrna/step05-merge/sce_rpf_merged.txt \
 --rpf /mnt/t64/rensc/serp/GSE67387/ribo/ribo/step05-merge/sce_rpf_merged.txt \
 --stop \
 -m 50 \
 -f 0 \
 -s E \
 --tis 10 \
 --tts 5 \
 -o E_ncr


