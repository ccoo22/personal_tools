wget https://www.metaboanalyst.ca/resources/libs/compound_db.rds -O metaboanalyst.compound_db.rds
wget https://www.metaboanalyst.ca/resources/libs/syn_nms.rds -O metaboanalyst.syn_nms.rds

for species in hsa mmu rno bta gga dre dme cel sce osa ath smm pfa tbr eco bsu ppu sau tma syf mlo
do
 wget https://www.metaboanalyst.ca/resources/libs/kegg/2018/metpa/$species.rda -O metaboanalyst.kegg.$species.rda
done


