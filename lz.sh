# 18-10-2018 JHZ

# extraction of SNP positions from LocusZoom 1.4 databases

export hg19=/rds-d4/user/jhz22/hpc-work/locuszoom_1.4/data/database/locuszoom_hg19.db
sqlite3 $hg19 <<END
.tables
.indices
.separator "\t"
.output snp_pos.tsv
select * from snp_pos;
END
gzip -f snp_pos.tsv
