R -q --no-save <<END

library(openxlsx)
library(dplyr)

require(openxlsx)
xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
# Supplementary Table 4
ST4 <- read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=FALSE, cols=c(6:8,11:13,23:25), rows=6:1986)
names(ST4) <- c("SNP","Chr","Pos","A1","A2","EAF","b","se","p")
write.table(ST4, file="plasmaproteins", row.names=FALSE, col.names=FALSE, quote=FALSE)

END

