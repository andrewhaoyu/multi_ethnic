TestTable <-
  matrix(c(5, 0, 1365, 1635),
         nrow = 2,
         dimnames = list(Case = c("Cases", "Control"),
                         Carrier = c("Carrier", "Not Carrier")))
fisher.test(TestTable)


TestTable <-
  matrix(c(2, 0, 305, 1635),
         nrow = 2,
         dimnames = list(Case = c("Cases", "Control"),
                         Carrier = c("Carrier", "Not Carrier")))
fisher.test(TestTable)
library(data.table)
pheno <- fread("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/EUR/pheno_summary_out_GA/phenotypes_rho1_1.phen")
