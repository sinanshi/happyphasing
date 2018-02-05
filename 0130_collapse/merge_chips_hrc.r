library(data.table)
library(stringr)
if (!exists("h"))
    load("h_884983_64940.RData")

if (!exists("std_"))
    std_ <- fread("humanomni1-quad_v1-0_h-b37.strand")

if (!exists("leg"))
    leg <- fread("legend")

std <- std_[V2 == 20]

dic <- list()
dic[["A"]] <- "T"
dic[["T"]] <- "A"
dic[["C"]] <- "G"
dic[["G"]] <- "C"

splt <- apply(std, 1, function(x){
                     s <- str_split_fixed(x[["V6"]], "", 2)
                     if (x[["V5"]] == "-") s <- c(dic[[s[1]]], dic[[s[2]]])
                     s
})


std <- data.table(std, ref=splt[1, ], minor=splt[2, ])
strand <- std[,c("V3", "ref", "minor")]
names(strand) <- c("position", "a0", "a1")

aa = data.table(strand$a0, strand$a1)
aa1 <- apply(aa, 1, sort)
aa = data.table(leg$a0, leg$a1)
aa2 <- apply(aa, 1, sort)


strand$a0 <- aa1[1, ]
strand$a1 <- aa1[2, ]
leg$a0 <- aa2[1, ]
leg$a1 <- aa2[2, ]

m <- merge(strand, leg, by=c("position", "a0", "a1"))
index <- leg$position %in% m$position & !duplicated(leg$position)

h_sub <- h[index, ]
row.names(h_sub) <- leg$position[index]
