library(Seurat)
library(SeuratDisk)
library(SeuratObject)

setwd("/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/pancreas_split")

split1_seed317 <- LoadH5Seurat("./pancreas_seed317_split1_seurat.h5Seurat")
split2_seed317 <- LoadH5Seurat("./pancreas_seed317_split2_seurat.h5Seurat")
split1_seed320 <- LoadH5Seurat("./pancreas_seed320_split1_seurat.h5Seurat")
split2_seed320 <- LoadH5Seurat("./pancreas_seed320_split2_seurat.h5Seurat")

spliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="spliced"),
                        SeuratObject::LayerData(split2_seed317,assay="spliced"))
unspliced_seed317 <- list(SeuratObject::LayerData(split1_seed317,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed317,assay="unspliced"))
spliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="spliced"),
                        SeuratObject::LayerData(split2_seed320,assay="spliced"))
unspliced_seed320 <- list(SeuratObject::LayerData(split1_seed320,assay="unspliced"),
                          SeuratObject::LayerData(split2_seed320,assay="unspliced"))

### produce scatterplots of correlations
total = spliced_seed317[[1]] + spliced_seed317[[2]]
p <- nrow(total)
cor_spliced_seed317 <- sapply(1:nrow(spliced_seed317[[1]]),function(i){ 
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(spliced_seed317[[1]][i,]+1), log10(spliced_seed317[[2]][i,]+1) ) } )
# sum(is.na(cor_spliced_seed317)) = 3
nacor_s317 <- which(is.na(cor_spliced_seed317))

cor_unspliced_seed317 <- sapply(1:nrow(unspliced_seed317[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(unspliced_seed317[[1]][i,]+1), log10(unspliced_seed317[[2]][i,]+1) ) } )
# sum(is.na(cor_unspliced_seed317)) = 22
nacor_u317 <- which(is.na(cor_unspliced_seed317))

cor_spliced_seed320 <- sapply(1:nrow(spliced_seed320[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(spliced_seed320[[1]][i,]+1), log10(spliced_seed320[[2]][i,]+1) ) } )
# sum(is.na(cor_spliced_seed320)) = 2
nacor_s320 <- which(is.na(cor_spliced_seed320))

cor_unspliced_seed320 <- sapply(1:nrow(unspliced_seed320[[1]]),function(i){
  if(i %% floor(p/10) == 0) cat('*')
  cor( log10(unspliced_seed320[[1]][i,]+1), log10(unspliced_seed320[[2]][i,]+1) ) } )
# sum(is.na(cor_unspliced_seed320)) = 23
nacor_u320 <- which(is.na(cor_unspliced_seed320))


gene_ind_317 <- c(3, 5,    9,   10,   13,   17,   18,   20,   21,   24,   25,
         27,   31,   32,   35,   36,   38,   41,   42,   43,   48,   49,
         50,   55,   56,   57,   59,   70,   72,   75,   78,   79,   80,
         81,   83,   89,  101,  105,  107,  110,  112,  113,  114,  122,
        126,  129,  131,  133,  135,  136,  139,  140,  142,  144,  145,
        150,  152,  153,  154,  156,  158,  159,  160,  162,  163,  164,
        165,  168,  169,  170,  176,  177,  187,  188,  189,  191,  192,
        195,  197,  201,  206,  207,  215,  217,  218,  219,  220,  226,
        228,  232,  235,  239,  241,  242,  243,  244,  245,  251,  252,
        256,  262,  265,  266,  267,  268,  273,  280,  281,  282,  285,
        287,  289,  294,  297,  309,  310,  317,  318,  319,  321,  325,
        335,  338,  341,  346,  355,  356,  372,  378,  380,  381,  382,
        383,  389,  392,  394,  396,  400,  401,  405,  406,  412,  413,
        417,  425,  429,  432,  437,  441,  444,  446,  447,  455,  457,
        460,  462,  469,  474,  476,  478,  482,  484,  485,  487,  492,
        494,  496,  497,  498,  500,  502,  503,  507,  511,  512,  516,
        518,  526,  527,  530,  533,  534,  541,  542,  543,  560,  564,
        568,  569,  571,  572,  573,  574,  579,  582,  589,  590,  596,
        597,  600,  602,  605,  606,  613,  618,  621,  622,  623,  626,
        630,  632,  634,  637,  639,  641,  648,  649,  653,  655,  661,
        662,  664,  668,  671,  674,  675,  676,  678,  683,  690,  692,
        696,  697,  698,  701,  711,  715,  720,  723,  730,  732,  735,
        738,  740,  743,  747,  748,  750,  754,  756,  758,  759,  767,
        769,  781,  784,  794,  796,  802,  806,  808,  810,  813,  814,
        826,  827,  831,  832,  839,  841,  842,  856,  860,  864,  868,
        872,  883,  886,  891,  893,  897,  899,  901,  904,  906,  907,
        912,  914,  919,  923,  924,  927,  929,  930,  934,  936,  939,
        940,  944,  945,  947,  949,  950,  951,  953,  956,  957,  958,
        959,  965,  966,  971,  973,  974,  975,  977,  990,  993,  995,
        996,  997, 1000, 1001, 1006, 1010, 1015, 1016, 1022, 1028, 1030,
       1035, 1036, 1037, 1038, 1039, 1043, 1044, 1049, 1050, 1053, 1063,
       1065, 1067, 1068, 1072, 1073, 1074, 1076, 1079, 1081, 1086, 1088,
       1089, 1090, 1097, 1098, 1099, 1106, 1131, 1133, 1136, 1137, 1142,
       1145, 1146, 1153, 1156, 1157, 1163, 1165, 1167, 1168, 1170, 1171,
       1172, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1185, 1186, 1188,
       1191, 1195, 1196, 1198, 1200, 1202, 1203, 1207, 1211, 1215, 1216,
       1218, 1219, 1222, 1223, 1226, 1228, 1245, 1246, 1247, 1248, 1254,
       1255, 1260, 1263, 1265, 1268, 1272, 1273, 1274, 1275, 1277, 1279,
       1281, 1287, 1292, 1297, 1298, 1302, 1303, 1305, 1308, 1309, 1312,
       1313, 1322, 1323, 1324, 1327, 1332, 1336, 1337, 1338, 1339, 1346,
       1347, 1358, 1364, 1366, 1368, 1369, 1373, 1374, 1376, 1377, 1378,
       1383, 1384, 1385, 1386, 1387, 1389, 1391, 1395, 1397, 1407, 1408,
       1409, 1410, 1412, 1421, 1424, 1434, 1438, 1439, 1442, 1446, 1447,
       1448, 1453, 1457, 1461, 1464, 1465, 1467, 1469, 1471, 1475, 1481,
       1482, 1486, 1488, 1492, 1497, 1503, 1504, 1505, 1508, 1518, 1520,
       1523, 1530, 1535, 1537, 1541, 1556, 1557, 1560, 1564, 1565, 1572,
       1575, 1576, 1587, 1588, 1589, 1591, 1596, 1598, 1601, 1602, 1604,
       1608, 1611, 1613, 1615, 1618, 1619, 1622, 1625, 1627, 1634, 1637,
       1641, 1642, 1643, 1644, 1646, 1649, 1655, 1656, 1657, 1658, 1660,
       1662, 1665, 1666, 1668, 1670, 1673, 1678, 1679, 1682, 1687, 1688,
       1694, 1696, 1697, 1698, 1701, 1702, 1710, 1711, 1713, 1716, 1717,
       1719, 1721, 1726, 1727, 1728, 1733, 1734, 1736, 1738, 1741, 1744,
       1752, 1757, 1758, 1760, 1762, 1766, 1768, 1771, 1773, 1774, 1775,
       1781, 1790, 1791, 1794, 1799, 1801, 1804, 1806, 1807, 1808, 1811,
       1814, 1816, 1818, 1819, 1820, 1824, 1828, 1830, 1832, 1835, 1836,
       1838, 1839, 1841, 1849, 1850, 1853, 1859, 1860, 1861, 1878, 1879,
       1880, 1885, 1887, 1894, 1895, 1896, 1898, 1899, 1903, 1904, 1906,
       1907, 1913, 1916, 1918, 1921, 1922, 1924, 1926, 1932, 1934, 1939,
       1941, 1942, 1947, 1949, 1953, 1958, 1963, 1964, 1972, 1974, 1981,
       1982, 1985, 1987, 1991, 1992, 1993, 1994)

color_values_317 <- rep("lightblue",2000)
color_values_317[gene_ind_317] <- "firebrick"


png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/corr_seed317.png"), height = 1000, width = 2000, units = "px", res = 300)
par(mfrow=c(1,2))
plot(cor_spliced_seed317, ylim = c(-1,1), pch = 16, col = color_values_317, 
     main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
plot(cor_unspliced_seed317, ylim = c(-1,1), pch = 16, col = color_values_317, 
     main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
graphics.off()

mean(cor_spliced_seed317,na.rm=T) # 0.03962642
mean(cor_spliced_seed317[gene_ind_317]) # 0.03746595
mean(cor_spliced_seed317[setdiff(1:2000,gene_ind_317)],na.rm=T) # 0.0406833

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/corr_seed317_hist_velogenes.png"), height = 1000, width = 2000, units = "px", res = 300)
par(mfrow=c(1,2))
hist(cor_spliced_seed317[setdiff(1:2000,c(gene_ind_317,nacor_s317))],col="lightblue",
    main="corr of unselected genes(pan,scv)",cex.main=.8,xlab="correlation")
abline(v=mean(cor_spliced_seed317[setdiff(1:2000,gene_ind_317)],na.rm=T),col="coral",lwd=1.2)
text(.1,850,label=paste0("mean=",round(mean(cor_spliced_seed317[setdiff(1:2000,gene_ind_317)],na.rm=T),4)),col="coral",cex=.5)
hist(cor_spliced_seed317[gene_ind_317],col="firebrick",main="corr of selected genes(pan,scv)",cex.main=.8,xlab="correlation")
abline(v=mean(cor_spliced_seed317[gene_ind_317]),col="coral",lwd=1.2)
text(.1,130,label=paste0("mean=",round(mean(cor_spliced_seed317[gene_ind_317]),4)),col="coral",cex=.5)
graphics.off()

## multiple genes (spliced counts)
gene_names <- rownames(LayerData(split1_seed317))
## plot 24 genes in two spliced count splits with greatest correlations (log10)
intersect(order(abs(cor_spliced_seed317),decreasing=T)[1:25],which(is.na(cor_spliced_seed317))) # integer(0)
intersect(order(abs(cor_spliced_seed317))[1:25],which(is.na(cor_spliced_seed317))) # integer(0)
corr_ordermax_317 <- order(abs(cor_spliced_seed317),decreasing=T)
corr_max24_317 <- corr_ordermax_317[1:24]
for (i_start in as.integer(seq(1,21,length.out=6))) {
    ind = as.integer(seq(i_start,i_start+3,length.out=4)) # 4 genes to be plotted, decreasing corr
    x = list()
    y = list()
    n = length(spliced_seed317[[1]][1,])
    for (j in 1:4) {
        i_gene = corr_max24_317[ind[j]]
        x[[j]] <- log10(spliced_seed317[[1]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
        y[[j]] <- log10(spliced_seed317[[2]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
    }
    png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/spliced_seed317_cor_gene_max",i_start,".png"),
        height = 800, width = 2800, units = "px", res = 300)
    par(mfrow=c(1,4))
    for (j in 1:4) {
        i_gene = corr_max24_317[ind[j]]
        plot(x[[j]], y[[j]], pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
             main=paste0("gene",i_gene," (",gene_names[i_gene],"), cor=",round(cor_spliced_seed317[i_gene],2)),cex=.7,cex.main=.7)
    }
    graphics.off()
}
 
### min corr (log10)
corr_ordermin_317 <- order(abs(cor_spliced_seed317),decreasing=F)
corr_min24_317 <- corr_ordermin_317[1:24]
for (i_start in as.integer(seq(1,21,length.out=6))) {
    ind = as.integer(seq(i_start,i_start+3,length.out=4)) # 4 genes to be plotted, decreasing corr
    x = list()
    y = list()
    n = length(spliced_seed317[[1]][1,])
    for (j in 1:4) {
        i_gene = corr_min24_317[ind[j]]
        x[[j]] <- log10(spliced_seed317[[1]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
        y[[j]] <- log10(spliced_seed317[[2]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
    }
    png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/spliced_seed317_cor_gene_min",i_start,".png"),
        height = 800, width = 2800, units = "px", res = 300)
    par(mfrow=c(1,4))
    for (j in 1:4) {
        i_gene = corr_min24_317[ind[j]]
        plot(x[[j]], y[[j]], pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
             main=paste0("gene",i_gene," (",gene_names[i_gene],"), cor=",round(cor_spliced_seed317[i_gene],2)),cex=.7,cex.main=.7)
    }
    graphics.off()
}

############################

gene_ind_320 <- c(1,    3,    7,    9,   10,   13,   17,   18,   19,   20,   21,
         22,   24,   25,   27,   31,   32,   35,   36,   38,   41,   42,
         43,   48,   50,   55,   56,   57,   59,   69,   70,   72,   75,
         78,   79,   80,   81,   83,   88,  101,  105,  107,  108,  110,
        111,  112,  113,  114,  122,  126,  129,  131,  135,  136,  140,
        142,  144,  145,  150,  153,  154,  156,  158,  159,  160,  162,
        163,  164,  165,  168,  169,  170,  176,  187,  189,  191,  195,
        197,  201,  206,  207,  214,  215,  217,  218,  219,  220,  226,
        228,  235,  239,  241,  242,  243,  244,  245,  251,  252,  256,
        262,  265,  266,  267,  268,  273,  280,  281,  282,  285,  287,
        289,  294,  309,  310,  311,  317,  318,  319,  321,  325,  335,
        343,  346,  355,  356,  362,  372,  378,  380,  382,  383,  384,
        389,  392,  396,  400,  405,  406,  412,  413,  417,  425,  429,
        432,  437,  441,  444,  446,  447,  450,  455,  457,  460,  462,
        463,  469,  474,  475,  476,  478,  482,  484,  485,  487,  492,
        494,  496,  497,  498,  500,  502,  503,  507,  511,  518,  526,
        530,  533,  534,  541,  542,  543,  560,  564,  568,  569,  571,
        572,  573,  574,  579,  589,  590,  592,  596,  597,  598,  600,
        602,  605,  606,  611,  618,  621,  622,  623,  626,  630,  632,
        634,  637,  639,  641,  648,  649,  653,  655,  661,  664,  668,
        671,  674,  675,  676,  678,  683,  690,  692,  696,  697,  698,
        701,  703,  708,  711,  715,  723,  730,  732,  733,  735,  738,
        740,  743,  747,  748,  750,  754,  758,  759,  767,  769,  781,
        784,  793,  796,  799,  807,  808,  810,  813,  814,  826,  827,
        831,  832,  839,  840,  841,  842,  856,  860,  864,  868,  872,
        874,  879,  883,  886,  891,  893,  897,  899,  900,  901,  906,
        907,  912,  914,  919,  923,  924,  927,  929,  930,  934,  936,
        939,  940,  944,  945,  947,  949,  950,  951,  953,  957,  958,
        959,  965,  966,  971,  973,  974,  975,  977,  983,  990,  993,
        995,  996, 1000, 1001, 1015, 1016, 1022, 1028, 1030, 1035, 1036,
       1037, 1038, 1039, 1046, 1048, 1050, 1054, 1063, 1065, 1067, 1068,
       1069, 1072, 1073, 1074, 1076, 1079, 1081, 1086, 1089, 1090, 1097,
       1099, 1106, 1107, 1131, 1133, 1136, 1142, 1145, 1146, 1153, 1156,
       1157, 1160, 1163, 1165, 1167, 1168, 1170, 1171, 1172, 1177, 1178,
       1179, 1180, 1181, 1182, 1183, 1185, 1186, 1188, 1191, 1195, 1196,
       1198, 1200, 1203, 1206, 1207, 1211, 1214, 1215, 1216, 1218, 1219,
       1222, 1223, 1226, 1228, 1245, 1246, 1247, 1248, 1252, 1254, 1255,
       1260, 1265, 1267, 1268, 1272, 1273, 1274, 1275, 1279, 1281, 1282,
       1284, 1287, 1291, 1292, 1297, 1298, 1303, 1305, 1306, 1308, 1309,
       1312, 1313, 1321, 1322, 1324, 1327, 1332, 1336, 1337, 1339, 1346,
       1347, 1358, 1364, 1366, 1368, 1369, 1371, 1372, 1373, 1374, 1376,
       1377, 1378, 1383, 1384, 1385, 1386, 1387, 1389, 1391, 1395, 1397,
       1404, 1407, 1408, 1409, 1410, 1412, 1421, 1424, 1429, 1438, 1439,
       1442, 1446, 1447, 1448, 1449, 1451, 1453, 1457, 1461, 1464, 1465,
       1467, 1469, 1471, 1475, 1481, 1482, 1488, 1492, 1497, 1503, 1504,
       1505, 1508, 1512, 1518, 1520, 1523, 1530, 1535, 1537, 1540, 1541,
       1556, 1557, 1560, 1564, 1565, 1572, 1575, 1576, 1587, 1588, 1589,
       1591, 1596, 1598, 1601, 1602, 1604, 1608, 1611, 1613, 1615, 1617,
       1618, 1622, 1625, 1627, 1634, 1637, 1641, 1642, 1643, 1644, 1646,
       1649, 1655, 1656, 1657, 1660, 1661, 1662, 1665, 1666, 1668, 1670,
       1673, 1678, 1679, 1682, 1688, 1694, 1697, 1698, 1701, 1702, 1710,
       1713, 1716, 1717, 1719, 1721, 1726, 1727, 1728, 1729, 1733, 1734,
       1736, 1738, 1741, 1744, 1752, 1758, 1760, 1762, 1768, 1771, 1773,
       1774, 1775, 1781, 1790, 1791, 1794, 1797, 1799, 1801, 1804, 1806,
       1807, 1808, 1810, 1811, 1814, 1818, 1819, 1820, 1824, 1828, 1832,
       1835, 1836, 1838, 1839, 1841, 1843, 1849, 1850, 1859, 1860, 1861,
       1878, 1879, 1880, 1885, 1887, 1894, 1895, 1896, 1898, 1899, 1903,
       1906, 1910, 1913, 1916, 1918, 1921, 1922, 1926, 1932, 1934, 1939,
       1941, 1942, 1949, 1953, 1958, 1963, 1964, 1972, 1974, 1981, 1982,
       1985, 1987, 1991, 1992, 1993, 1994)
color_values_320 <- rep("lightblue",2000)
color_values_320[gene_ind_320] <- "firebrick"

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/corr_seed320.png"), height = 1000, width = 2000, units = "px", res = 300)
par(mfrow=c(1,2))
plot(cor_spliced_seed320, ylim = c(-1,1), pch = 16, col = color_values_320, 
     main="Spliced",ylab="corr",cex=.7,cex.main=.8) 
plot(cor_unspliced_seed320, ylim = c(-1,1), pch = 16, col = color_values_320, 
     main="Unspliced",ylab="corr",cex=.7,cex.main=.8)
graphics.off()

mean(cor_spliced_seed320,na.rm=T) # 0.03916057
mean(cor_spliced_seed320[gene_ind_320]) # 0.03800043
mean(cor_spliced_seed320[setdiff(1:2000,gene_ind_320)],na.rm=T) # 0.03972638

png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/corr_seed320_hist_velogenes.png"), height = 1000, width = 2000, units = "px", res = 300)
par(mfrow=c(1,2))
hist(cor_spliced_seed320[setdiff(1:2000,c(gene_ind_320,nacor_s320))],col="lightblue",
    main="corr of unselected genes(pan,scv)",cex.main=.8,xlab="correlation")
abline(v=mean(cor_spliced_seed320[setdiff(1:2000,gene_ind_320)],na.rm=T),col="coral",lwd=1.2)
text(.1,850,label=paste0("mean=",round(mean(cor_spliced_seed320[setdiff(1:2000,gene_ind_320)],na.rm=T),4)),col="coral",cex=.5)
hist(cor_spliced_seed320[gene_ind_320],col="firebrick",main="corr of selected genes(pan,scv)",cex.main=.8,xlab="correlation")
abline(v=mean(cor_spliced_seed320[gene_ind_320]),col="coral",lwd=1.2)
text(.1,130,label=paste0("mean=",round(mean(cor_spliced_seed320[gene_ind_320]),4)),col="coral",cex=.5)
graphics.off()

## 24 genes with greatest correlations
corr_ordermax_320 <- order(abs(cor_spliced_seed320),decreasing=T)
corr_max24_320 <- corr_ordermax_320[1:24]
for (i_start in as.integer(seq(1,21,length.out=6))) {
    ind = as.integer(seq(i_start,i_start+3,length.out=4)) # 4 genes to be plotted, decreasing corr
    x = list()
    y = list()
    n = length(spliced_seed320[[1]][1,])
    for (j in 1:4) {
        i_gene = corr_max24_320[ind[j]]
        x[[j]] <- log10(spliced_seed320[[1]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
        y[[j]] <- log10(spliced_seed320[[2]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
    }
    png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/spliced_seed320_cor_gene_max",i_start,".png"),
        height = 800, width = 2800, units = "px", res = 300)
    par(mfrow=c(1,4))
    for (j in 1:4) {
        i_gene = corr_max24_320[ind[j]]
        plot(x[[j]], y[[j]], pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
             main=paste0("gene",i_gene," (",gene_names[i_gene],"), cor=",round(cor_spliced_seed320[i_gene],2)),cex=.7,cex.main=.7)
    }
    graphics.off()
}

### min corr
corr_ordermin_320 <- order(abs(cor_spliced_seed320),decreasing=F)
corr_min24_320 <- corr_ordermin_320[1:24]
for (i_start in as.integer(seq(1,21,length.out=6))) {
    ind = as.integer(seq(i_start,i_start+3,length.out=4)) # 4 genes to be plotted, increasing corr
    x = list()
    y = list()
    n = length(spliced_seed320[[1]][1,])
    for (j in 1:4) {
        i_gene = corr_min24_320[ind[j]]
        x[[j]] <- log10(spliced_seed320[[1]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
        y[[j]] <- log10(spliced_seed320[[2]][i_gene,]+1) + runif(n, min = -0.1, max = 0.1)
    }
    png(file = paste0("/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/pancreas/scvelo/Jun11_cor/spliced_seed320_cor_gene_min",i_start,".png"),
        height = 800, width = 2800, units = "px", res = 300)
    par(mfrow=c(1,4))
    for (j in 1:4) {
        i_gene = corr_min24_320[ind[j]]
        plot(x[[j]], y[[j]], pch=16,asp=T, col = rgb(0.5,0.5,0.5,0.1),xlab="split1",ylab="split2",
             main=paste0("gene",i_gene," (",gene_names[i_gene],"), cor=",round(cor_spliced_seed320[i_gene],2)),cex=.7,cex.main=.7)
    }
    graphics.off()
}



