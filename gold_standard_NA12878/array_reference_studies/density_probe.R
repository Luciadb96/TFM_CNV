library(ggplot2)

#Load files
header_bed = c("chr", "start", "end", "cnv_type", "number_cp", "variant_id", 
               "variant_region_acc", "variant_region_id", "validation", "num_map", 
               "count")

#study id **nstd64**
nstd46 <- read.table(file = "{path}/nstd46_NA12878_CNV_TR.bed", 
                     header = FALSE, sep = "\t", col.names = header_bed)
summary(nstd46)

#study id **estd195**
estd195 <- read.table(file = "{path}/estd195_NA12878_CNV_TR.bed",
                      header = FALSE, sep = "\t", col.names = header_bed)
summary(estd195)

#study id **estd20**
estd20 <- read.table(file = "{path}t/estd20_NA12878_CNV_TR.bed",
                     header = FALSE, sep = "\t", col.names = header_bed)
summary(estd20)

###############################################################################
###############################################################################

#Histograma x = CNVs caracterizadas para NA12878 en el estudio nstd46
#           y = Número de sondas por CNV caracterizada (bedtools intersect -c)

# (1) Arrays global
ggplot(nstd46, aes(x = variant_id, y = count)) +
  geom_bar(stat="identity") + theme(axis.line = element_line(size = 0.5), 
                                    axis.ticks = element_line(linetype = "dotted"), 
                                    axis.title = element_text(size = 15), 
                                    axis.text.x = element_blank(), 
                                    plot.title = element_text(size = 20, face = "bold",
                                                              hjust = 0.5)) +
  labs(title = "Histograma sondas/CNV (ID: nstd46)", 
       x = "CNVs caracterizadas para NA12878 (496 CNVs)", 
       y = "Número de sondas por CNV caracterizada") 

###############################################################################

#Histograma x = Número de sondas por CNV caracterizada (bedtools intersect -c)
#           y = Numero de CNVs (Frecuencia: número de observaciones)

# (1) Arrays global
table(nstd46$count)

#ggplot
ggplot(nstd46, aes(x=count)) + geom_bar(color = "black") +
  labs(title = "Histograma sondas/CNV (ID: nstd46)",
       x = "Número de sondas por CNV", y = "Frecuencia (Nº de CNVs)") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  geom_vline(aes(xintercept=mean(count)), color="red") +
  geom_vline(aes(xintercept=median(count)), color= "darkgreen") +
  scale_y_continuous(limits=c(0, 20))

###############################################################################
###############################################################################

#Histograma x = CNVs caracterizadas para NA12878 en el estudio estd195
#           y = Número de sondas por CNV caracterizada (bedtools intersect -c)

# (1) Arrays global
ggplot(estd195, aes(x = variant_id, y = count)) +
  geom_bar(stat="identity") + theme(axis.line = element_line(size = 0.5), 
                                    axis.ticks = element_line(linetype = "dotted"), 
                                    axis.title = element_text(size = 15), 
                                    axis.text.x = element_blank(), 
                                    plot.title = element_text(size = 20, face = "bold",
                                                              hjust = 0.5)) +
  labs(title = "Histograma sondas/CNV (ID: estd195)", 
       x = "CNVs caracterizadas para NA12878 (138 CNVs)", 
       y = "Número de sondas por CNV caracterizada") 

###############################################################################

#Histograma x = Número de sondas por CNV caracterizada (bedtools intersect -c)
#           y = Numero de CNVs (Frecuencia: número de observaciones)

# (1) Arrays global
table(estd195$count)

#ggplot
ggplot(estd195, aes(x=count)) + geom_bar(color = "black") +
  labs(title = "Histograma sondas/CNV (ID: estd195)",
       x = "Número de sondas por CNV", y = "Frecuencia (Nº de CNVs)") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  geom_vline(aes(xintercept=mean(count)), color="red") +
  geom_vline(aes(xintercept=median(count)), color= "darkgreen") +
  scale_y_continuous(limits=c(0, 20))

###############################################################################
###############################################################################

#Histograma x = CNVs caracterizadas para NA12878 en el estudio estd20
#           y = Número de sondas por CNV caracterizada (bedtools intersect -c)

# (1) Arrays global
ggplot(estd20, aes(x = variant_id, y = count)) +
  geom_bar(stat="identity") + theme(axis.line = element_line(size = 0.5), 
                                    axis.ticks = element_line(linetype = "dotted"), 
                                    axis.title = element_text(size = 15), 
                                    axis.text.x = element_blank(), 
                                    plot.title = element_text(size = 20, face = "bold",
                                                              hjust = 0.5)) +
  labs(title = "Histograma sondas/CNV (ID: estd20)", 
       x = "CNVs caracterizadas para NA12878 (954 CNVs)", 
       y = "Número de sondas por CNV caracterizada") 

###############################################################################

#Histograma x = Número de sondas por CNV caracterizada (bedtools intersect -c)
#           y = Numero de CNVs (Frecuencia: número de observaciones)

# (1) Arrays global
table(estd20$count)

#ggplot
ggplot(estd20, aes(x=count)) + geom_bar(color = "black") +
  labs(title = "Histograma sondas/CNV (ID: estd20)",
       x = "Número de sondas por CNV", y = "Frecuencia (Nº de CNVs)") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  geom_vline(aes(xintercept=mean(count)), color="red") +
  geom_vline(aes(xintercept=median(count)), color= "darkgreen") +
  scale_y_continuous(limits=c(0, 20))

