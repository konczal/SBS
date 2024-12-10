Fst = read.table("SBS_RNS_fst.25k.windowed.weir.fst", sep="\t", header=T)
Fst$ZFst = scale(Fst$WEIGHTED_FST, center = TRUE, scale = TRUE)
m<-mean(Fst$ZFst)
std<-sqrt(var(Fst$ZFst))
hist(Fst$ZFst, nclass = 500, xlim = c(-4,4), prob=TRUE, xlab = "ZFst", main="")
curve(dnorm(x, mean=0, sd=1), 
      col="red", lwd=2, add=TRUE, yaxt="n")
zfst2 <- Fst[Fst$ZFst > (2.5),]
write.table(zfst2, file = "SBS_RNS_zfst25.txt", sep="\t", row.names = FALSE, quote = FALSE)

#zfst3 <- Fst[Fst$ZFst > (3), ]
#write.table(zfst3, file = "SBS_RNS_zfst30.txt", sep="\t", row.names = FALSE, quote = FALSE)
