pdf("rs1426654_freq.pdf", width=8, height=8)
 
plot(freq$Allele_A, freq$lat, xlab="f(A) rs1426654", ylab="Latitude", 
     pch=16, cex=0.8, col=myColors[freq$superpop], xlim=c(0, 1))

legend('topleft', c('African', 'Admixed American', 'East Asia', 'European', 'South Asian'), 
     cex=0.8, col=c('red','blue','darkgreen','salmon','black'), pch=16, inset=0.02)

dev.off()


### World maps begin
pdf("WorldMap_rs1426654.pdf", width=12, height = 6)
map("worldHires", xlim=c(-120, 120), ylim=c(-12,72), col="gray", fill=FALSE)
#points(freq$long, freq$lat, pch=16, col="green")

for(m in 1:26)
{
  add.pie(z=c(freq$Allele_A[m], freq$Allele_G[m]), x=freq$long[m], y=freq$lat[m], radius=freq$N_CHR[m]/100, 
        col=c(alpha("orange", 0.6), alpha("purple", 0.6)), labels="")
  m=m+1
}

box()

dev.off()

