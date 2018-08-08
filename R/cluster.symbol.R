cluster.symbol <- function(K)
{
symbol.matrix <- cbind(paste("pch=",seq(0:20),sep=""),
c("square","circle","triangle point up","plus","cross","diamond","triangle point down","square cross","star","diamond plus","circle plus","triangles up and down","square plus","circle cross","square and triangle down","filled square","filled circle","filled triangle point-up","filled diamond","solid circle","bullet"))
colnames(symbol.matrix) <- c("pch","symbol")
symbol.name <- symbol.matrix[1:K,2]
return(symbol.name)
}
