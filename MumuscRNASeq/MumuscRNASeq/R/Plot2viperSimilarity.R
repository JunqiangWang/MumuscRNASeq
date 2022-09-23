




#' This function compares the viper similarity ofphenotype datasets
#' @param pa.dat.query The query dataset
#' @param pa.dat.ref The reference dataset
#' @param nn The param in function viperSimilarity
#' @param method The param in function viperSimilarity
#' @author Junqiang Wang
#' @export
Plot2viperSimilarity<-function(pa.dat.query, pa.dat.ref, nn=50, method=c("two.sided", "greater", "less")){

message("if the rownames are different, the rows with common rownames will be used!")

pa.common<-intersect(rownames(pa.dat.ref), rownames(pa.dat.query))

pa.dat.ref<-pa.dat.ref[match(pa.common, rownames(pa.dat.ref)), ]
pa.dat.query<-pa.dat.query[match(pa.common, rownames(pa.dat.query)), ]


pa.dat<-cbind(pa.dat.ref, pa.dat.query)


vp.similarity<-viperSimilarity(pa.dat, nn=nn, method=method)

vp.similarity<-vp.similarity[-1:-ncol(pa.dat.ref), 1:ncol(pa.dat.ref)]


#plot the data use red and blue

set.seed(77)

require(RColorBrewer)
require(ggplot2)
require(reshape2)


# convert to long format
df <- melt(vp.similarity)
summary(df)

df$viper_similarity<-df$value


p1 <- ggplot(df, aes(x = Var1, y = Var2, fill = viper_similarity)) +
  geom_raster() +
  theme_minimal(base_size = 16)

limit <- max(abs(df$viper_similarity)) * c(-1, 1)

p2<-p1 +
  scale_fill_distiller(palette = 'RdBu', limit = limit)

p2

}

















