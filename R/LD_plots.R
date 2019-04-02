# correlation plots (LD triangles)

make_ld_plot <- function(X, meas) {

    stopifnot(meas %in% c("r", "D'"))
   
    require(LDheatmap)
    require(chopsticks)

    colnames(X)<- paste(1:ncol(X), "  ", sep="")
    gX <- as(X, "snp.matrix")

    cat("LD plot display:\n")
    ld <- LDheatmap(gX, flip=TRUE, name="", title=NULL, LDmeasure = meas,
                    add.map= T, geneMapLocation = 0.01, geneMapLabelX=1000)
}