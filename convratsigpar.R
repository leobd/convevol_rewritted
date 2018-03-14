#########################################################################################################
##           Function from Stayton's convevol R package (2014),                                        ##
##           modified for use in Botton-Divet et al (2016)                                             ##
#########################################################################################################

convratsigpar<-function (phyl, phendata, convtips, nsim) 
{
    if (class(phyl) != "phylo") 
        stop("your tree must be class 'phylo.'")
    if (nrow(phendata) != length(phyl$tip)) 
        stop("your data matrix must have the same number of rows as tips in the tree.")
    if (length(convtips)==2){obs <- convratpar(phyl, phendata, convtips)}
    if (length(convtips)>2){obs <- convratpar(phyl, phendata, convtips)$mean}
    
    C <- vcv.phylo(phyl)
    vcv <- phyl.vcv(phendata, C, 0)
     simphendata <- sim.char(phyl, vcv$R, nsim, model = c("BM"),root = 0)
simresults<-matrix(data =unlist(foreach(i=1:nsim) %dopar% {convratpar(phyl, simphendata[, , i], convtips)$mean}),nrow = nsim, ncol = 4, byrow = TRUE)
    C1greater <- 0
    C2greater <- 0
    C3greater <- 0
    C4greater <- 0
    for (i in 1:nsim) {
        if (simresults[i,1] >=  obs[1]) {
            C1greater <- C1greater + 1
        }
        if (simresults[i,2] >= obs[2]) {
            C2greater <- C2greater + 1
        }
        if (simresults[i,3] >= obs[3]) {
            C3greater <- C3greater + 1
        }
        if (simresults[i,4] >= obs[4]) {
            C4greater <- C4greater + 1
        }
    }
    answer <- c(C1greater/(nsim + 1), C2greater/(nsim + 1), C3greater/(nsim + 
        1), C4greater/(nsim + 1))
    names(answer) <- c("C1", "C2", "C3", "C4")
    C1cutoff <- sort(simresults[,1])[round(nsim * 0.95)]
    C2cutoff <- sort(simresults[,2])[round(nsim * 0.95)]
    C3cutoff <- sort(simresults[,3])[round(nsim * 0.95)]
    C4cutoff <- sort(simresults[,4])[round(nsim * 0.95)]
    Cutoffs <- c(C1cutoff, C2cutoff, C3cutoff, C4cutoff)
    output<-list(Cobs=obs,frequ=answer,Csimul=simresults,Cutoffs=Cutoffs)
    output
}
