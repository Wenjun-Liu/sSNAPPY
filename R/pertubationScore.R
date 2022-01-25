rel<-c("activation","compound","binding/association","expression","inhibition",
       "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
       "inhibition_dephosphorylation","dissociation","dephosphorylation",
       "activation_dephosphorylation","state change","activation_indirect effect",
       "inhibition_ubiquination","ubiquination", "expression_indirect effect",
       "inhibition_indirect effect","repression","dissociation_phosphorylation",
       "indirect effect_phosphorylation","activation_binding/association",
       "indirect effect","activation_compound","activation_ubiquination")

.weightedAdjMatrix <- function(list, beta, ){

    if(is.null(beta)){
        beta=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
        names(beta)<-rel
    }else{
        if(!all(names(beta) %in% rel) | length(names(beta))!=length(rel)){
            stop(paste("beta must be a numeric vector of length",length(rel), "with the following names:", "\n", paste(rel,collapse=",")))
        }
    }


    test <- sapply(names(datpT), function(x){
        g2gInteraction <- sapply(rel, function(y){
            datpT[[x]][[y]] * beta[y]

        }, simplify = FALSE)
        g2gInteraction <-   Reduce('+', g2gInteraction)

        numDownstream <- sapply(rel, function(y){
            datpT[[x]][[y]] * abs(sign(beta[y]))
        }, simplify = FALSE)
        numDownstream <- Reduce('+', numDownstream)
        numDownstream <- apply(numDownstream, 2, sum)
        numDownstream[numDownstream==0]<-1
        B <- t(t(g2gInteraction)/numDownstream)
        diag(B) <- diag(B)-1
        B

    }, simplify = FALSE)


}

perturbationScore <- function()
