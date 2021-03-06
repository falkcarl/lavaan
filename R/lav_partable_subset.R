# YR 11 feb 2017: initial version

# given a parameter table (PT), extract a part of the model:
# eg.:
# - only the measurement model (with saturated latent variables
# - only the stuctural part
# - a single measurement model (1 factor only)
# ...


# FIXME: if we have more than 1 factor, we remove the structural
# part, but should we add ALL correlations among the latent variables?
# (NO for now)
lav_partable_subset_measurement_model <- function(PT = NULL,
                                                  lavpta = NULL,
                                                  lv.names = NULL,
                                                  idx.only = FALSE) {

    # PT
    PT <- as.data.frame(PT, stringsAsFactors = FALSE)

    # lavpta
    if(is.null(lavpta)) {
        lavpta <- lav_partable_attributes(PT)
    }

    # ngroups
    ngroups <- lavpta$ngroups

    # lv.names: list with element per group
    if(is.null(lv.names)) {
        lv.names <- lavpta$vnames$lv.regular
    } else if(!is.list(lv.names)) {
        lv.names <- list(lv.names)
    }
    
    # keep rows idx
    keep.idx <- integer(0L)

    # remove not-needed measurement models
    for(g in 1:ngroups) {
        # indicators for latent variables we keep
        IND.idx <- which(  PT$op == "=~"              &
                           PT$lhs %in% lv.names[[g]]  &
                           PT$group == g )
        IND <- PT$rhs[ IND.idx ]

        # keep =~
        keep.idx <- c(keep.idx, IND.idx)

        # keep ~~
        OV.VAR.idx <- which( PT$op == "~~"   &
                             PT$lhs %in% IND &
                             PT$rhs %in% IND &
                             PT$group == g )
        keep.idx <- c(keep.idx, OV.VAR.idx)

        LV.VAR.idx <- which( PT$op == "~~"             &
                             PT$lhs %in% lv.names[[g]] &
                             PT$rhs %in% lv.names[[g]] &
                             PT$group == g )
        keep.idx <- c(keep.idx, LV.VAR.idx)

        # intercepts indicators
        OV.INT.idx <- which( PT$op == "~1"    &
                             PT$lhs %in% IND  &
                             PT$group == g )
        keep.idx <- c(keep.idx, OV.INT.idx)

        # intercepts latent variables
        LV.INT.idx <- which( PT$op == "~1"                 &
                             PT$lhs %in% lv.names[[g]]  &
                             PT$group == g )
        keep.idx <- c(keep.idx, LV.INT.idx)

        # thresholds
        TH.idx <- which( PT$op == "|"    &
                         PT$lhs %in% IND &
                         PT$group == g )
        keep.idx <- c(keep.idx, TH.idx)

        # scaling factors
        SC.idx <- which( PT$op == "~*~"  &
                         PT$lhs %in% IND &
                         PT$group == g )
        keep.idx <- c(keep.idx, SC.idx)

        # FIXME: ==, :=, <, >, == involving IND...
    }

    if(idx.only) {
        return(keep.idx)
    }
    
    PT <- PT[keep.idx,,drop = FALSE]

    # clean up
    PT <- lav_partable_complete(PT)

    # check if we have enough indicators?
    # TODO
    
    PT
}

lav_partable_subset_structural_model <- function(PT = NULL,
                                                 lavpta = NULL,
                                                 idx.only = FALSE) {

    # PT
    PT <- as.data.frame(PT, stringsAsFactors = FALSE)

    # lavpta
    if(is.null(lavpta)) {
        lavpta <- lav_partable_attributes(PT)
    }

    # ngroups
    ngroups <- lavpta$ngroups

    # eqs.names
    eqs.x.names <- lavpta$vnames$eqs.x
    eqs.y.names <- lavpta$vnames$eqs.y

    # keep rows idx
    keep.idx <- integer(0L)

    # remove not-needed measurement models
    for(g in 1:ngroups) {

        # eqs.names
        eqs.names <- unique( c(lavpta$vnames$eqs.x[[g]],
                               lavpta$vnames$eqs.y[[g]]) )

        # regressions
        reg.idx <- which(PT$op == "~" & PT$group == g &
                         PT$lhs %in% eqs.names &
                         PT$rhs %in% eqs.names)

        # the variances
        var.idx <- which(PT$op == "~~" & PT$group == g &
                         PT$lhs %in% eqs.names &
                         PT$rhs %in% eqs.names &
                         PT$lhs == PT$rhs)

        # optionally covariances (exo!)
        cov.idx <- which(PT$op == "~~" & PT$group == g &
                         PT$lhs %in% eqs.names &
                         PT$rhs %in% eqs.names &
                         PT$lhs != PT$rhs)

        # means/intercepts
        int.idx <- which(PT$op == "~1" & PT$group == g &
                         PT$lhs %in% eqs.names)

        keep.idx <- c(keep.idx, reg.idx, var.idx, cov.idx, int.idx)
    }

    if(idx.only) {
        return(keep.idx)
    }

    PT <- PT[keep.idx, , drop = FALSE]

    # clean up
    PT <- lav_partable_complete(PT)

    PT
}
