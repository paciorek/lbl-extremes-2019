calc_riskRatio_pot <- function(returnValue, y1, y2, x1 = NULL, x2 = x1,
                               threshold1, threshold2 = threshold1, locationFun1 = NULL,
                               locationFun2 = locationFun1, scaleFun1 = NULL, scaleFun2 = scaleFun1,
                               shapeFun1 = NULL, shapeFun2 = shapeFun1, nBlocks1 = nrow(x1),
                               nBlocks2 = nrow(x2), blockIndex1 = NULL, blockIndex2 = NULL, firstBlock1 = 1,
                               firstBlock2 = 1, index1 = NULL, index2 = NULL, nReplicates1 = 1,
                               nReplicates2 = 1, replicateIndex1 = NULL, replicateIndex2 = NULL,
                               weights1 = NULL, weights2 = NULL, proportionMissing1 = NULL,
                               proportionMissing2 = NULL, xNew1 = NULL, xNew2 = NULL, declustering = NULL,
                               upperTail = TRUE, scaling1 = 1, scaling2 = 1, ciLevel = 0.90, ciType,
                               bootSE, bootControl = NULL, lrtControl = NULL, 
                               optimArgs = NULL, initial1 = NULL,
                               initial2 = NULL, logScale1 = NULL, logScale2 = NULL,
                               getReturnCalcs = FALSE, getParams = FALSE, getFit = FALSE  ) {

    if(is.null(returnValue) || is.null(y1) || is.null(y2))  # for Python interface
        stop("calc_riskRatio_binom: argument 'returnValue' or 'y1' or 'y2' is missing, with no default")

    if(missing(bootSE) || is.null(bootSE))
        bootSE <- !missing(ciType) && sum(ciType %in% names(climextRemes:::bootTypes))
    if(missing(returnValue))
        stop("calc_riskRatio_pot: 'returnValue' must be provided.")
    if(length(returnValue) > 1)
        stop("calc_riskRatio_pot: please provide a single value for 'returnValue'.")
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE)
    ciLabels <- as.character(c(( 1 - ciLevel ) / 2, 1 - ( 1 - ciLevel ) / 2))
    if(!is.null(xNew1)) {
        xNew1tmp <- try(as.data.frame(xNew1))
        if(class(xNew1tmp) == 'try-error') stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
        m <- nrow(xNew1tmp)
    } else {
        if(is.null(x1))
            m <- 1 else {
                       x1tmp <- try(as.data.frame(x1))
                       if(class(x1tmp) == 'try-error') stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
                       m <- nrow(x1tmp)
                   }
    }
    results <- list()

    if(missing(ciType)) {
        ciType <- ""
    } else {
        wh <- setdiff(ciType, c('delta', 'lrt', 'koopman', names(climextRemes:::bootTypes)))
        if(length(wh))
            stop("calc_riskRatio_pot: ", paste(ciType[wh], collapse = ','), " in 'ciType' are not valid types.")
    }

    bootTypesUse <- ciType[ciType %in% names(climextRemes:::bootTypes)]
    if(length(bootTypesUse)) {
        bControl <- list(seed = 1, n = 250, by = "block", getSample = FALSE)
        bControl[names(bootControl)] <- bootControl
        bootControl <- bControl
        bootControlTmp <- bootControl
        bootControlTmp$getSample <- TRUE
        bootControlTmp$getSampleSE <- TRUE
    } 
    
    fit1 <- fit_pot(y1, x = x1, threshold = threshold1, locationFun = locationFun1,
                    scaleFun = scaleFun1, shapeFun = shapeFun1, nBlocks = nBlocks1,
                    blockIndex = blockIndex1, firstBlock = firstBlock1, index = index1,
                    nReplicates = nReplicates1, replicateIndex = replicateIndex1,
                    weights = weights1, proportionMissing = proportionMissing1, returnValue = returnValue,
                    xNew = xNew1, declustering = declustering, upperTail = upperTail,
                    scaling = scaling1, bootSE = bootSE, bootControl = bootControlTmp,
                    optimArgs = optimArgs, initial = initial1, logScale = logScale1,
                    getFit = TRUE, getParams = getParams, .getInputs = TRUE) 
    fit2 <- fit_pot(y2, x = x2, threshold = threshold2, locationFun = locationFun2,
                    scaleFun = scaleFun2, shapeFun = shapeFun2, nBlocks = nBlocks2,
                    blockIndex = blockIndex2, firstBlock = firstBlock2, index = index2,
                    nReplicates = nReplicates2, replicateIndex = replicateIndex2,
                    weights = weights2, proportionMissing = proportionMissing2, returnValue = returnValue,
                    xNew = xNew2, declustering = declustering, upperTail = upperTail,
                    scaling = scaling2, bootSE = bootSE, bootControl = bootControlTmp,
                    optimArgs = optimArgs, initial = initial2, logScale = logScale2,
                    getFit = TRUE, getParams = getParams, .getInputs = TRUE)
    if(fit1$info$failure || fit2$info$failure) {
        warning("calc_riskRatio_pot: fitting failed for one of two datasets.")
        results$logRiskRatio <- results$se_logRiskRatio <- results$riskRatio <- rep(NA, m)
        for(type in ciType)
            results[[paste0('ci_riskRatio_', type)]] <- drop(matrix(NA, m, 2))
        if(bootSE) 
            results$se_logRiskRatio_boot <- rep(NA, m)
    } else {
        if(length(fit1$logReturnProb) != length(fit2$logReturnProb))
            stop("calc_riskRatio_pot: number of return probabilities calculated for each model fit must be the same; for nonstationary models this is determined by the number of covariate set inputs (provided in 'x' or 'xNew').")
        results$logRiskRatio <- fit1$logReturnProb - fit2$logReturnProb
        
        results$se_logRiskRatio <- sqrt(fit1$se_logReturnProb^2 + fit2$se_logReturnProb^2)
        results$riskRatio <- exp(results$logRiskRatio)
        if('delta' %in% ciType) {
            results$ci_riskRatio_delta <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio,
                                                    results$logRiskRatio + z_alpha*results$se_logRiskRatio))
            colnames(results$ci_riskRatio_delta) <- ciLabels
            results$ci_riskRatio_delta <- drop(results$ci_riskRatio_delta)
        }
        
        if(length(bootTypesUse)) {
            bootData <- matrix(drop(fit1$logReturnProb_boot - fit2$logReturnProb_boot), ncol = m)
            if(bootSE)
                results$se_logRiskRatio_boot <- apply(bootData, 2, sd)
            
            fake_data <- cbind(sample(c(0,1), size = 5, replace = TRUE),
                   sample(c(0,1), size = 5, replace = TRUE))
            logRRfun <- function(dat, ind) {
                log(sum(dat[ind,1])/sum(dat[ind,2]))
            }
            bootInput <- boot::boot(fake_data, logRRfun, R = bootControl$n)

            for(type in bootTypesUse) {
                typeName <- paste0('ci_riskRatio_', type)
                results[[typeName]] <- matrix(NA, m, 2)
                colnames(results[[typeName]]) <- ciLabels
            }
            for(i in seq_len(m)) {
                bootInput$t0 = c(results$logRiskRatio[i], results$se_logRiskRatio[i]^2)
                logRRvals <- bootData[ , i]
                logRRvals[logRRvals == Inf] = 1e6 # boot.ci can't handle Inf
                bootInput$t <- cbind(logRRvals,
                                     fit1$logReturnProb_boot_se[ , i, 1]^2 + fit2$logReturnProb_boot_se[ , i, 1]^2)
                bootResults <- try(boot::boot.ci(bootInput, conf = ciLevel, type = unlist(bootTypes[bootTypesUse])))
                if(!is(bootResults, 'try-error')) {
                    for(type in bootTypesUse) {
                        typeName <- paste0('ci_riskRatio_', type)
                        if(type == 'boot_norm') {
                            results[[typeName]][i, ] <- exp(bootResults[[bootTypesCols[[type]]]][2:3])
                        } else results[[typeName]][i, ] <- exp(bootResults[[bootTypesCols[[type]]]][4:5])
                        if(i == m) results[[typeName]] <- drop(results[[typeName]])
                    }
                }
            }
            if(bootControl$getSample) results$riskRatio_boot <- drop(exp(bootData))
        }
        
        if('lrt' %in% ciType) {
            lControl <- list(bounds = c(.01, 100))
            lControl[names(lrtControl)] <- lrtControl
            lrtControl <- lControl
            oArgs <- list(method = "Nelder-Mead", lower = -Inf, upper = Inf, control = list())
            oArgs[names(optimArgs)] <- optimArgs
            if(!upperTail) returnValue <- -returnValue
            results$ci_riskRatio_lrt <- climextRemes:::calc_riskRatio_lrt(fit1, fit2, returnValue, ciLevel = ciLevel, bounds = lrtControl$bounds, type = "PP", optimArgs = oArgs)
            if(is.null(dim(results$ci_riskRatio_lrt))) {
                names(results$ci_riskRatio_lrt) <- ciLabels
            } else colnames(results$ci_riskRatio_lrt) <- ciLabels
        }
    }
    if(getFit || getParams || getReturnCalcs) {
        fit1$inputs <- fit2$inputs <- NULL
        if(!getFit) {
            fit1$fit <- fit2$fit <- fit1$info <- fit2$info <- NULL
        }
        if(!getReturnCalcs) {
            fit1$logReturnProb <- fit1$se_logReturnProb <- fit1$logReturnPeriod <- fit1$se_logReturnPeriod <- NULL
            fit2$logReturnProb <- fit2$se_logReturnProb <- fit2$logReturnPeriod <- fit2$se_logReturnPeriod <- NULL
            if(bootSE) {
                fit1$se_logReturnProb_boot <- fit1$se_logReturnPeriod <- NULL
                fit2$se_logReturnProb_boot <- fit2$se_logReturnPeriod <- NULL
            }
        }
        results$fit1 <- fit1
        results$fit2 <- fit2
    }
    return(results)
} # end calc_riskRatio_pot
