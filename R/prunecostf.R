
prunecostf <- function(edges, data, coly, colx, method=1,ind_col) 
{
  sswt <- sswf(data, unique(as.integer(edges)), coly, colx, method,ind_col)
  cores <- get.coresOption()
  if (is.null(cores)) {
    parallel <- "no"
  }
  else {
    parallel <- ifelse(get.mcOption(), "multicore", "snow")
  }
  ncpus <- ifelse(is.null(cores), 1L, cores)
  cl <- NULL
  if (parallel == "snow") {
    cl <- get.ClusterOption()
    if (is.null(cl)) {
      parallel <- "no"
      warning("no cluster in ClusterOption, parallel set to no")
    }
  }
  if (nrow(edges) < 300)
    parallel <- "no"
  if (parallel == "snow") {
    if (requireNamespace("parallel", quietly = TRUE)) {
      sI <- parallel::splitIndices(nrow(edges), length(cl))
      sswp <- do.call("c", parallel::parLapply(cl, sI,
                                               sapply, function(i) {
                                                 pruned.ids <- prunemst(rbind(edges[i, ], edges[-i,
                                                                                                ]), only.nodes = TRUE)
                                                 sum(sapply(pruned.ids, function(j) sswf(data,
                                                                                        j, coly, colx, method,ind_col)))
                                               }))
    }
    else {
      stop("parallel not available")
    }
  }
  else if (parallel == "multicore") {
    if (requireNamespace("parallel", quietly = TRUE)) {
      sI <- parallel::splitIndices(nrow(edges), ncpus)
      out <- parallel::mclapply(sI, sapply, function(i) {
        pruned.ids <- prunemst(rbind(edges[i, ], edges[-i,
                                                       ]), only.nodes = TRUE)
        sum(sapply(pruned.ids, function(j) sswf(data,
                                               j, coly, colx, method,ind_col)))
      }, mc.cores = ncpus)
      sswp <- do.call("c", out)
    }
    else {
      stop("parallel not available")
    }
  }
  else {
    sswp <- sapply(1:nrow(edges), function(i) {
      pruned.ids <- prunemst(rbind(edges[i, ], edges[-i,
                                                     ]), only.nodes = TRUE)
      
        sum(sapply(pruned.ids, function(j) sswf(data, j, coly, colx, method,ind_col))) #ORIGINALE

      # sumparz = matrix(NA, nrow = length(pruned.ids), ncol=1)      
      # for (k in 1:length(pruned.ids)) {
      #   if (lengths(pruned.ids)[k]<4) {
      #     sumparz[k,1] = sswf(data, pruned.ids[[k]], coly, colx, method=1)
      #   }
      #   else {
      #     sumparz[k,1] = sswf(data, pruned.ids[[k]], coly, colx, method)
      #        }
      # }
      # sum(sumparz[,1]) 
      # 
      
      
      
      
    })
    }
      return(sswt - sswp)
}

