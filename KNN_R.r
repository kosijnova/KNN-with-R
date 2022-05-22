#KNN
#Nomralizacja - MinMax oraz ZScore
MinMax <- function(x, new_min = 0, new_max = 1) {
  return(((x - min(x)) / (max(x) - min(x))) * (new_max - new_min) + new_min)
}

ZScore <- function(x) {
  return ((x - mean(x))/sd(x))  
}

KNNtrain <- function(X, y_tar, k = 2, XminNew = 0, XmaxNew = 1) {
  if (k <= 0) {
    stop("K mniejsze bądź równe 0!")
  }
  else if (!is.matrix(X) & !is.data.frame(X)) {
    stop("X nie jest macierzą bądź ramką danych!")
  }
  else if (anyNA(X) | anyNA(y_tar)) {
    stop("Brakujące dane w X lub y_tar!")
  }
  else if (nrow(X) != length(y_tar)) {
    stop("Wrong size of X and Y!")
  }
  else {
    df <- data.frame(matrix(0,nrow = nrow(X),ncol = ncol(X)))
    minOrg <- c()
    maxOrg <- c()
    minmaxNew <- c()
    for (i in 1:length(X)) {
      if (all(sapply(X[i], class) == "numeric")) {
        minOrg <- append(minOrg,min(X[i]))
        maxOrg <- append(maxOrg,max(X[i]))
        X[i] <- as.vector(unlist(MinMax(X[i],XminNew,XmaxNew)))
        minmax <- c(c(min(X[i]),max(X[i])))
        names(minmax) <- c(paste(colnames(X[i]),"min",sep="_"),paste(colnames(X[i]),"max",sep="_"))
        minmaxNew <- append(minmaxNew,minmax)
      }
      else if (all(sapply(X[i], class) == "factor")) {
        maxOrg <- append(maxOrg,max(as.vector(as.numeric(unlist(X[i])))))
        minOrg <- append(minOrg,min(as.vector(as.numeric(unlist(X[i])))))
        X[i] <- factor(as.vector(unlist(X[i])))
      }
      else {
        maxOrg <- append(maxOrg,max(as.vector(unlist(X[i]))))
        minOrg <- append(minOrg,min(as.vector(unlist(X[i]))))
        X[i] <- (X[i])
      }
    }
    
    names(maxOrg) <- colnames(X)
    names(minOrg) <- colnames(X)
    attr(X,paste("minOrg")) <- minOrg
    attr(X,paste("maxOrg")) <- maxOrg
    attr(X,paste("minmaxNew")) <- minmaxNew
    knn <- list("X" = X, "y_tar" = y_tar, "k" = k)
    
    return (knn)
  }
}

#Miary niepodobieństwa
#zmienne ilorazowe (numeryczne)
d_euklides <- function(x_i, x_n) {
  return (sqrt(sum((x_i - x_n)^2)))
}

#zmienne nominalne
d_hamming <- function(x_i, x_n) {
  return (sum(x_i != x_n)/length(x_i))
}

#zmienne przedziałowe 
d_interval <- function(moc, x_i, x_n) {
  return (sum(abs(x_i - x_n)/(moc - 1)))
}

d_czebyszew <- function(x_i, x_n) {
  return (max(abs(x_i - x_n)))
}

#zmienne mieszane
d_gower <- function(dane,x_i,x_n) {
  if (length(x_i) != length(x_n)) {
    stop("\nwrong vectors length!")
  }
  else {
    out <- 0
    for (i in 1:length(x_i)) {
      if (all(sapply(x_i[i],class) == "numeric") | any(sapply(x_i[i],class) == "ordered")) {
        dist <- abs(as.numeric(paste(x_i[i])) - as.numeric(paste(x_n[i])))/(max(as.numeric(paste(dane[[names(x_i[i])]]))) - as.numeric(paste(min(dane[[names(x_n[i])]]))))
        out <- out + dist
      }
      else {
        dist <- sum(x_i[i] != x_n[i])
        out <- out + dist
      }
    }
    return (as.numeric(out)/length(x_i))
  }
}

#funkcje pomocnicze
#moda - dominanta
getMode <- function(x) {
  distinct <- unique(x)
  distinct_max <- distinct[which.max(tabulate(match(x,distinct)))]
  return (distinct_max)
}

#wartości unikalne
getDistinct <- function(x) {
  return (unique(x))
}

#Predykcje
KNNpred <-function(KNNmodel, X, Ncores = 1) {
  colCheck <- colnames(KNNmodel) == colnames(X)
  if (anyNA(X)) {
    stop("Brakujące dane w X!")
  }
  else if (all(colCheck) == FALSE) {
    stop("Złe kolumny!")
  }
  else {
    colTypes <- unique(unlist(lapply(X, class)))
    
    cores <- parallel::detectCores()
    if (Ncores > cores) {
      stop("\nToo much cores!")
    }
    
    cl <- parallel::makeCluster(Ncores, outfile = "")
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl,"d_euklides")
    parallel::clusterExport(cl,"getMode")
    parallel::clusterExport(cl,"getDistinct")
    parallel::clusterExport(cl,"d_hamming")
    parallel::clusterExport(cl,"d_gower")
    parallel::clusterExport(cl,"d_czebyszew")
    parallel::clusterExport(cl,"d_interval")
    parallel::clusterExport(cl,"KNNmodel",envir = environment())
    parallel::clusterExport(cl,"X",envir = environment())
    
    if (isTRUE(all.equal(c("ordered","factor"),colTypes))) {
      #skala porządkowa
      KNNmodel$X <- sapply(KNNmodel$X,as.numeric)
      X <- sapply(X,as.numeric)
      y <- KNNmodel$y
      if (is.factor(y)) { 
        print("Skala porządkowa - Klasyfikacja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_czebyszew(x, X[j,])))
        
        pred <- double(nPred)
        out <- c()
        parallel::clusterExport(cl,"pred",envir = environment())
        parallel::clusterExport(cl,"out",envir = environment())
        
        out <- (parSapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred[i] <- getMode(KNNmodel$y[kNaj])
          tab <- as.data.frame(table(KNNmodel$y[kNaj]))
          p <- data.frame(tab, "p" = tab$Freq/sum(tab$Freq))
          pr <- rbind(p$p)
          out <- rbind(out,as.vector(pr))
        }))
        
        out <- as.data.frame(t(out))
        out <- out[1:length(unique(KNNmodel$y))]
        colnames(out) <- lapply(1:length(unique(KNNmodel$y)), function(i) paste(sort(unique(KNNmodel$y))[i],sep = "_"))
        out$max <- apply(out, 1, max)
        pred <- colnames(out)[max.col(out, ties.method = "first")]
        drops <- c("max")
        out <- out[ , !(names(out) %in% drops)]
        
        return (data.frame("y_tar" = KNNmodel$y,out,y_hat = pred))
      } 
      else {
        print("Skala porządkowa - Regresja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_czebyszew(x, X[j,])))
        
        out <- (parLapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred <- mean(KNNmodel$y[kNaj])
        }))
        
        return (unlist(out))
      } 
    } #porządkowa
    else if (isTRUE(all.equal("numeric",colTypes))) {
      #skala ilorazowa -> odległość Euklidesa
      y <- KNNmodel$y
      
      parallel::clusterExport(cl,"y",envir = environment())
      
      if (is.factor(y)) {
        print("Skala ilorazowa - Klasyfikacja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_euklides(x, X[j,])))
        
        pred <- double(nPred)
        out <- c()
        parallel::clusterExport(cl,"pred",envir = environment())
        parallel::clusterExport(cl,"out",envir = environment())
        
        out <- (parSapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred[i] <- getMode(KNNmodel$y[kNaj])
          tab <- as.data.frame(table(KNNmodel$y[kNaj]))
          p <- data.frame(tab, "p" = tab$Freq/sum(tab$Freq))
          pr <- rbind(p$p)
          out <- rbind(out,as.vector(pr))
        }))
        
        out <- as.data.frame(t(out))
        out <- out[1:length(unique(KNNmodel$y))]
        colnames(out) <- lapply(1:length(unique(KNNmodel$y)), function(i) paste(sort(unique(KNNmodel$y))[i],sep = "_"))
        out$max <- apply(out, 1, max)
        pred <- colnames(out)[max.col(out, ties.method = "first")]
        drops <- c("max")
        out <- out[ , !(names(out) %in% drops)]
        
        return (data.frame("y_tar" = KNNmodel$y,out,y_hat = pred))
      }
      else {
        print("Skala ilorazowa - Regresja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_euklides(x, X[j,])))
        
        out <- (parLapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred <- mean(KNNmodel$y[kNaj])
        }))
        
        return (unlist(out))
      }
    } #ilorazowa
    else if (isTRUE(all.equal("factor",colTypes))) {
      #skala nominalna -> odległość Hamminga
      y <- KNNmodel$y
      if (is.factor(y)) {
        print("Skala nominalna - Klasyfikacja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_hamming(x, X[j,])))
        
        pred <- double(nPred)
        out <- c()
        parallel::clusterExport(cl,"pred",envir = environment())
        parallel::clusterExport(cl,"out",envir = environment())
        
        out <- (parSapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred[i] <- getMode(KNNmodel$y[kNaj])
          tab <- as.data.frame(table(KNNmodel$y[kNaj]))
          p <- data.frame(tab, "p" = tab$Freq/sum(tab$Freq))
          pr <- rbind(p$p)
          out <- rbind(out,as.vector(pr))
        }))
        
        out <- as.data.frame(t(out))
        out <- out[1:length(unique(KNNmodel$y))]
        colnames(out) <- lapply(1:length(unique(KNNmodel$y)), function(i) paste(sort(unique(KNNmodel$y))[i],sep = "_"))
        out$max <- apply(out, 1, max)
        pred <- colnames(out)[max.col(out, ties.method = "first")]
        drops <- c("max")
        out <- out[ , !(names(out) %in% drops)]
        
        return (data.frame("y_tar" = KNNmodel$y,out,y_hat = pred))
      }
      else {
        print("Skala nominalna - Regresja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_hamming(x, X[j,])))
        
        out <- (parLapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred <- mean(KNNmodel$y[kNaj])
        }))
        
        return (unlist(out))
      }
    } #nominalna
    else {
      #skala mieszana -> odległość Gowera
      y <- KNNmodel$y
      if (is.factor(y)) {
        print("Skala mieszana - Klasyfikacja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        parallel::clusterExport(cl,"y",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j)
            d_gower(KNNmodel$X, x, X[j,])))
        
        pred <- double(nPred)
        out <- c()
        parallel::clusterExport(cl,"pred",envir = environment())
        parallel::clusterExport(cl,"out",envir = environment())
        
        out <- (parSapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred[i] <- getMode(KNNmodel$y[kNaj])
          tab <- as.data.frame(table(KNNmodel$y[kNaj]))
          p <- data.frame(tab, "p" = tab$Freq/sum(tab$Freq))
          pr <- rbind(p$p)
          out <- rbind(out,as.vector(pr))
        }))
        
        out <- as.data.frame(t(out))
        out <- out[1:length(unique(KNNmodel$y))]
        colnames(out) <- lapply(1:length(unique(KNNmodel$y)), function(i) paste(sort(unique(KNNmodel$y))[i],sep = "_"))
        out$max <- apply(out, 1, max)
        pred <- colnames(out)[max.col(out, ties.method = "first")]
        drops <- c("max")
        out <- out[ , !(names(out) %in% drops)]
        
        return (data.frame("y_tar" = KNNmodel$y,out,y_hat = pred))
      }
      else {
        print("Skala mieszana - Regresja")
        nTrain <- nrow(KNNmodel$X)
        nPred <- nrow(X)
        
        parallel::clusterExport(cl,"nTrain",envir = environment())
        parallel::clusterExport(cl,"nPred",envir = environment())
        
        odl <- parApply(cl, KNNmodel$X, 1, function(x) 
          sapply(1:nPred, function(j) d_gower(KNNmodel$X, x, X[j,])))
        
        out <- (parLapply(cl, 1:nPred, function(i) {
          kNaj <- order(odl[,i])
          kNaj <- kNaj[1:KNNmodel$k]
          pred <- mean(KNNmodel$y[kNaj])
        }))
        
        return (unlist(out))
      }
    } #mieszana
  }
  
  parallel:stopCluster(cl)
  
}
