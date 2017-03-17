# function from gamlr::corr to calculate correlation between simple triplet matrix
# and matrix, to avoid filling memory with a giant regular matrix
cor_slam <- function (x, y)
{
      if (!inherits(x, "simple_triplet_matrix")) {
            return(cor(x, y))
      }
      n <- nrow(x)
      v <- t(normalize(y))
      r <- tcrossprod_simple_triplet_matrix(t(x)/sdev(x), v)/(nrow(x) - 
                                                                    1)
      dimnames(r) <- list(dimnames(x)[[2]], dimnames(y)[[2]])
      return(r)
}  
# this one also comes from gamlr...
normalize <-  function (x, m = NULL, s = sdev(x)) 
{
      if (!is.null(ncol(x))) 
            if (length(s) != ncol(x)) 
                  stop("length(s)!=ncol(x)")
      s[s == 0] <- 1
      if (is.simple_triplet_matrix(x)) {
            x$v <- x$v/s[x$j]
            return(x)
      }
      x <- as.matrix(x)
      if (is.null(m)) 
            m <- col_means(x)
      return(t((t(x) - m)/s))
}

# and this one also...
sdev <- function (x) 
{
      if (!inherits(x, "simple_triplet_matrix")) {
            return(apply(as.matrix(x), 2, sd))
      }
      n <- nrow(x)
      sqrt(col_sums(x^2)/(n - 1) - col_sums(x)^2/(n^2 - n))
      return(sqrt(col_sums(x^2)/(n - 1) - col_sums(x)^2/(n^2 - 
                                                               n)))
}

# need a function from this one also... can't be bothered copy-pasting...


# here they are all wrapped up: calculate the correlation, subset by corlimit and put in order
findAssocsBig <- function(u, word, corlimit){
      suppressWarnings(x.cor <-  cor_slam(          u[ ,!u$dimnames$Terms == word],        
                                                    as.matrix(u[  ,u$dimnames$Terms == word ])  )  )  
      x <- sort(round(x.cor[ (x.cor[ ,word ] > corlimit), ], 3), decreasing = TRUE)
      return(x)
}

