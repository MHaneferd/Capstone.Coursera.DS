
## Cluster
```{r, cache=FALSE,message=FALSE}
####################
#####3
####### Cluster
#######
#################

# This sparse percentage denotes the proportion of empty elements. A sparse parameter of 0.7 means that we select from those terms which are less than 70% empty. We set the sparsity at 0.95 so that terms with at least that sparse percentage will be removed. So the terms we accept can be very empty- at most 95% empty. Then we can coerce the TD matrix into a regular matrix.
tdm <- TermDocumentMatrix(docs, control=list(wordLengths=c(3, 15), bounds = list(global = c(1000,10000))))

tdms <- removeSparseTerms(tdm, 0.99)
dim(tdms)
# Delete rows that are zero:
col_freqr <-  col_sums(tdms)
tdms <- tdms[,row_freqr>0]
dim(tdms)

set.seed(1000)
sam <- sample(seq(1, 12056), 12000)
tdms <- tdms[,sam]

memory.limit(size=60000)
m <- as.matrix(tdms)

memory.limit(size=60000)
#rownames(m) <- 1:nrow(m)

## hierarchical clustering
library(proxy)
#compute distance between document vectors
d <- dist(m, method="cosine")

hc <- hclust(d, method="average")

plot(hc, hang=-1)

cl <- cutree(hc, 15)

table(cl)

findFreqTerms(tdms[cl==1], 5)


#run hierarchical clustering using Ward’s method
groups <- hclust(d,method="ward.D")
#plot dendogram, use hang to ensure that labels fall below tree
plot(groups, hang=-1)

rect.hclust(groups, k=5)




# dtm <- DocumentTermMatrix(docs, control=list(wordLengths=c(4, 15)))

# Here we have told R to include only those words that occur in  500 to 2000 documents. We have also enforced  lower and upper limit to length of the words included (between 4 and 20 characters).


# Inspect dtm
dtm

freqr <-  col_sums(dtm)
length(freqr)
#create sort order (asc)
ordr <- order(freqr,decreasing=TRUE)
#inspect most frequently occurring terms
freqr[head(ordr)]
# The output:
#    one   will   like    can   time   just 
# 136263 116009 111796 109971 108426 100612 
# Dous not tell much, lets see what the least output is:
#inspect least frequently occurring terms
freqr[tail(ordr)]

# keep only  words that occur >1000 times in all docs
#dtm2 <- dtm[,which( freqr > 30000)]
#dtm2
#freqr <-  col_sums(dtm2)

# Delete rows that are zero:
row_freqr <-  row_sums(dtm)
dtm <- dtm[row_freqr>0,]

# I will investigate the words freedom and sea a little bit further. To see the associations
# We can also find associations with a word, specifying a correlation limit

assosiations <- findAssocs(dtm, c("film","game"), corlimit=0.2)


d <- data.frame(word = colnames(dtm),freq=freqr)
dd <- head(d[order(d$freq,decreasing=TRUE),], 20)

barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word, col ="lightblue", main ="Most frequent words", ylab = "Word frequencies")


##########################################
#############
#############    Clustering
#############
##########################################
## do tfxidf
dtm_tfxidf <- weightTfIdf(dtm)
inspect(dtm_tfxidf[1:10, 5001:5010])

## do document clustering

### k-means (this uses euclidean distance)
m <- as.matrix(dtm_tfxidf)
rownames(m) <- 1:nrow(m)

### don't forget to normalize the vectors so Euclidean makes sense
norm_eucl <- function(m) m/apply(m, MARGIN=1, FUN=function(x) sum(x^2)^.5)
m_norm <- norm_eucl(m)


### Remove NA's
m_norm[is.na(m_norm)] <- 0

### cluster into 10 clusters
cl <- kmeans(m_norm, 10)
cl

table(cl$cluster)

### show clusters using the first 2 principal components
plot(prcomp(m_norm)$x, col=cl$cl)

findFreqTerms(dtm_tfxidf[cl$cluster==1], 50)
inspect(reuters[which(cl$cluster==1)])


### this is going to take 4-ever (O(n^2))
memory.limit(size=40000)
d <- dist(m, method="cosine")
hc <- hclust(d, method="average")
plot(hc, hang=-1)

cl <- cutree(hc, 5)
table(cl)

findFreqTerms(dtm[cl==1], 50)




###################################
# Quantitive analyzis of text
words <- dtms %>% as.matrix %>% colnames %>% (function(x) x[nchar(x) < 20])

length(words)

head(words, 15)

summary(nchar(words))

table(nchar(words))

dist_tab(nchar(words))






###################################

freq <- rollup(tdm, 2, na.rm=TRUE, FUN = sum)
class(freq)

# Explore the document term matrix
sumfreq <- rowSums(as.matrix(freq))
length(sumfreq)







# Clustering

# clus: a skmeans object
# dtm: a Document Term Matrix
# first: eg. 10 most frequent words per cluster
# unique: if FALSE all words of the DTM will be used
#         if TRUE only cluster specific words will be used 



# result: List with words and frequency of words 
#         If unique = TRUE, only cluster specific words will be considered.
#         Words which occur in more than one cluster will be ignored.



mfrq_words_per_cluster <- function(clus, dtm, first = 10, unique = TRUE){
      if(!any(class(clus) == "skmeans")) return("clus must be an skmeans object")
      
      dtm <- as.simple_triplet_matrix(dtm)
      indM <- table(names(clus$cluster), clus$cluster) == 1 # generate bool matrix
      
      hfun <- function(ind, dtm){ # help function, summing up words
            if(is.null(dtm[ind, ]))  dtm[ind, ] else  col_sums(dtm[ind, ])
      }
      frqM <- apply(indM, 2, hfun, dtm = dtm)
      
      if(unique){
            # eliminate word which occur in several clusters
            frqM <- frqM[rowSums(frqM > 0) == 1, ] 
      }
      # export to list, order and take first x elements 
      res <- lapply(1:ncol(frqM), function(i, mat, first)
            head(sort(mat[, i], decreasing = TRUE), first),
            mat = frqM, first = first)
      
      names(res) <- paste0("CLUSTER_", 1:ncol(frqM))
      return(res)
}

# dtm3 has no rows with zero.

LDA(dtm, 30)

rownames(dtm3) <- paste0("Doc_", 1:dim(dtm3)[1])
clus <- skmeans(dtm3, 3)
# weightBin(x)
# Download CLUTO at:
# http://glaros.dtc.umn.edu/gkhome/cluto/cluto/download

file.create("cluto.file.test")
con <- file("cluto.file.test", open = "rw") 
writeCM(dtm3, con)
readCM(con, clabel = NULL)

close(con)

cluto <- doc2mat(dtm3)

clus <- skmeans(dtm3, method = "CLUTO", k = 4, control = list(nruns = 20))

mfrq_words_per_cluster(clus, dtm3)
mfrq_words_per_cluster(clus, dtm3, unique = FALSE)


```