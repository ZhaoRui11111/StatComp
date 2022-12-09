#' @title A illustration dataset
#' @name document
#' @description A dataset used to illustrate the performance of \code{vaccR} and \code{vaccC}.
#' @examples
#' \dontrun{
#' document
#' }
NULL

#' Title a function helps for LDA model visualization.
#'
#' @param data a document-word list
#'
#' @return no return
#' @importFrom "LDAvis" "createJSON"
#' @importFrom "LDAvis" "serVis"
#' @importFrom "lda" "lda.collapsed.gibbs.sampler"
#' @importFrom "stats" "rmultinom"
#' 
#'@export
LDAVIS <- function(data){
  doc.list <- strsplit(data, "[[:space:]]+")
  term.table <- table(unlist(doc.list))
  term.table <- sort(term.table, decreasing = TRUE)

  del <- names(term.table) %in%  term.table < 5
  term.table <- term.table[!del]
  vocab <- names(term.table)
  get.terms <- function(x){
    return(get_termsC(x,vocab))
  }
  documents <- lapply(doc.list, get.terms)
  doc.length <- sapply(documents, function(x) sum(x[2, ]))  
  term.frequency <- as.integer(term.table) 

  K <- 13
  G <- 3000
  alpha <- 0.02
  eta <- 0.02
  
  fit <- lda.collapsed.gibbs.sampler(documents = documents, K = K, vocab = vocab, num.iterations = G, alpha = alpha, eta = eta, initial = NULL, burnin = 0,compute.log.likelihood = TRUE)
  
  theta <- t(apply(fit$document_sums + alpha, 2, function(x) x/sum(x)))
  
  phi <- t(apply(t(fit$topics) + eta, 2, function(x) x/sum(x)))
  
  MovieReviews <- list(phi = phi,theta = theta,doc.length = doc.length,vocab = vocab,term.frequency = term.frequency)
  
  json <- createJSON(phi = MovieReviews$phi, 
                     theta = MovieReviews$theta, 
                     doc.length = MovieReviews$doc.length, 
                     vocab = MovieReviews$vocab, 
                     term.frequency = MovieReviews$term.frequency)
  
  serVis(json, out.dir = 'vis', open.browser = TRUE)
}

#' Title a function helps for LDA model
#'
#' @param data a document-word list
#' @param num_topic  the number of topics
#' @param num_word   the number of the top words to return
#'
#' @return  a list with top num_word words of every topic
#' 
#'@export
LDAR <- function(data,num_topic,num_word){
  maxdoclen=0
  for (i in 1:length(data)) {
    line <- strsplit(data[i], "[[:space:]]+")
    if (length(line[[1]])>maxdoclen)
      maxdoclen = length(line[[1]])
  }
  tokens <- strsplit(data, "[[:space:]]+")
  word_dict <- table(unlist(tokens))
  word_dict <- sort(word_dict, decreasing = TRUE)
  for (i in 1:length(word_dict)) {
    word_dict[i]=i
  }
  
  vocasize=length(word_dict)
  docnum=length(tokens)  
  n_features = num_word
  iter_times = 200
  
  K=num_topic
  beta = 0.1
  alpha = 0.5
  p <-  numeric(K)
  nw <- matrix(0,nrow = vocasize,ncol = K)
  nwsum <-  numeric(K)
  nd <-  matrix(0,nrow = docnum,ncol = K)
  ndsum <-  numeric(docnum)
  z1 <- z2 <- matrix(0,nrow = docnum,ncol = maxdoclen)
  theta <-  matrix(0,nrow = docnum,ncol = K)
  phi1 <- phi2 <- matrix(0,nrow = K,ncol = vocasize)
  

  for (i in 1:docnum) {
    tokensi <- tokens[[i]]
    ndsum[i] <- length(tokensi)
    for (j in 1:length(tokensi)) {
      word <- tokensi[j]
      topic_index <- sample(1:K, 1)
      word_id = as.numeric(word_dict[word]) 
      z1[i,j] <- topic_index
      z2[i,j] <- word_id
      nw[word_id,topic_index] =nw[word_id,topic_index]+1   
      nd[i,topic_index] = nd[i,topic_index]+1              
      nwsum[topic_index] = nwsum[topic_index]+1          
    }
  }

  Vbeta = vocasize * beta
  Kalpha = K * alpha
  for (iter in 1:iter_times) {
    for (i in 1:docnum) {
      for (j in 1:ndsum[i]) {
        word_id = z2[i,j]
        topic = z1[i,j]
        nw[word_id,topic] =nw[word_id,topic]- 1
        nd[i,topic] =nd[i,topic]- 1
        nwsum[topic] = nwsum[topic]-1
        for (k in 1:K) {
          p[k] = (nw[word_id,k] + beta) / (nwsum[k] + Vbeta) * (nd[i,k] + alpha) / (ndsum[i] + Kalpha)
        }
        s <- rmultinom(1,1,p/sum(p))
        t <- which(s==max(s))
        nw[word_id,t] = nw[word_id,t]+1
        nwsum[t] = nwsum[t]+1
        nd[i,t] = nd[i,t]+1
        z1[i,j] = t
        z2[i,j] = word_id
      }
    }
  }
  
  for (i in 1:docnum) {
    for (j in 1:K) {
      theta[i,j] = (nd[i,j] + alpha) / (ndsum[i] + Kalpha)
    }
  }
  for (i in 1:K) {
    for (j in 1:vocasize) {
      phi1[i,j]=(nw[j,i]+beta)/(ndsum[i]+Vbeta)
      phi2[i,j]=j
    }
  }
  result <- list()
  for (i in 1:K) {
    indexphi <- order(phi1[i,],decreasing = TRUE)
    hotword<- character(length = n_features)
    for (j in 1:n_features) {
      topicword = word_dict[phi2[i,indexphi[j]]]
      hotword[j] <- names(topicword)
    }
    result[[i]] <- hotword
  }
  return(result)
}