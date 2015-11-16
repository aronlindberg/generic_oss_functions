# Function for reading sequence data WITH timestamps
read_seqdata <- function(data, startdate, stopdate){
  data <- read.table(data, sep = ",", header = TRUE)
  data <- subset(data, select = c("pull_req_id", "action", "created_at"))
  colnames(data) <- c("id", "event", "time")
  data <- sqldf(paste0("SELECT * FROM data WHERE strftime('%Y-%m-%d', time,
                       'unixepoch', 'localtime') >= '",startdate,"' AND strftime('%Y-%m-%d', time,
                       'unixepoch', 'localtime') <= '",stopdate,"'"))
  data$end <- data$time
  data <- data[with(data, order(time)), ]
  data$time <- match(data$time, unique(data$time))
  data$end <- match(data$end, unique(data$end))
  (data)
}

# Function for reading sequence data withOUT timestamps
# Also replaces "synchronized" and "subscribed" with "NA"
read_seqdata_notime <- function(data, startdate, stopdate){
  data <- read.table(data, sep = ",", header = TRUE)
  data <- subset(data, select = c("pull_req_id", "action", "created_at"))
  data <- subset(data, action!="synchronize")
  # data <- subset(data, action!="assigned")
  data <- subset(data, action!="subscribed")
  data <- subset(data, action!="unsubscribed")
  # data <- subset(data, action!="merged")
  # data <- subset(data, action!="mentioned")
  # data <- subset(data, action!="referenced")
  data <- subset(data, action!="head_ref_cleaned")
  data <- subset(data, action!="head_ref_deleted")
  data <- subset(data, action!="head_ref_restored")
  colnames(data) <- c("id", "event", "time")
  data <- sqldf(paste0("SELECT * FROM data WHERE strftime('%Y-%m-%d', time,
                       'unixepoch', 'localtime') >= '",startdate,"' AND strftime('%Y-%m-%d', time,
                       'unixepoch', 'localtime') <= '",stopdate,"'"))
  data.split <- split(data$event, data$id)
  list.to.df <- function(arg.list) {
    max.len  <- max(sapply(arg.list, length))
    arg.list <- lapply(arg.list, `length<-`, max.len)
    as.data.frame(arg.list)
  }
  data <- list.to.df(data.split)
  data <- t(data)
  (data)
}

# Function for opening network data for processinr
read_seqdata_for_network <- function(data, startdate, stopdate){
  data <- read.table(data, sep = ",", header = TRUE)
  data <- subset(data, select = c("pull_req_id", "user", "action", "created_at"))
  data <- subset(data, action!="synchronize")
  # data <- subset(data, action!="assigned")
  data <- subset(data, action!="subscribed")
  data <- subset(data, action!="unsubscribed")
  # data <- subset(data, action!="mentioned")
  # data <- subset(data, action!="referenced")
  # data <- subset(data, action!="merged")
  data <- subset(data, action!="head_ref_cleaned")
  data <- subset(data, action!="head_ref_deleted")
  data <- subset(data, action!="head_ref_restored")
  colnames(data) <- c("id", "actor", "event", "time")
  data <- sqldf(paste0("SELECT * FROM data WHERE strftime('%Y-%m-%d', time,
                       'unixepoch', 'localtime') >= '",startdate,"' AND strftime('%Y-%m-%d', time,
                       'unixepoch', 'localtime') <= '",stopdate,"'"))
  data$end <- data$time
  data <- data[with(data, order(time)), ]
  data$time <- match( data$time , unique( data$time ) )
  data$end <- match( data$end , unique( data$end ) )
  slmax <- max(data$time)
  (data)
}

# Function for counting the number of events
event_count <- function(data){
  sequences.sts <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  event.count <- seqstatf(sequences.sts)
  (sum(event.count$Freq))
}

# Function for storing sequence length in a variable
sequence_length <- function(data){
  # slmax <- max(data$time)
  # sequences.seqe <- seqecreate(data)
  # sequences.sts <- seqformat(data, from="SPELL", to="DSS", begin="time", end="end", id="id", status="event", limit=slmax)
  # sequences.sts <- seqdef(data, informat="SPELL", var=c("id", "time", "end", "event"), process=FALSE, right = "DEL", left = "DEL", gaps = "DEL")
  sequences.sts <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  sequences.length <- seqlength(sequences.sts)
  (sequences.length)
}

# Function for calculating entropies
entropies <- function(data){
  # slmax <- max(data$time)
  # sequences.sts <- seqformat(data, from="SPELL", to="STS", begin="time", end="end", id="id", status="event", limit=slmax)
  # sequences.sts <- seqdef(sequences.sts)
  sequences.sts <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  sequences.ent <- seqient(sequences.sts, norm = TRUE) # This stores the entropies
  (sequences.ent)
}

# Function for calculating subsequences
subsequences <- function(data){
  # slmax <- max(data$time)
  # sequences.seqe <- seqecreate(data)
  # sequences.sts <- seqformat(data, from="SPELL", to="DSS", begin="time", end="end", id="id", status="event", limit=slmax)
  sequences.sts <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  # sequences.sts <- seqdef(sequences.sts, right = "DEL", left = "DEL", gaps = "DEL")
  sub.sequences <- seqsubsn(sequences.sts, DSS = FALSE)
  (sub.sequences)
}

# Function for calculating turbulence
turbulence <- function(data){
  # slmax <- max(data$time)
  # sequences.sts <- seqformat(data, from="SPELL", to="STS", begin="time", end="end", id="id", status="event", limit=slmax)
  sequences.sts <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  # sequences.sts <- seqdef(sequences.sts)
  sequences.turb <- seqST(sequences.sts) # This stores the entropies
  (sequences.turb)
}

# Function for generating a dissimilarity value
dissimilarity <- function(data, label){
  # slmax <- max(data$time)
  # sequences.seqe <- seqecreate(data)
  # sequences.seqe <- seqformat(data, from="SPELL", to="STS", begin="time", end="end", id="id", status="event", limit=slmax)
  # sequences.sts <- seqdef(sequences.seqe, left = "DEL", right = "DEL", gaps = "DEL")
  sequences.sts <- seqdef(data, left = "DEL", gaps = "DEL")
  ccost <- seqsubm(sequences.sts, method = "CONSTANT", cval = 2, with.missing=TRUE)
  sequences.OM <- seqdist(sequences.sts, method = "OM", norm = TRUE, sm = ccost, with.missing=TRUE)
  clusterward <- agnes(sequences.OM, diss = TRUE, method = "ward")
  plot(clusterward, which.plots = 2)
  # (clusterward$ac)
  (sequences.OM)
}

# Function for least-squares normalization of subseq by length
normalization <- function(order, length){
  data <- as.data.frame(cbind(order, length))
  model <- lm(order ~ length, data)
  order_normalized <- order-(model$coefficients[1]+(model$coefficients[2]*length))
  (order_normalized)
}

# Function for normalizing subsequences by max_subseq_possible
normalization2 <- function(order, length){
  order_normalized <- order/2^max(length)
  (order_normalized*1000000)
}

# Function for clustering the data
clustering <- function(data, k, c){
  data <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  ccost <- seqsubm(data, method = "CONSTANT", cval = 2, with.missing=TRUE)
  data.om <- seqdist(data, method = "OM", norm = TRUE, sm = ccost, with.missing=TRUE)
  clusterward <- agnes(data.om, diss = TRUE, method = "ward")
  (clusterward)
}

# Function for creating a distance matrix (also needed for clustering)
distance_matrix <- function(data){
  data <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  ccost <- seqsubm(data, method = "CONSTANT", cval = 2, with.missing=TRUE)
  data.om <- seqdist(data, method = "OM", norm = TRUE, sm = ccost, with.missing=TRUE)
  (data.om)
}

# Cutting clusters
cluster_cut <- function(data, clusterward, n_clusters, name_clusters){
  data <- seqdef(data, left = "DEL", gaps = "DEL", right = "DEL")
  cluster4 <- cutree(clusterward, k = n_clusters)
  cluster4 <- factor(cluster4, labels = paste("Type", as.character(seq(from = 1, to = n_clusters, by = 1))))
  (data[cluster4==name_clusters,])
}

# This function makes sure I get the pagination right
digest_header_links <- function(x) {
  y <- x$headers$link
  if(is.null(y)) {
    # message("No links found in header.")
    m <- matrix(0, ncol = 3, nrow = 4)
    links <- as.data.frame(m)
    names(links) <- c("rel", "per_page", "page")
    return(links)
  }
  y %>%
    str_split(", ") %>% unlist %>%  # split into e.g. next, last, first, prev
    str_split_fixed("; ", 2) %>%    # separate URL from the relation
    plyr::alply(2) %>%              # workaround: make into a list
    as.data.frame() %>%        # convert to data.frame, no factors!
    setNames(c("URL", "rel")) %>%   # sane names
    dplyr::mutate_(rel = ~ str_match(rel, "next|last|first|prev"),
                   per_page = ~ str_match(URL, "per_page=([0-9]+)") %>%
                     `[`( , 2) %>% as.integer,
                   page = ~ str_match(URL, "&page=([0-9]+)") %>%
                     `[`( , 2) %>% as.integer,
                   URL = ~ str_replace_all(URL, "<|>", ""))
}

download_pulls_query <- function(owner, repo){
  # This function pulls down data on all the pull requests.
  pull <- function(i){
    commits <- get.pull.request.commits(owner = owner, repo = repo, id = i, ctx = get.github.context(), per_page=100)
    links <- digest_header_links(commits)
    number_of_pages <- links[2,]$page
    if (number_of_pages != 0)
      try_default(for (n in 1:number_of_pages){
        if (as.integer(commits$headers$`x-ratelimit-remaining`) < 5)
          Sys.sleep(as.integer(commits$headers$`x-ratelimit-reset`)-as.POSIXct(Sys.time()) %>% as.integer())
        else
          get.pull.request.commits(owner = owner, repo = repo, id = i, ctx = get.github.context(), per_page=100, page = n)
      }, default = NULL)
    else 
      return(commits)
  }
  
  list <- pr_id_list
  
  pull_lists <- lapply(list, pull)
  return(pull_lists)
}

download_commits_query <- function(owner, repo){
  # This is a function for getting all the correct SHAs (ignores parent and tree SHAs)
  all_shas <- lapply(
    pull_lists, 
    function(x) unlist(x)[grep("sha$",names(unlist(x)))]
  )
  sha_list <- lapply(all_shas, function(x) x[attr(x, "names")=="content.sha"])
  # this removes all the NULL values
  # sha_list_clean <- sha_list[ ! sapply(sha_list, is.null) ]
  get_commits <- function(sha){
    commits <- get.commit(git = NULL, ctx = get.github.context(), owner = owner, repo = repo, sha = sha, per_page = 100)
    links <- digest_header_links(commits)
    number_of_pages <- links[2,]$page
    if (number_of_pages != 0)
      try_default(for (n in 1:number_of_pages){
        if (as.integer(commits$headers$`x-ratelimit-remaining`) < 5)
          Sys.sleep(as.integer(commits$headers$`x-ratelimit-reset`)-as.POSIXct(Sys.time()) %>% as.integer())
        else
          get.commit(git = NULL, ctx = get.github.context(), owner = owner, repo = repo, sha = sha, per_page = 100, page = n)
      }, default = NULL)
    else 
      return(commits)
  }
  commit_lists0 <- lapply(sha_list, function(x) lapply(x, get_commits))
  return(commit_lists0)
}

# turn the data into event-graphs
graph_transformation <- function(dat){
  dat2 <- ddply(dat, .(id), function(d){
    data.frame(
      event = d$event[-1],
      from = d$event[-NROW(d)],
      to = d$event[-1],
      time = paste(d$time[-NROW(d)], d$time[-1], sep = "-")
    )
  })
  rc1_edge_list <- cbind(dat2["from"], dat2["to"])
  network_data_object <- try_default(as.network(rc1_edge_list[, 1:2], hyper=FALSE, loops=TRUE, multiple=TRUE), default = NULL)
  (network_data_object)
}

# This function extracts the degrees for each event (i.e. the number of sequential transitions it has with other event types)
extract_degrees <- function(x){
  list_rc <- split(x, x$id)
  sna_list <- lapply(list_rc, graph_transformation)
  
  edgecount_list <- list()
  for (i in 1:length(sna_list)){
    if(is.null(sna_list[[i]])) { 
      edgecount_list[i] <- list(NULL) 
    } else {
      edgecount_list[i] <- network.edgecount(sna_list[[i]])
    }
  }
  for (i in 1:length(edgecount_list)){
    if(is.null(edgecount_list[[i]])) { 
      edgecount_list[[i]] <- 0 
    } else {
      edgecount_list[[i]] <- edgecount_list[[i]]
    }
  }
  
  (unlist(edgecount_list))
}