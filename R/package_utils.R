#'        Quiet version of the function "require"
#'
#'  Loads packages quietly
#'
#' @export
Require<- function(pkg,...) {
  suppressWarnings(suppressMessages(require(pkg, character.only=TRUE, quietly = TRUE, ...)))
}

#'     Get the number of R processes running
#'
#'  Uses 'tasklist' to count the number of `Rscript.exe` instances
#'
#' Modified from https://stackoverflow.com/questions/15935931/detect-number-of-running-r-instances-in-windows-within-r
#' @export
getNumRprocess <- function(){

  targets=c("Rscript.exe","rsession.exe","rstudio.exe","Rgui.exe", "R")

  switch(Sys.info()[['sysname']],
      Windows = {
        progs <- system("tasklist", intern = TRUE)
        progs <- vapply(stringr::str_split(progs, "[[:space:]]"), "[[", "", i = 1)
        sum(progs %in% targets)
      },
      Linux   =  {
        user_name <- Sys.info()[["user"]]
        res <- system(paste0("pgrep -c -u ",user_name," '",paste(targets,collapse="|"),"'"), intern = TRUE)
        if(!is.null(attr(res,"status")))
          0
        else
          as.numeric(res)
      }
  )
}

#'     Get R processes CPU usage
#'
#'  Uses 'tasklist' to count the number of `Rscript.exe` instances
#'
#' Modified from https://stackoverflow.com/questions/47318401/r-how-to-check-how-many-cores-cpu-usage-available
#' @export
getCpuUsageRprocs <- function(targets=c("Rscript.exe","rsession.exe","rstudio.exe","Rgui.exe", "R")){

  switch(Sys.info()[['sysname']],
         Windows = {
           procs <- system("wmic path Win32_PerfFormattedData_PerfProc_Process get Name,PercentProcessorTime", intern = TRUE)
           procs %<>%
             strsplit(" ") %>%
             lapply(function(x) {x <- x[x != ""];data.frame(process = x[1], cpu = x[2])}) %>%
             do.call(rbind,.) %>%
             dplyr::filter(process %in% targets)

           if(length(procs)>0L)
             procs %>% dplyr::pull(cpu) %>% sum()
           else
             0
         },
         Linux   =  {
           user_name <- Sys.info()[["user"]]
           as.numeric(system(paste0("pgrep -c -u ",user_name," '",paste(targets,collapse="|"),"'"), intern = TRUE))
         }
  )
}

#'     Get the free available memory
#'
#'  Extract the amount of memory available on the `Windows` and `Linux` machines
#'
#' Slightly modified from https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
#' @export
getFreeMemory <- function(units='KB') {

  switch(Sys.info()[["sysname"]],
    Windows={
      x <- system2("wmic", args = "OS get FreePhysicalMemory /Value", stdout = TRUE)
      x <- x[grepl("FreePhysicalMemory", x)]
      x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
      x <- gsub("\r", "", x, fixed = TRUE)
      switch(units, KB = as.integer(x), MB=as.integer(x)/1024, GB = as.integer(x)/1024/1024)
    },
    Linux= {
      x <- system2('free', args='-k', stdout=TRUE)
      x <- strsplit(x[2], " +")[[1]][4]
      switch(units, KB = as.integer(x), MB=as.integer(x)/1024, GB = as.integer(x)/1024/1024)
    }
  )
}

getMemUsageRprocess <- function(units='KB', targets=c("Rscript.exe","rsession.exe","rstudio.exe","Rgui.exe", "R")){

  switch(Sys.info()[['sysname']],
         Windows = {
           memtasktab <- system("tasklist", intern = TRUE) %>%
             stringr::str_split("[[:space:]]") %>%
             `[`(-c(1:3)) %>%
             sapply(function(x) `[`(x, c(1,length(x)-1)), simplify = FALSE) %>%
             do.call(rbind,.) %>%
             as.data.frame(row.names=NULL) %>%
             setNames(c("Procs","MemUse")) %>%
             dplyr::mutate(MemUse=as.numeric(gsub(",","",MemUse)))

           mem_usage <- memtasktab %>%
             dplyr::filter(Procs %in% targets) %>%
             dplyr::pull(MemUse) %>%
             sum()
           switch(units, KB=mem_usage, MB=mem_usage/1024, GB=mem_usage/1024/1024)
         },
         Linux= {
           user_name <- Sys.info()[["user"]]
           memtasktab <- system(
             paste0("ps -u ",user_name," -o user,%cpu,%mem,stat,cmd | grep '", paste(gsub("\\.exe","",targets),collapse="|"),
                    "' | grep -v grep"),
             intern = TRUE)
           if(length(memtasktab)==0L)
             mem_usage = 0
           else
             mem_usage <- memtasktab %>%
             strsplit(" ") %>%
               `[`(-c(1)) %>%
             lapply(function(x) data.frame())
           switch(units, KB = as.integer(x), MB=as.integer(x)/1024, GB = as.integer(x)/1024/1024)
         }
  )
}


#' Unregister a parallel backend
#'
#' Provides a mechanism to unregister a foreach backend
#'
#' https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
#' @export
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


#'    Count the number of decimal
#' Modified from https://stat.ethz.ch/pipermail/r-help/2012-July/317676.html
#' @param x A character string or numeric vector.
#' @export
decimalnumcount<-Vectorize(function(x){
  if(inherits(x,"numeric"))
    x<- as.character(x)

  if(!grepl(".",x,fixed=TRUE))
    return(0)

  x<-gsub("(.*)?(?>\\.)|([0]*$)","",x,perl=T)
  nchar(x, type="width")
})

# A helper function that erases all of y from x:
# st_erase = function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y)))

#' utility function: find the minimum position for each row of a matrix, breaking ties at random.
#' @export
min.col <- function(m, ...) max.col(-m,...)


#' helper function that split a vector x in n chunks of approximatively equal size
#' @export
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


#' utility function: map function that transforms input by applying a function to each element and returning
#' a data frame with a number of columns equal to the number of element in the input.
#' @export
map_lfd <- function(x, fun, ...){
  num_index <- length(x)
  ret <- lapply(1:num_index, FUN=function(i) fun(`[`(x,i)), ...)
  res <- as.data.frame(do.call(cbind, ret))
  names(res) <- names(x)
  res
}


#' utility function: strip filename extension.
#' @export
strip_extension <- function(fn){
  stringr::str_extract(fn, ".+?(?=\\.[a-z]{3}$)")
}


#'              Most frequent value(s)
#'
#' Find the most frequent n value(s) from a vector
#'
#' @param x A vector
#' @param n A numeric integer specifying the number of values to return. Default is 1, i.e. the most frequent value.
#' @param ... Additional parameters to be passed to the function \code{table}.
#' @seealso \code{table}
#' @export
modal <- function(x, n=1, ...){
  sort(table(x,...),decreasing=TRUE)[1:n]
}



#' Remove temporary raster files
#' @export
removeTMPFiles <- function(h=24) {

  # remove files in the temp folder that are > h hours old
  warnopt <- getOption('warn')
  on.exit(options('warn'= warnopt))

  tmpdir <- raster::tmpDir(create=FALSE)
  if (!is.na(tmpdir)) {

    d <- raster:::.removeTrailingSlash(tmpdir)
    f <- list.files(path=d, pattern='r_tmp*', full.names=TRUE, include.dirs=TRUE)
    #		f <- list.files(path=d, pattern='[.]gr[di]', full.names=TRUE, include.dirs=TRUE)
    fin <- file.info(f)
    dif <- Sys.time() - fin$mtime
    dif <- as.numeric(dif, units="hours")

    f <- f[which(dif > h)]
    unlink(f, recursive=TRUE)
  }
  options('warn'=warnopt)
}

