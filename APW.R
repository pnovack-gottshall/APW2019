## Analytical Paleobiology Workshop GitHub repository
## June 28-30, 2019, Riverside, CA


## GitHub recs:

# From http://happygitwithr.com/rstudio-git-github.html

# Now is the perfect time to make sure you don't need to keep providing username
# and password on each push. Pick one:
#    Credential caching for HTTPS access, chapter 11 at http://happygitwithr.com/credential-caching.html#credential-caching
#    Set up SSH keys, chapter 12. at http://happygitwithr.com/ssh-keys.html#ssh-keys



## Phil Novack-Gottshall <pnovack-gottshall@ben.edu> ---------------------------

## Function to quickly download base, midpoint, and top ages for PBDB intervals.
## Useful for plotting time series.
#  scale_level: eons = 1, eras = 2, periods = 3, subperiods = 4, & epochs = 5
PBDB_interval_ages <- function(scale_level = 4) {
  strat_names <-
    read.csv("https://www.paleobiodb.org/data1.2/intervals/list.csv?all_records&vocab=pbdb")
  intervals <-
    strat_names[which(strat_names$scale_level == scale_level), ]
  mids <- apply(intervals[, 9:10], 1, mean)
  return(data.frame(interval = intervals$interval_name, base = intervals$max_ma,
                    midpoint = mids, top = intervals$min_ma))
}

# Examples:
PBDB_interval_ages()
PBDB_interval_ages(scale_level = 1)

## Not run:
require(geoscale)
require(paleoTS)
ages <- PBDB_interval_ages(scale_level = 5)
sim <- sim.GRW(ns = length(ages$midpoint), ms = -0.1, vs = 0.5)
geoscalePlot(ages$midpoint, sim$mm, units = c("Epoch", "Period"), tick.scale = "Period", boxes = "Epoch", cex.age = 0.65, cex.ts=0.75, cex.pt = 1,
             age.lim = c(540, 0), ts.col = TRUE, vers = "ICS2015", type = "l", abbrev = "Period")



## Function to add higher taxonomic names (phylum, class, order, etc.) for PBDB
## genus (and subgenus) names
# gen.seq = Vector (sequence) of number of genus names to process
# gen.names = Vector of PBDB genus (and subgenus) names
# pbdb = data frame of all PBDB occurrences (downloaded from PBDB)
#
# output is a list, with each item the taxonomy for a single genus
prep.PBDB <- function(gen.seq = 1, gen.names, pbdb) {
  scales <- c("phylum", "subphylum", "class", "subclass", "order", "suborder",
              "superfamily", "family", "subfamily", "genus", "subgenus")
  out <- data.frame(Phylum = character(1), Subphylum = character(1),
                    Class = character(1), Subclass = character(1),
                    Order = character(1), Suborder = character(1),
                    Superfamily = character(1), Family = character(1),
                    Subfamily = character(1), Genus = character(1),
                    Subgenus = character(1), Species = "sp.",
                    stringsAsFactors = FALSE)
  out$Genus <- as.character(gen.names[g])
  wh <- which(pbdb$accepted_name == out$Genus & (pbdb$taxon_rank == "genus" |
                                                   pbdb$taxon_rank == "subgenus"))[1]
  out$early_age <- as.numeric(pbdb$firstapp_max_ma[wh])
  out$late_age <- as.numeric(pbdb$lastapp_min_ma[wh])
  if (pbdb$accepted_rank[wh] == "subgenus")
    out$Subgenus <- as.character(out$Genus)
  parent <- pbdb[which(pbdb$accepted_no == pbdb$parent_no[wh]), ][1, ]
  repeat {
    if (parent$accepted_rank %in% scales)
      out[1, which(scales == parent$accepted_rank)] <-
        as.character(parent$accepted_name)
    parent <-
      pbdb[which(pbdb$accepted_no == parent$parent_no),][1,]
    if (all(is.na(parent)))
      break
  }
  return(out)
}

## Function to unpack the list into a compact data frame
# Note that with a large PBDB object that this still takes some time to
# complete. Consider tweaking to speed up.
unpack.PBDB <- function(prep) {
  Seq <- seq(length(prep))
  out <- prep[[1]]
  for (i in 2:length(prep)) {
    out[i, ] <- prep[[i]]
  }
  return(out)
}

# Examples:
## Not run:
# Following for all PBDB animals:
# pbdb <- read.csv("www.paleobiodb.org/data1.2/taxa/list.csv?base_name=Metazoa&show=app&vocab=pbdb")
# Faster example, all trilobites:
pbdb <- read.csv("www.paleobiodb.org/data1.2/taxa/list.csv?base_name=Metazoa&show=app&vocab=pbdb")
head(pbdb)
# If above read.csv() calls do not work, paste URL above into web browser to
# download, then save file to your working directory and read directly.

# Examples here run through the first 1000 genera, using:
num.gen <- 1000 # all would be nrow(pbdb)

# Serial version
(start.t <- Sys.time())
# Use genera and subgenera
which.gsg <- which(pbdb$accepted_rank=="genus" | pbdb$accepted_rank=="subgenus")
n.gen <- length(unique(pbdb$accepted_name[which.gsg]))
x <- data.frame(Phylum=character(n.gen), Subphylum=character(n.gen),
                Class=character(n.gen), Subclass=character(n.gen), Order=character(n.gen),
                Suborder=character(n.gen), Superfamily=character(n.gen), Family=character(n.gen),
                Subfamily=character(n.gen), Genus=sort(unique(pbdb$accepted_name[which.gsg])),
                Subgenus=character(n.gen), Species=rep("sp.", n.gen), early_age=numeric(n.gen),
                late_age=numeric(n.gen), stringsAsFactors=FALSE)
scales <- c("phylum", "subphylum", "class", "subclass", "order", "suborder",
            "superfamily", "family", "subfamily", "genus", "subgenus")
for(g in 1:num.gen) {
  wh <- which(pbdb$accepted_name == x$Genus[g] & (pbdb$taxon_rank == "genus" |
                                                    pbdb$taxon_rank == "subgenus"))[1]
  x$early_age[g] <- pbdb$firstapp_max_ma[wh]
  x$late_age[g] <- pbdb$lastapp_min_ma[wh]
  if(pbdb$accepted_rank[wh] == "subgenus") x$Subgenus[g] <- as.character(x$Genus[g])
  parent <- pbdb[which(pbdb$accepted_no == pbdb$parent_no[wh]), ][1, ]
  repeat {
    if(parent$accepted_rank %in% scales)
      x[g, which(scales == parent$accepted_rank)] <- as.character(parent$accepted_name)
    parent <- pbdb[which(pbdb$accepted_no == parent$parent_no), ][1, ]
    if(all(is.na(parent))) break
  }
}
(Sys.time() - start.t)
head(x)

# Serial version takes 232.3 seconds for 1000 (or ~4.7 hours for all)

# Re-do using lapply (only 1 CPU)
(start.t <- Sys.time())
which.gsg <- which(pbdb$accepted_rank=="genus" | pbdb$accepted_rank=="subgenus")
gen.names <- sort(unique(pbdb$accepted_name[which.gsg]))
# gen.seq <- seq_along(gen.names)
gen.seq <- 1:num.gen  # for all, use 1:nrow(pbdb)
prep <- lapply(X = gen.seq, FUN = prep.PBDB, gen.names = gen.names, pbdb = pbdb)
output <- unpack.PBDB(prep)
(Sys.time() - start.t)
# Serial lapply version takes 41.2 seconds for 1000 (or ~49 mins for all)


# Now compare using parallel computing:
library(snowfall)
(t.start0 <- Sys.time())
# Initialize
which.gsg <- which(pbdb$accepted_rank=="genus" | pbdb$accepted_rank=="subgenus")
gen.names <- sort(unique(pbdb$accepted_name[which.gsg]))
# gen.seq <- seq_along(gen.names)
gen.seq <- 1:num.gen     # for all, use seq_along(gen.names)
# Set up computer cluster
cpus <- 8					# Number of CPUs to cluster together
# sfSetMaxCPUs(cpus)			# Use if plan more than 32 CPUs
sfInit(parallel=T, cpus=cpus, slaveOutfile="initfile") # Initialize cluster
stopifnot(sfCpus()==cpus)		# Confirm set up CPUs properly
stopifnot(sfParallel()==TRUE)		# Confirm now running in parallel
sfExportAll()				# Export all libraries, files, & objects
# Execute the function
(t.start1 <- Sys.time())
prep <- NA
prep <- sfLapply(x = gen.seq, fun = prep.PBDB, gen.names = gen.names, pbdb = pbdb) # Version without load-balancing
# prep <- sfClusterApplyLB(x = gen.seq, fun = prep.PBDB, gen.names = gen.names, pbdb = pbdb) 	# Version using load-balancing
sfStop()
(Sys.time() - t.start1)
output2 <- unpack.PBDB(prep)
(Sys.time() - t.start0)
(Sys.time() - t.start1)

# Parallel lapply version takes 60.2 secs for 1000 (21 s for the sfLapply) (or 27.8 min for all)
# Parallel lapply version (w/o export time) takes 19.7 secs for 1000 (or 23.7 mins for all)
# Parallel lapplyLB version takes 642 secs for 1000 (or 12.8 hrs for all)
# Parallel lapplyLB version (w/o export time) takes 602 secs for 1000 (or 12.0 hrs for all)

# Compare output
rbind(x[1000, ], output[1000, ], output2[1000, ])





## Jordan Claytor ----------------------------------------------------------


## YiYing Deng -------------------------------------------------------------


## Erin Dillon -------------------------------------------------------------

#Hi Angelina! How's it going? I wrote you a piece of code.
library(car)
data("mtcars")
head(mtcars)
boxplot(mpg~cyl,data=mtcars)

## Emma Dunne --------------------------------------------------------------


## Angelina Ivkic ----------------------------------------------------------


## Mark Juhn ---------------------------------------------------------------


## Bridget Kelly -----------------------------------------------------------


## Miranta Kouvari ---------------------------------------------------------


## Florencia Leone -------------------------------------------------------------


## Ekaterina Larina --------------------------------------------------------


## Pablo S. Milla Carmona --------------------------------------------------


## Selina Robson -----------------------------------------------------------


## Anwesha (Mimi) Saha -----------------------------------------------------


## Jaime Andres Villafaña Navea --------------------------------------------


## Tom Womack --------------------------------------------------------------


## Michelle Zill -----------------------------------------------------------


