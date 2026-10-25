# Lecture et validation des données originales
path <- "data/leukemia.csv"
if (!file.exists(path)) stop("Placer leukemia.csv dans data/ (usage local uniquement).")
raw <- readLines(path, warn=FALSE)
raw <- sub('^"', '', raw)
raw <- sub('"$', '', raw)
raw <- gsub('""', '"', raw, fixed=TRUE)
d <- read.csv(text=paste(raw,collapse="\n"), na.strings="NA", stringsAsFactors=FALSE)
if ("X" %in% names(d)) d$X <- NULL
needed <- c("patid","srv","srvstat","rel","relstat","yrel","age","score")
stopifnot(all(needed %in% names(d)), nrow(d)==1977L, !anyDuplicated(d$patid))
stopifnot(all(d$srvstat %in% c(0,1)), all(d$relstat %in% c(0,1)))
stopifnot(all(d$srv >= 0), all(d$rel >= 0))
stopifnot(all(d$rel[d$relstat==1] <= d$srv[d$relstat==1]))
# yrel est structurellement manquant pour les personnes sans rechute.
adults <- subset(d, age >= 18 & complete.cases(d[,c("srv","srvstat","rel","relstat","age","score")]))
stopifnot(nrow(adults)==1835L)
# Priorité à la rechute observée ; sinon décès observé ; sinon censure.
# Dans ces données, rel==srv pour les décès sans rechute.
make_first <- function(x) {
 x$time <- ifelse(x$relstat==1, x$rel, x$srv)
 x$status <- ifelse(x$relstat==1, 1L, ifelse(x$srvstat==1,2L,0L))
 x$type <- factor(c("event_free","relapse","death")[x$status+1L],levels=c("event_free","relapse","death"))
 x
}
d <- make_first(d); adults <- make_first(adults)
print(table(d$type)); print(table(adults$type)); print(colSums(is.na(d)))
stopifnot(identical(as.integer(table(d$type)), c(836L,456L,685L)))
saveRDS(d,"data/leucemie.rds"); saveRDS(adults,"data/leucemie2.rds")
