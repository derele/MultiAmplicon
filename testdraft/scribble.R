## with hyena data
foo <- MA[, which(colnames(MA)%in%c("2018_22_hyena_-S1045-P1-FLD0007",
                                    "2018_22_hyena_-S1074-P1-FLD0036",
                                    "2018_22_hyena_-S1066-P1-FLD0028"))]

readsFL <- lapply(seq_along(foo@PairedReadFileSet), function(i) {
    rawFiles <- foo@PairedReadFileSet@readsF[[i]]
    rawFiles <- rawFiles[file.exists(rawFiles)]
    readFastq(rawFiles)
})

snames <- colnames(foo)
names(readsFL) <- snames

readsF <- readsFL[["2018_22_hyena_-S1045-P1-FLD0007"]]

strat <- lapply(foo@stratifiedFiles, function(x) {
    grep("2018_22_hyena_-S1045-P1-FLD0007", x@readsF, value=TRUE)
})

readsF_stratified <- readFastq(unlist(strat))

length(readsF_stratified[!ShortRead::id(readsF_stratified) %in%
                         ShortRead::id(readsF)])


ShortRead::id(readsF_stratified)
ShortRead::id(readsF_stratified)[duplicated(ShortRead::id(readsF_stratified))]

sread(readsF_stratified)[duplicated(ShortRead::id(readsF_stratified))]


pastry <- paste(ShortRead::id(readsF_stratified), sread(readsF_stratified), collaps="_")

### OKAY really both read and id are duplicated
table(duplicated(pastry))
