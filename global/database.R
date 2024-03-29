# Study class
setClass("Study1", representation(
  studyName="character",
  species="character",
  MCMC="matrix"))

setGeneric("meta", function(object) {
  standardGeneric("meta")
})

setMethod("meta", signature("Study1"), function(object) {
  data.frame(species=object@species, studyName=object@studyName)
})


# Database class
setClass("Database",
         representation(
           name="character",
           dir="character",
           meta.file="character",
           working.path="character"
         ),
         prototype(
           dir="",
           meta.file="",
           working.path=""
         ),
         validity = function(object) {
           errors <- character()
           if (length(errors) == 0) TRUE else errors
         }
)

# custom contructor to set meta
setMethod("initialize", "Database",
          function(.Object, name) {
            .Object <- callNextMethod()
            .Object@dir      <- paste(DB.dir, name, sep="/")
            .Object@meta.file <- paste(".Database", name, "meta/meta", sep="/")
            .Object@working.path <- paste(".Database", name, "working/working", sep="/")
            if (!file.exists(.Object@dir)) dir.create(.Object@dir, recursive=T)
            dir.path <- paste(".Database", name, "meta",  sep="/")
            if (!file.exists(dir.path)) dir.create(dir.path, recursive=T)
            dir.path <- paste(".Database", name, "working",  sep="/")
            if (!file.exists(dir.path)) dir.create(dir.path, recursive=T)
            if (!file.exists(.Object@working.path))
              file.create(.Object@working.path, overwrite=F)
            studies <- c()
            studies <- DB.load(.Object, list.files(path=.Object@dir))
            print(paste("number of studies saved: " ,length(studies), sep=""))
            db.meta <- data.frame(
              "species"=character(0),
              "studyNames"=character(0)
            )
            if(length(studies) > 0) {
              db.meta <- lapply(studies, function(study_use) meta(study_use))
              db.meta <- do.call(rbind, db.meta)
            }
            DB.sync(.Object, db.meta)
            .Object
          }
)

# Return Database meta information as data.frame
setMethod("meta", signature("Database"), function(object) {
  readRDS(db@meta.file)
})

# write database meta data to file, should be called everytime when
# database is modified
DB.sync <- function(db, db.meta) {
  saveRDS(db.meta, file=db@meta.file)
}

# save x to db as file
DB.save <- function(db, study_use) {
  #  if(class(study) != "Study") stop("study must be Study")
  saveRDS(study_use, file=paste(db@dir, study_use@studyName, sep="/"))
  studies <- DB.load(db, list.files(path=db@dir))
  db.meta <- lapply(studies, function(study_use) meta(study_use))
  db.meta <- do.call(rbind, db.meta)
  DB.sync(db, db.meta)
}

# load file from db
DB.load <- function(db, studies) {
  res <- c()
  for(study in studies) {
    res <- c(res, readRDS(paste(db@dir, study, sep='/')))
  }
  res
}

MergedDB.load <- function(db){
  if(file.exists(paste(DB.load.working.dir(db), "MergedDB.rds", sep="/"))){
    readRDS(paste(DB.load.working.dir(db), "MergedDB.rds", sep="/"))
  }else{
    return(NA)
  }
}

MergedSpecies.load <- function(db){
  if(file.exists(paste(DB.load.working.dir(db), 
                       "MergedSpecies.rds", sep="/"))){
    readRDS(paste(DB.load.working.dir(db), 
                  "MergedSpecies.rds", sep="/"))
  }else{
    return(NA)
  }
}

MergedStudyNames.load <- function(db){
  if(file.exists(paste(DB.load.working.dir(db), 
                       "MergedStudyNames.rds", sep="/"))){
    readRDS(paste(DB.load.working.dir(db), 
                  "MergedStudyNames.rds", sep="/"))
  }else{
    return(NA)
  }
}

MergedPM.load <- function(db){
  if(file.exists(paste(DB.load.working.dir(db), 
                       "MergedPM.rds", sep="/"))){
    readRDS(paste(DB.load.working.dir(db), 
                  "MergedPM.rds", sep="/"))
  }else{
    return(NA)
  }
}

# delete file from db
DB.delete <- function(db, studies) {
  file.remove(paste(db@dir, list.files(path=db@dir)[as.numeric(studies)], sep="/"))
  db.meta <- meta(db)
  db.meta <- db.meta[!(rownames(db.meta) %in% studies),]
  if(nrow(db.meta)>0){
    rownames(db.meta) = 1:nrow(db.meta)
  }
  DB.sync(db, db.meta)
}

# list all files in db
DB.ls <- function(db) {
  meta(db)$studyName
}

DB.set.working.dir <- function(db, path){
  file.con <- file(db@working.path)
  writeLines(path, file.con)
  close(file.con)
}

DB.load.working.dir <- function(db){
  file.con <- file(db@working.path)
  path <- readLines(file.con)
  close(file.con)
  if (length(path) == 0)
    stop(MSG.no.working.dir)
  else
    path
}

