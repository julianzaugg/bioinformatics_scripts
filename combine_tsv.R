# Combine all tsv files in a directory into a single excel spreadsheet

library(openxlsx)


file.path2 = function(..., fsep = .Platform$file.sep){
  
  out <- gsub("//", "/", file.path(..., fsep = fsep))
}

options <- commandArgs(trailingOnly = TRUE)


in_dir <- options[1]
out_dir <- options[2]
out_name <- options[3]

if (!file.exists(in_dir)){
  stop("Invalid input directory")
}

setwd(in_dir)
in_dir <- getwd()

if (is.na(out_dir)){
  out_dir <- normalizePath(".")
}

if (!file.exists(out_dir)){
  stop("Invalid output directory.")
}
out_dir <- normalizePath(out_dir)

if (is.na(out_name)){
  out_name <- "combined"
}

#print(out_dir)

tsv_filenames <- as.character(list.files(in_dir,"\\.tsv"))
wb <- createWorkbook()
# sheet_num <- 1
for (filename in tsv_filenames){
  file_path <- paste(in_dir,filename, sep ="/")
  read_try <- try(read.table(file_path, sep = "\t", header = T), silent = T)
  if(class(read_try) == "try-error"){
    print(paste("Could not read", file_path))
    next
  }else
    mydata <- read.table(file_path, sep = "\t", header = T)
    reduced_name <- gsub(".tsv", "", filename)
    addWorksheet(wb, sheetName = reduced_name)
    writeData(wb, sheet = reduced_name, x = mydata)
    # sheet_num <- sheet_num + 1
}

out_spreadsheet_path <- paste(out_dir,"/",out_name,".xlsx", sep = "")
saveWorkbook(wb, file = out_spreadsheet_path)
