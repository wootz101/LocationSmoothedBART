# Aug15th Cervix functional Dataset with correct dimensions

library(EBImage)
library(oro.dicom)
library(tidyverse)
library(glcm) 
source("~/Desktop/shapeHist2_July15.R")
setwd("/Users/Wootz/Desktop/ORGAN DATA")
# Test new shape hist function

#folderName <- "~/Desktop/pyfun/Good_Kidney2/"
#g_Files <- list.files("~/Desktop/pyfun/Good_Kidney2/")
#folderName_oneDrive <- "~/Desktop/cervix_dicoms/"




folderName <- "~/Desktop/pyfun/Good_Bladder/"
g_Files <- list.files("~/Desktop/pyfun/Good_Bladder/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

test1_output = shape_function_metadata()


test_output = shapeHist(g_Files, folderName, folderName_oneDrive)

shapeHist(g_Files[17], folderName, folderName_oneDrive)

one_drive_Files <- list.files(folderName_oneDrive)
subname = substr(g_Files[17], start = nchar(g_Files[k]) - 20, stop = nchar(g_Files[k]) - 4)
matches <- sapply(one_drive_Files, grepl, pattern = subname)
dicom_file = one_drive_Files[matches]
print(dicom_file)
metadata_files <- list.files(paste0(folderName_oneDrive, dicom_file))
print(metadata_files[2])
#print(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[2] ))
metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[1] ))



# Get all the metadata needed
metadata <- NULL
metadata_return = NULL
for (i in 1:36) {
  one_drive_Files <- list.files(folderName_oneDrive)
  dicom_file = one_drive_Files[i]
  print(dicom_file)
  
  subname = substr(dicom_file, start = 0, stop= nchar(dicom_file) -11)
  metadata_files <- list.files(paste0(folderName_oneDrive, dicom_file))
  metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[1] ))
  #print(metadata$hdr$name)
  sliceThickness = as.numeric(metadata$hdr$value[metadata$hdr$name == "SliceThickness"])
  ps <- metadata$hdr$value[metadata$hdr$name == "PixelSpacing"]
  pixelSpace = as.numeric(strsplit(ps, " ")[[1]][1])
  
  print(paste("Thickness: ", sliceThickness, "Pixel Spacing: ", pixelSpace ))
  metadata_return = rbind(metadata_return, c(subname, sliceThickness,  pixelSpace))
  rm(metadata)
}

metadata_return1_36 = metadata_return

metadata <- NULL
metadata_return = NULL
for (i in 38:87) {
  one_drive_Files <- list.files(folderName_oneDrive)
  dicom_file = one_drive_Files[i]
  print(dicom_file)
  subname = substr(dicom_file, start = 0, stop= nchar(dicom_file) -11)
  metadata_files2 <- list.files(paste0(folderName_oneDrive, dicom_file))
  metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[1] ))
  #print(metadata$hdr$name)
  sliceThickness = as.numeric(metadata$hdr$value[metadata$hdr$name == "SliceThickness"])
  ps <- metadata$hdr$value[metadata$hdr$name == "PixelSpacing"]
  pixelSpace = as.numeric(strsplit(ps, " ")[[1]][1])
  
  print(paste("Thickness: ", sliceThickness, "Pixel Spacing: ", pixelSpace ))
  metadata_return = rbind(metadata_return, c(subname, sliceThickness,  pixelSpace))
  rm(metadata)
}
metadata_return38_87 = metadata_return


metadata <- NULL
metadata_return = NULL
for (i in 89:110) {
  one_drive_Files <- list.files(folderName_oneDrive)
  dicom_file = one_drive_Files[i]
  print(dicom_file)
  subname = substr(dicom_file, start = 0, stop= nchar(dicom_file) -11)
  metadata_files2 <- list.files(paste0(folderName_oneDrive, dicom_file))
  metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[1] ))
  #print(metadata$hdr$name)
  sliceThickness = as.numeric(metadata$hdr$value[metadata$hdr$name == "SliceThickness"])
  ps <- metadata$hdr$value[metadata$hdr$name == "PixelSpacing"]
  pixelSpace = as.numeric(strsplit(ps, " ")[[1]][1])
  
  print(paste("Thickness: ", sliceThickness, "Pixel Spacing: ", pixelSpace ))
  metadata_return = rbind(metadata_return, c(subname, sliceThickness,  pixelSpace))
  rm(metadata)
}
metadata_return89_110 = metadata_return


metadata <- NULL
metadata_return = NULL
for (i in 112:140) {
  one_drive_Files <- list.files(folderName_oneDrive)
  dicom_file = one_drive_Files[i]
  print(dicom_file)
  subname = substr(dicom_file, start = 0, stop= nchar(dicom_file) -11)
  metadata_files2 <- list.files(paste0(folderName_oneDrive, dicom_file))
  metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[1] ))
  #print(metadata$hdr$name)
  sliceThickness = as.numeric(metadata$hdr$value[metadata$hdr$name == "SliceThickness"])
  ps <- metadata$hdr$value[metadata$hdr$name == "PixelSpacing"]
  pixelSpace = as.numeric(strsplit(ps, " ")[[1]][1])
  
  print(paste("Thickness: ", sliceThickness, "Pixel Spacing: ", pixelSpace ))
  metadata_return = rbind(metadata_return, c(subname, sliceThickness,  pixelSpace))
  rm(metadata)
}
metadata_return112_140 = metadata_return



metadata <- NULL
metadata_return = NULL
for (i in c(37,88,111) ) {
  one_drive_Files <- list.files(folderName_oneDrive)
  dicom_file = one_drive_Files[i]
  print(dicom_file)
  subname = substr(dicom_file, start = 0, stop= nchar(dicom_file) -11)
  metadata_files2 <- list.files(paste0(folderName_oneDrive, dicom_file))
  metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[3] ))
  #print(metadata$hdr$name)
  sliceThickness = as.numeric(metadata$hdr$value[metadata$hdr$name == "SliceThickness"])
  ps <- metadata$hdr$value[metadata$hdr$name == "PixelSpacing"]
  pixelSpace = as.numeric(strsplit(ps, " ")[[1]][1])
  
  print(paste("Thickness: ", sliceThickness, "Pixel Spacing: ", pixelSpace ))
  metadata_return = rbind(metadata_return, c(subname, sliceThickness,  pixelSpace))
  rm(metadata)
}
metadata_return_rand = metadata_return

metadata_return1_36
metadata_return_rand[1,]
metadata_return38_87
metadata_return_rand[2,]
metadata_return89_110
metadata_return_rand[3,]
metadata_return112_140

rbind(metadata_return1_36,
      metadata_return_rand[1,],
      metadata_return38_87,
      metadata_return_rand[2,],
      metadata_return89_110,
      metadata_return_rand[3,],
      metadata_return112_140)


all_meta=as.data.frame(rbind(metadata_return1_36,
                    metadata_return_rand[1,],
                    metadata_return38_87,
                    metadata_return_rand[2,],
                    metadata_return89_110,
                    metadata_return_rand[3,],
                    metadata_return112_140))
all_meta$V2 = as.numeric(all_meta$V2)
all_meta$V3 = as.numeric(all_meta$V3)
names(all_meta) = c("name", "sliceThickness", "pixelSpacing")

write.csv(all_meta, "Cervix_Metadata.csv", row.names = FALSE)

one_drive_Files <- list.files(folderName_oneDrive)
subname = substr(g_Files[k], start = nchar(g_Files[k]) - 20, stop = nchar(g_Files[k]) - 4)
matches <- sapply(one_drive_Files, grepl, pattern = subname)
#print(one_drive_Files[matches])
dicom_file = one_drive_Files[matches]
print(dicom_file)

#############################################################################################
#### Get the metadata from the chosen file
#paste0(folderName_oneDrive, dicom_file)
#metadata <- readDICOM(paste0(folderName_oneDrive, dicom_file), verbose = TRUE,
#                 recursive = FALSE)

metadata_files <- list.files(paste0(folderName_oneDrive, dicom_file))
#print(metadata_files[2])
#print(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[2] ))
metadata <- readDICOMFile(paste0(folderName_oneDrive, dicom_file,'/', metadata_files[2] ))










folderName <- "~/Desktop/pyfun/Good_Bladder/"
g_Files <- list.files("~/Desktop/pyfun/Good_Bladder/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

cervix_metadata <- read_csv("/Users/Wootz/Desktop/ORGAN DATA/Cervix_metadata.csv")

test_output = shapeHist_metadata(g_Files, folderName, cervix_metadata)

test_output_moreFunctions = shape_function_metadata(g_Files[1], folderName, cervix_metadata)


#test_output = shapeHist(g_Files, folderName, folderName_oneDrive)
test_output

dim(test_output[[2]])

# Get the numeric values 
as.numeric(test_output[[2]][ ,3:13])

matrix(as.numeric(test_output[[2]]),
       nrow = dim(test_output[[2]])[1],
       ncol = dim(test_output[[2]])[2])

matrix(as.numeric(test_output[[2]][ ,3:13]),
       nrow = dim(test_output[[2]][ ,3:13])[1],
       ncol = dim(test_output[[2]][ ,3:13])[2])



test_output[[1]]
dim(test_output[[1]])

test_output[[2]]
test_output[[3]]


write.csv(test_output[[1]], "GoodBladder_hist_Sept4.csv", row.names = FALSE)
write.csv(test_output[[2]], "GoodBladder_2Dfunc_Sept4.csv", row.names = FALSE)
write.csv(test_output[[3]], "GoodBladder_metadata_Sept4.csv", row.names = FALSE)
##############################################################################



folderName <- "~/Desktop/pyfun/Bad_Bladder/"
g_Files <- list.files("~/Desktop/pyfun/Bad_Bladder/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

bad_output = shapeHist_metadata(g_Files, folderName, cervix_metadata)



bad_output


bad_output[[1]]
bad_output[[2]]
bad_output[[3]]


write.csv(bad_output[[1]], "BadBladder_hist_Sept4.csv", row.names = FALSE)
write.csv(bad_output[[2]], "BadBladder_2Dfunc_Sept4.csv", row.names = FALSE)
write.csv(bad_output[[3]], "BadBladder_metadata_Sept4.csv", row.names = FALSE)
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
#Rectum

folderName <- "~/Desktop/pyfun/good_rectum/"
g_Files <- list.files("~/Desktop/pyfun/good_rectum/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

test_output = shapeHist_metadata(g_Files, folderName, cervix_metadata)




dim(test_output[[2]])



test_output[[1]]
dim(test_output[[1]])

test_output[[2]]
test_output[[3]]


write.csv(test_output[[1]], "good_rectum_hist_Sept4.csv", row.names = FALSE)
write.csv(test_output[[2]], "good_rectum_2Dfunc_Sept4.csv", row.names = FALSE)
write.csv(test_output[[3]], "good_rectum_metadata_Sept4.csv", row.names = FALSE)
##############################################################################



folderName <- "~/Desktop/pyfun/bad_rectum/"
g_Files <- list.files("~/Desktop/pyfun/bad_rectum/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

bad_output = shapeHist_metadata(g_Files, folderName, cervix_metadata)



bad_output


bad_output[[1]]
bad_output[[2]]
bad_output[[3]]


write.csv(bad_output[[1]], "bad_rectum_hist_Sept4.csv", row.names = FALSE)
write.csv(bad_output[[2]], "bad_rectum_2Dfunc_Sept4.csv", row.names = FALSE)
write.csv(bad_output[[3]], "bad_rectum_metadata_Sept4.csv", row.names = FALSE)



############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
#PAN



folderName <- "~/Desktop/pyfun/good_PAN/"
g_Files <- list.files("~/Desktop/pyfun/good_PAN/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

test_output = shapeHist_metadata(g_Files, folderName, cervix_metadata)




dim(test_output[[2]])



test_output[[1]]
dim(test_output[[1]])

test_output[[2]]
test_output[[3]]


write.csv(test_output[[1]], "good_PAN_hist_Sept4.csv", row.names = FALSE)
write.csv(test_output[[2]], "good_PAN_2Dfunc_Sept4.csv", row.names = FALSE)
write.csv(test_output[[3]], "good_PAN_metadata_Sept4.csv", row.names = FALSE)
##############################################################################



folderName <- "~/Desktop/pyfun/bad_PAN/"
g_Files <- list.files("~/Desktop/pyfun/bad_PAN/")
folderName_oneDrive <- "/Users/Wootz/Library/CloudStorage/OneDrive-InsideMDAnderson/cervixfinalst/"
one_drive_Files <- list.files(folderName_oneDrive)
one_drive_Files

bad_output = shapeHist_metadata(g_Files, folderName, cervix_metadata)



bad_output


bad_output[[1]]
bad_output[[2]]
bad_output[[3]]


write.csv(bad_output[[1]], "bad_PAN_hist_Sept4.csv", row.names = FALSE)
write.csv(bad_output[[2]], "bad_PAN_2Dfunc_Sept4.csv", row.names = FALSE)
write.csv(bad_output[[3]], "bad_PAN_metadata_Sept4.csv", row.names = FALSE)






