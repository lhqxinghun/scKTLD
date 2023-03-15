################## GM12878 50kb
rm(list = ls())

library(Matrix)
source("/media/biology/datadisk/bkup_d1/liuerhu/hicda_liuerhu/seq-depth/funDownsample.R")
source("/media/biology/datadisk/bkup_d1/liuerhu/hicda_liuerhu/seq-depth/sparse2matrix.R")

resolution = 50000
## hg19
chrlen =c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 
          159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 
          115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 
          59128983, 63025520, 48129895, 51304566, 155270560)
  
## mm9
# chrlen = c(197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 
#            152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 
#            120284312, 125194864, 103494974, 98319150, 95272651, 90772031,
#            61342430, 166650296)

contact_target = 350000
chrname = paste0("chr",c(1:22, "X"))
contact = c(118858632, 192651361, 312979693, 100002514, 192751614, 640355881, 247935436, 332214022, 387488865,  66361476,  64148484,  75352536,  61665490,  82120761)
#for (c in 1:length(chrname))
for (c in 1:1)
{
  #c = 1
  print(paste0("chr",c))
  filelist = list.files(paste0("/media/biology/2b52abd3-7feb-4bc9-8c84-a8408d93f7ec/dump/GM12878/", chrname[c], "/50K"), full.names = TRUE)
  
  for (f in 6:length(filelist))
  {
    #f = 1
    print(paste0("file",f))
    file = filelist[f]
    sparsemat = read.table(file)
    
    n_sample = contact[f]/contact_target
    hicmat = sparse2matrix(sparsemat, floor(chrlen[c]/resolution)+1, resolution)
    hicmat_ds = funDownsample(hicmat, n_sample)
    
    sparsedf <- as.data.frame(summary(Matrix(hicmat_ds)))
    sparsedf[,1] <- sparsedf[,1] - 1
    sparsedf[,2] <- sparsedf[,2] - 1
    
    write.table(sparsedf, file = paste0("/media/biology/2b52abd3-7feb-4bc9-8c84-a8408d93f7ec/dump/GM12878/", chrname[c], "/50K_ds350K/", basename(file)), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}



################## GM12878 25kb
rm(list = ls())

library(Matrix)
source("/media/biology/datadisk/bkup_d1/liuerhu/hicda_liuerhu/seq-depth/funDownsample.R")
source("/media/biology/datadisk/bkup_d1/liuerhu/hicda_liuerhu/seq-depth/sparse2matrix.R")

resolution = 25000
reference = "hg19"

chrlen = NULL
if (reference == "hg19")
{
  chrlen =c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 
            159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 
            115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 
            59128983, 63025520, 48129895, 51304566, 155270560)
}
if (reference == "mm9")
{
  chrlen = c(197195432, 181748087, 159599783, 155630120, 152537259, 149517037,
               152524553, 131738871, 124076172, 129993255, 121843856, 121257530,
               120284312, 125194864, 103494974, 98319150, 95272651, 90772031,
               61342430, 166650296)
}

if(is.null(reference))
{
  stop("No correct referece genome selected")
}

contact_target = 350000
chrname = paste0("chr",c(1:22, "X"))
contact = c(118858632, 192651361, 312979693, 100002514, 192751614, 640355881, 247935436, 332214022, 387488865,  66361476,  64148484,  75352536,  61665490,  82120761)
#for (c in 1:length(chrname))
for (c in 1:1)
{
  #c = 1
  print(paste0("chr",c))
  filelist = list.files(paste0("/media/biology/2b52abd3-7feb-4bc9-8c84-a8408d93f7ec/dump/GM12878/", chrname[c], "/25K"), full.names = TRUE)
  
  for (f in 6:length(filelist))
  {
    #f = 1
    print(paste0("file",f))
    file = filelist[f]
    sparsemat = read.table(file)
    
    n_sample = contact[f]/contact_target
    hicmat = sparse2matrix(sparsemat, floor(chrlen[c]/resolution)+1, resolution)
    hicmat_ds = funDownsample(hicmat, n_sample)
    
    sparsedf <- as.data.frame(summary(Matrix(hicmat_ds)))
    sparsedf[,1] <- sparsedf[,1] - 1
    sparsedf[,2] <- sparsedf[,2] - 1
    
    write.table(sparsedf, file = paste0("/media/biology/2b52abd3-7feb-4bc9-8c84-a8408d93f7ec/dump/GM12878/", chrname[c], "/25K_ds350K/", basename(file)), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}

