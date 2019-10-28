# Goal: merge 1000 genome map, hap and legend file into one
#map file
# the file header are position, COMBINED_rate(cm/Mb), Genetic_Map (cM)
# since all the files share the same header
# we will use head and tail function in linux to merge them together into one file
# use head function
setwd("/dcl01/chatterj/data/hzhang1/1000GP_Phase3")
#head command will take the first row of the file
take.head = paste0("head -1 genetic_map_chr1_combined_b37.txt >> map_all.txt")
system(take.head)
#tail -n +2 will take the file from 2 to the end
#-q tell it not to print the header with the file name
merge.code <- paste0("tail -n +2 -q")
for(i in 1:22){
  #take out all the files with chri with grep command
  #fixed option = T only take out the exact match
  
  
  merge.code <- paste0(merge.code,
                       " genetic_map_chr",
                       i,
                       "_combined_b37.txt")  
  
  
  
  
}
#>> will add to the file, not overwrite it as >
merge.code <- paste0(merge.code,
                     " >> map_all.txt")
system(merge.code)

