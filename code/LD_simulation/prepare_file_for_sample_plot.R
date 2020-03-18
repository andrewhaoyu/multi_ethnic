#we are going to use the gene MTHFR to show the LD difference (chr1 pos 11785723 to 11806103)
#get the position number of the simulated genotypes from the gen file
#select the particular regions for different ethnic groups
#cd /data/zhangh24/multi_ethnic/result/LD_simulation/EUR
#awk -F " " '($3 >= 11785723) && ($3 <= 11806103)' chr1_1.controls.gen >> chr1_MTHFR.gen
#do similar things for all other ethnic groups
#create the tag file for 1KG to use
#awk '{print $3}' chr1_MTHFR_EUR.gen >> MTHFR.tag 
#mv MTHFR.tag 
