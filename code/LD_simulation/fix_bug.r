#the EUR chr2 data have some corruputed part
#take those out
tag <- as.data.frame(fread("/gpfs/gsfs11/users/zhangh24/multi_ethnic/result/LD_simulation/EUR/chr2.tag.info.txt",header=F))
idx <- which(tag$V3==239637061)

#head -700000  chr2.combined.tag.gen	>  chr2.combined.tag.gen.temp	    
