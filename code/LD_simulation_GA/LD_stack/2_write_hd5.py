# -*- coding:utf-8 -*-
#!/usr/bin/python





"""


Read LD blocks and write as hdf5





"""


import sys


import scipy as sp


import h5py





i = int(sys.argv[1])

eth = ["eur", "afr", "amr", "eas", "sas"]






OUT_DIR = '/data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_' + eth[i]

BLK_DIR = '/data/zhangh24/software/PRScsx/1KGLD_MEGA/ldblk_1kg_'+ eth[i] + '/block_ld'




with open(BLK_DIR + '/snplist_ldblk/blk_chr') as ff:


    blk_chr = [int(line.strip()) for line in ff]





with open(BLK_DIR + '/snplist_ldblk/blk_size') as ff:


    blk_size = [int(line.strip()) for line in ff]





n_blk = len(blk_chr)


n_chr = max(blk_chr)




chrom = 1




for chrom in range(1,n_chr+1):


    print('... parse chomosome %d ...' % chrom)





    chr_name = OUT_DIR + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'


    hdf_chr = h5py.File(chr_name, 'w')


    blk_cnt = 0


    for blk in range(n_blk):


        if blk_chr[blk] == chrom:


            if blk_size[blk] > 0:


                blk_file = BLK_DIR + '/ldblock_' + str(blk+1) + '.ld'
                

                with open(blk_file) as ff:


                    ld = [[float(val) for val in (line.strip()).split()] for line in ff]


                print('blk %d size %s' % (blk+1, sp.shape(ld)))





                snp_file = BLK_DIR + '/snplist_ldblk/snplist_blk' + str(blk+1)


                with open(snp_file) as ff:


                    snplist = [line.strip() for line in ff]


            else:


                ld = []; snplist = []





            blk_cnt += 1


            hdf_blk = hdf_chr.create_group('blk_%d' % blk_cnt)


            hdf_blk.create_dataset('ldblk', data=sp.array(ld), compression="gzip", compression_opts=9)


            hdf_blk.create_dataset('snplist', data=snplist, compression="gzip", compression_opts=9)











