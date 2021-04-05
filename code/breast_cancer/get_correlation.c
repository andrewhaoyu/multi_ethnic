#include <stdlib.h>
#include <math.h>
#include <stdio.h>
// plinkBed, the file name of plink bed file
// NumSample, Number of Samples in Plink
// NumSNP, Number of SNPs in prs_model
// idx, SNP's (in prs_model) index in plink bim file
// chr, SNP's (in prs_model) chromosome
// beta_tau, beta * tau in prs model
// adj_cor, output summation of beta*tau*cor
// D, dosage dataMatrix, each column is a SNP,  each row is a sample, missing value must be set to -9

double correlation(int *Xs, int *Ys, int *N) {
  double sumX = 0;
  double sumX2 = 0;
  double sumY = 0;
  double sumY2 = 0;
  double sumXY = 0;
  int x, y;
  int n = 0;
  for(int i = 0; i < *N; i++) {
    x = Xs[i];
    y = Ys[i];
    if(abs(x) > 2 || abs(y) > 2)
      continue;
    sumX += x;
    sumX2 += x * x;
    sumY += y;
    sumY2 += y * y;
    sumXY += x * y;
    n++;
  }
  
  // double stdX = sumX2 - sumX * sumX / n;
  // double stdY = sumY2 - sumY * sumY / n;
  //
    // if(n == 0 || stdX <= 0 || stdY <= 0 || fpclassify(stdX) == FP_ZERO || fpclassify(stdY) == FP_ZERO)
      // 	return 0;
  //
    // stdX = sqrt(stdX);
    // stdY = sqrt(stdY);
    
    double stdX = sqrt(sumX2 - sumX * sumX / n);
    double stdY = sqrt(sumY2 - sumY * sumY / n);
    double cov = sumXY - sumX * sumY / n;
    
    double cor = cov / stdX / stdY;
    // printf( "Correlation is : %f\n", cor);
    return cor;
}

void getAdjCorrelation (const char **plinkBed, int *NumSample, int *NumSNP, int *idx, int *chrs, int *pos, int *pos_thr, double *beta_tau, double *adj_cor){
  
  int NSam = *NumSample;
  int NSNP = *NumSNP;
  int sampleOffSet = NSam % 4;
  int snpBulkSize = NSam/4;
  if(sampleOffSet > 0)
    snpBulkSize++;
  int *Gmap = (int *) malloc (4 * 256 * sizeof(int));
  int *G = (int *) malloc (NSNP * NSam * sizeof(int));
  if(G == NULL){
    printf("Out of memory!\n");
    exit(1);
  }
  printf("Memory size of genotype matrix G : %fM \n", NSNP * NSam * sizeof(int) / pow(2, 20));
  unsigned char buffer; // note: 1 byte
  int tmp;
  int i, j, k;
  long int snp_idx;
  
  // creat a map for a single byte to 4 genotype in the same order with plink.fam file
  // counts are the number of Allele A1 in bim file
  for(i = 0; i < 256; i++){
    buffer = (unsigned char) i;
    k = i * 4;
    for(j = 0; j < 4; j++){
      tmp = (buffer >> (j * 2)) & 3;
      switch(tmp){
        case 0:
          Gmap[k + j] = 2;
          break;
          case 2:
            Gmap[k + j] = 1;
            break;
            case 3:
              Gmap[k + j] = 0;
              break;
              default:
                Gmap[k + j] = -9;
      }
    }
  }
  
  
  FILE *ptr_plink = fopen(*plinkBed, "rb");
  if(!ptr_plink){
    printf("Unable to open file!");
    exit(1);
  }
  // fseek(ptr_plink, 3, SEEK_SET);
  
  for(j = 0; j < NSNP; j++){
    snp_idx = abs(idx[j]) - 1;
    fseek(ptr_plink, 3 + snpBulkSize * snp_idx, SEEK_SET);
    if(idx[j] < 0){
      for(i = 0; i <= NSam - 4; i=i+4){
        fread(&buffer, 1, 1, ptr_plink);
        // tmp = tmp * 4
        tmp = (((int) buffer) << 2);
        for(k = 0; k < 4; k++){
          G[i+k+j*NSam] = 2 - Gmap[tmp + k];
        }
      }
      if(sampleOffSet > 0){
        fread(&buffer, 1, 1, ptr_plink);
        // tmp = tmp * 4
        tmp = (((int) buffer) << 2);
        for(k = 0; k < sampleOffSet; k++){
          G[i+k+j*NSam] = 2 - Gmap[tmp + k];
        }
      }
      
    }else{
      for(i = 0; i <= NSam - 4; i=i+4){
        fread(&buffer, 1, 1, ptr_plink);
        // tmp = tmp * 4
        tmp = (((int) buffer) << 2);
        for(k = 0; k < 4; k++){
          G[i+k+j*NSam] = Gmap[tmp + k];
        }
      }
      if(sampleOffSet > 0){
        fread(&buffer, 1, 1, ptr_plink);
        // tmp = tmp * 4
        tmp = (((int) buffer) << 2);
        for(k = 0; k < sampleOffSet; k++){
          G[i+k+j*NSam] = Gmap[tmp + k];
        }
      }
      
    }
  }
  // printf("adj_cor[0] : %f \n", adj_cor[0]);
  adj_cor[0] = 0;
  // printf("adj_cor[0] : %f \n", adj_cor[0]);
  for(i = 0; i < NSNP - 1; i++){
    int *Gi = &G[i*NSam];
    for(j = i + 1; j < NSNP; j++){
      if(chrs[i] != chrs[j])
        continue;
      if(abs(pos[i] - pos[j]) > *pos_thr)
        continue;
      int *Gj = &G[j*NSam];
      double cur_cor = correlation(Gi, Gj, &NSam);
      if(isnormal(cur_cor))
        adj_cor[0] = adj_cor[0] + cur_cor * beta_tau[i] * beta_tau[j];
      // else
        // 	printf("i: %d\tj: %d\tcur_cor: %f\tadj_cor: %f\n", i, j, cur_cor, adj_cor[0]);
    }
  }
  printf("adj_cor : %f \n", adj_cor[0]);
  fclose(ptr_plink);
  free(G);
  free(Gmap);
}


double correlationDosage(double *Xs, double *Ys, int *N) {
  double sumX = 0;
  double sumX2 = 0;
  double sumY = 0;
  double sumY2 = 0;
  double sumXY = 0;
  double x, y;
  int n = 0;
  for(int i = 0; i < *N; i++) {
    x = Xs[i];
    y = Ys[i];
    if(fabs(x) > 2 || fabs(y) > 2)
      continue;
    sumX += x;
    sumX2 += x * x;
    sumY += y;
    sumY2 += y * y;
    sumXY += x * y;
    n++;
  }
  
  // double stdX = sumX2 - sumX * sumX / n;
  // double stdY = sumY2 - sumY * sumY / n;
  //
    // if(n == 0 || stdX <= 0 || stdY <= 0 || fpclassify(stdX) == FP_ZERO || fpclassify(stdY) == FP_ZERO)
      // 	return 0;
  //
    // stdX = sqrt(stdX);
    // stdY = sqrt(stdY);
    
    double stdX = sqrt(sumX2 - sumX * sumX / n);
    double stdY = sqrt(sumY2 - sumY * sumY / n);
    double cov = sumXY - sumX * sumY / n;
    
    double cor = cov / stdX / stdY;
    // if(isnormal(cor))
      // 	printf( "Correlation is : %f\n", cor);
    // else
      //	printf( "Correlation is not available NA!\n");
    return cor;
}

void getAdjCorrelationDosage (double *D, int *NumSample, int *NumSNP, int *chrs, int *pos, int *pos_thr, double *beta_tau, double *adj_cor){
  
  int NSam = *NumSample;
  int NSNP = *NumSNP;
  
  int i, j;
  
  adj_cor[0] = 0;
  // printf("adj_cor[0] : %f \n", adj_cor[0]);
  for(i = 0; i < NSNP - 1; i++){
    double *Di = &D[i*NSam];
    for(j = i + 1; j < NSNP; j++){
      if(chrs[i] != chrs[j])
        continue;
      if(abs(pos[i] - pos[j]) > *pos_thr)
        continue;
      double *Dj = &D[j*NSam];
      double cur_cor = correlationDosage(Di, Dj, &NSam);
      if(isnormal(cur_cor))
        adj_cor[0] = adj_cor[0] + cur_cor * beta_tau[i] * beta_tau[j];
      // else
        // 	printf("i: %d\tj: %d\tcur_cor: %f\tadj_cor: %f\n", i, j, cur_cor, adj_cor[0]);
    }
  }
  printf("adj_cor : %f \n", adj_cor[0]);
  
}