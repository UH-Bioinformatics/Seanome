#include <zlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


/*https://github.com/lh3/readfq*/
int main(int argc, char * argv[]){
  if(argc < 2 || argc > 3){
    fprintf(stderr, "USAGE: %s <fastq file> [prefix]\n",argv[0]);
    return 0;
  }
  gzFile fp;
  kseq_t *seq;
  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);
  size_t idx = 0;
  if(argc == 3){
    while (kseq_read(seq) >= 0){
      fprintf(stdout, "@%"PRId64"_%s\n%s\n+%s_%"PRId64"\n%s\n", idx, argv[2], seq->seq.s, argv[2] , idx, seq->qual.s);
      ++ idx;
    }
  } else {
    while (kseq_read(seq) >= 0){
      fprintf(stdout, "@%"PRId64"\n%s\n+%"PRId64"\n%s\n", idx, seq->seq.s, idx, seq->qual.s);
      ++ idx;
    }
  } 
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}
