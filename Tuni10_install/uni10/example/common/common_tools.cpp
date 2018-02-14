#include "common_tools.h"

void progressbar(long long int progress, long long int begin, long long int end, bool isflush){

  long long int total_prgress = end - begin;
  int barWidth = 50;
  
  fprintf(stdout, "[");
  int pos = barWidth * (float)progress/(float)total_prgress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) 
      fprintf(stdout, "=");
    else if (i == pos && progress != end) 
      fprintf(stdout, ">") ;
    else if (i == pos && progress == end) 
      fprintf(stdout, "=") ;
    else 
      fprintf(stdout, " ") ;
  }



  if(isflush){
    fprintf(stdout, "] %.2f%s", (float)progress / (float)total_prgress * 100., "%\r");
    fflush(stdout);
  }
  else
    fprintf(stdout, "] %.2f%s", (float)progress / (float)total_prgress * 100., "%");

}
