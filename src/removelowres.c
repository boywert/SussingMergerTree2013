#include "allvars.h"

void SubtreeFilter();
void resetHaloTable();

char FilePrefix[MAXSTRING];
char HBTFilePrefix[MAXSTRING];
char FolderName[MAXSTRING];
char HBTFolderName[MAXSTRING];
char TreeFolder[MAXSTRING];
char ppFolder[MAXSTRING];
char inputFile[MAXSTRING];
char gadgetfolder[MAXSTRING];
char *trimwhitespace(char *str);

int main(int argc, char **argv)
{
  int isTree;
  MyIDtype i;
  char buffer[MAXSTRING];
  FILE* fp;

  //FilePrefix = calloc(MAXSTRING,sizeof(char));
  //FolderName = calloc(MAXSTRING,sizeof(char));
  //  strcpy(FolderName,"/export/virgo/Boyd/Millgas/62.5_ALEX");
  //  strcpy(FilePrefix,"62.5_dm_");

  sprintf(FolderName,"/gpfs/data/Millgas/cs390/SUSSING2013/datasetIII_unfiltered");
  sprintf(FilePrefix,"datasetIII_");
  sprintf(SnapTimeFile,"/gpfs/data/Millgas/cs390/SUSSING2013/datasetIII_unfiltered/a4_data_snaplist.txt");
  sprintf(gadgetfolder,"/gpfs/data/aquarius/halo_data/Aq-A/4");
  getSnapTime();

  read_singlesnap(61);

  return 0;
}
