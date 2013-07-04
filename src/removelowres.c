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
char gadgetPrefix[MAXSTRING];
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
#ifdef AQUARIUS
  sprintf(FolderName,"/gpfs/data/Millgas/cs390/SUSSING2013/datasetIII_unfiltered");
  sprintf(FilePrefix,"datasetIII_");
  sprintf(SnapTimeFile,"/gpfs/data/Millgas/cs390/SUSSING2013/datasetIII_unfiltered/a4_data_snaplist.txt");
  sprintf(gadgetfolder,"/gpfs/data/aquarius/halo_data/Aq-A/4");
  sprintf(gadgetPrefix,"snap_C02_400_");

#else
  sprintf(FolderName,"/scratch/cs390/SUSSING2013/datasetI");
  sprintf(FilePrefix,"62.5_dm_");
  sprintf(SnapTimeFile,"/scratch/cs390/SUSSING2013/datasetI/data_snaplist.txt");
  sprintf(gadgetfolder,"/scratch/cs390/SUSSING2013/62.5_ori");
  sprintf(gadgetPrefix,"62.5_dm_");
#endif
  getSnapTime();

  read_singlesnap(500);

  return 0;
}
