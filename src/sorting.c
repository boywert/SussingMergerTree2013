#include "allvars.h"
int compareHaloID(const void *v1, const void *v2)
{
    const struct progtable *u1 = v1;
    const struct progtable *u2 = v2;
    if(u1->haloID < u2->haloID)
      return -1;
    else if(u1->haloID > u2->haloID)
      return 1;
    else if(u1->haloID == u2->haloID)
      return 0;
}

int compareID(const void *v1, const void *v2)
{
    const struct HALOPROPS *u1 = v1;
    const struct HALOPROPS *u2 = v2;
    if(u1->ID < u2->ID)
      return -1;
    else if(u1->ID > u2->ID)
      return 1;
    else if(u1->ID == u2->ID)
      return 0;
}

int compareMyIDType(const void *v1, const void *v2)
{
  if (*(MyIDtype *)v1 < *(MyIDtype *)v2)
    return -1;
  else if (*(MyIDtype *)v1 > *(MyIDtype *)v2)
    return 1;
  else if (*(MyIDtype *)v1 == *(MyIDtype *)v2)
    return 0;
}

int compareParticleEnergy(const void *v1, const void *v2)
{
    const struct particle_data *u1 = v1;
    const struct particle_data *u2 = v2;
    if(u1->ParticleEnergy < u2->ParticleEnergy)
      return -1;
    else if(u1->ParticleEnergy > u2->ParticleEnergy)
      return 1;
    else if(u1->ParticleEnergy == u2->ParticleEnergy)
      return 0;
}

MyIDtype IDsearch( MyIDtype searchKey )
{
  MyIDtype middle,low,high;
  low = 0;
  high = TotNhalos-1;

  if(searchKey < IDmap[0] || searchKey > IDmap[TotNhalos-1])
    {
      return NULLPOINT;
    }
  while ( low <= high) {
    middle = ( low + high ) / 2;
 
    if ( searchKey == IDmap[ middle ] ) 
      {
	return middle;
      } 
    else if ( searchKey < IDmap[ middle ] ) 
      {
	high = middle - 1;
      } 
    else 
      {
	low = middle + 1;
      }
   }
  return NULLPOINT;
}
MyIDtype Generalsearch( MyIDtype searchKey,MyIDtype n_array ,const void *Array  )
{
  MyIDtype middle,low,high,i;
  MyIDtype *pool = (MyIDtype*) Array;
  //printf("start search\n");
  /*
  for(i =0; i< n_array; i++)
    {
      printf("%llu\n",pool[i]);
    }
  */
  low = 0;
  high = n_array-1;

  if(searchKey < pool[0] || searchKey > pool[n_array-1])
    {
      return NULLPOINT;
    }
  while ( low <= high) {
    middle = ( low + high ) / 2;
 
    if ( searchKey == pool[ middle ] ) 
      {
	return middle;
      } 
    else if ( searchKey < pool[ middle ] ) 
      {
	high = middle - 1;
      } 
    else 
      {
	low = middle + 1;
      }
   }
  return NULLPOINT;
}
