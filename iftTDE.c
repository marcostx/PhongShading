#include "ift.h"

iftSet *ObjectBorders(iftImage *label)
{
  iftAdjRel *A;
  iftSet    *B=NULL;

    A = iftSpheric(1.0);

  for (int p=0; p < label->n; p++) {
    if (label->val[p] != 0){ // p is an object voxel
      iftVoxel u = iftGetVoxelCoord(label,p);

      for (int i=1; i < A->n; i++) {
      	iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      	if (iftValidVoxel(label,v)){
      	  int q = iftGetVoxelIndex(label,v);
      	  if (label->val[q]!=label->val[p]){ // q belongs to another
      					     // object/background,
      					     // then p is border
      	    iftInsertSet(&B,p);
      	    break;
  	      }
  	     }else{ // p is at the border of the image, then it is border
	         iftInsertSet(&B,p);
	       break;
	     }
      }
    }
  }

  iftDestroyAdjRel(&A);
  return B;
}


iftImage *iftTDE(iftImage *label)
{
  iftSet    *B=ObjectBorders(label);
  iftImage  *tde, *root;
  iftGQueue *Q;
  iftAdjRel *A;

    A = iftSpheric(3.0);

  tde    = iftCreateImage(label->xsize, label->ysize, label->zsize);
  root   = iftCreateImage(label->xsize, label->ysize, label->zsize);
  Q      = iftCreateGQueue(IFT_QSIZE, tde->n, tde->val);

  for (int p=0; p < label->n; p++)
    if (label->val[p]!=0)
      tde->val[p] = IFT_INFINITY_INT;

  while (B != NULL)
  {
    int p           = iftRemoveSet(&B);
    root->val[p]    = p;
    tde->val[p]     = 0;
    iftInsertGQueue(&Q, p);
  }

  // Image Foresting Transform

  while (!iftEmptyGQueue(Q))
  {
    int p       = iftRemoveGQueue(Q);
    iftVoxel u  = iftGetVoxelCoord(label, p);
    iftVoxel rp = iftGetVoxelCoord(label,root->val[p]);

    for (int i = 1; i < A->n; i++)
    {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(label, v))
      {
        int q = iftGetVoxelIndex(label, v);

      	if ((label->val[q]==label->val[p])&&(tde->val[q]>tde->val[p])){
      	  int tmp = iftSquaredVoxelDistance(v,rp);
      	  if (tmp < tde->val[q]) {
      	    //if (tde->val[q] != IFT_INFINITY_INT)
      	    if (Q->L.elem[q].color == IFT_GRAY){
      		  iftRemoveGQueueElem(Q,q);
      	    }
      	    root->val[q]     = root->val[p];
      	    tde->val[q]      = tmp;
      	    iftInsertGQueue(&Q, q);
      	  }
      	}
      }
    }
  }

  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);
  iftDestroyImage(&root);

  return(tde);
}

int main(int argc, char *argv[])
{
  iftImage       *label=NULL, *tde=NULL;
  timer          *tstart=NULL;
  int             MemDinInicial, MemDinFinal;

  MemDinInicial = iftMemoryUsed(1);

  if (argc!=3)
    iftError("Usage: iftTDE <input-label.[*]> <output-TDE.[*]>","main");


  tstart = iftTic();

  label   = iftReadImageByExt(argv[1]);
  tde     = iftTDE(label);

  iftWriteImageByExt(tde,argv[2]);

  iftDestroyImage(&label);
  iftDestroyImage(&tde);

  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));

  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);

  return(0);
}
