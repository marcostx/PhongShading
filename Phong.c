#include <stdio.h>
#include "ift.h"

#define  NUM_NORMALS    65161

#define GetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define GetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define GetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))


#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])

int RADIUS_THRESHOLD=2;

typedef struct phong_model {
  float      ka;     
  float      kd;     
  float      ks;     
  float      ns;     
  int        ndists; 
  iftVector *normal; 
  float     *Idist;  
} PhongModel; 

typedef struct plane {
  iftPoint  pos;    
  iftVector normal; 
} Plane;

typedef struct _surface_rendering_buffers {
  float  depth;         /* distance to the viewing plane, whose center
         is initially positioned at
         (diag/2,diag/2,-diag/2) */ 
  float  opacity;       /* accumulated opacity of the voxels projected
         onto the viewing plane */
  int    voxel;         /* voxels in the scene that have been projected
         onto the viewing plane */
  int    object;        /* object in the scene whose voxel has been
         projected onto the viewing plane */
} SRBuffers;

typedef struct object_attr {
  float opacity;          
  float green;
  float red;
  float blue; 
  char visibility;        
} ObjectAttributes; 


typedef struct gc {
  ObjectAttributes *object;       
  int               numberOfObjects;
  float             overall_opac; 
  char              proj_mode;    
  PhongModel       *phong;        
  //ViewDir          *viewdir;      
  SRBuffers        *surf_render;  
  iftImage         *scene;        
  iftImage         *label;        
  Plane            *face;         
  iftImage         *normal;       
  iftImage         *opacity;      
} GraphicalContext;




int VoxelValue(iftImage *img, iftVoxel v)
{
    return img->val[GetVoxelIndex(img, v)];
}

int sign( int x ){
    if(x >= 0)
        return 1;
    return -1;
}

typedef struct volume
{
    iftMatrix* orthogonal;
    iftMatrix* center;
}iftVolumeFaces;

int isValidPoint(iftImage *img, iftVoxel u)
{
    if ((u.x >= 0) && (u.x < img->xsize) &&
        (u.y >= 0) && (u.y < img->ysize) &&
        (u.z >= 0) && (u.z < img->zsize)){
        return 1;
    }
    else{
        return 0;
    }
}

int diagonalSize(iftImage *img)
{
    return ROUND(sqrt((double) (img->xsize * img->xsize) + (img->ysize * img->ysize) + (img->zsize * img->zsize)));
}


float MatrixInnerProduct(iftMatrix *A, iftMatrix *B)
{
  float result = 0;

  for (int i = 0; i < A->ncols; i++)
    result += (A->val[i] * B->val[i]);

  return result;
}

iftVoxel GetVoxelCoord(iftImage *img, int p)
{
    iftVoxel u;

    u.x = GetXCoord(img, p);
    u.y = GetYCoord(img, p);
    u.z = GetZCoord(img, p);

    return u;
}

/* this function creas the rotation/translation matrix for the given theta */
iftMatrix *createTransformationMatrix(iftImage *img, int xtheta, int ytheta)
{
    iftMatrix *resMatrix = NULL;

    iftVector v1 = {.x = (float)img->xsize / 2.0, .y = (float)img->ysize / 2.0, .z = (float)img->zsize / 2.0};
    iftMatrix *transMatrix1 = iftTranslationMatrix(v1);

    iftMatrix *xRotMatrix = iftRotationMatrix(IFT_AXIS_X, -xtheta);
    iftMatrix *yRotMatrix = iftRotationMatrix(IFT_AXIS_Y, -ytheta);

    float D = sqrt(img->xsize*img->xsize + img->ysize*img->ysize);
    iftVector v2 = {.x = -(D / 2.0), .y = -(D / 2.0), .z = -(D / 2.0)};
    iftMatrix *transMatrix2 = iftTranslationMatrix(v2);


    resMatrix = iftMultMatricesChain(4, transMatrix1, xRotMatrix,yRotMatrix, transMatrix2);

    return resMatrix;
}

int DDA(iftImage *img, iftVoxel p1, iftVoxel pn)
{
    int n, k;
    //iftVoxel p;
    iftVoxel p;
    float J=0, max=0;
    int Dx,Dy,Dz;
    float dx=0,dy=0,dz=0;
    //iftCreateMatrix* J;


    if (p1.x == pn.x && p1.y == pn.y && p1.z == pn.z)
        n=1;
    else{
        Dx=pn.x - p1.x;
        Dy=pn.y - p1.y;
        Dz=pn.z - p1.z;

        if( abs(Dx) >= abs(Dy) && abs(Dx) >= abs(Dz) ){
            n = abs(Dx)+1;
            dx = sign(Dx);
            dy = (dx * Dy)/Dx;
            dz = (dx * Dz)/Dx;
        }
        else if ( abs(Dy) >= abs(Dx) && abs(Dy) >= abs(Dz)){
            n = abs(Dy)+1;
            dy = sign(Dy);
            dx = (dy * Dx)/Dy;
            dz = (dy * Dz)/Dy;
        }
        else{
            n = abs(Dz)+1;
            dz = sign(Dz);
            dx = (dz * Dx)/Dz;
            dy = (dz * Dy)/Dz;
        }
    }

    p.x = p1.x;
    p.y = p1.y;
    p.z = p1.z;

    //int validP=0;
    for (k = 1; k < n; k++)
    {
        if (isValidPoint(img,p)){

        }

        p.x = p.x + dx;
        p.y = p.y + dy;
        p.z = p.z + dz;
    }

    return (int)max;
}

iftVector *createNormalTable()
{
    int        i, gamma, alpha;
    float      gamma_rad, alpha_rad;
    iftVector *normaltable;

    /* creates normal look-up table */

    normaltable = (iftVector *)calloc(NUM_OF_NORMALS, sizeof(iftVector));

    normaltable[0].x = 0.0;
    normaltable[0].y = 0.0;
    normaltable[0].z = 0.0;

    i = 1;
    for (gamma = -90; gamma <= 90; gamma++)
    {
        gamma_rad = (PI * gamma) / 180.0;
        for (alpha = 0; alpha < 360; alpha++)
        {
            alpha_rad = (PI * alpha) / 180.0;
            normaltable[i].x = cos(gamma_rad) * cos(alpha_rad);
            normaltable[i].y = cos(gamma_rad) * sin(alpha_rad);
            normaltable[i].z = sin(gamma_rad);
            i++;
        }
    }

    return (normaltable);
}

PhongModel *createPhongModel(iftImage *scene)
{
    PhongModel *phong = (PhongModel *) malloc(sizeof(PhongModel) * 6);

    phong->ka     = 0.1;
    phong->kd     = 0.7;
    phong->ks     = 0.2;
    phong->ns     = 5.0;
    phong->normal = createNormalTable();
    phong->ndists = (int)(2.0 * diagonalSize(scene) + 1);
    phong->Idist  = (float *) malloc(phong->ndists*sizeof(float));



    for (int d = 0; d < phong->ndists; d++)
        phong->Idist[d] = (float)0.8 * (phong->ndists - d) / (float)phong->ndists + 0.2;

    return (phong);
}

int MaximumValue(iftImage *img)
{
    int p, Imax = -999999999;

    for (p = 0; p < img->n; p++)
        if (img->val[p] > Imax)
            Imax = img->val[p];

    return Imax;
}

ObjectAttributes *createObjectAttr(iftImage *label, int *numberOfObjects)
{
    ObjectAttributes *object;
    *numberOfObjects = MaximumValue(label);

    object = (ObjectAttributes *)calloc(*numberOfObjects + 1, sizeof(ObjectAttributes));

    object[0].opacity    = 0;
    object[0].red        = 0;
    object[0].green      = 0;
    object[0].blue       = 0;
    object[0].visibility = 0;

    /* default for objects */

    for (int i = 1; i <= *numberOfObjects; i++)
    {
        object[i].opacity    = 1;
        object[i].red        = 1;
        object[i].green      = 1;
        object[i].blue       = 1;
        object[i].visibility = 1;
    }

    return (object);
}

void setSceneFaces(GraphicalContext *gc)
{
    iftImage *scene = gc->scene;

    gc->face[0].pos.x = scene->xsize / 2.0;
    gc->face[0].pos.y = scene->ysize / 2.0;
    gc->face[0].pos.z = 0;
    gc->face[1].pos.x = scene->xsize / 2.0;
    gc->face[1].pos.y = scene->ysize / 2.0;
    gc->face[1].pos.z = scene->zsize - 1;
    gc->face[2].pos.x = scene->xsize / 2.0;
    gc->face[2].pos.y = 0;
    gc->face[2].pos.z = scene->zsize / 2.0;
    gc->face[3].pos.x = scene->xsize / 2.0;
    gc->face[3].pos.y = scene->ysize - 1;
    gc->face[3].pos.z = scene->zsize / 2.0;
    gc->face[4].pos.x = 0;
    gc->face[4].pos.y = scene->ysize / 2.0;
    gc->face[4].pos.z = scene->zsize / 2.0;
    gc->face[5].pos.x = scene->xsize - 1;
    gc->face[5].pos.y = scene->ysize / 2.0;
    gc->face[5].pos.z = scene->zsize / 2.0;

    gc->face[0].normal.x =  0;
    gc->face[0].normal.y =  0;
    gc->face[0].normal.z = -1;
    gc->face[1].normal.x =  0;
    gc->face[1].normal.y =  0;
    gc->face[1].normal.z =  1;
    gc->face[2].normal.x =  0;
    gc->face[2].normal.y = -1;
    gc->face[2].normal.z =  0;
    gc->face[3].normal.x =  0;
    gc->face[3].normal.y =  1;
    gc->face[3].normal.z =  0;
    gc->face[4].normal.x = -1;
    gc->face[4].normal.y =  0;
    gc->face[4].normal.z =  0;
    gc->face[5].normal.x =  1;
    gc->face[5].normal.y =  0;
    gc->face[5].normal.z =  0;

}

GraphicalContext *createGC(iftImage *scene, iftImage *imageLabel)
{
    GraphicalContext *gc;
    
    gc = (GraphicalContext *) calloc(1, sizeof(GraphicalContext));
    
    gc->scene          = iftCopyImage(scene);
    gc->phong          = createPhongModel(scene);
    //gc->viewdir        = NULL;
    gc->numberOfObjects= 0;
    gc->overall_opac   = 1.0;
    gc->face           = (Plane *) malloc(sizeof(Plane) * 6);
    //setViewDir(gc, 0, 0);
    setSceneFaces(gc);
    
    gc->label       = iftCopyImage(imageLabel);
    gc->object      = createObjectAttr(imageLabel, &gc->numberOfObjects);
    //gc->surf_render = createSRBuffers(imageLabel);
    
    gc->opacity       = NULL;
    gc->normal        = NULL;
    
    return (gc);
}


iftVector columnVectorMatrixToVector(iftMatrix *m)
{
    iftVector v = {.x = iftMatrixElem(m, 0, 0), .y = iftMatrixElem(m, 0, 1), .z = iftMatrixElem(m, 0, 2)};
    v.x = iftAlmostZero(v.x) ? 0.0 : v.x;
    v.y = iftAlmostZero(v.y) ? 0.0 : v.y;
    v.z = iftAlmostZero(v.z) ? 0.0 : v.z;

    return v;
}


int ComputeIntersection(iftMatrix *Tpo, iftImage *img, iftMatrix *Tn, iftVolumeFaces *vf, iftVoxel *p1, iftVoxel *pn)
{

    float lambda[6] = { -1};
    float max=-9999999.9, min=9999999.9;
    int i;
    p1->x=pn->x=p1->y=pn->y,p1->z=pn->z=-1;
    iftMatrix *Nj = iftCreateMatrix(1, 3);
    iftMatrix *DiffCandP0 = iftCreateMatrix(1, 3);
    float NdotNj = 0, DiffShiftDotNj = 0;
    iftVector v1, v2, v3;
    iftVoxel v;


    for (i = 0; i < 6; i++) {
      iftMatrixElem(Nj, 0, 0) = vf[i].orthogonal->val[0];
      iftMatrixElem(Nj, 0, 1) = vf[i].orthogonal->val[1];
      iftMatrixElem(Nj, 0, 2) = vf[i].orthogonal->val[2];
      v1 = columnVectorMatrixToVector(Nj);
      v2 = columnVectorMatrixToVector(Tn);

      NdotNj = iftVectorInnerProduct(v1,v2);

      //NdotNj = MatrixInnerProduct(Tn,Nj);

      if (NdotNj != 0){

        iftMatrixElem(DiffCandP0, 0, 0) = vf[i].center->val[0] - Tpo->val[0];
        iftMatrixElem(DiffCandP0, 0, 1) = vf[i].center->val[1] - Tpo->val[1];
        iftMatrixElem(DiffCandP0, 0, 2) = vf[i].center->val[2] - Tpo->val[2];

        v3 = columnVectorMatrixToVector(DiffCandP0);

        DiffShiftDotNj = iftVectorInnerProduct(v1,v3);
        //DiffShiftDotNj = MatrixInnerProduct(Nj,DiffCandP0);

        //if(NdotNj == 0.000000){
        //  continue;
        //}

        lambda[i]=(float) DiffShiftDotNj / NdotNj;


        v.x = Tpo->val[0] + lambda[i] * Tn->val[0];
        v.y = Tpo->val[1] + lambda[i] * Tn->val[1];
        v.z = Tpo->val[2] + lambda[i] * Tn->val[2];

        if (isValidPoint(img, v))
        {
          if (lambda[i] < min){
              p1->x = v.x;
              p1->y = v.y;
              p1->z = v.z;
              min = lambda[i];
          }
          if (lambda[i] > max) {
              pn->x = v.x;
              pn->y = v.y;
              pn->z = v.z;
              max = lambda[i];
          }
        }
      }
    }

    iftDestroyMatrix(&Nj);
    iftDestroyMatrix(&DiffCandP0);


    if ((p1->x != -1) && (pn->x != -1))
      return 1;
    else
      return 0;
}



int LinearInterpolationValue(iftImage *img, float x, float y, float z)
{
    iftVoxel u[8];
    float dx = 1.0;
    float dy = 1.0;
    float dz = 1.0;
    float  P12, P34, P56, P78;
    float aux1, aux2;
    int Pi;

    if ((int) (x + 1.0) == img->xsize)
        dx = 0.0;
    if ((int) (y + 1.0) == img->ysize)
        dy = 0.0;
    if ((int) (z + 1.0) == img->zsize)
        dz = 0.0;

    //closest neighbour in each direction
    u[0].x = (int)x;        u[0].y = (int)y;          u[0].z = (int)z;
    u[1].x = (int)(x + dx); u[1].y = (int)y;          u[1].z = (int)z;
    u[2].x = (int)x;        u[2].y = (int)(y + dy);   u[2].z = (int)z;
    u[3].x = (int)(x + dx); u[3].y = (int)(y + dy);   u[3].z = (int)z;
    u[4].x = (int)x;        u[4].y = (int)y;          u[4].z = (int)(z + dz);
    u[5].x = (int)(x + dx); u[5].y = (int)y;          u[5].z = (int)(z + dz);
    u[6].x = (int)x;        u[6].y = (int)(y + dy);   u[6].z = (int)(z + dz);
    u[7].x = (int)(x + dx); u[7].y = (int)(y + dy);   u[7].z = (int)(z + dz);


    P12 = (float)iftImgVal2D(img,u[1].x,u[1].y) * (x - u[0].x) + (float)iftImgVal2D(img,u[0].x,u[0].y) * (u[1].x - x);
    P34 = (float)iftImgVal2D(img,u[3].x,u[3].y) * (x - u[2].x) + (float)iftImgVal2D(img,u[2].x,u[2].y) * (u[3].x - x);
    P56 = (float)iftImgVal2D(img,u[5].x,u[5].y) * (x - u[4].x) + (float)iftImgVal2D(img,u[4].x,u[4].y) * (u[5].x - x);
    P78 = (float)iftImgVal2D(img,u[7].x,u[7].y) * (x - u[6].x) + (float)iftImgVal2D(img,u[6].x,u[6].y) * (u[7].x - x);
    aux1 = P34 *  (y - u[0].y) + P12 * (u[2].y - y);
    aux2 = P56 * (y - u[0].y) + P78 * (u[2].y - y);
    Pi  = (int)aux2 * (z - u[0].z) + aux1 * (u[4].z - z);

    return Pi;
}



void DestroyVF(iftVolumeFaces *vf)
{
    int i;

    for (i = 0; i < 6; i++)
    {
        iftDestroyMatrix(&vf[i].orthogonal);
        iftDestroyMatrix(&vf[i].center);
    }
    free(vf);
}


iftMatrix *voxelToMatrix(iftVoxel v)
{
    iftMatrix *voxMat = iftCreateMatrix(1, 4);
    iftMatrixElem(voxMat, 0, 0) = v.x;
    iftMatrixElem(voxMat, 0, 1) = v.y;
    iftMatrixElem(voxMat, 0, 2) = v.z;
    iftMatrixElem(voxMat, 0, 3) = 1;

    return voxMat;
}

iftMatrix *imagePixelToMatrix(iftImage *img, int p)
{
    iftMatrix *pixMat = iftCreateMatrix(1, 4);
    iftMatrixElem(pixMat, 0, 0) = p % img->xsize;
    iftMatrixElem(pixMat, 0, 1) = p / img->xsize;
    iftMatrixElem(pixMat, 0, 2) = 0;
    iftMatrixElem(pixMat, 0, 3) = 1;

    return pixMat;
}




int main(int argc, char *argv[])
{
    //if (argc != 6)
    //    Error("Run: ./main <filename> <output> <xtheta> <ytheta> , "main");

    //char buffer[512];

    //float tx, ty;
    GraphicalContext *gc;
    char *imgFileName = iftCopyString(argv[1]);
    char *imgLabelFileName = iftCopyString(argv[2]);
    iftImage *img = iftReadImageByExt(imgFileName);
    iftImage *imgLabel = iftReadImageByExt(imgLabelFileName);
    //tx = atof(argv[4]);
    //ty = atof(argv[5]);
    
    //iftImage *output = NULL;

    gc = createGC(img, imgLabel);
    printf("Done!\n");
    // sprintf(buffer, "data/%.1f%.1f%s", tx, ty, argv[2]);
    // iftImage *normalizedImage= iftNormalize(output,0,255);

    // iftWriteImageByExt(normalizedImage, buffer);
    // iftDestroyImage(&img);
    // iftDestroyImage(&output);
    return 0;
}
