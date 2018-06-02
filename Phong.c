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
  iftMatrix        *Tinv;
  iftVector        vDir;         
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

iftVoxel GetVoxelCoord(iftImage *img, int p)
{
    iftVoxel u;

    u.x = GetXCoord(img, p);
    u.y = GetYCoord(img, p);
    u.z = GetZCoord(img, p);

    return u;
}

float PhongShading(GraphicalContext *gc, int p, iftVector N, float dist)
{
    float cos;
    float cos2, aux, phong_val = 0.0;
    float thresold=1E-05;

    cos  = iftVectorInnerProduct(gc->vDir, N);
    if (cos > 1.0)
        cos = 1.0;

    if (cos > thresold)   // |angle| <= 90° 
    {

        cos2 = 2 * cos * cos - 1;

        if (cos2 <= thresold)   //|angle| >= 45°
            aux = 0.;
        else
        {
            aux = 1.;
            for (int k = 0; k < gc->phong->ns; k++)
                aux = aux * cos2;
        }
        phong_val = gc->phong->ka + gc->phong->Idist[(int)dist] * (gc->phong->kd * cos + gc->phong->ks * aux);
    }

    return phong_val;
}

void DDA(GraphicalContext* gc, iftMatrix* Tp0, iftVoxel p1, iftVoxel pn, float* red, float* green, float* blue)
{
    int n, k, idx;
    iftVoxel v;
    iftVoxel p;
    float max=0;
    int Dx,Dy,Dz;
    float dx=0, dy=0, dz=0, phong_val, dist;
    iftVector N;


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
        v.x=ROUND(p.x);
        v.y=ROUND(p.y);
        v.z=ROUND(p.z);
        idx = iftGetVoxelIndex(gc->scene, v);

        
        if (gc->label->val[idx] != 0)
        {
            if (gc->object[gc->label->val[idx]].visibility != 0)
            {
                dist =sqrtf((p.x-Tp0->val[0])*(p.x-Tp0->val[0])+(p.y-Tp0->val[1])*(p.y-Tp0->val[1])+(p.z-Tp0->val[2])*(p.z-Tp0->val[2]));
                
                N.x  = gc->phong->normal[gc->normal->val[idx]].x;
                N.y  = gc->phong->normal[gc->normal->val[idx]].y;
                N.z  = gc->phong->normal[gc->normal->val[idx]].z;
                phong_val = PhongShading(gc, idx, N, dist);
                *red   += phong_val * gc->object[gc->label->val[idx]].red;
                *green += phong_val * gc->object[gc->label->val[idx]].green;
                *blue  += phong_val * gc->object[gc->label->val[idx]].blue;
            }
        }
    

        p.x = p.x + dx;
        p.y = p.y + dy;
        p.z = p.z + dz;
    }

    if (*red < 0) *red = 0;
    if (*green < 0) *green = 0;
    if (*blue < 0) *blue = 0;
    if (*red > 1) *red = 1;
    if (*green > 1) *green = 1;
    if (*blue > 1) *blue = 1;

    //return (int)max;
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

    // Phong constantss
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

    //object = (ObjectAttributes *)calloc(*numberOfObjects + 1, sizeof(ObjectAttributes));
    object = (ObjectAttributes *)calloc(2 + 1, sizeof(ObjectAttributes));

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

    printf("setting\n");
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





GraphicalContext *createGC(iftImage *scene, iftImage *imageLabel, float tilt, float spin)
{
    GraphicalContext *gc;

    gc = (GraphicalContext *) calloc(1, sizeof(GraphicalContext));

    gc->scene          = iftCopyImage(scene);
    gc->phong          = createPhongModel(scene);
    
    gc->numberOfObjects= 0;
    gc->overall_opac   = 1.0;
    gc->face           = (Plane *) malloc(sizeof(Plane) * 6);
    setSceneFaces(gc);

    // computing transformations
    iftVector v1 = {.x = (float)scene->xsize / 2.0, .y = (float)scene->ysize / 2.0, .z = (float)scene->zsize / 2.0};
    iftMatrix *transMatrix1 = iftTranslationMatrix(v1);

    iftMatrix *xRotMatrix = iftRotationMatrix(IFT_AXIS_X, -tilt);
    iftMatrix *yRotMatrix = iftRotationMatrix(IFT_AXIS_Y, -spin);

    float D = sqrt(scene->xsize*scene->xsize + scene->ysize*scene->ysize);
    iftVector v2 = {.x = -(D / 2.0), .y = -(D / 2.0), .z = -(D / 2.0)};
    iftMatrix *transMatrix2 = iftTranslationMatrix(v2);

    gc->Tinv = iftMultMatricesChain(3, transMatrix1, xRotMatrix,yRotMatrix, transMatrix2);

    gc->label       = iftCopyImage(imageLabel);
    gc->object      = createObjectAttr(imageLabel, &gc->numberOfObjects);

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


iftVector columnVectorVoxelToVector(iftVoxel m)
{
    iftVector v = {.x = m.x, .y = m.y, .z = m.z};
    v.x = iftAlmostZero(v.x) ? 0.0 : v.x;
    v.y = iftAlmostZero(v.y) ? 0.0 : v.y;
    v.z = iftAlmostZero(v.z) ? 0.0 : v.z;

    return v;
}


int computeIntersection(GraphicalContext* gc, iftMatrix *Tpo, iftMatrix *Tn, iftVoxel *p1, iftVoxel *pn)
{

    //float max=, min=9999999.9;
    int i;
    p1->x=pn->x=p1->y=pn->y,p1->z=pn->z=0;
    float lambda, l1=9999999.9, ln=-9999999.9;
    float innerP = 0, innerPV=0;
    iftVector v1, v2, v3;
    iftVoxel V, P; 

    for (i = 0; i < 6; i++) {
        v1 = gc->face[i].normal;

        v2 = columnVectorMatrixToVector(Tn);
        innerP = iftVectorInnerProduct(v1,v2);

        if (fabs(innerP) > 1E-04)
        {

            V.x     = gc->face[i].pos.x - Tpo->val[0];
            V.y     = gc->face[i].pos.y - Tpo->val[1];
            V.z     = gc->face[i].pos.z - Tpo->val[2];

            v3 = columnVectorVoxelToVector(V);
            innerPV = iftVectorInnerProduct(v1, v3);
            lambda= innerPV / innerP;
            P.x = ROUND(Tpo->val[0] + lambda * Tn->val[0]);
            P.y = ROUND(Tpo->val[1] + lambda * Tn->val[1]);
            P.z = ROUND(Tpo->val[2] + lambda * Tn->val[2]);
            if (isValidPoint(gc->scene, P))
            {
                if (lambda < l1)
                    l1 = lambda;
                if (lambda > ln)
                    ln = lambda;
            }   
        }
    }

    if (l1 < ln)
    {
        printf("entrou\n");
        p1->x = ROUND(Tpo->val[0] + l1 * Tn->val[0]);
        p1->y = ROUND(Tpo->val[1] + l1 * Tn->val[1]);
        p1->z = ROUND(Tpo->val[2] + l1 * Tn->val[2]);

        pn->x = ROUND(Tpo->val[0] + ln * Tn->val[0]);
        pn->y = ROUND(Tpo->val[1] + ln * Tn->val[1]);
        pn->z = ROUND(Tpo->val[2] + ln * Tn->val[2]);

        return 1;
    }
    return 0;

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




iftImage* phongRender(GraphicalContext *gc)
{

    float diagonal = 0;
    iftImage *outputImage = NULL;
    
    iftMatrix *Mtemp,*Tpo;
    iftVoxel p1, pn;
    iftColor  RGB, YCbCr;
    float     r = 0.0, g = 0.0, b = 0.0;

    iftMatrix* vec = iftCreateMatrix(1, 4);

    iftMatrixElem(vec, 0, 0) = 0;
    iftMatrixElem(vec, 0, 1) = 0;
    iftMatrixElem(vec, 0, 2) = 1;
    iftMatrixElem(vec, 0, 3) = 0;

    iftMatrix* tVec = iftMultMatrices(gc->Tinv, vec);
    diagonal = sqrt((gc->scene->xsize * gc->scene->xsize) + (gc->scene->ysize * gc->scene->ysize) + (gc->scene->zsize * gc->scene->zsize));
    outputImage = iftCreateImage(diagonal, diagonal, 1);


    for (int p = 0; p < outputImage->n; p++)
    {
        
        Mtemp = imagePixelToMatrix(outputImage,p);
        iftMatrixElem(Mtemp, 0, 2) = diagonal/2;

        Tpo =  iftMultMatrices(gc->Tinv, Mtemp);

        if (computeIntersection(gc, Tpo, tVec, &p1, &pn))
        {
            printf("passou\n");

            DDA(gc, Tpo, p1, pn, &r, &g, &b);

            RGB.val[0]     = (int)(255.0 * r);
            RGB.val[1]     = (int)(255.0 * g);
            RGB.val[2]     = (int)(255.0 * b);
            YCbCr          = iftRGBtoYCbCr(RGB,1);

            iftSetYCbCr(outputImage, p, YCbCr);

        }
    }

    return outputImage;
}


int main(int argc, char *argv[])
{
    //char buffer[512];

    float tilt, spin;
    GraphicalContext *gc;
    char *imgFileName = iftCopyString(argv[1]);
    char *imgLabelFileName = iftCopyString(argv[2]);
    iftImage *output = NULL;


    iftImage *img = iftReadImageByExt(imgFileName);
    iftImage *imgLabel = iftReadImageByExt(imgLabelFileName);

    tilt = atof(argv[3]);
    spin = atof(argv[4]);

    gc = createGC(img, imgLabel, tilt, spin);
    output   = phongRender(gc);
    printf("Done\n");
    
    // sprintf(buffer, "data/%.1f%.1f%s", tx, ty, argv[2]);


    // iftImage *normalizedImage= iftNormalize(output,0,255);
    //
    // iftWriteImageByExt(normalizedImage, buffer);
    // iftDestroyImage(&img);
    // iftDestroyImage(&output);
    return 0;
}
