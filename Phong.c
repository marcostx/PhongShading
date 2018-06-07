#include <stdio.h>
#include <math.h>
#include "ift.h"

#define  NUM_NORMALS    65161

#define GetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define GetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define GetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))


#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])

int RADIUS_THRESHOLD=2;
int maxDist=164;

typedef struct phong_model {
  float      ka;
  float      kd;
  float      ks;
  float      ns;
  int        ndists;
  iftVector *normal;
  float     *Idist;
} PhongModel;


typedef struct volume
{
    iftMatrix* normal;
    iftMatrix* center;
}iftVolumeFaces;

typedef struct object_attr {
  float opacity;
  float green;
  float red;
  float blue;
  char visibility;
} ObjectAttributes;


typedef struct _surface_rendering_buffers {
  float  depth;
  float  opacity;
  int    voxel;
  int    object;
} SRBuffers;

typedef struct gc {
  ObjectAttributes *object;
  int               numberOfObjects;
  float             overall_opac;
  char              proj_mode;
  PhongModel       *phong;
  iftMatrix        *Tinv;
  iftVector        vDir;
  SRBuffers        *surf_render;
  iftImage         *scene;
  iftImage         *label;
  iftVolumeFaces   *faces;
  iftImage         *normal;
  iftImage         *opacity;
} GraphicalContext;


int sign( int x ){
    if(x >= 0)
        return 1;
    return -1;
}

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

float PhongShading(GraphicalContext *gc, int p, iftVector N, float dist)
{
    float cos_,arcos;
    float cos_2, aux, phong_val = 0.0;
    float threshold=1E-05;


    cos_  =   iftVectorInnerProduct(gc->vDir, N);
    //if (cos_ > 1.0)
    //    cos_ = 1.0;

    arcos = acos(cos_);
    if (arcos <= PI/2 && arcos > 0)
    {
        cos_2 = 2 * cos_ * cos_ - 1;

        if (arcos <= PI/4 && arcos > 0){
          aux = 1.;
          for (int k = 0; k < gc->phong->ns; k++)
              aux = aux * cos_2;

        }

        else
        {
          aux = 0.;
        }

        phong_val = gc->phong->ka + gc->phong->Idist[(int)dist] * (gc->phong->kd * cos_ + gc->phong->ks * aux);
        //phong_val = gc->phong->ka+ (1-(dist/maxDist)) * (gc->phong->kd * cos_ + gc->phong->ks * aux);

        //printf("%f\n", phong_val);
    }

    return phong_val;
}




float DDA(GraphicalContext* gc, iftMatrix* Tp0, iftVector p1, iftVector pn)
{
    int n, k, idx;
    iftVoxel v;
    iftVector p;
    float max=0,J=0.0;
    int Dx,Dy,Dz;
    float dx=0, dy=0, dz=0, phong_val, dist;
    iftVector N;

    //*red = *green = *blue = 0.0;

    if (p1.x == pn.x && p1.y == pn.y && p1.z == pn.z){
        n=1;
    }
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
        v.x=(int)p.x;
        v.y=(int)p.y;
        v.z=(int)p.z;

        idx = iftGetVoxelIndex(gc->scene, v);

        if (gc->label->val[idx] != 0)
        {
            if (gc->object[gc->label->val[idx]].visibility != 0)
            {
                //return 1.;
                dist =sqrtf((p.x-Tp0->val[0])*(p.x-Tp0->val[0])+(p.y-Tp0->val[1])*(p.y-Tp0->val[1])+(p.z-Tp0->val[2])*(p.z-Tp0->val[2]));

                N.x  = -gc->phong->normal[gc->normal->val[idx]].x;
                N.y  = -gc->phong->normal[gc->normal->val[idx]].y;
                N.z  = -gc->phong->normal[gc->normal->val[idx]].z;
                phong_val = PhongShading(gc, idx, N, dist);
                J += (float) phong_val;
                break;
                // *red   += phong_val * gc->object[gc->label->val[idx]].red;
                // *green += phong_val * gc->object[gc->label->val[idx]].green;
                // *blue  += phong_val * gc->object[gc->label->val[idx]].blue;
            }
        }

        p.x = p.x + dx;
        p.y = p.y + dy;
        p.z = p.z + dz;
    }
    if (J > 1) J = 1;
    if (J < 0) J = 0;
    // if (*red < 0) *red = 0;
    // if (*green < 0) *green = 0;
    // if (*blue < 0) *blue = 0;
    // if (*red > 1) *red = 1;
    // if (*green > 1) *green = 1;
    // if (*blue > 1) *blue = 1;

    return J;
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
    //phong->ndists = (int)(2 * diagonalSize(scene) + 1);
    //phong->Idist  = (float *) malloc(phong->ndists*sizeof(float));
    phong->ndists = (int)(maxDist);
    phong->Idist  = (float *) malloc(phong->ndists*sizeof(float));
    for (int d = 0; d < phong->ndists; d++){
        phong->Idist[d] = (float) (phong->ndists - d) / (float)phong->ndists;
    }

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

iftVolumeFaces* createVF(GraphicalContext* gc){
    iftImage *scene = gc->scene;
    int Nx = scene->xsize;
    int Ny = scene->ysize;
    int Nz = scene->zsize;

    iftVolumeFaces *vf = (iftVolumeFaces *) malloc(sizeof(iftVolumeFaces) * 6);

    for (int i = 0; i < 6; i++)
    {
        vf[i].normal = iftCreateMatrix(1, 4);
        vf[i].center = iftCreateMatrix(1, 4);
    }

     // Face of Plane XY
  vf[0].normal->val[0] = 0;
  vf[0].normal->val[1] = 0;
  vf[0].normal->val[2] = -1;
  vf[0].normal->val[3] = 1;

  vf[0].center->val[0] = Nx / 2;
  vf[0].center->val[1] = Ny / 2;
  vf[0].center->val[2] = 0;
  vf[0].center->val[3] = 1;

  // Face of Plane XZ
  vf[1].normal->val[0] = 0;
  vf[1].normal->val[1] = -1;
  vf[1].normal->val[2] = 0;
  vf[1].normal->val[3] = 1;

  vf[1].center->val[0] = Nx / 2;
  vf[1].center->val[1] = 0;
  vf[1].center->val[2] = Nz / 2;
  vf[1].center->val[3] = 1;

  // Face of Plane YZ
  vf[2].normal->val[0] = -1;
  vf[2].normal->val[1] = 0;
  vf[2].normal->val[2] = 0;
  vf[2].normal->val[3] = 1;

  vf[2].center->val[0] = 0;
  vf[2].center->val[1] = Ny / 2;
  vf[2].center->val[2] = Nz / 2;
  vf[2].center->val[3] = 1;

  // Face of Opposite Plane XY
  vf[3].normal->val[0] = 0;
  vf[3].normal->val[1] = 0;
  vf[3].normal->val[2] = 1;
  vf[3].normal->val[3] = 1;

  vf[3].center->val[0] = Nx / 2;
  vf[3].center->val[1] = Ny / 2;
  vf[3].center->val[2] = Nz - 1;
  vf[3].center->val[3] = 1;

  // Face of Opposite Plane XZ
  vf[4].normal->val[0] = 0;
  vf[4].normal->val[1] = 1;
  vf[4].normal->val[2] = 0;
  vf[4].normal->val[3] = 1;

  vf[4].center->val[0] = Nx / 2;
  vf[4].center->val[1] = Ny - 1;
  vf[4].center->val[2] = Nz / 2;
  vf[4].center->val[3] = 1;

  // Face of Opposite Plane YZ
  vf[5].normal->val[0] = 1;
  vf[5].normal->val[1] = 0;
  vf[5].normal->val[2] = 0;
  vf[5].normal->val[3] = 1;

  vf[5].center->val[0] = Nx-1;
  vf[5].center->val[1] = Ny / 2;
  vf[5].center->val[2] = Nz / 2;
  vf[5].center->val[3] = 1;

  return vf;

}

int GetNormalIndex(iftVector N)
{

    int gamma, alpha, idx;

    if ((N.x == 0.0) && (N.y == 0.0) && (N.z == 0.0))
    {
        return (0);
    }

    gamma = (int)(asin(N.z) * 180.0 / PI);
    alpha = (int)(atan2(N.y, N.x) * 180.0 / PI);
    if (alpha < 0)
        alpha += 360;
    idx = ((gamma + 90) * 360) + alpha + 1;

    return idx;
}

void computeNormals(GraphicalContext* gc)
{
    iftImage   *borders;
    iftAdjRel  *A   = iftSpheric(3.0);
    float      *mag = (float *) malloc(A->n*sizeof(float));
    float      diff;
    int        i, p, q, idx;
    iftVoxel   u, v;
    iftVector  N;

    // img borders
    borders    = iftObjectBorders(gc->label, A);

    gc->normal = iftCreateImage(gc->label->xsize, gc->label->ysize, gc->label->zsize);

    A   = iftSpheric(3.0);
    for (i = 0; i < A->n; i++)
        mag[i] = sqrtf(A->dx[i] * A->dx[i] + A->dy[i] * A->dy[i] + A->dz[i] * A->dz[i]);

    for (int p=0; p < borders->n; p++)
    {

        if (borders->val[p] != 0){
            u = iftGetVoxelCoord(gc->label, p);
            N.x = N.y = N.z = 0.0;

            for (i = 1; i < A->n; i++)
            {
                v = iftGetAdjacentVoxel(A, u ,i);

                if (isValidPoint(gc->label, v))
                {
                    q = iftGetVoxelIndex(gc->label, v);

                    diff = gc->label->val[q] - gc->label->val[p];

                    N.x  += diff * A->dx[i] / mag[i];
                    N.y  += diff * A->dy[i] / mag[i];
                    N.z  += diff * A->dz[i] / mag[i];
                }
            }
            N = iftNormalizeVector(N);
            // v.x = ROUND(u.x + N.x);
            // v.y = ROUND(u.y + N.y);
            // v.z = ROUND(u.z + N.z);
            // if (isValidPoint(gc->label, v)) {
            //     q = iftGetVoxelIndex(gc->label, v);
            //
            //     if (gc->label->val[q] != 0) {
            //        N.x = -N.x;
            //        N.y = -N.y;
            //        N.z = -N.z;
            //     }
            // }
            N.x = -N.x; N.y = -N.y; N.z = -N.z;
            gc->normal->val[p] = GetNormalIndex(N);
        }
    }


    free(mag);
    iftDestroyAdjRel(&A);
}




GraphicalContext *createGC(iftImage *scene, iftImage *imageLabel, float tilt, float spin)
{
    GraphicalContext *gc;

    gc = (GraphicalContext *) calloc(1, sizeof(GraphicalContext));

    gc->scene          = iftCopyImage(scene);
    gc->phong          = createPhongModel(scene);

    gc->numberOfObjects = 0;
    gc->overall_opac    = 1.0;
    gc->faces           = createVF(gc);

    // computing transformations
    iftVector v1 = {.x = (float)scene->xsize / 2.0, .y = (float)scene->ysize / 2.0, .z = (float)scene->zsize / 2.0};
    iftMatrix *transMatrix1 = iftTranslationMatrix(v1);

    iftMatrix *xRotMatrix = iftRotationMatrix(IFT_AXIS_X, -tilt);
    iftMatrix *yRotMatrix = iftRotationMatrix(IFT_AXIS_Y, -spin);

    float D = sqrt(scene->xsize*scene->xsize + scene->ysize*scene->ysize + scene->zsize*scene->zsize);
    iftVector v2 = {.x = -(D / 2.0), .y = -(D / 2.0), .z = -(D / 2.0)};
    iftMatrix *transMatrix2 = iftTranslationMatrix(v2);

    gc->Tinv = iftMultMatricesChain(4, transMatrix1, xRotMatrix,yRotMatrix, transMatrix2);
    gc->vDir.x = 0; gc->vDir.y = 0; gc->vDir.z = -1;
    iftMatrix* vec = iftCreateMatrix(1, 4);

    iftMatrixElem(vec, 0, 0) =  0;
    iftMatrixElem(vec, 0, 1) =  0;
    iftMatrixElem(vec, 0, 2) = -1;
    iftMatrixElem(vec, 0, 3) =  0;

    vec   = iftMultMatrices(gc->Tinv, vec);

    gc->vDir.x = iftMatrixElem(vec, 0, 0); gc->vDir.y = iftMatrixElem(vec, 0, 1); gc->vDir.z = iftMatrixElem(vec, 0, 2);

    gc->label       = iftCopyImage(imageLabel);
    gc->object      = createObjectAttr(imageLabel, &gc->numberOfObjects);

    gc->opacity     = NULL; 

    computeNormals(gc);

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


int computeIntersection(GraphicalContext* gc, iftMatrix *Tpo, iftMatrix *Tn, iftVector *p1, iftVector *pn)
{

    float max=-9999999.9, min=9999999.9;
    int i;
    p1->x=pn->x=p1->y=pn->y,p1->z=pn->z=-1.;
    float lambda, l1=9999999.9, ln=-9999999.9;
    float innerP = 0, innerPV=0;
    iftMatrix* Nj = iftCreateMatrix(1,3);
    iftVector v1, v2, v3;
    iftMatrix *V = iftCreateMatrix(1, 3);
    iftVector P;

    for (i = 0; i < 6; i++) {
        iftMatrixElem(Nj, 0, 0) = gc->faces[i].normal->val[0];
        iftMatrixElem(Nj, 0, 1) = gc->faces[i].normal->val[1];
        iftMatrixElem(Nj, 0, 2) = gc->faces[i].normal->val[2];
        v1 = columnVectorMatrixToVector(Nj);
        v2 = columnVectorMatrixToVector(Tn);
        innerP = iftVectorInnerProduct(v1,v2);

        if (innerP != 0)
        {

            iftMatrixElem(V, 0, 0)     = gc->faces[i].center->val[0] - Tpo->val[0];
            iftMatrixElem(V, 0, 1)     = gc->faces[i].center->val[1] - Tpo->val[1];
            iftMatrixElem(V, 0, 2)     = gc->faces[i].center->val[2] - Tpo->val[2];


            v3 = columnVectorMatrixToVector(V);
            innerPV = iftVectorInnerProduct(v1, v3);
            lambda= (float) innerPV / innerP;
            //printf("%f\n", lambda);

            P.x = Tpo->val[0] + lambda * Tn->val[0];
            P.y = Tpo->val[1] + lambda * Tn->val[1];
            P.z = Tpo->val[2] + lambda * Tn->val[2];
            iftVoxel test;
            test.x = (int)P.x;
            test.y = (int)P.y;
            test.z = (int)P.z;
            if (isValidPoint(gc->scene, test))
            {
                if (lambda < min){
                    p1->x = P.x;
                    p1->y = P.y;
                    p1->z = P.z;
                    min = lambda;
                }
                if (lambda > max) {
                    pn->x = P.x;
                    pn->y = P.y;
                    pn->z = P.z;
                    max = lambda;
                }
            }
        }
    }
    iftDestroyMatrix(&Nj);
    iftDestroyMatrix(&V);

    if ((p1->x != -1) && (pn->x != -1) && (min<max)){
        return 1;
    }
    else
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
    float intensity;

    iftMatrix *Mtemp,*Tpo;
    iftVector p1, pn;
    iftColor  RGB, YCbCr;


    iftMatrix* vec = iftCreateMatrix(1, 4);

    iftMatrixElem(vec, 0, 0) =  0;
    iftMatrixElem(vec, 0, 1) =  0;
    iftMatrixElem(vec, 0, 2) = -1;
    iftMatrixElem(vec, 0, 3) = 0;

    iftMatrix* tVec = iftMultMatrices(gc->Tinv, vec);
    diagonal = sqrt((gc->scene->xsize * gc->scene->xsize) + (gc->scene->ysize * gc->scene->ysize) + (gc->scene->zsize * gc->scene->zsize));
    outputImage = iftCreateImage(diagonal, diagonal, 1);
    //iftSetCbCr(outputImage, 128);

    for (int p = 0; p < outputImage->n; p++)
    {
        //printf("Step : %d\n", p);
        float     r = 0.0, g = 0.0, b = 0.0;
        Mtemp = imagePixelToMatrix(outputImage,p);
        iftMatrixElem(Mtemp, 0, 2) = diagonal/2;

        Tpo =  iftMultMatrices(gc->Tinv, Mtemp);

        if (computeIntersection(gc, Tpo, tVec, &p1, &pn))
        {
            intensity = DDA(gc, Tpo, p1, pn);
            outputImage->val[p] = (int)(255.0 * intensity);
            //printf("%d %d %d\n", r,g,b);

            // RGB.val[0]     = (int)(255.0*r);
            // RGB.val[1]     = (int)(255.0*g);
            // RGB.val[2]     = (int)(255.0*b);

            // YCbCr          = RGBtoYCbCr(RGB);

            // //iftSetYCbCr(outputImage,p,YCbCr);

            // outputImage->val[p] = YCbCr.val[0];
            // outputImage->Cb[p]  = YCbCr.val[1];
            // outputImage->Cr[p]  = YCbCr.val[2];
        }
    }

    return outputImage;
}


int main(int argc, char *argv[])
{
    char buffer[512];

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

    sprintf(buffer, "data/test3.png");

    iftImage *normalizedImage= iftNormalize(output,0,255);

    iftWriteImageByExt(normalizedImage, buffer);
    iftDestroyImage(&img);
    iftDestroyImage(&output);
    return 0;
}
