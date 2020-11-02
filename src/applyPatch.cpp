
#include <malloc.h>
#include <stdio.h>


typedef struct patch_param {
    long label;
    int xmin,xmax ;
    int ymin,ymax ;
    int zmin,zmax ;
    long size ;
    int reg_level ;
    long indexInfile;
} patch_param;

typedef struct scatter3d {
    short x ;
    short y ;
    short z ;
    short amp ;
} scatter3d;


int headerGetParam(FILE *ph, long ind, patch_param *praram)
{
    char ligne[4096] ;
    fseek(ph, 0, SEEK_SET);    
    fgets(ligne, 4096, ph);
    fgets(ligne, 4096, ph);
    int cont = 1;

    while ( cont )
    {
        int read = fscanf(ph, "%ld %d %d %d %d %d %d %ld %d %ld", &praram->label,
                  &praram->xmin, &praram->xmax,
                  &praram->ymin, &praram->ymax,
                  &praram->zmin, &praram->zmax,
                  &praram->size,
                  &praram->reg_level,
                  &praram->indexInfile);
        if ( read != 10 ) return 0;
        if ( praram->label == ind ) return 1;
    }
}

void dataCreate(patch_param param, short *data)
{
    char *dataFileName = "/data/PLI/DIR_PROJET/UMC-NK/DATA/3D/UMC_small/DATA/SEISMIC/seismic3d.HR_NEAR.cnx";

    scatter3d *scatter;
    scatter = (scatter3d*)calloc(param.size, sizeof(scatter3d));
    FILE *pd = fopen(dataFileName, "r");
    fseek(pd, param.indexInfile, SEEK_SET);
    fread(scatter, sizeof(scatter3d), param.size, pd);

    for (int i=0; i<700*1500; i++)data[i] = -1;
    for (int i=0; i<param.size; i++)
    {
        data[scatter[i].z*1500+scatter[i].y] = scatter[i].x;
    }
    fclose(pd);
}


int main(int argc, char **argv)
{
    char *headerFileName = "/data/PLI/DIR_PROJET/UMC-NK/DATA/3D/UMC_small/DATA/SEISMIC/seismic3d.HR_NEAR.H_cnx";
    char *dataFileName = "/data/PLI/DIR_PROJET/UMC-NK/DATA/3D/UMC_small/DATA/SEISMIC/seismic3d.HR_NEAR.cnx";

    int ind = 2249;
    patch_param param;
    FILE *ph = fopen(headerFileName, "r");
    // FILE *pd = fopen(dataFileName, "r");
    short *data = (short*)calloc(1500*700, sizeof(data));
    char outFilename[1000];

    for (long i=0; i<751078; i++)
    {
        fprintf(stderr, "idx: %d %d\n", i, 1180);
        ind = i;
        int ret = headerGetParam(ph, ind, &param);
        fprintf(stderr, "ret: %d\n %ld %d %d %d %d %d %d %ld %d %ld\n", 
        ret,
        param.label, param.xmin, param.xmax,
        param.ymin, param.ymax,
        param.zmin, param.zmax,
        param.size, param.reg_level, param.indexInfile);
        if ( ret == 0 ) continue;   
        dataCreate(param, data);
        sprintf(outFilename, "/data/PLI/jacques/cnx/surface_%.5ld.raw", i);
        FILE *p0 = fopen(outFilename, "w");
        fwrite(data, sizeof(short), 1500*700, p0);
        fclose(p0);
    }
    fclose(ph);
}


















/*
extern "C" {
#define Linux 1
#include "image.h"
#include "comOCR.h"
}
extern "C" {
#define Linux 1
    int c_wrcmim_OCR(struct image* im,struct com_OCR* comOCR);
    int c_ircmim_OCR(struct image* im,struct com_OCR* comOCR);
    struct image *image_(char *nom,char *mode,char *verif,struct nf_fmt *gfmt);
    //int lireopt(char opt,char fopt,char varopt,char f1,char v1,char f2,char v2) ;
    int inr_xargc = 0;
    char** inr_xargv = NULL;
}
void  getPatch(short* rgt, long dimx, long dimy, long dimz, short* patch,short *outPatch);
struct com_OCR  comOCR ;
struct com_OCR  *pcomOCR ;


typedef struct RgtSeed {
    long x, y, z;
    double rgtValue;
} RgtSeed;
template<typename DataType, typename PatchType>
std::vector<RgtSeed> getSeedsFromPatch(const DataType* rgt, long dimx, long dimy, long dimz, const PatchType* patchTab, PatchType nonVal, bool filterSeeds=true, long seedMaxOutNumber=1000) ;
typedef struct patch_param {
    long label;
    int xmin,xmax ;
    int ymin,ymax ;
    int zmin,zmax ;
    long size ;
    int reg_level ;
    long indexInfile;
} patch_param  ;
typedef struct scatter3d {
    short x ;
    short y ;
    short z ;
    short amp ;
} scatter3d ;
//using namespace std;

main(int argc,char *argv[])
{
    char *rgt_in, *rgt_out, *cnx_hd, *cnx ;
    struct image *rgtIn, *rgtOut;
    struct nf_fmt *ifmt;
    char ligne[4096] ;
    ifmt = (struct nf_fmt *) malloc(sizeof(struct nf_fmt))  ;
    int CodeRetour ;
    pcomOCR = &comOCR;
    rgt_in = argv[1] ;
    cnx_hd = argv[2] ;
    cnx    = argv[3] ;
    rgt_out = argv[4] ;
    struct patch_param *cnxL = NULL ;
    int m_cnx = 100 ;
    cnxL = (struct patch_param *) realloc(cnxL, m_cnx*sizeof(patch_param)) ;
    long nb_cnx=0 ;
    FILE *pcnx_hd_File =fopen(cnx_hd,"r") ;
    FILE *pcnx_File =fopen(cnx,"r") ;
    FILE *pdump_File =fopen("dump.raw","w") ;
    fgets(ligne, 4096, pcnx_hd_File);
    fgets(ligne, 4096, pcnx_hd_File);

    while( fscanf(pcnx_hd_File, "%ld %d %d %d %d %d %d %ld %d %ld",&(cnxL[nb_cnx].label),
                  &(cnxL[nb_cnx].xmin), &(cnxL[nb_cnx].xmax),
                  &(cnxL[nb_cnx].ymin), &(cnxL[nb_cnx].ymax),
                  &(cnxL[nb_cnx].zmin), &(cnxL[nb_cnx].zmax),
                  &(cnxL[nb_cnx].size),
                  &(cnxL[nb_cnx].reg_level),
                  &(cnxL[nb_cnx].indexInfile)) == 10) {
        if(nb_cnx >= m_cnx-1) {
            m_cnx += 100 ;
            cnxL = (struct patch_param *) realloc(cnxL, m_cnx*sizeof(patch_param)) ;
        }
        nb_cnx ++ ;
    }
    rgtIn = image_(rgt_in,"e"," ",ifmt);
    CodeRetour=c_ircmim_OCR(rgtIn,pcomOCR);
    rgtOut = image_(rgt_out,"c","a",ifmt);
    CodeRetour=c_wrcmim_OCR(rgtOut,pcomOCR);
    long dimx = ifmt->NDIMX ;
    long dimy = ifmt->NDIMY ;
    long dimz = ifmt->NDIMZ ;

    short *rgt1, *rgt2;
    rgt1 = (short *) malloc(dimx*dimy*dimz*sizeof(short)) ;
    rgt2 = (short *) malloc(dimx*dimy*dimz*sizeof(short)) ;
    memset(rgt2,0,dimx*dimy*dimz*sizeof(short));
    short *patch_I = NULL, *patch_O = NULL ;
    patch_I = (short *) malloc(dimy*dimz*sizeof(short)) ;
    patch_O = (short *) malloc(dimy*dimz*sizeof(short)) ;

    for(int iz=0; iz < dimz ; iz ++) {
        c_lect(rgtIn,dimy,&(rgt1[iz*dimx*dimy])) ;
    }
    struct scatter3d *tabcnx =NULL ;
    long size_tabcnx = 1000 ;
    tabcnx = (struct scatter3d *) realloc(tabcnx, size_tabcnx*sizeof(scatter3d)) ;
    for(long label = 0; label < nb_cnx ; label ++) {
	    if(cnxL[label].size > 5000) {
        printf(" label %ld / %ld \n", label + 1, nb_cnx) ;
        long ind_seek =  cnxL[label].indexInfile;
        if(cnxL[label].size > size_tabcnx) {
            size_tabcnx = cnxL[label].size ;
            tabcnx = (struct scatter3d *) realloc(tabcnx, size_tabcnx*sizeof(scatter3d)) ;
        }
        fseek(pcnx_File,ind_seek,SEEK_SET);
        fread(tabcnx,sizeof(scatter3d),cnxL[label].size,pcnx_File) ;
        memset(patch_I,0,dimy*dimz*sizeof(short)) ;
        memset(patch_O,0,dimy*dimz*sizeof(short)) ;
        for(int i=0; i < cnxL[label].size; i++) {
            long ind = tabcnx[i].z*dimy + tabcnx[i].y;
            patch_I[ind] = tabcnx[i].x ;
        }
        getPatch(rgt1, dimx, dimy, dimz,  patch_I, patch_O) ;
	//fwrite(patch_I, sizeof(short), dimy*dimz,pdump_File) ;
	//fwrite(patch_O, sizeof(short), dimy*dimz,pdump_File) ;
	//for(int j=0; j < dimy*dimz; j++) {
//		int ix =  patch_O[j] ;
//		patch_O[j] = rgt1[j*dimx + ix] ;
//	}
//	fwrite(patch_O, sizeof(short), dimy*dimz,pdump_File) ;
        for(int j=0; j < dimy*dimz; j++) {
            int x=patch_O[j] ;
            rgt2[j*dimx + x] += 0.001*cnxL[label].size ;
        }
    }
    }
    double *buff ;
    buff = (double *) malloc(dimx*sizeof(double)) ;
    for(int j=0; j < dimy*dimz; j++) {
        memset(buff,0, dimx*sizeof(double)) ;
        double som = 0 ;
        buff[0] = rgt2[j*dimx] ;
        for(int k=1; k < dimx; k++) {
            buff[k] = buff[k-1] + rgt2[j*dimx + k] ;
        }
        if(buff[dimx-1] > 0.0) {
            for(int k=0; k < dimx; k++) {
                buff[k] = 32000*buff[k]/buff[dimx-1] ;
            }
        }
        for(int k=0; k < dimx; k++) {
            rgt2[j*dimx + k] = (short) buff[k] ;
        }
    }
    for(int iz=0; iz < dimz ; iz ++) {
        c_ecr(rgtOut,dimy,&(rgt2[iz*dimx*dimy])) ;
    }
    c_fermnf(rgtIn) ;
    c_fermnf(rgtOut) ;
    fclose(pcnx_hd_File) ;
    fclose(pcnx_File) ;
    free(rgt1);
    free(rgt2) ;
    free(buff) ;
    free(patch_I) ;
    free(patch_O) ;
}

template<typename DataType, typename PatchType>
std::vector<RgtSeed> getSeedsFromPatch(const DataType* rgt, long dimx, long dimy, long dimz, const PatchType* patchTab, PatchType nonVal, bool filterSeeds, long seedMaxOutNumber) {
    int IsEdge;
    int nbseeds = 0 ;
    std::vector<RgtSeed> newSeeds;
    long iz, iy, ix, vz, vy;
    for(iz = 1; iz < dimz-1; iz +=1) {
        for(iy=1; iy < dimy-1; iy +=1) {
            if(patchTab[iz*dimy+iy] != nonVal) {
                IsEdge = 0 ;
                for (vz=-1; vz <= 1 ; vz ++) {
                    for(vy=-1; vy <=1; vy++) {
                        if(patchTab[(iz+vz)*dimy + iy+vy] == nonVal) {
                            IsEdge = 1 ;
                            nbseeds ++ ;
                            RgtSeed newSeed;
                            ix = patchTab[iz*dimy + iy] ;
                            if (ix<0) {
                                ix = 0;
                            } else if (ix>=dimx) {
                                ix = dimx-1;
                            }
                            newSeed.x = ix;
                            newSeed.y = iy;
                            newSeed.z = iz;
                            newSeed.rgtValue = rgt[newSeed.x + newSeed.y*dimx +newSeed.z*dimy*dimx];

                            //printf(" seed %d ix %d iy %d iz %d newSeed.rgtValue %lf \n",nbseeds, ix , iy,iz,newSeed.rgtValue)  ;
                            newSeeds.push_back(newSeed);
                        }
                        if(IsEdge == 1) break ;
                    }
                    if(IsEdge == 1) break ;
                }
            }
        }
    }
    int Nfiltered = seedMaxOutNumber;
    if (newSeeds.size()<Nfiltered || !filterSeeds) {
	std::sort(newSeeds.begin(), newSeeds.end(), [](const RgtSeed& a, const RgtSeed& b){
		return a.rgtValue<b.rgtValue;
	});
        return newSeeds;
    } else {
        //int maxSeedBySlice = 20;
        std::vector<RgtSeed> filteredSeeds;
        filteredSeeds.resize(Nfiltered);
        unsigned timeSeed = std::chrono::system_clock::now().time_since_epoch().count();

        std::shuffle (newSeeds.begin(), newSeeds.end(), std::default_random_engine(timeSeed));
        for (std::size_t i=0; i<Nfiltered; i++) {
            filteredSeeds[i] = newSeeds[i];
        }

	std::sort(filteredSeeds.begin(), filteredSeeds.end(), [](const RgtSeed& a, const RgtSeed& b){
		return a.rgtValue<b.rgtValue;
	});
        return filteredSeeds;
    }
}

template<typename DataType, typename PatchType>
std::pair<PatchType, PatchType> getSeedsFromPatch(DataType* rgt, long dimx, long dimy, long dimz, PatchType* inPatchTab, PatchType* outPatchTab, const std::vector<RgtSeed>& seeds, long distancePower, PatchType nonVal) {
    long Nseeds = seeds.size();
    #pragma omp parallel for schedule (dynamic)
    for (long iz = 0; iz < dimz; iz++ ) {
        std::vector<double> dist;
        dist.resize(Nseeds);
        for (long iy = 0; iy < dimy; iy++ ) {
            PatchType value;
            double somCoef = 0;
            double weightedIso = 0;
            bool seedFound = false;
            if (inPatchTab[iz*dimy + iy]!=nonVal) {
                value = inPatchTab[iz*dimy + iy];
                seedFound = true;
            } else {
                for(long i=0; (i < Nseeds) && !seedFound; i++) {
                    long distance = ((iy - seeds[i].y)*(iy-seeds[i].y) + (iz -seeds[i].z)*(iz-seeds[i].z));
                    if (distance!=0) {
                        dist[i] = 1.0 / std::pow(distance,distancePower) ;
                        somCoef += dist[i] ;
                    } else {
                        value = seeds[i].x;
                        seedFound = true;
                    }
                }
            }
            if (!seedFound) {
                long ix=0;
                for(int i=0; i < Nseeds ; i++) {
                    while ( (rgt[iz*dimy*dimx + iy*dimx + ix])  < seeds[i].rgtValue  && ix<dimx ) {
                        ix++;
                    }
                    if (ix>=dimx) {
                        ix = dimx - 1;
                    }
                    weightedIso += ix*dist[i] ;
                }
                if (somCoef==0.0f) {
                    somCoef = 1;
                }
                value = weightedIso / somCoef;

            }
            outPatchTab[iz*dimy + iy] = value;

        }
    }
}

void  getPatch(short* rgt, long dimx, long dimy, long dimz, short* patch,short *outPatch) {
    std::vector<RgtSeed> seeds =  getSeedsFromPatch<short, short>(rgt, dimx, dimy, dimz, patch, 0);
    getSeedsFromPatch<short, short>(rgt, dimx, dimy, dimz, patch, outPatch, seeds, 5, 0);
}
*/