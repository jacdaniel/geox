/******************************************************************************

	Fichier  an_pacet1

	Role	 Recherche des anomalies dans une tranche

******************************************************************************/


#include <math.h>
#include <stdio.h>
#define POSITIF 0
#define NEGATIF 0xffffffff
#define RAIDE  (-99)
typedef struct {
    short x ;
    short y ;
    short z ;
    short amp ;
} scatter3d ;
typedef struct {
    short y ;
    short z ;
} scatter2d ;


static void proximite3d() ;

extern int debug_;
extern int xargc;
extern char **xargv;
void fifo_in3d(scatter3d *Cand3D,short x,short y,short z,short amp,long *indC ) ;

void MinMax_cp ( short tabsis[], char *tabpol, char *tabext, int dimx, int dimy ) ;
void MinMax ( short tabsis[], char *tabpol, char *tabext, int dimx, int dimy ) ;
void compact(short *a_in,short *a_out,short *i_out,int n_in,int *n_out) ;
void kcnx3d(char *input,int vmin,int sizemin,int nbiter) ;
static void proximite3d(char *input, scatter3d *Cand3D,short  *image_don,short *image_label,int dimx,int dimy,int dimz,int vminG,int sizemin) ;
int xargc;
char **xargv;

int main(int argc, char *argv[])
{

    int ihelp=0,vmin=20,vmax=255;
    char in[80],out[80],out_cha[80];
    int i,j,k,taille=1,sizemin = 1000, nbiter = 3 ;
    xargc = argc;
    xargv = argv;

//    debug_ = lireopt("-D","",0,"",0,"",0);
//    ihelp  = lireopt("-h","",0,"",0,"",0);
//    lireopt("","%s",in,"",0,"",0);
//    lireopt("-vmin","",0,"%d",&vmin,"",0);
//    lireopt("-sizemin","",0,"%d",&sizemin,"",0);
//    lireopt("-iter","",0,"%d",&nbiter,"",0);
    if(ihelp == 0)
        kcnx3d(in,vmin,sizemin,nbiter);

    else
    {
        printf(" usage : cnx3d input output1 output2 -vi (int)  \n");
        printf("    input1  : image1  8  bits       input  file\n");
        printf("    output1 : env_sup  16 bits       output  file\n");
        printf("    output2 : env_inf  16 bits       output  file\n");
        printf("    -vi (Def = 1000)  \n");
    }

}
void kcnx3d(char *input,int vmin,int sizemin,int nbiter)
{
    short *image_don=NULL,*buf=NULL ;
    short *image_label=NULL  ;
    scatter3d *Cand3D=NULL ;
    long  taille =0,ind,i;
    struct image *nf,*nf1,*nf2,*nfo,*nfo1,*image_();
    struct nf_fmt gfmt;
    int *ifmt,ifmt1[9],ifmt2[9];
    int dimx,dimy,dimz ;
    int dimx1,dimy1,dimz1 ;
    int ii,jj,l,CodeRetour;
    long ix,iy,iz;
    int ix1,iy1,iz1;
    struct fc_cont *nfcha ;

    nf = image_(input,"e"," ",&gfmt) ;
    ifmt = gfmt.lfmt ;
    pcomOCR = &comOCR;
    CodeRetour=c_ircmim_OCR(nf,pcomOCR);
    dimx = ifmt[4] ;
    dimy = ifmt[5] ;
    dimz = ifmt[7] ;

    taille = (long) dimx*dimy*dimz  ;
    image_don = (short * )realloc(image_don,taille*sizeof(short)) ;
    image_label = (short *) realloc(image_label,taille*sizeof(short )) ;
    buf = (short *) malloc(dimx*dimy*sizeof(short)) ;
    Cand3D = (scatter3d *) realloc(Cand3D,dimz*dimy*sizeof(scatter3d));
    for(iz=0; iz < dimz; iz ++) {

        ind = (long) (iz*dimx*dimy) ;
        /*printf(" read iz %ld: %d ind %ld \n", iz,dimz,ind) ;*/
        c_lect(nf,dimy,&(image_don[ind])) ;
    }
    printf(" fin lecture bloc \n") ;
    short  *tabc, *index_c,*tabp, *tabn;
    tabn =(short *) malloc(dimx*sizeof(short)) ;
    tabp =(short *) malloc(dimx*sizeof(short)) ;
    tabc =(short *) malloc(dimx*sizeof(short)) ;
    index_c =(short *) malloc(dimx*sizeof(short)) ;
    char *tabpol = (char *) malloc(dimx*sizeof(char)) ;
    char *tabext = (char *) malloc(dimx*sizeof(char)) ;
    int n_out,iter,k;
    for(iz=0; iz < dimz*dimy; iz++) {
        /* if(iz % dimy == 0) printf("processing iz %ld / %d \n",iz / dimy,dimz) ; */
        memset(tabpol,0,dimx ) ;
        memset(tabext,0,dimx ) ;
        ind = (long) (iz*dimx) ;
        MinMax_cp ( &(image_don[ind]),tabpol,tabext,dimx,1 ) ;
        memset(tabn,0,2*dimx) ;
        memset(tabp,0,2*dimx) ;
        short v_deb = image_don[ind] ;
        short v_fin = image_don[ind + dimx-1] ;
        for(k=0; k < dimx ; k++) {
            if(tabext[k] == -1) {
                tabn[k] = image_don[ind+k] ;
                image_don[ind+k] = -1;
            } else if(tabext[k] == 1) {
                tabp[k] = image_don[ind+k] ;
                image_don[ind+k] = 1;
            } else {
                image_don[ind+k] = 0 ;
            }
        }


        for (iter = 0; iter < nbiter; iter++) {
            compact(tabp,tabc,index_c,dimx,&n_out) ;
            memset(tabpol,0,dimx ) ;
            memset(tabext,0,dimx ) ;
            if(n_out < 3) break ;
            MinMax_cp ( tabc,tabpol,tabext,n_out,1 ) ;

            memset(tabp,0,2*dimx) ;
            for(i=1; i < n_out  -1 ; i++) {
                if(tabext[i] == 1) {
                    tabp[index_c[i]] = tabc[i] ;
                    image_don[ind+index_c[i]] = iter + 1 ;
                }

            }
            tabp[0] = v_deb ;
            tabp[dimx-1] = v_fin ;
        }
        /* on traite les minima locaux */
        for (iter = 0; iter < nbiter; iter++) {
            compact(tabn,tabc,index_c,dimx,&n_out) ;
            memset(tabpol,0,dimx ) ;
            memset(tabext,0,dimx ) ;
            if(n_out < 3) break ;
            MinMax_cp ( tabc,tabpol,tabext,n_out,1 ) ;
            memset(tabn,0,2*dimx) ;
            for(i=1; i < n_out  -1 ; i++) {
                if(tabext[i] == -1) {
                    tabn[index_c[i]] = tabc[i] ;
                    image_don[ind+index_c[i]] = -iter - 1 ;
                }
            }
            tabn[0] = v_deb ;
            tabn[dimx-1] = v_fin ;
        }
    }
    printf(" fin Regression \n") ;
    memset(image_label,0,2*taille) ;

    proximite3d(input, Cand3D,image_don,image_label,
                dimx,dimy,dimz,nbiter,sizemin);
    c_fermnf(nf) ;
}

static void proximite3d(char *input, scatter3d *Cand3D,short  *image_don,short *image_label,int dimx,int dimy,int dimz,int nbiter,int sizemin)
{
    long x,y,z,ix,iy,iz,k,l,indC0=0 ;
    long longue,posx,posy,posz;
    int test,label,label_pores;
    int iok,i;
    int xmin,xmax,ymin,ymax,zmin,zmax ;
    int isolatedLabel,nc ;
    int ipar[10],nbr=1 ;
    size_t IndexInFile=0;
    short amp ;
    long ind, indCour,indC, taille=((long) dimx) *dimy*dimz ;
    double GlobalSize = 0;
    char *present ;
    present = (char *) malloc(dimz*dimy*sizeof(char)) ;
    int   *ichn=NULL ;
    ichn = (int * ) realloc (ichn, (dimz+1) *sizeof(int));
    int ln  ;
    char *out=NULL, *out1=NULL,*in_temp = NULL ;
    out = (char *) realloc (out, 2048*sizeof(char)) ;
    out1 = (char *) realloc (out1, 2048*sizeof(char)) ;
    in_temp = (char *) realloc (in_temp, 2048*sizeof(char)) ;
    short *patch = NULL ;
    patch = (short *) realloc(patch, dimy*dimz*sizeof(short)) ;
    int nb_patch = 0 ;

    ln = strlen(input) ;
    printf( " input %s ln %d\n",input,ln) ;
    strncpy(in_temp,input,ln -3) ;
    sprintf(out,"%s.cnx",in_temp) ;
    printf(" out %s\n",out) ;
    sprintf(out1,"%s.H_cnx",in_temp) ;
    printf(" out1 %s\n",out1) ;
    FILE *fcnx, *fH_cnx ;
    FILE *pdump_File =fopen("dump.raw","w") ;
    fcnx = fopen(out,"w") ;
    fH_cnx = fopen(out1,"w") ;
    fprintf(fH_cnx, "label  xmin xmax ymin ymax zmin zmax size  Reg_Level indexInFile  \n") ;
    fprintf(fH_cnx, "dimx %d dimy %d dimz %d \n",dimx,dimy,dimz) ;
    indC0 = dimy*dimz ;
    Cand3D = (scatter3d * ) realloc(Cand3D,indC0*sizeof(scatter3d)) ;



    label=0 ;
    for(i=0; i < dimx*dimy*dimz; i++) {
        image_label[i]=0;
    }
    GlobalSize = 0;
    isolatedLabel=0;
    label=1;
    int vmin ;
    for(vmin = nbiter + 1  ; vmin >= 2 ; vmin --) {
        printf(" regres level %d\n", vmin) ;
        for(z=0; z < dimz; z++) {
            for(y=0; y < dimy; y++) {
                for(x=0; x<dimx ; x++) {
                    ind = x+y*dimx+z*dimx*dimy ;
                    if (abs(image_don[ind]) ==  vmin && image_label[ind] == 0 ) {
                        int signe =  1 ;
                        if(image_don[ind] < 0) signe = -1 ;
                        label ++ ;
                        longue = 0;
                        image_label[ind]= label ;
                        amp = image_don[ind] ;
                        indC = 0 ;
                        fifo_in3d(Cand3D,x,y,z,amp,&indC);
                        memset(present, 0, dimz*dimy) ;
                        present[z*dimy+y]=1 ;
                        iok = 1;
                        indCour = 0;
                        xmin = xmax=x;
                        ymin=ymax=y;
                        zmin=zmax=z;
                        while(iok == 1 )
                        {
                            if(indC > indC0 - 27 ) {
                                indC0 += (long) dimy*dimz ;
                                Cand3D = (scatter3d * ) realloc(Cand3D,indC0*sizeof(scatter3d)) ;
                            }
                            posx=Cand3D[indCour].x ; /* init postion de recherche en x (echantillon temps) */
                            posy=Cand3D[indCour].y ; /* init postion de recherche en y (numero de trace) */
                            posz=Cand3D[indCour].z ; /* init postion de recherche en z (indice profil) */

                            for(iz = posz-1; iz <=posz+1; iz ++ ) {
                                for (iy=posy-1; iy<=posy+1; iy++) {
                                    for (ix=posx-1; ix<=posx+1; ix++) {
                                        ind = ix + iy*dimx + iz*dimx*dimy ;
                                        if (iz >= 0 && iz < dimz && iy>=0 && iy<dimy && ix>=0 && ix<dimx &&
                                                image_label[ind] == 0 && signe*image_don[ind] ==  vmin && present[iz*dimy+iy] == 0 ) {
                                            longue++;
                                            image_label[ind]=label ;
                                            amp = image_don[ind] ;
                                            xmin = MIN(xmin,ix) ;
                                            xmax=MAX(xmax,ix);
                                            ymin = MIN(ymin,iy) ;
                                            ymax=MAX(ymax,iy);
                                            zmin = MIN(zmin,iz) ;
                                            zmax=MAX(zmax,iz);
                                            fifo_in3d(Cand3D,ix,iy,iz,amp,&indC);
                                            present[iz*dimy+iy] = 1 ;
                                        }
                                    }
                                }
                            }
                            if(indCour == (indC+1) ) iok = 0; /* plus de candidat pour continuer on arrete */
                            indCour ++;		/* on va chercher le candidat suivant qui sert de germe pour continuer la propagation */
                        }
                        if(xmin > 0 && xmax < dimx-1 &&  ymin > 0 && ymin < dimy-1 && zmin > 0 && zmax <  dimz-1) isolatedLabel ++ ;

                        if(longue >  sizemin ) {
                            printf(" label %d type %d vmin %d longue %ld indC %ld xmin %d ymin %d zmin %d xmax %d ymax %d zmax %d ",label, signe, vmin,longue,indC,xmin,ymin,zmin,xmax,ymax,zmax) ;
                            /* on filtre les trous de la composante connexe */
                            proximite2d (present, dimy,dimz,ymin,ymax,zmin,zmax) ;
                            iok = 1 ;
                            indCour = 0 ;
                            while(iok == 1 )
                            {
                                if(indC > indC0 - 27 ) {
                                    indC0 += (long) dimy*dimz ;
                                    Cand3D = (scatter3d * ) realloc(Cand3D,indC0*sizeof(scatter3d)) ;
                                }
                                posx=Cand3D[indCour].x ; /* init postion de recherche en x (echantillon temps) */
                                posy=Cand3D[indCour].y ; /* init postion de recherche en y (numero de trace) */
                                posz=Cand3D[indCour].z ; /* init postion de recherche en z (indice profil) */

                                for(iz = posz-1; iz <=posz+1; iz ++ ) {
                                    for (iy=posy-1; iy<=posy+1; iy++) {
                                        for (ix=posx-1; ix<=posx+1; ix++) {
                                            ind = ix + iy*dimx + iz*dimx*dimy ;
                                            if (iz >= 0 && iz < dimz && iy>=0 && iy<dimy && ix>=0 && ix<dimx &&
                                                    image_label[ind] == 0 && signe*image_don[ind] > 0  && present[iz*dimy+iy] == 2 ) {
                                                longue++;
                                                image_label[ind]=label ;
                                                amp = image_don[ind] ;
                                                xmin = MIN(xmin,ix) ;
                                                xmax=MAX(xmax,ix);
                                                ymin = MIN(ymin,iy) ;
                                                ymax=MAX(ymax,iy);
                                                zmin = MIN(zmin,iz) ;
                                                zmax=MAX(zmax,iz);
                                                fifo_in3d(Cand3D,ix,iy,iz,amp,&indC);
                                                present[iz*dimy+iy] = 1 ;
                                            }
                                        }
                                    }
                                }
                                if(indCour == (indC+1) ) iok = 0; /* plus de candidat pour continuer on arrete */
                                indCour ++;		/* on va chercher le candidat suivant qui sert de germe pour continuer la propagation */
                            }
                            printf(" new longue %ld \n",longue) ;
                            /* ecriture dans un contour */
                            ipar[0] = label ;
                            ipar[1] = xmin ;
                            ipar[2] = xmax ;
                            ipar[3] = ymin ;
                            ipar[4] = ymax ;
                            ipar[5] = zmin ;
                            ipar[6] = zmax ;
                            ipar[7] = indC ;
                            ipar[8] = vmin ;
                            nbr++ ;
                            fprintf(fH_cnx, "%d \t %d \t %d \t %d \t %d \t %d \t %d \t %ld \t %d \t %ld \n",label,xmin,xmax,ymin,ymax,zmin,zmax,indC,vmin,IndexInFile) ;
                            IndexInFile += indC*sizeof(scatter3d) ;
                            fwrite(Cand3D,sizeof(scatter3d), indC, fcnx) ;
                            if(nb_patch < 10) {
                                memset(patch,0,dimy*dimz*sizeof(short)) ;
                                int l ;
                                for(l=0; l < indC ; l ++) {
                                    long ind = Cand3D[l].z*dimy + Cand3D[l].y ;
                                    patch[ind] = Cand3D[l].x ;
                                }
                                fwrite(patch, sizeof(short), dimy*dimz,pdump_File) ;
                                nb_patch ++ ;
                            }
                        } else {
                            for(i=0; i < indC; i++) {
                                long cx=Cand3D[i].x ;
                                long cy=Cand3D[i].y ;
                                long cz=Cand3D[i].z ;
                                ind = cx  + cy*dimx + cz*dimx*dimy ;
                                image_label[ind]= 1 ;
                            }
                            label -- ;
                        }

                    }
                }
            }
        }
    }
    free(present) ;
    close(fcnx) ;
    close(fH_cnx) ;
}


void fifo_in3d(scatter3d *Cand3D,short x,short y,short z,short amp,long *indC )
{
    /* on met dans une file de type fifo les coordonnes bloc x,y,z d un point candidat a la propagation */
    Cand3D[*indC].x = x;
    Cand3D[*indC].y = y;
    Cand3D[*indC].z = z;
    Cand3D[*indC].amp = amp;
    (*indC)  ++ ;
}

void fifo_in2d( scatter2d *Cand2D,short y,short z,long *indC )
{
    /* on met dans une file de type fifo les coordonnes bloc x,y,z d un point candidat a la propagation */
    Cand2D[*indC].y = y;
    Cand2D[*indC].z = z;
    (*indC)  ++ ;
}
void MinMax_cp ( short tabsis[], char *tabpol, char *tabext, int dimx, int dimy )

{
    /* --------------DECLARATIONS VARIABLES LOCALES--------------------     */

    int lig;        /* Indice de parcourt de ligne entre 1 et dimy          */
    int ij; /* Indice 1D du point dans les tableaux 2D              */
    int ideb;       /* Indice de debut de plage de meme sens                */
    int idebp;      /* Indice de debut de plateau                           */
    int id1,id2; /* Difference entre deux points consecutifs                */
    int imax;       /* sens de la plage ---> 1 croissant .... 0 decroissant */
    int ifin;       /* indice fin de la plage                               */
    int ifinlig; /* indice de fin de ligne                          */
    int k;  /* indice pour les while                                */
    /* ----------------------------------------------------------------     */

    for (lig = 0; lig < dimy; lig++) {
        ideb = ij = lig * dimx;
        ifinlig = ideb + dimx-1;

        /*   saut du premier plateau */
        while (ij < ifinlig && tabsis[ij+1] == tabsis[ij]) ij+=1;

        /*  si le premier plateau n'est pas une trace raide  */
        if (ij < ifinlig - 1) {
            id1= tabsis[ij+1] - tabsis[ij];

            imax = (id1 > 0) ? 1 : 0;

            while(ij < ifinlig) {
                ij +=1;

                /*   poursuite de la pente */
                while (ij < ifinlig &&
                        id1 * (tabsis[ij+1] - tabsis[ij]) > 0) ij += 1;

                /*   fin de la pente   */
                if (ij <= ifinlig) {

                    /*  debut du plateau */
                    idebp = ij;
                    while (ij < ifinlig && tabsis[ij+1] == tabsis[ij]) ij+=1;

                    if (ij < ifinlig) {
                        id2 = tabsis[ij+1] - tabsis[ij];
                        if (id1 * id2 < 0) {

                            /*  un extrema trouve  */
                            if (id1 <  0) {
                                ifin = ij - (ij - idebp) / 2;
                                /* bset_(tabext, &ideb);*/
                                tabext[ideb] = 1 ;
                                /* bset_(tabext, &ifin); */
                                tabext[ifin] = -1 ;
                                for (k = ideb; k < ifin; k++) tabpol[k] = 1 ; /*bset_(tabpol, &k); */
                                imax = 1;
                            } else {
                                ideb = ij - (ij - idebp) / 2;
                                imax = 0;
                            }
                            id1 = id2;
                        }
                    } else {
                        if (id1  < 0) {
                            ifin = ij;
                            /*bset_(tabext, &ideb); */
                            tabext[ideb] = 1 ;
                            for (k = ideb; k <= ifin; k++) tabpol[k] = 1 ; /*bset_(tabpol, &k); */
                        }
                    }
                } else {
                    if (id1 < 0) {
                        ifin = ij;
                        /*bset_(tabext, &ideb);*/
                        tabext[ideb] = 1 ;
                        for (k = ideb; k < ifin; k++) tabpol[k] = 1 ; /*bset_(tabpol, &k); */
                    }
                }
            }
        }
    }
}
void compact(short *a_in,short *a_out,short *i_out,int n_in,int *n_out)


{   int k,n_o;

    n_o = 0;
    for(k=0; k<n_in; k++) {
        if(a_in[k] != 0.) {
            a_out[n_o] = a_in[k];
            i_out[n_o] = k;
            n_o++;
        }
    }
    *n_out = n_o;
}
void proximite2d (char *present,int  dimy,int  dimz,int yminG,int ymaxG,int zminG,int zmaxG) {
    int iy,iz,i; 
    long  indC, indCour,indC0= 0 ;
    int ymin,ymax,zmin,zmax ;
    int izG, iyG ;
    int posy,posz ;
    long ind , longue;
    int iok ;
    scatter2d  * Cand2D = NULL ;
    indC0 = (ymaxG - yminG +1)*(zmaxG-zminG + 1) ;
    Cand2D = (scatter2d  *) malloc(indC0*sizeof(scatter2d ));
    for(izG = zminG ; izG <=zmaxG ; izG ++) {
        for (iyG = yminG ; iyG <= ymaxG ; iyG ++) {
            if(present[izG*dimy + iyG] == 0) {

                longue = 0;
                indC = 0 ;
                fifo_in2d(Cand2D,iyG,izG,&indC);
                present[izG*dimy+iyG]=1 ;
                iok = 1;
                indCour = 0;
                ymin=ymax=iyG;
                zmin=zmax=izG;
                while(iok == 1 )
                {
/*
                    if(indC > indC0 - 27 ) {
                        indC0 += 1000 ;
                        Cand2D = (scatter2d * ) realloc(Cand2D,indC0*sizeof(scatter2d)) ;
                    }
*/
                    posy=Cand2D[indCour].y ; /* init postion de recherche en y (numero de trace) */
                    posz=Cand2D[indCour].z ; /* init postion de recherche en z (indice profil) */

                    for(iz = posz-1; iz <=posz+1; iz ++ ) {
                        for (iy=posy-1; iy<=posy+1; iy++) {
                            ind = iy + iz*dimy ;
                            if (iz >= zminG && iz < zmaxG && iy>=yminG && iy<ymaxG  && present[iz*dimy+iy] == 0 ) {
                                longue++;
                                ymin = MIN(ymin,iy) ;
                                ymax=MAX(ymax,iy);
                                zmin = MIN(zmin,iz) ;
                                zmax=MAX(zmax,iz);
                                fifo_in2d(Cand2D,iy,iz,&indC);
                                present[iz*dimy+iy] = 1 ;
                            }
                        }
                    }
                if(indCour == (indC+1) ) iok = 0; /* plus de candidat pour continuer on arrete */
                indCour ++;		/* on va chercher le candidat suivant qui sert de germe pour continuer la propagation */
                }
            }
            if(ymin > yminG && ymax < ymaxG && zmin > zminG && zmax <  zmaxG ) {
                for(i=0; i < indC; i++) {
                    long cy=Cand2D[i].y ;
                    long cz=Cand2D[i].z ;
                    ind = cy + cz*dimy ;
                    present[ind]= 2 ;
                }
            }
        }
    }
    free(Cand2D) ;
}
