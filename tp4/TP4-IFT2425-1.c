//------------------------------------------------------
// Module  : TP4-IFT2425-1.c
// Author  : Etienne Collin & Justin Villeneuve
// Date    : 2025-04-05
// Version : 1.0
// Language: C++
// Note    :
//------------------------------------------------------

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmath>
#include <new>

/************************************************************************/
/* WINDOWS                                                              */
/************************************************************************/
#include <X11/Xutil.h>

Display* display;
int screen_num;
int depth;
Window root;
Visual* visual;
GC gc;

//------------------------------------------------
// DEFINITIONS -----------------------------------
//------------------------------------------------
#define CARRE(X) ((X) * (X))

#define OUTPUT_FILE "Tp4-Img-I.pgm"
#define VIEW_PGM "open"
#define DEBUG 0

// Constants Model
#define X_1 0.0
#define Y_1 1.0
#define X_2 -1.0 / sqrt(2.0)
#define Y_2 -1.0 / 2.0
#define X_3 1.0 / sqrt(2.0)
#define Y_3 -1.0 / 2.0
#define C 0.25
#define R 0.1
#define D 0.3

#define X_1_INI 0.2   // Initial x
#define X_2_INI 0.0   // Initial x'
#define X_3_INI -1.6  // Initial y
#define X_4_INI 0.0   // Initial y'

// Constants Runge-Kutta
#define H 0.0001  // Time step size
#define T_0 0.0   // Initial time
#define T_F 30.0  // Final time
#define NB_INTERV (T_F - T_0) / H

// Constants Image
#define WIDTH 512
#define HEIGHT 512
#define MAX_X 4.0
#define MAX_Y 4.0
#define EVOL_GRAPH 3000

#define WHITE 255
#define GREYWHITE 230
#define GREY 200
#define GREYDARK 120
#define BLACK 0

//------------------------------------------------
// GLOBAL CST ------------------------------------
//------------------------------------------------
float Xmin = -(MAX_X / 2.0);
float Xmax = MAX_X / 2.0;
float Ymin = -(MAX_Y / 2.0);
float Ymax = MAX_Y / 2.0;

float xx_1 = ((WIDTH / MAX_X) * X_1) + (WIDTH / 2);
float yy_1 = (-(HEIGHT / MAX_Y) * Y_1) + (HEIGHT / 2);
float xx_2 = ((WIDTH / MAX_X) * X_2) + (WIDTH / 2);
float yy_2 = (-(HEIGHT / MAX_Y) * Y_2) + (HEIGHT / 2);
float xx_3 = ((WIDTH / MAX_X) * X_3) + (WIDTH / 2);
float yy_3 = (-(HEIGHT / MAX_Y) * Y_3) + (HEIGHT / 2);

/************************************************************************/
/* OPEN_DISPLAY()                                                       */
/************************************************************************/
int open_display() {
    if ((display = XOpenDisplay(NULL)) == NULL) {
        printf("Connection impossible\n");
        return (-1);
    } else {
        screen_num = DefaultScreen(display);
        visual = DefaultVisual(display, screen_num);
        depth = DefaultDepth(display, screen_num);
        root = RootWindow(display, screen_num);
        return 0;
    }
}

/************************************************************************/
/* FABRIQUE_WINDOW()                                                    */
/* Cette fonction crée une fenetre X et l'affiche à l'écran.            */
/************************************************************************/
Window fabrique_window(char* nom_fen, int x, int y, int width, int height, int zoom) {
    Window win;
    XSizeHints size_hints;
    XWMHints wm_hints;
    XClassHint class_hints;
    XTextProperty windowName, iconName;

    char* name = nom_fen;

    if (zoom < 0) {
        width /= -zoom;
        height /= -zoom;
    }
    if (zoom > 0) {
        width *= zoom;
        height *= zoom;
    }

    win = XCreateSimpleWindow(display, root, x, y, width, height, 1, 0, 255);

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.min_width = width;
    size_hints.min_height = height;

    XStringListToTextProperty(&name, 1, &windowName);
    XStringListToTextProperty(&name, 1, &iconName);
    wm_hints.initial_state = NormalState;
    wm_hints.input = True;
    wm_hints.flags = StateHint | InputHint;
    class_hints.res_name = nom_fen;
    class_hints.res_class = nom_fen;

    XSetWMProperties(display, win, &windowName, &iconName, NULL, 0, &size_hints, &wm_hints, &class_hints);

    gc = XCreateGC(display, win, 0, NULL);

    XSelectInput(display, win,
                 ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask | ButtonMotionMask |
                     PointerMotionHintMask | StructureNotifyMask);

    XMapWindow(display, win);
    return (win);
}

/****************************************************************************/
/* CREE_XIMAGE()                                                            */
/* Crée une XImage à partir d'un tableau de float                           */
/* L'image peut subir un zoom.                                              */
/****************************************************************************/
XImage* cree_Ximage(float** mat, int z, int length, int width) {
    int lgth, wdth, lig, col, zoom_col, zoom_lig;
    float somme;
    unsigned char pix;
    unsigned char* dat;
    XImage* imageX;

    // Zoom positif et négatif
    if (z > 0) {
        lgth = length * z;
        wdth = width * z;

        dat = (unsigned char*)malloc(lgth * (wdth * 4) * sizeof(unsigned char));
        if (dat == NULL) {
            printf("Impossible d'allouer de la memoire.");
            exit(-1);
        }

        for (lig = 0; lig < lgth; lig = lig + z) {
            for (col = 0; col < wdth; col = col + z) {
                pix = (unsigned char)mat[lig / z][col / z];
                for (zoom_lig = 0; zoom_lig < z; zoom_lig++) {
                    for (zoom_col = 0; zoom_col < z; zoom_col++) {
                        dat[((lig + zoom_lig) * wdth * 4) + ((4 * (col + zoom_col)) + 0)] = pix;
                        dat[((lig + zoom_lig) * wdth * 4) + ((4 * (col + zoom_col)) + 1)] = pix;
                        dat[((lig + zoom_lig) * wdth * 4) + ((4 * (col + zoom_col)) + 2)] = pix;
                        dat[((lig + zoom_lig) * wdth * 4) + ((4 * (col + zoom_col)) + 3)] = pix;
                    }
                }
            }
        }
    } else {
        z = -z;
        lgth = (length / z);
        wdth = (width / z);

        dat = (unsigned char*)malloc(lgth * (wdth * 4) * sizeof(unsigned char));
        if (dat == NULL) {
            printf("Impossible d'allouer de la memoire.");
            exit(-1);
        }

        for (lig = 0; lig < (lgth * z); lig = lig + z) {
            for (col = 0; col < (wdth * z); col = col + z) {
                somme = 0.0;
                for (zoom_lig = 0; zoom_lig < z; zoom_lig++) {
                    for (zoom_col = 0; zoom_col < z; zoom_col++) {
                        somme += mat[lig + zoom_lig][col + zoom_col];
                    }
                }

                somme /= (z * z);
                dat[((lig / z) * wdth * 4) + ((4 * (col / z)) + 0)] = (unsigned char)somme;
                dat[((lig / z) * wdth * 4) + ((4 * (col / z)) + 1)] = (unsigned char)somme;
                dat[((lig / z) * wdth * 4) + ((4 * (col / z)) + 2)] = (unsigned char)somme;
                dat[((lig / z) * wdth * 4) + ((4 * (col / z)) + 3)] = (unsigned char)somme;
            }
        }
    }

    imageX = XCreateImage(display, visual, depth, ZPixmap, 0, (char*)dat, wdth, lgth, 16, wdth * 4);
    return (imageX);
}

//------------------------------------------------
// FUNCTIONS -------------------------------------
//------------------------------------------------
//-------------------------//
//-- Matrice de Double ----//
//-------------------------//
//---------------------------------------------------------
// Alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* dmatrix_allocate_1d(int hsize) {
    float* matrix;
    matrix = new float[hsize];
    return matrix;
}

//----------------------------------------------------------
// Alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** dmatrix_allocate_2d(int vsize, int hsize) {
    float** matrix;
    float* imptr;

    matrix = new float*[vsize];
    imptr = new float[(hsize) * (vsize)];
    for (int i = 0; i < vsize; i++, imptr += hsize) {
        matrix[i] = imptr;
    }
    return matrix;
}

//----------------------------------------------------------
// Libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_dmatrix_1d(float* pmat) { delete[] pmat; }

//----------------------------------------------------------
// Libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_dmatrix_2d(float** pmat) {
    delete[] (pmat[0]);
    delete[] pmat;
}

//----------------------------------------------------------
// SaveImagePgm
//----------------------------------------------------------
void SaveImagePgm(char* name, float** mat, int lgth, int wdth) {
    int i, j;
    char buff[300];
    FILE* fic;

    // Extension
    strcpy(buff, name);

    // Ouverture fichier
    fic = fopen(buff, "wb");
    if (fic == NULL) {
        printf("Probleme dans la sauvegarde de %s", buff);
        exit(-1);
    }
    printf("\n Sauvegarde de %s au format pgm\n", buff);

    // Sauvegarde de l'entete
    fprintf(fic, "P5");
    fprintf(fic, "\n# IMG Module");
    fprintf(fic, "\n%d %d", wdth, lgth);
    fprintf(fic, "\n255\n");

    // Enregistrement
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) fprintf(fic, "%c", (char)mat[i][j]);
    }

    // Fermeture fichier
    fclose(fic);
}

//------------------------------------------------------------------------
// Plot_point
//
// Affiche entre x dans [-MAX_X/2  MAX_X/2]
//               y dans [-MAX_Y/2  MAX_Y/2]
//------------------------------------------------------------------------
void plot_point(float** MatPts, float** MatPict, int NbPts) {
    int x_co, y_co;
    int i, j, k;

    // Init
    for (i = 0; i < HEIGHT; i++) {
        for (j = 0; j < WIDTH; j++) {
            MatPict[i][j] = GREYWHITE;
        }
    }

    for (i = 0; i < HEIGHT; i++) {
        for (j = 0; j < WIDTH; j++) {
            if ((fabs(i - yy_1) + fabs(j - xx_1)) < 10) MatPict[i][j] = GREYDARK;
            if ((fabs(i - yy_2) + fabs(j - xx_2)) < 10) MatPict[i][j] = GREYDARK;
            if ((fabs(i - yy_3) + fabs(j - xx_3)) < 10) MatPict[i][j] = GREYDARK;
        }
    }

    // Loop
    for (k = 0; k < NbPts; k++) {
        x_co = (int)((WIDTH / MAX_X) * MatPts[k][0]);
        y_co = -(int)((HEIGHT / MAX_Y) * MatPts[k][1]);
        y_co += (HEIGHT / 2);
        x_co += (WIDTH / 2);
        if (DEBUG) {
            printf("[%d::%d]", x_co, y_co);
        }
        if ((x_co < WIDTH) && (y_co < HEIGHT) && (x_co > 0) && (y_co > 0)) {
            MatPict[y_co][x_co] = BLACK;
        }
    }
}

//------------------------------------------------------------------------
// Fill_Pict
//------------------------------------------------------------------------
void Fill_Pict(float** MatPts, float** MatPict, int PtsNumber, int NbPts) {
    int i, j;
    int x_co, y_co;
    int k, k_Init, k_End;

    // Init
    for (i = 0; i < HEIGHT; i++) {
        for (j = 0; j < WIDTH; j++) {
            if (MatPict[i][j] != GREYWHITE) {
                MatPict[i][j] = GREY;
            }
            if ((fabs(i - yy_1) + fabs(j - xx_1)) < 10) {
                MatPict[i][j] = GREYDARK;
            }
            if ((fabs(i - yy_2) + fabs(j - xx_2)) < 10) {
                MatPict[i][j] = GREYDARK;
            }
            if ((fabs(i - yy_3) + fabs(j - xx_3)) < 10) {
                MatPict[i][j] = GREYDARK;
            }
        }
    }

    // Loop
    k_Init = PtsNumber;
    k_End = (k_Init + EVOL_GRAPH) % NbPts;
    for (k = k_Init; k < k_End; k++) {
        k = (k % NbPts);
        x_co = (int)((WIDTH / MAX_X) * MatPts[k][0]);
        y_co = -(int)((HEIGHT / MAX_Y) * MatPts[k][1]);
        y_co += (HEIGHT / 2);
        x_co += (WIDTH / 2);
        if ((x_co < WIDTH) && (y_co < HEIGHT) && (x_co > 0) && (y_co > 0)) {
            MatPict[y_co][x_co] = BLACK;
        }
    }
}

//------------------------------------------------
// FONCTIONS TPs ---------------------------------
//------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// Question 1
///////////////////////////////////////////////////////////////////////////////

void f(float t, float u[4], float out[4]) {
    float x[3];
    x[0] = X_1;
    x[1] = X_2;
    x[2] = X_3;

    float y[3];
    y[0] = Y_1;
    y[1] = Y_2;
    y[2] = Y_3;

    float sum1 = 0, sum2 = 0;

    for (int i = 0; i < 3; i++) {
        float denom = powf(x[i] - u[0], 2) + powf(y[i] - u[2], 2) + powf(D, 2);  // The denominator in the sum

        sum1 += (x[i] - u[0]) / powf(denom, 3.0 / 2.0);
        sum2 += (y[i] - u[2]) / powf(denom, 3.0 / 2.0);
    }

    float f0 = u[1];
    float f1 = sum1 - R * u[1] - C * u[0];
    float f2 = u[3];
    float f3 = sum2 - R * u[3] - C * u[2];

    out[0] = f0;
    out[1] = f1;
    out[2] = f2;
    out[3] = f3;
}

// Uses the Runge Kutta Fehlberg method
// Updates the vector u
// Out is related to the f function
void calculate_next_rk_value(float t, float u[4], float out[4]) {
    // ================================================
    // Calculate k1
    // ================================================
    f(t, u, out);

    float k1[4];
    for (int i = 0; i < 4; i++) {
        k1[i] = H * out[i];
    }

    // ================================================
    // Calculate k2
    // ================================================
    float temp[4];
    for (int i = 0; i < 4; i++) {
        temp[i] = u[i] + (k1[i] / 4.0);
    }
    f(t + (H / 4.0), temp, out);

    // Compute k2
    float k2[4];
    for (int i = 0; i < 4; i++) {
        k2[i] = H * out[i];
    }

    // ================================================
    // Calculate k3
    // ================================================
    for (int i = 0; i < 4; i++) {
        temp[i] = u[i] + (3.0 / 32.0) * k1[i] + (9.0 / 32.0) * k2[i];
    }
    f(t + (3.0 / 8.0) * H, temp, out);

    float k3[4];
    for (int i = 0; i < 4; i++) {
        k3[i] = H * out[i];
    }

    // ================================================
    // Calculate k4
    // ================================================
    for (int i = 0; i < 4; i++) {
        temp[i] = u[i] + (1932.0 / 2197.0) * k1[i] - (7200.0 / 2197.0) * k2[i] + (7296.0 / 2197.0) * k3[i];
    }
    f(t + (12.0 / 13.0) * H, temp, out);

    float k4[4];
    for (int i = 0; i < 4; i++) {
        k4[i] = H * out[i];
    }

    // ================================================
    // Calculate k5
    // ================================================
    for (int i = 0; i < 4; i++) {
        temp[i] = u[i] + (439.0 / 216.0) * k1[i] - 8 * k2[i] + (3680.0 / 513.0) * k3[i] - (845.0 / 4104.0) * k4[i];
    }
    f(t + H, temp, out);

    float k5[4];
    for (int i = 0; i < 4; i++) {
        k5[i] = H * out[i];
    }

    // ================================================
    // Calculate k6
    // ================================================
    for (int i = 0; i < 4; i++) {
        temp[i] = u[i] - (8.0 / 27.0) * k1[i] + 2.0 * k2[i] - (3544.0 / 2565.0) * k3[i] + 1859.0 / 4104.0 * k4[i] -
                  (11.0 / 40.0) * k5[i];
    }

    f(t + (H / 2.0), temp, out);

    float k6[4];
    for (int i = 0; i < 4; i++) {
        k6[i] = H * out[i];
    }

    // ================================================
    // Update u
    // ================================================
    for (int i = 0; i < 4; i++) {
        u[i] = u[i] + (16.0 / 135.0) * k1[i] + (6656.0 / 12825.0) * k3[i] + (28561.0 / 56430.0) * k4[i] -
               (9.0 / 50.0) * k5[i] + (2.0 / 55.0) * k6[i];
    }
}

void fill_Mat_Pts(float** MatPts) {
    // Initial conditions of ODE
    // u[0] = x, u[1] = x', u[2] = y, u[3] = y'
    float u[4] = {X_1_INI, X_2_INI, X_3_INI, X_4_INI};

    // For f function
    float out[4] = {0.0, 0.0, 0.0, 0.0};

    // Fill MatPts matrix
    for (int k = T_0; k < (int)(NB_INTERV); k++) {
        MatPts[k][0] = u[0];  // u[0] = x
        MatPts[k][1] = u[2];  // u[2] = y

        // Updates U
        calculate_next_rk_value(k * H, u, out);
    }
}

//----------------------------------------------------------
// MAIN
//----------------------------------------------------------
int main(int argc, char** argv) {
    int i, j, k;
    int flag_graph;
    int zoom;

    XEvent ev;
    Window win_ppicture;
    XImage* x_ppicture;
    char nomfen_ppicture[100];
    char BufSystVisu[100];

    // AllocMemory
    float** MatPict = dmatrix_allocate_2d(HEIGHT, WIDTH);
    float** MatPts = dmatrix_allocate_2d((int)(NB_INTERV), 2);

    // Init
    for (i = 0; i < HEIGHT; i++) {
        for (j = 0; j < WIDTH; j++) {
            MatPict[i][j] = GREYWHITE;
        }
    }
    for (i = 0; i < 2; i++) {
        for (j = 0; j < (int)(NB_INTERV); j++) {
            MatPts[i][j] = 0.0;
        }
    }
    flag_graph = 1;
    zoom = 1;

    //---------------------------------------------------------------------
    // Question 1
    //---------------------------------------------------------------------
    fill_Mat_Pts(MatPts);

    //---------------------------------------------------------------------
    // Fin Question 1
    //---------------------------------------------------------------------

    // Affichage des Points dans MatPict
    plot_point(MatPts, MatPict, (int)(NB_INTERV));

    // Save & Visu de MatPict
    SaveImagePgm((char*)OUTPUT_FILE, MatPict, HEIGHT, WIDTH);
    strcpy(BufSystVisu, VIEW_PGM);
    strcat(BufSystVisu, " ");
    strcat(BufSystVisu, OUTPUT_FILE);
    strcat(BufSystVisu, " &");
    system(BufSystVisu);

    // Affiche Statistique
    printf("\n\n Stat:  Xmin=[%.2f] Xmax=[%.2f] Ymin=[%.2f] Ymax=[%.2f]\n", Xmin, Xmax, Ymin, Ymax);

    //---------------------------------------------------------------------
    //-------------- visu sous XWINDOW de l'évolution de MatPts -----------
    //---------------------------------------------------------------------
    if (flag_graph) {
        // Uuverture Session Graphique
        if (open_display() < 0) {
            printf(" Impossible d'ouvrir une session graphique");
        }
        sprintf(nomfen_ppicture, "Évolution du Graphe");
        win_ppicture = fabrique_window(nomfen_ppicture, 10, 10, HEIGHT, WIDTH, zoom);
        x_ppicture = cree_Ximage(MatPict, zoom, HEIGHT, WIDTH);

        printf("\n\n Pour quitter,appuyer sur la barre d'espace");
        fflush(stdout);

        // Boucle Evolution
        for (i = 0; i < HEIGHT; i++) {
            for (j = 0; j < WIDTH; j++) {
                MatPict[i][j] = GREYWHITE;
            }
        }
        for (k = 0;;) {
            k = ((k + EVOL_GRAPH) % (int)(NB_INTERV));
            Fill_Pict(MatPts, MatPict, k, (int)(NB_INTERV));
            XDestroyImage(x_ppicture);
            x_ppicture = cree_Ximage(MatPict, zoom, HEIGHT, WIDTH);
            XPutImage(display, win_ppicture, gc, x_ppicture, 0, 0, 0, 0, x_ppicture->width, x_ppicture->height);
            usleep(10000);  // Si votre machine est lente mettre un nombre plus petit
        }
    }

    // Retour
    printf("\n Fini... \n\n\n");
    return 0;
}
