//------------------------------------------------------
// Module  : Tp-IFT2425-II.1.c
// Author  : Etienne Collin & Justin Villeneuve
// Date    : 2025-01-23
// Version : 1.0
// Language: C++
// Note    :
//------------------------------------------------------

//------------------------------------------------
// FICHIERS INCLUS
//------------------------------------------------
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cstdio>
#include <iostream>
#include <new>

//------------------------------------------------
// DEFINITIONS
//------------------------------------------------
#define CARRE(X) ((X) * (X))
#define CUBE(X) ((X) * (X) * (X))

//-------------------------
// Windows
//-------------------------
// #include <X11/X.h>
// #include <X11/Xlib.h>
#include <X11/Xutil.h>
// #include <X11/Xos.h>
// #include <X11/Xatom.h>
// #include <X11/cursorfont.h>

Display *display;
int screen_num;
int depth;
Window root;
Visual *visual;
GC gc;

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
Window fabrique_window(char *nom_fen, int x, int y, int width, int height, int zoom) {
    Window win;
    XSizeHints size_hints;
    XWMHints wm_hints;
    XClassHint class_hints;
    XTextProperty windowName, iconName;

    char *name = nom_fen;

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
    wm_hints.input = true;
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
XImage *cree_Ximage(float **mat, int z, int length, int width) {
    int lgth, wdth, lig, col, zoom_col, zoom_lig;
    float somme;
    unsigned char pix;
    unsigned char *dat;
    XImage *imageX;

    // Zoom positiv else Zoom negatifv
    if (z > 0) {
        lgth = length * z;
        wdth = width * z;

        dat = (unsigned char *)malloc(lgth * (wdth * 4) * sizeof(unsigned char));
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

        dat = (unsigned char *)malloc(lgth * (wdth * 4) * sizeof(unsigned char));
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

    imageX = XCreateImage(display, visual, depth, ZPixmap, 0, (char *)dat, wdth, lgth, 16, wdth * 4);
    return (imageX);
}

//-------------------------//
//-- Matrice de Flottant --//
//-------------------------//
//----------------------------------------------------------
// Alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float *fmatrix_allocate_1d(int hsize) {
    float *matrix;
    matrix = new float[hsize];
    return matrix;
}

//----------------------------------------------------------
// Alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float **fmatrix_allocate_2d(int vsize, int hsize) {
    float **matrix;
    float *imptr;

    matrix = new float *[vsize];
    imptr = new float[(hsize) * (vsize)];
    for (int i = 0; i < vsize; i++, imptr += hsize) {
        matrix[i] = imptr;
    }
    return matrix;
}
//----------------------------------------------------------
// Libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_fmatrix_1d(float *pmat) { delete[] pmat; }

//----------------------------------------------------------
// Libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_fmatrix_2d(float **pmat) {
    delete[] (pmat[0]);
    delete[] pmat;
}

//----------------------------------------------------------
// Sauvegarde de l'image de nom <name> au format pgm
//----------------------------------------------------------
void SaveImagePgm(char *bruit, char *name, float **mat, int lgth, int wdth) {
    int i, j;
    char buff[300];
    FILE *fic;

    // Extension
    strcpy(buff, bruit);
    strcat(buff, name);
    strcat(buff, ".pgm");

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
        for (j = 0; j < wdth; j++) {
            fprintf(fic, "%c", (char)mat[i][j]);
        }
    }

    // Fermeture fichier
    fclose(fic);
}

//----------------------------------------------------------
// Recal
//----------------------------------------------------------
void Recal(float **mat, int lgth, int wdth) {
    int i, j;
    float max, min, tmp;

    // Recherche du min
    min = mat[0][0];
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) {
            if (mat[i][j] < min) {
                min = mat[i][j];
            }
        }
    }

    // Plus min
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) {
            mat[i][j] -= min;
        }
    }

    // Recherche du max
    max = mat[0][0];
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) {
            if (mat[i][j] > max) max = mat[i][j];
        }
    }

    // Recalibre la matrice
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) {
            mat[i][j] *= (255 / max);
        }
    }
}

//----------------------------------------------------------
// Egalisation Histogramme
//----------------------------------------------------------
void Egalise(float **img, int lgth, int wdth, int thresh) {
    int i, j;
    float tmp;
    float nb;
    float HistoNg[256];
    float FnctRept[256];

    // Calcul Histogramme Ng
    for (i = 0; i < 256; i++) {
        HistoNg[i] = 0.0;
    }

    nb = 0;
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) {
            tmp = img[i][j];
            if (tmp > thresh) {
                HistoNg[(int)(tmp)]++;
                nb++;
            }
        }
    }

    for (i = 0; i < 256; i++) {
        HistoNg[i] /= (float)(nb);
    }

    // Calcul Fnct Repartition
    for (i = 0; i < 256; i++) {
        FnctRept[i] = 0.0;
    }

    for (i = 0; i < 256; i++) {
        if (i > 0)
            FnctRept[i] = FnctRept[i - 1] + HistoNg[i];
        else
            FnctRept[i] = FnctRept[i];
    }

    for (i = 0; i < 256; i++) {
        FnctRept[i] = (int)((FnctRept[i] * 255) + 0.5);
    }

    // Egalise
    for (i = 0; i < lgth; i++) {
        for (j = 0; j < wdth; j++) {
            img[i][j] = FnctRept[(int)(img[i][j])];
        }
    }
}

//----------------------------------------------------------
//----------------------------------------------------------
// Auxilary functions --------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
// Takes a pixel as an argument. Returns the real part of the number c
// Associated with this pixel
double find_c_real_part(int k, int width) { return 2.0 * (k - width / 1.35) / (width - 1.0); }

// Takes a pixel as an argument. Returns the imaginary part of the number c
// Associated with this pixel
double find_c_imaginary_part(int l, int length) { return 2.0 * (l - length / 2.0) / (length - 1.0); }

double calculate_modulus(double real_part, double imaginary_part) {
    return sqrt(pow(real_part, 2) + pow(imaginary_part, 2));
}

// Function used at each iteration to test if the sequence diverges, using
// The modulus test
bool is_divergent(double real_part, double imaginary_part) {
    double modulus = calculate_modulus(real_part, imaginary_part);
    // The sequence diverges to infinity if the modulus of the number is greater than 2
    // Else cannot conclude that the sequence diverges
    return modulus > 2.0;
}

// Calculates the next real number in the sequence
// Z[n+1] = z[n]^2 + c  =>  x[n+1] = x[n]^2 - y[n]^2 + Re(c)
// X represent the real part of x_n
// Y represents the imaginary part of x_n
double calculate_next_real(double c_real, double real, double imaginary) {
    return pow(real, 2) - pow(imaginary, 2) + c_real;
}

// Calculates the next imaginary number in the sequence
// Z[n+1] = z[n]^2 + c  =>  y[n+1] = 2 * x[n] * y[n] + Im(c)
// X represent the real part of x_n
// Y represents the imaginary part of x_n
double calculate_next_imaginary(double c_imaginary, double real, double imaginary) {
    return 2.0 * real * imaginary + c_imaginary;
}

// Determines if a pixel belongs to Mandlebrot's set
bool is_in_mandelbrot(int k, int l, int length, int width, int num_of_iterations) {
    // Compute the real and imaginary parts of the number c associated with the pixel
    double c_real = find_c_real_part(k, width);
    double c_imaginary = find_c_imaginary_part(l, length);

    // Initialize the first number in the sequence
    double real = 0;
    double imaginary = 0;

    for (int _ = 0; _ < num_of_iterations; _++) {
        // Compute next number in the sequence
        double new_real = calculate_next_real(c_real, real, imaginary);
        double new_imaginary = calculate_next_imaginary(c_imaginary, real, imaginary);

        // Update the current number in the sequence
        real = new_real;
        imaginary = new_imaginary;

        // If the sequence diverges, the pixel does not belong to Mandlebrot's set
        if (is_divergent(real, imaginary)) {
            return false;
        }
    }

    // We have not found that the sequence diverges, so the pixel belongs to Mandlebrot's set
    return true;
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc, char **argv) {
    int i, j, k;
    bool flag_graph;
    int zoom;

    // Pour Xwindow
    //------------
    XEvent ev;
    Window win_ppicture;
    XImage *x_ppicture;
    char nomfen_ppicture[100];
    int length, width;

    length = width = 512;
    float **Graph2D = fmatrix_allocate_2d(length, width);
    flag_graph = true;
    zoom = 1;

    // Init
    for (i = 0; i < length; i++) {
        for (j = 0; j < width; j++) {
            Graph2D[i][j] = 0.0;
        }
    }

    //--------------------------------------------------------------------------------
    // PROGRAMME
    // ---------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    // Display mandelbrot set
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
            if (is_in_mandelbrot(j, i, length, width, 200)) {
                Graph2D[i][j] = 0.0;
            } else {
                Graph2D[i][j] = 255.0;
            }
        }
    }

    //--------------------------------------------------------------------------------
    //---------------- visu sous XWINDOW
    //---------------------------------------------
    //--------------------------------------------------------------------------------

    // Recalage-Egalise le graph
    Recal(Graph2D, length, width);
    Egalise(Graph2D, length, width, 0.0);

    if (flag_graph) {
        // Ouverture session graphique
        if (open_display() < 0) printf(" Impossible d'ouvrir une session graphique");
        sprintf(nomfen_ppicture, "Graphe : ");
        win_ppicture = fabrique_window(nomfen_ppicture, 10, 10, width, length, zoom);
        x_ppicture = cree_Ximage(Graph2D, zoom, length, width);

        // Sauvegarde
        SaveImagePgm((char *)"", (char *)"FractalMandelbrot_QII.1", Graph2D, length, width);
        printf("\n\n Pour quitter,appuyer sur la barre d'espace");
        fflush(stdout);

        // Boucle d'evenements
        while (true) {
            XNextEvent(display, &ev);
            switch (ev.type) {
                case Expose:
                    XPutImage(display, win_ppicture, gc, x_ppicture, 0, 0, 0, 0, x_ppicture->width, x_ppicture->height);
                    break;
                case KeyPress:
                    XDestroyImage(x_ppicture);
                    XFreeGC(display, gc);
                    XCloseDisplay(display);
                    flag_graph = 0;
                    break;
            }
            if (!flag_graph) break;
        }
    }

    // Retour sans probleme
    printf("\n Fini... \n\n\n");
    return 0;
}
