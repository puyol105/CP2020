#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define TIME_RESOLUTION 1000000 // time measuring resolution (us)

int **matrix;
int **img;
double initial_time;
double clearcache[30000000];

void clearCache();
void start(void);
double stop_timer();
void prepareMatrixs(char *path, int nlines, int ncols);

int neighboursMatrix(int i, int j) {
  int sum = 0;

  // Soma todos os vizinhos de P1;
  sum = matrix[i - 1][j] + matrix[i - 1][j + 1] + matrix[i][j + 1] +
        matrix[i + 1][j + 1] + matrix[i + 1][j] + matrix[i + 1][j - 1] +
        matrix[i][j - 1] + matrix[i - 1][j - 1];

  // Verifica a condição do Algoritmo;
  if (2 <= sum && sum <= 6) {
    return 1;
  }
  return 0;
}

int neighboursImg(int i, int j) {
  int sum = 0;

  // Soma todos os vizinhos de P1;
  sum = img[i - 1][j] + img[i - 1][j + 1] + img[i][j + 1] + img[i + 1][j + 1] +
        img[i + 1][j] + img[i + 1][j - 1] + img[i][j - 1] + img[i - 1][j - 1];

  // Verifica a condição do Algoritmo;
  if (2 <= sum && sum <= 6) {
    return 1;
  }
  return 0;
}

int transitionsMatrix(int i, int j) {

  int count = 0, arr[9], n = 9;

  // Guarda no array os vizinhos de P1.
  arr[0] = matrix[i - 1][j];     // p2
  arr[1] = matrix[i - 1][j + 1]; // p3
  arr[2] = matrix[i][j + 1];
  arr[3] = matrix[i + 1][j + 1];
  arr[4] = matrix[i + 1][j];
  arr[5] = matrix[i + 1][j - 1];
  arr[6] = matrix[i][j - 1];
  arr[7] = matrix[i - 1][j - 1];
  arr[8] = matrix[i - 1][j]; // p2

  // Verifica as transições que existem.
  for (i = 0; i < n - 1; i++) {
    if (arr[i] == 0 && arr[i + 1] == 1)
      count++;
  }

  // Verifica a condição do Algoritmo;
  if (count == 1)
    return 1;

  return 0;
}

int transitionsImg(int i, int j) {

  int count = 0, arr[9], n = 9;

  // Guarda no array os vizinhos de P1.
  arr[0] = img[i - 1][j];     // p2
  arr[1] = img[i - 1][j + 1]; // p3
  arr[2] = img[i][j + 1];
  arr[3] = img[i + 1][j + 1];
  arr[4] = img[i + 1][j];
  arr[5] = img[i + 1][j - 1];
  arr[6] = img[i][j - 1];
  arr[7] = img[i - 1][j - 1];
  arr[8] = img[i - 1][j]; // p2

  // Verifica as transições que existem.
  for (i = 0; i < n - 1; i++) {
    if (arr[i] == 0 && arr[i + 1] == 1)
      count++;
  }

  // Verifica a condição do Algoritmo;
  if (count == 1)
    return 1;

  return 0;
}

// Função referente à iteração na matriz para gerar o skeleton;
void makeSkeleton(char *path, int nlines, int ncols) {

  int i, j, count, iter = 1;

  prepareMatrixs(path, nlines, ncols);

  // limpar cache
  clearCache();
  start();

  int endlines = nlines - 2;
  int endcols = ncols - 2;

  while (iter > 0) {
    count = 0;

#pragma omp parallel for collapse(2)
    for (i = 1; i <= endlines; i++) {
      for (j = 1; j <= endcols; j++) {

        if (matrix[i][j] == 1) {

          // Vai a cada posição da matriz e testa as condições do algoritmo,
          // caso se verifique, coloca essa posição a 0;
          //(!p4 && p6 && p2) && !(p4 && p6 && p8)

          if (!(matrix[i][j + 1] && matrix[i + 1][j] && matrix[i - 1][j]) &&
              !(matrix[i][j + 1] && matrix[i + 1][j] && matrix[i][j - 1]) &&
              neighboursMatrix(i, j) && transitionsMatrix(i, j)) {

            img[i][j] = 0;
            count++;

          } else
            img[i][j] = 1;
        } else
          img[i][j] = 0;
      }
    }

    iter = count;
    count = 0;

#pragma omp barrier
#pragma omp parallel for collapse(2)
    for (i = 1; i <= endlines; i++) {
      for (j = 1; j <= endcols; j++) {

        if (img[i][j] == 1) {

          // Vai a cada posição da matriz e testa as condições do algoritmo,
          // caso se verifique, coloca essa posição a 0;
          //!(p4 && p2 && p8) && !(p2 && p6 && p8)

          if (!(img[i][j + 1] && img[i - 1][j] && img[i][j - 1]) &&
              !(img[i - 1][j] && img[i + 1][j] && img[i][j - 1]) &&
              neighboursImg(i, j) && transitionsImg(i, j)) {
            matrix[i][j] = 0;
            count++;
          }

          else
            matrix[i][j] = 1;
        }

        else
          matrix[i][j] = 0;
      }
    }

    if (iter == 0)
      iter = count;

#pragma omp barrier
  }

  stop_timer();

  // Ficheiro para escrever o resultado final;
  FILE *f = fopen("result.txt", "w");

  // Escrever no ficheiro e dar free à matriz;
  for (i = 0; i < nlines; i++) {
    for (j = 0; j < ncols; j++) {
      fprintf(f, "%d", matrix[i][j]);
    }
    fprintf(f, "\n");
    free(matrix[i]);
    free(img[i]);
  }

  fclose(f);
}

void prepareMatrixs(char *path, int nlines, int ncols) {
  int j, i;
  int var;
  char c;

  FILE *file = fopen(path, "r");

  // Alocar memória para as linhas da matriz;
  img = malloc(nlines * sizeof(int *));
  matrix = malloc(nlines * sizeof(int *));

  // Retirar os elementos do .txt e preencher a matriz;
  char aux[ncols];
  for (i = 0; i < nlines; i++) {
    matrix[i] = malloc(ncols * sizeof(int));
    img[i] = malloc(ncols * sizeof(int));

    fscanf(file, "%s", aux);
    // printf("%s\n", aux);

    for (j = 0; j < ncols; j++) {
      c = aux[j];
      var = atoi(&c);
      img[i][j] = var;
      matrix[i][j] = var;
      //           printf("%d", img[i][j]);
    }
    //    printf("\n");
  }
  // printf("\n");

  // Fechar o ficheiro;
  fclose(file);
}

int main(int argc, char **argv) {
  int nlines = 0, ncols = 0;
  char c;

  FILE *file;

  // Abrir o ficheiro;
  file = fopen(argv[1], "r");

  // Percorre uma linha no ficheiro para descobrir o número de colunas;
  for (c = getc(file); c != '\n'; c = getc(file)) {
    ncols++;
  }
  // Percorre o ficheiro para descobrir o número de linhas;
  for (; c != EOF; c = getc(file)) {
    if (c == '\n') {
      nlines = nlines + 1;
    }
  }

  // Voltar para o inicio do ficheiro;
  fclose(file);

  // Função que aplica o algoritmo;
  makeSkeleton(argv[1], nlines, ncols);

  return 0;
}

void clearCache(void) {
  int i;
  for (i = 0; i < 30000000; ++i)
    clearcache[i] = i;
}

void start(void) {
  double time = omp_get_wtime();
  initial_time = time * TIME_RESOLUTION;
}

double stop_timer() {
  double time = omp_get_wtime();
  double final_time = time * TIME_RESOLUTION;

  printf("time: %f\n", final_time - initial_time);
  return final_time - initial_time;
}
